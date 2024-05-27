#!/usr/bin/env perl
# Copyright [2021] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases features_overlap clone_Transcript);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(clone_Exon);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Getopt::Long qw(:config no_ignore_case);

use FindBin;
use lib "$FindBin::Bin/../conf";

use LoggerConfig;
use Utils qw(check_gtf_file setup_fasta set_slice create_analysis_object set_transcript set_exon set_translation set_canonical_transcript set_attributes \
build_translation);

my $genome_file;
my $gtf_file;
my $slice_name;
my $compute_translations;

my $analysis_name = "ensembl";
my $module_name = "Anno";

GetOptions(
            'gtf_file=s' => \$gtf_file,
            'genome_file=s' => \$genome_file,
            'slice_name=s' => \$slice_name,
            'compute_translations!' => \$compute_translations);

my $analysis = create_analysis_object($analysis_name, $module_name);
# Validate GTF file
check_gtf_file($gtf_file);
# Setup fasta
setup_fasta($genome_file);

# Read fasta file, store seq name and sequence, filter for the slice name if available 
my $slice_hash = { map { $_->seq_region_name => $_ } @{ fetch_slices_from_fasta($genome_file, $slice_name) } };
my %strand_conversion = ('+' => '1', '-' => '-1', '.' => undef, '?' => undef);

say "Reading GTF";
my $transcripts_by_gene_id = read_gtf($gtf_file, \%strand_conversion, $slice_hash);

my $genes = [];
foreach my $gene_id (keys(%$transcripts_by_gene_id)) {
    my $transcripts = $transcripts_by_gene_id->{$gene_id};
    my $new_gene = set_gene($transcripts)
    set_canonical_transcript($new_gene);
    push(@$genes,$new_gene);
}
# Write output
my $transcript_seq_file = $gtf_file.".cdna.fa";
my $transcript_canonical_seq_file = $gtf_file.".canonical.cdna.fa";
my $protein_seq_file = $gtf_file.".prot.fa";
my $protein_canonical_seq_file = $gtf_file.".canonical.prot.fa";

open(OUT1,">".$transcript_seq_file);
open(OUT2,">".$transcript_canonical_seq_file);
open(OUT3,">".$protein_seq_file);
open(OUT4,">".$protein_canonical_seq_file);
foreach my $gene (@$genes) {
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
        my $transcript_seq = $transcript->seq->seq();
        say OUT1 ">".$transcript->stable_id."\n".$transcript_seq;
        if($transcript->is_canonical()) {
        say OUT2 ">".$transcript->stable_id."\n".$transcript_seq;
        }

        if($transcript->translation()) {
        my $translation_seq = $transcript->translation->seq();
        say OUT3 ">".$transcript->stable_id."\n".$translation_seq;
        if($transcript->is_canonical()) {
            say OUT4 ">".$transcript->stable_id."\n".$translation_seq;
        }
        }
    } 
}
close OUT1;
close OUT2;
close OUT3;
close OUT4;

exit;


=head2 fetch_slices_from_fasta

Read fasta file, store seq name and sequence, filter for the slice name if available 

=cut

sub fetch_slices_from_fasta {
    my ($genome_file, $filter_slice_name) = @_;

    my @slices;
    my $current_slice_name = "";
    my $current_seq = "";

    open IN, '<', $genome_file or die "Cannot open $genome_file: $!";
    while (<IN>) {
        my $line = $_;
        chomp($line);
        if ($line =~ /^>(.+)/) {
            if ($current_slice_name) {
                
                if (!$filter_slice_name || $current_slice_name eq $filter_slice_name) {
                    push @slices,set_slice ($current_slice_name, 1, length($current_seq),1);
                }
                $current_seq = "";
            }
            $current_slice_name = $1;
        } else {
            $current_seq .= $_;
        }
    }
    close IN;
    #add the last sequence
    if ($current_seq && (!$filter_slice_name || $current_slice_name eq $filter_slice_name)) {
        push @slices,set_slice ($current_slice_name, 1, length($current_seq),1);
    }

    return \@slices;
}

=head2 read_gtf

Process gtf file and get list of transcripts

=cut

sub read_gtf {
    my ($gtf_file, $strand_conversion, $slice_hash) = @_;
    my ($current_gene_id, $current_transcript_id, $current_translation_coords, \
    $current_biotype, $exons, $transcripts, $transcripts_by_gene_id) = ("", "", "", "", [], [], {});

    open (IN, $gtf_file) or die "Cannot open $gtf_file: $!";
    while (<IN>) {
        my $line = $_;
        my @eles = split("\t", $line);
        next unless @eles == 9;
        my ($gtf_region, $gtf_type, $gtf_start, $gtf_end, $gtf_strand) = @eles[0, 2, 3, 4, 6];
        
        $gtf_strand = $strand_conversion{$gtf_strand} || die "Issue with parsing strand";
        next unless $slice_hash->{$gtf_region};
        my $gtf_slice = $slice_hash->{$gtf_region};
        my $attributes = set_attributes($eles[8]);        
        my $gtf_gene_id         = $attributes->{'gene_id'};
        my $gtf_transcript_id   = $attributes->{'transcript_id'};
        my $gtf_biotype         = $attributes->{'biotype'};
        my $gtf_translation_coords = $attributes->{'translation_coords'};
        #$transcripts_by_gene_id->{$gtf_gene_id} //= [];

        if ($gtf_type eq 'transcript' && !$current_transcript_id) {
            # This is for the very first transcript
            ($current_transcript_id, $current_gene_id) = ($gtf_transcript_id, $gtf_gene_id);
            if ($gtf_translation_coords){
                $current_translation_coords = $gtf_translation_coords;
                $gtf_biotype         = $gtf_biotype ? || 'protein_coding';     
            }  
            $current_biotype = $gtf_biotype ? || 'not_set';
        } elsif ($gtf_type eq 'transcript') {
            process_transcript($exons, $current_transcript_id, $current_biotype,$transcripts_by_gene_id,$current_translation_coords);
            ($current_transcript_id, $current_gene_id, $exons) = \
            ($gtf_transcript_id, $gtf_gene_id, []);
            if ($gtf_translation_coords){
                $gtf_biotype         = $gtf_biotype ? || 'protein_coding';     
            }  
            $current_biotype = $gtf_biotype ? || 'not_set';
            $current_translation_coords = $gtf_translation_coords ? $gtf_translation_coords : "";
        } elsif ($gtf_type eq 'exon') {
            my $exon = set_exon($gtf_start, $gtf_end, $gtf_strand, $gtf_slice)
            push @$exons, $exon;
        } else {
            die "Found an unexpected type in the GTF file, expected transcript or exon, found: $gtf_type";
        }
    }
    # Process the last transcript
    process_transcript($exons, $current_transcript_id, $current_biotype,$transcripts_by_gene_id);

    say "Finished reading GTF";
    close IN;
    return $transcripts_by_gene_id;
}

=head2 process_transcript

Process exons and create new transcript

=cut

sub process_transcript {
    my ($exons, $current_transcript_id, $current_biotype,$transcripts_by_gene_id,$current_translation_coords) = @_;
    return unless scalar(@$exons);

    @$exons = sort { $a->strand == 1 ? $a->start <=> $b->start : $b->start <=> $a->start } @$exons;
    my $transcript = set_transcript ($exons, $current_transcript_id, $current_biotype, $exons->[0]->slice)
    $transcript->analysis($analysis);
    if($current_translation_coords) {
        build_translation($transcript,$current_translation_coords);
    } elsif($compute_translations) {
        compute_translation($transcript);
    }

    push @{$transcripts_by_gene_id->{$current_gene_id}}, $transcript;

    return $transcripts_by_gene_id;
}




=head2 set_gene

Define gene object.

=cut
sub set_gene {
    my ($transcripts) = @_;
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->strand(${$transcripts}[0]->strand());
    $gene->biotype(${$transcripts}[0]->biotype());
    $gene->slice(${$transcripts}[0]->slice());
    $gene->analysis(${$transcripts}[0]->analysis());
    foreach my $transcript (@$transcripts) {
        $gene->add_Transcript($transcript);
    }
    return $gene;
}