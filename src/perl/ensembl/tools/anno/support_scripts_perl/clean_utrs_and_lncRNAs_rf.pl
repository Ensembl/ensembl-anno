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
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw (attach_Slice_to_Transcript attach_Analysis_to_Transcript calculate_exon_phases);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Getopt::Long qw(:config no_ignore_case);

use FindBin;
use lib "$FindBin::Bin/../conf";

use LoggerConfig;
use Utils qw(check_gtf_file setup_fasta set_slice create_analysis_object set_transcript set_exon set_translation set_canonical_transcript \
parse_region_details set_attributes build_translation);

my $analysis_name = "ensembl";
my $module_name = "Anno";

my $genome_file;
my $input_gtf_file;
my $output_gtf_file;
my $region_details;
my $compute_translations;

GetOptions(
            'input_gtf_file=s'      => \$input_gtf_file,
            'output_gtf_file=s'     => \$output_gtf_file,
            'genome_file=s'         => \$genome_file,
            'region_details=s'      => \$region_details,
            'compute_translations!' => \$compute_translations,
            'analysis_name=s'       => \$analysis_name,
            );


my $analysis = create_analysis_object($analysis_name, $module_name);
# Validate GTF file
check_gtf_file($gtf_file);
# Setup fasta
setup_fasta($genome_file);

# Parse region details
my ($region_name, $region_start, $region_end) = parse_region_details($region_details);
my $slice = set_slice($region_name, $region_start, $region_end, 1);

my %strand_conversion = ('+' => '1', '-' => '-1', '.' => undef, '?' => undef);

say "Reading GTF";
my ($genes_to_process,$genes_to_copy) =read_gtf($input_gtf_file,$slice,$analysis,%strand_conversion);
say "Genes loaded, have ".scalar(@$genes_to_process)." coding genes to process. A further ".scalar(@$genes_to_copy)." non-coding genes will be copied";

my $cleaned_genes = clean_genes($genes_to_process);
say "Finished cleaning genes. Have ".scalar(@$cleaned_genes)." cleaned genes";


=head2 read_gtf

Process gtf file and get list of transcripts

=cut

sub read_gtf {
    my ($input_gtf_file,$slice,$analysis,%strand_conversion) = @_;
    my ($genes_to_process, $genes_to_copy, $current_gene_id, $current_transcript_id, $current_translation_coords, \
        $current_biotype, $exons, $transcripts, $transcripts_by_gene_id) = ([], [], "", "", "", "", [], [], {});

    open (IN, $gtf_file) or die "Cannot open $gtf_file: $!";
    while (<IN>) {
        my $line = $_;
        my @eles = split("\t", $line);
        next unless @eles == 9;
        my ($gtf_region, $gtf_type, $gtf_start, $gtf_end, $gtf_strand) = @eles[0, 2, 3, 4, 6];
        #skip if the seq region name is different from the input one
        next unless $gtf_region eq $slice->seq_region_name();
        $gtf_strand = $strand_conversion{$gtf_strand} || die "Issue with parsing strand";
        #next unless $slice_hash->{$gtf_region};
        my $gtf_slice = $slice_hash->{$gtf_region};
        my $attributes = set_attributes($eles[8]);        
        my $gtf_gene_id         = $attributes->{'gene_id'};
        my $gtf_transcript_id   = $attributes->{'transcript_id'};
        my $gtf_biotype         = $attributes->{'biotype'};
        my $gtf_translation_coords = $attributes->{'translation_coords'} || '';
        $transcripts_by_gene_id->{$gtf_gene_id} //= [];

        if ($gtf_type eq 'transcript' && !$current_transcript_id) {
            # This is for the very first transcript
            ($current_transcript_id, $current_gene_id) = ($gtf_transcript_id, $gtf_gene_id);
            if ($gtf_translation_coords){
                $current_translation_coords = $gtf_translation_coords;
                $gtf_biotype         = $gtf_biotype ? || 'protein_coding';     
            }  
            $current_biotype = $gtf_biotype ? || 'not_set';
            next;
        } elsif ($gtf_type eq 'transcript') {
            process_transcript($exons, $current_transcript_id, $current_biotype,$transcripts_by_gene_id, $current_translation_coords);
            ($current_transcript_id, $current_gene_id, $exons) = \
            ($gtf_transcript_id, $gtf_gene_id, []);
            if ($gtf_translation_coords){
                $gtf_biotype         = $gtf_biotype ? || 'protein_coding';     
            }  
            $current_biotype = $gtf_biotype ? || 'not_set';
            $current_translation_coords = $gtf_translation_coords ? $gtf_translation_coords : "";
        } elsif ($gtf_type eq 'exon') {
            my $exon = set_exon($gtf_start, $gtf_end, $gtf_strand, $gtf_slice)
            $exon->{'region_name'} = $gtf_region;
            push @$exons, $exon;
        } else {
            die "Found an unexpected type in the GTF file, expected transcript or exon, found: $gtf_type";
        }
    }
    # Process the last transcript
    process_transcript($exons, $current_transcript_id, $current_biotype,$transcripts_by_gene_id);

    say "Finished reading GTF";
    close IN;

    my $genes_by_region = {};
    foreach my $gene_id (keys(%$transcripts_by_gene_id)) {
        my $transcripts = $transcripts_by_gene_id->{$gene_id};
        my $exons = ${$transcripts}[0]->get_all_Exons();
        my $region_name = ${$exons}[0]->{'region_name'};
        $genes_by_region->{$region_name}->{$gene_id} = $transcripts;
    }

    foreach my $gene_id (keys %{$genes_by_region->{$region_name}}) {
    my $gene_is_coding = 0;
    my $transcripts = $genes_by_region->{$region_name}->{$gene_id};

    foreach my $transcript (@$transcripts) {
        attach_Slice_to_Transcript($transcript, $slice);
        attach_Analysis_to_Transcript($transcript, $analysis);
        
        if ($transcript->{'translation_coords'} || $compute_translations) {
            ($transcript->{'translation_coords'}) ? build_translation($transcript, $transcript->{'translation_coords'}) : compute_translation($transcript);
            $gene_is_coding = 1;
        }
    }
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis($analysis);
    $gene->stable_id($gene_id);
    foreach my $transcript (@$transcripts) {
        $gene->add_Transcript($transcript);
    }
    push (@{ ($gene_is_coding) ? $genes_to_process : $genes_to_copy }, $gene);
    }

    return($genes_to_process,$genes_to_copy);
}

=head2 process_transcript

Process exons and create new transcript

=cut

sub process_transcript {
    my ($exons, $current_transcript_id, $current_biotype,$transcripts_by_gene_id, $current_translation_coords) = @_;
    return unless scalar(@$exons);

    @$exons = sort { $a->strand == 1 ? $a->start <=> $b->start : $b->start <=> $a->start } @$exons;
    my $transcript = set_transcript ($exons, $current_transcript_id, $current_biotype, $exons->[0]->slice)
    $transcript->analysis($analysis);
    $transcript->{'translation_coords'} = $current_translation_coords;
    #push @{ $transcripts_by_gene_id->{$current_gene_id} //= [] }, $transcript;
    if($transcripts_by_gene_id->{$current_gene_id}) {
        push(@{$transcripts_by_gene_id->{$current_gene_id}},$transcript);
    } else {
        $transcripts_by_gene_id->{$current_gene_id} = [$transcript];
    }
    return $transcripts_by_gene_id;
}



