# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
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

=head1 NAME

Src::Perl::EnsEMBL::Tools::Anno::Support_scripts_perl::Utils - utilities for genset finalisation 

=head1 SYNOPSIS

  use Src::Perl::EnsEMBL::Tools::Anno::Support_scripts_perl::Utils  qw(sub);

  or 

  use Src::Perl::EnsEMBL::Tools::Anno::Support_scripts_perl::Utils 

  to get all methods

=head1 DESCRIPTION

=cut


package Src::Perl::EnsEMBL::Tools::Anno::Support_scripts_perl::Utils; #or just Utils??
#Make sure that the MyUtils.pm module is in the same directory or is in a directory included in your Perl library path (@INC).
use strict;
use warnings;
use feature 'say';
use Exporter;

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

use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
  check_gtf_file
  set_slice  
  create_analysis_object        
  set_transcript
  set_exon
  set_translation
  set_canonical_transcript
  parse_region_details
  set_attributes
  build_translation
  write_to_gtf_file
  build_gtf_record
  create_single_transcript_genes
  sort_features_by_slice
  run_gene_builder_for_slice
            );



=head2 check_gtf_file

Validate GTF file.

=cut

sub check_gtf_file {
    my ($gtf_file) = @_;
    die "Could not open the GTF file, path used: $gtf_file" unless -e $gtf_file;
}



=head2 set_slice

Define slice object.

=cut
sub set_slice {
    my ($region_name, $region_start, $region_end, $strand) = @_;
    my $region_length = $region_end - $region_start + 1;
    my $slice = Bio::EnsEMBL::Slice->new(
        -start             => $region_start,
        -end               => $region_end,
        -strand            => $strand,
        -seq_region_name   => $region_name,
        -seq_region_length => $region_length
    );
    return $slice;
}

=head2 create_analysis_object

Define analysis object.

=cut
sub create_analysis_object {
    my ($analysis_name, $module_name) = @_;

    my $analysis = Bio::EnsEMBL::Analysis->new(
        -logic_name => $analysis_name,
        -module     => $module_name
    );

    return $analysis;
}

=head2 set_transcript

Define transcript object.

=cut
sub set_transcript {
    my ($exons, $transcript_id, $biotype, $slice) = @_;
    my $transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $exons);
    $transcript->stable_id($transcript_id);
    $transcript->biotype($biotype);
    $transcript->slice($slice);
    return $transcript;
}

=head2 set_exon

Define exon object.

=cut
sub set_exon {
    my ($start, $end, $strand, $slice) = @_;
    my $exon = Bio::EnsEMBL::Exon->new(
        -START     => $start,
        -END       => $end,
        -STRAND    => $strand,
        -SLICE     => $slice,
        -PHASE     => -1,
        -END_PHASE => -1);
    return $exon;
}

=head2 set_translation

Define transcript object.

=cut
sub set_translation {
    my ($start_exon, $end_exon) = @_;
    my $translation = Bio::EnsEMBL::Translation->new();
    $translation->start_Exon($start_exon);
    $translation->start(1);
    $translation->end_Exon($end_exon);
    $translation->end($end_exon->length());
    return $translation;
}


=head2 set_canonical_transcript

    This function determines and sets the canonical transcript for a given gene.
    The canonical transcript is selected based on the following criteria:
    
    1. Prefer transcripts with a translation (coding transcripts).
    2. Among coding transcripts, choose the one with the longest translation.
    3. If translations are of equal length, choose the transcript with the longest total length.
    4. If no coding transcripts are available, choose the longest non-coding transcript.

    The chosen transcript is then set as the canonical transcript for the gene.

=cut

sub set_canonical_transcript {
    my ($gene) = @_;

    my $transcripts = $gene->get_all_Transcripts();
    my $current_canonical = pop(@$transcripts);

    # Iterate through each transcript and set the canonical transcript based on certain criteria
    foreach my $transcript (@$transcripts) {
        my $current_translation = $current_canonical->translation();
        my $new_translation = $transcript->translation();
        # Update the canonical transcript if:
        # - The current transcript has no translation, and the new one does
        # - The new translation is longer than the current canonical translation
        # - If translations are of equal length, but the new transcript is longer overall
        # - If both transcripts are non-coding, choose the longest one
        if (!$current_translation || ($current_translation->length() < $new_translation->length()) ||
            ($current_translation->length() == $new_translation->length() && $current_canonical->length() < $transcript->length()) ||
            !($current_canonical->translation()) && $current_canonical->length() < $transcript->length()) {
            $current_canonical = $transcript;
        }
    }

    # Set the canonical transcript for the gene
    $gene->canonical_transcript($current_canonical);
    return $gene;
}

=head2 set_attributes

Get attributes.

=cut

sub set_attributes {
    my ($attribute_string) = @_;
    my $attribute_pairs = {};

    my @attribute_array = split(";",$attribute_string);
    foreach my $attribute (@attribute_array) {
        my @pairs = split(" ",$attribute);
        if(scalar(@pairs) == 2) {
        $pairs[1] =~ s/\"//g;
        $attribute_pairs->{$pairs[0]} = $pairs[1];
        } else {
        $attribute =~ s/ //g;
        $attribute_pairs->{$attribute} = 1;
        }
    }

    return($attribute_pairs);
}

=head2 create_single_transcript_genes

Create one gene per transcript.

=cut

sub create_single_transcript_genes {
    my ($transcripts) = @_;
    # Creates single transcript genes for use with things like genebuilder
    my $single_transcript_genes = [];
    foreach my $transcript (@$transcripts) {
        my $gene = Bio::EnsEMBL::Gene->new();
        $gene->stable_id($transcript->stable_id());
        $gene->biotype($transcript->biotype);
        $gene->slice($transcript->slice());
        $gene->analysis($analysis);
        $gene->add_Transcript($transcript);
        push(@$single_transcript_genes,$gene);
    }
    return($single_transcript_genes);
}

=head2 sort_features_by_slice

    Group features by slice name

=cut

sub sort_features_by_slice {
    my ($features) = @_;

    # Create a hash to store features grouped by slice name
    my $features_by_slice = {};

    # Group features by slice name 
    #//= operator to initialize the array reference only if it doesn't already exist for a given slice name.
    push @{$features_by_slice->{$_->slice->name} //= []}, $_ for @$features;

    return $features_by_slice;
}

=head2 parse_region_details

Parse region details.

=cut
sub parse_region_details {
    my ($region_details) = @_;
    unless ($region_details =~ /(.+)\.rs(\d+)\.re(\d+)/) {
        die "Issues with the seq region details";
    }
    return ($1, $2, $3);
}


=head2 sort_transcripts_by_slice

NOT used

=cut
sub sort_transcripts_by_slice {
    my ($transcripts) = @_;
    my %sorted_transcript_hash;
    #checks if there's already an arrayref associated with the current transcript's seq_region_name.
    #If it's not defined (or evaluates to false), it assigns an empty arrayref [] to it.
    #then pushes the current transcript onto the array corresponding to its seq_region_name.
    push @{$sorted_transcript_hash{$_->seq_region_name} //= []}, $_ foreach @$transcripts;

    return \%sorted_transcript_hash;
}


=head2 build_translation

Process translation corrds to build translation and calculate exon phases.

=cut
sub build_translation {
    my ($transcript, $translation_coords) = @_;

    unless ($translation_coords =~ /(\d+):(\d+):(\d+):(\d+):(\d+):(\d+)/) {
        throw("Issue parsing translation coords, coords string:\n$translation_coords");
    }
    my ($start_exon_start, $start_exon_end, $start_exon_offset, $end_exon_start, $end_exon_end, $end_exon_offset) = ($1, $2, $3, $4, $5, $6);
    my $exons = $transcript->get_all_Exons();
    my ($start_exon, $end_exon);
    foreach my $exon (@$exons) {
        if ($exon->start() == $start_exon_start && $exon->end() == $start_exon_end) {
            $start_exon = $exon;
        }
        if ($exon->start() == $end_exon_start && $exon->end() == $end_exon_end) {
            $end_exon = $exon;
        }
    unless ($start_exon && $end_exon) {
        throw("Could not find matching start/end exon from the translation coords, translation coords:\n$translation_coords");
    }

    my $translation = set_translation($start_exon,  $end_exon);
    $translation->start($start_exon_offset);
    $translation->end($end_exon_offset);
    $transcript->translation($translation);
    calculate_exon_phases($transcript, 0);
}
}

=head2 write_to_gtf_file

    Write geneset into gtf, transcript sequences in cdna file and protein sequences in protein files.

=cut


sub write_to_gtf_file {
    my ($genes_to_write, $output_gtf_file, $analysis_name, $write_cdna_prot, $final_biotype) = @_;

    # Prepare additional output files for transcript and protein 
    my $output_transcript_seq_file = "$output_gtf_file.cdna";
    my $output_transcript_prot_file = "$output_gtf_file.prot";

    # Open output files
    open(my $out_gtf, ">", $output_gtf_file);
    #open(my $out_cdna, ">", $output_transcript_seq_file);
    #open(my $out_prot, ">", $output_transcript_prot_file);
    if ($write_cdna_prot) {
        open($out_cdna, ">", "$output_gtf_file.cdna") or die "Cannot open cDNA file: $!";
        open($out_prot, ">", "$output_gtf_file.prot") or die "Cannot open protein file: $!";
    }

    my $gene_count = 1;
    my $transcript_count = 1;

    # Iterate through each gene and its transcripts
    foreach my $gene (@$genes_to_write) {
        my $gene_id = "${analysis_name}_${gene_count}";
        my $transcripts = $gene->get_all_Transcripts();

        foreach my $transcript (@$transcripts) {
            my $transcript_id = "${analysis_name}_${transcript_count}";
            my $exons = $transcript->get_all_Exons();
            my $translation = $transcript->translation();
            my $biotype = $final_biotype // $transcript->biotype(); # Use provided biotype or default to $transcript->biotype()

            # Build and print GTF record
            my $gtf_record = build_gtf_record($transcript, $exons, $analysis_name, $gene_id, $transcript_id, $translation, $biotype);
            say $out_gtf join("\n", @$gtf_record);

            # Print transcript sequence
            #say $out_cdna ">${transcript_id}\n" . $transcript->seq->seq();

            # Print protein sequence if available
            #say $out_prot ">${transcript_id}\n" . ($translation ? $translation->seq() : "");
            # If requested, write cDNA and protein sequences
            if ($write_cdna_prot) {
                # Print transcript sequence
                say $out_cdna ">${transcript_id}\n" . $transcript->seq->seq();
                # Print protein sequence if available
                say $out_prot ">${transcript_id}\n" . ($translation ? $translation->seq() : "");
            }
            $transcript_count++;
        }

        $gene_count++;
    }

    # Close output files
    close $out_gtf;
    close $out_cdna if $write_cdna_prot;
    close $out_prot if $write_cdna_prot;
}

=head2 build_gtf_record

    Define transcript and exon records for the GTF file

=cut

sub build_gtf_record {
    my ($transcript, $exons, $analysis_name, $gene_id, $transcript_id, $translation, $transcript_biotype) = @_;

    my $record = [];
    my $strand = ($transcript->strand() == -1) ? "-" : "+";

    # Build attributes for the transcript
    my $transcript_attribs = qq{gene_id "${gene_id}"; transcript_id "${transcript_id}"; biotype "${transcript_biotype}";};
    if ($translation) {
        # Encode translation information
        my $start_exon = $translation->start_Exon();
        my $end_exon = $translation->end_Exon();
        my $translation_coords = qq{ translation_coords \"${start_exon->start()}:${start_exon->end()}:${translation->start()}:${end_exon->start()}:${end_exon->end()}:${translation->end()}\";};
        $transcript_attribs .= $translation_coords;
    }
    if($transcript->is_canonical()) {
        $transcript_attribs .= " canonical_transcript;"
    }
    # Construct GTF line for the transcript
    my @transcript_cols = ($transcript->slice->seq_region_name(), $analysis_name, 'transcript', $transcript->start(), $transcript->end(), '.', $strand, '.', $transcript_attribs);
    my $transcript_line = join("\t", @transcript_cols);
    push(@$record, $transcript_line);

    # Build attributes for exons
    my $exon_attribs_generic = qq{gene_id "${gene_id}"; transcript_id "${transcript_id}";};
    my $exon_rank = 1;

    # Construct GTF lines for each exon
    foreach my $exon (@$exons) {
        my $exon_attribs = qq{${exon_attribs_generic} exon_number "${exon_rank}";};
        my @exon_cols = ($transcript->slice->seq_region_name(), $analysis_name, 'exon', $exon->start(), $exon->end(), '.', $strand, '.', $exon_attribs);
        my $exon_line = join("\t", @exon_cols);
        push(@$record, $exon_line);
        $exon_rank++;
    }

    return $record;
}

=head2 run_gene_builder_for_slice

    Run gene builder for a slice

=cut

sub run_gene_builder_for_slice {
    my ($slice, $slice_genes, $analysis, $output_biotype) = @_;

    # Create a new GeneBuilder Runnable
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
        -query                    => $slice,
        -analysis                 => $analysis,
        -genes                    => $slice_genes,
        -output_biotype           => $output_biotype // 'gbuild',
        -max_transcripts_per_cluster => 100,
        -min_short_intron_len     => 7,
        -max_short_intron_len     => 15,
        -blessed_biotypes         => {},
        -skip_readthrough_check   => 1,
        -max_exon_length          => 50000,
        -coding_only              => 1,
    );

    # Run the GeneBuilder
    $runnable->run();

    # Return the collapsed genes for the slice
    return $runnable->output();
}



