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
  setup_fasta
  set_slice  
  create_analysis_object        
  set_transcript
  set_exon
  set_translation
  set_canonical_transcript
  parse_region_details
  set_attributes
  build_translation
            );



=head2 check_gtf_file

Validate GTF file.

=cut

sub check_gtf_file {
    my ($gtf_file) = @_;
    die "Could not open the GTF file, path used: $gtf_file" unless -e $gtf_file;
}

=head2 setup_fasta

Set up fasta file

=cut
sub setup_fasta {
    my $genome_file = @_;

    if ($genome_file && -e $genome_file) {
        setup_fasta(-FASTA => $genome_file);
    }

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

        # This array will store exactly what to write. By default, we would write the gene as constructed to this point.
        # However, if the run is just to process a particular type of transcripts into a cleaned set, we may want to
        # pull those back out into single transcript genes for further processing later. In that case, the array will
        # be replaced with single transcript genes

=cut

sub set_canonical_transcript {
    my ($gene) = @_;

    my $transcripts = $gene->get_all_Transcripts();
    my $current_canonical = pop(@$transcripts);

    # Iterate through each transcript and set the canonical transcript based on certain criteria
    foreach my $transcript (@$transcripts) {
        my $current_translation = $current_canonical->translation();
        my $new_translation = $transcript->translation();

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