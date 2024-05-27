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
use List::Util qw(min max);

use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(clone_Exon);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils
  qw(calculate_exon_phases clone_Transcript features_overlap);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils
  qw(contains_internal_stops compute_translation);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;

use FindBin;
use lib "$FindBin::Bin/../conf";

use LoggerConfig;
use Utils qw(check_gtf_file setup_fasta set_slice create_analysis_object set_transcript set_exon set_translation set_canonical_transcript parse_region_details);

my ($host, $port, $user, $pass, $dbname, $dna_host, $dna_port, $dna_user, $dna_dbname);
my $specify_strand;

my $analysis_name   = "anno";
my $module_name     = "anno";
my $region_details  = '';
my $good_biotype    = 'transcriptomic';
my $bad_biotype     = 'transcriptomic_flagged';
my $all_cds_exons   = 0;
my $join_transcripts = 0;
my $clean_transcripts = 0;
my $consider_cds = 0;
my $add_utr = 1;
my $write_single_transcript_genes = 1;
my $set_canonical    = 1;
my $skip_db_write    = 0;
my $genome_file      = '';
my $input_gtf_file   = '';
my $output_gtf_file  = '';
my $final_biotype    = 'not_set';

my $min_multi_exon_cds_length = 300;
my $min_single_exon_cds_length = 450;
my $min_intron_length = 75;
my $stats = {
    cds_exon_length      => 150,
    genomic_span         => 20000,
    cds_exons            => 5,
    cds_length           => 300,
    cds_ratio            => 0.85,
    small_orf_length     => 300,
    small_genomic_span   => 1000,
    single_exon_cds_length => 450,
};
GetOptions(
    'region_details=s'       => \$region_details,
    'specify_strand=s'       => \$specify_strand,
    'cds_search!'            => \$consider_cds,
    'skip_db_write!'         => \$skip_db_write,
    'input_gtf_file=s'       => \$input_gtf_file,
    'output_gtf_file=s'      => \$output_gtf_file,
    'genome_file=s'          => \$genome_file,
    'final_biotype=s'        => \$final_biotype,
    'join_transcripts!'      => \$join_transcripts,
    'clean_transcripts!'     => \$clean_transcripts,
    'all_cds_exons!'         => \$all_cds_exons
);

#say "$input_gtf_file\n$output_gtf_file\n$region_details\n";
$LoggerConfig::logger->info("$input_gtf_file\n$output_gtf_file\n$region_details\n");
# Validate GTF file
check_gtf_file($input_gtf_file);
# Setup fasta
setup_fasta($genome_file);

my $analysis = create_analysis_object($analysis_name, $module_name);
# Parse region details
my ($region_name, $region_start, $region_end) = parse_region_details($region_details);
my $slice = set_slice($region_name, $region_start, $region_end, 1);

# Process GTF file and get transcripts in the slice
my $transcripts = process_gtf_file($input_gtf_file, $region_name, $slice, $specify_strand, $slice), $analysis;

#my $transcripts_by_slice = sort_transcripts_by_slice($transcripts);
#foreach my $slice_name (keys %$transcripts_by_slice) {
#say "Processing region: $slice_name";

#my $initial_transcripts = $transcripts_by_slice->{$slice_name};
my $sorted_transcripts = remove_overlapping_exons($transcripts);

#if the condition is true builds a new array with the modified transcripts (ALWAYS DISABLED NOW)
my $cloned_transcripts = $join_transcripts ? join_transcripts([map { clone_Transcript($_->biotype('orig')) } @$sorted_transcripts]) : [];
# if the condition is true combine cloned transcripts with the original ones
my $joined_transcripts = $join_transcripts ? [@$cloned_transcripts, @$sorted_transcripts] : $sorted_transcripts;

#say "Transcript count after joining: " . scalar(@$joined_transcripts);
$LoggerConfig::logger->info("Transcript count after joining: " . scalar(@$joined_transcripts));
#say "Computing translations";
$LoggerConfig::logger->info("Computing translations");
foreach my $transcript (@$joined_transcripts) {
    if ($all_cds_exons) {
        my $exons      = $transcript->get_all_Exons();
        my $start_exon = $exons->[0];
        my $end_exon   = $exons->[-1];

        my $translation = set_translation($start_exon, $end_exon);

        $transcript->translation($translation);
        calculate_exon_phases($transcript, 0);
    } else {
        compute_translation($transcript);
    }
}

my $cleaned_transcripts;
if ($clean_transcripts) {
    $LoggerConfig::logger->info("Cleaning transcripts");
    #say "Cleaning transcripts";
    $cleaned_transcripts = clean_initial_transcripts($joined_transcripts);
} else {
    $cleaned_transcripts = $joined_transcripts;
}
#say "Transcripts post initial cleaning: " . scalar(@$cleaned_transcripts);
#say "Processing transcripts to look for additional ORFs in UTR regions";
$LoggerConfig::logger->info("Transcripts post initial cleaning: " . scalar(@$cleaned_transcripts));
$LoggerConfig::logger->info("Processing transcripts to look for additional ORFs in UTR regions");
#Process UTR or CDS
my $processed_transcripts = [];
if ($consider_cds) {
    my $processed_exon_strings = {};
    while (scalar(@$cleaned_transcripts)) {
        my $transcript = pop(@$cleaned_transcripts);
        process_transcript($transcript, $cleaned_transcripts, $processed_transcripts, $processed_exon_strings);
    }
    
    say "Transcripts after scanning UTRs for potential ORFs: " . scalar(@$processed_transcripts);
} else {
    $processed_transcripts = $cleaned_transcripts;
} 

my $initial_genes = create_single_transcript_genes($processed_transcripts);
#say "Initial genes: " . scalar(@$initial_genes);
$LoggerConfig::logger->info("Initial genes: " . scalar(@$initial_genes));
my $select_genes = select_geneset($initial_genes);
#say "Select genes: " . scalar(@$select_genes);
$LoggerConfig::logger->info("Select genes: " . scalar(@$select_genes));
my $final_genes = build_final_geneset($select_genes);
#say "Final genes: " . scalar(@$final_genes);
$LoggerConfig::logger->info("Final genes: " . scalar(@$final_genes));
my $genes_to_write = prep_genes_for_writing($final_genes, $write_single_transcript_genes, $analysis);

#say "Number of prepped genes: " . scalar(@$genes_to_write);
$LoggerConfig::logger->info("Number of prepped genes: " . scalar(@$genes_to_write));

if ($set_canonical) {
    # Iterate through each gene and set its canonical transcript
    set_canonical_transcript($_) foreach @$genes_to_write;
}

write_to_gtf_file($genes_to_write,$output_gtf_file,$analysis_name,$final_biotype);

exit;



=head1 METHODS





=head2 process_gtf_file

Read GTF file.

=cut
sub process_gtf_file {
    my ($input_gtf_file, $region_name, $region_start, $region_end, $specify_strand, $slice, $analysis) = @_;
    my %strand_conversion = ('+' => '1', '-' => '-1', '.' => undef, '?' => undef);
    #say "Reading in GTF";
    $LoggerConfig::logger->info("Reading in GTF");
    open(IN,$input_gtf_file) or die "Could not open GTF file: $input_gtf_file $!";
    my $exons = [];
    my $transcripts = [];

    my $current_transcript_id = "";
    my $current_gtf_region_name = "";
    my $gtf_region_name;
    my $new_record = 0;

    while (<IN>) {
        #chomp;
        my $current_line = $_;
        next if $current_line =~ /^\#/;
        my @columns = split("\t", $current_line);
        $gtf_region_name = $columns[0];
        #skip if the seq region name is different from the input one
        next if $gtf_region_name ne $region_name;
        my $type = $columns[2];
        if ( $type eq 'transcript' && $current_transcript_id ) {
        $new_record = 1;
        next;
        }
        next if $type ne 'exon';
        my $start  = $columns[3];
        my $end    = $columns[4];
        my $strand = $strand_conversion{$columns[6]} // 1;  
        if ( $specify_strand && $specify_strand != $strand ) {
        next;
        }
        my $phase = ( $columns[7] =~ /\./ ) ? undef : $columns[7];
        my $attributes = set_attributes( $columns[8] );
        my $gene_id       = $attributes->{'gene_id'};
        my $transcript_id = $attributes->{'transcript_id'};
        if ( $attributes->{'biotype'} ) {
        $biotype = $attributes->{'biotype'} || 'not set';
        }
        $exon = set_exon($start, $end, $strand, $slice);
        # This is weak if the transcript id is not unique
        if ($new_record) {
        $new_record = 0;
        $exons = [ sort { $a->start <=> $b->start } @{$exons} ] if $exons->[0]->strand == 1;
        $exons = [ sort { $b->start <=> $a->start } @{$exons} ] if $exons->[0]->strand != 1;
        my $transcript = set_transcript($exons, $current_transcript_id, $biotype, $slice)
        compute_translation($transcript);
        push( @{$transcripts}, $transcript );
        $exons = [];
        #it is used to deine strand for the next transcript, all the overlapping exons will be later cleaned
        push( @{$exons}, $exon ); 
        }
        else {
        push( @{$exons}, $exon );
        }
        $current_transcript_id = $transcript_id;
        $current_gtf_region_name = $gtf_region_name;
        }
        # Process remaining exons
        if ( scalar(@$exons) ) {
        $exons = [ sort { $a->start <=> $b->start } @{$exons} ] if $exons->[0]->strand == 1;
        $exons = [ sort { $b->start <=> $a->start } @{$exons} ] if $exons->[0]->strand != 1;
        my $transcript = set_transcript($exons, $current_transcript_id, $biotype, $slice)
        push( @{$transcripts}, $transcript );      
    }
    return $transcripts;
    #say "Finished reading GTF";
    $LoggerConfig::logger->info("Finished reading GTF");
}

=head2 set_attributes

Set attribute dictionary

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
        }
    }

    return($attribute_pairs);
}

=head2 remove_overlapping_exons

Identify and filter any overlapping exon.

=cut
sub remove_overlapping_exons {
    my ($transcripts) = @_;
    my $filtered_transcripts = [];

    for my $transcript (@$transcripts) {
        my $exons = $transcript->get_all_Exons();
        unless (exons_overlap($exons)) {
            push(@$filtered_transcripts, $transcript);
    }
    return $filtered_transcripts;
}
}

=head2 exons_overlap

Identy overlapping exon.

=cut
sub exons_overlap {
    my ($exons) = @_;
    for my $i (my $i=0; $i<scalar(@$exons)-1; $i++) {
        my $exon1 = $exons->[$i];
        my $exon2 = $exons->[$i + 1];
        if (features_overlap($exon1, $exon2)) {
            return 1;    # Overlapping exons found
        }
    }
    return 0;            # No overlapping exons
}

=head2 join_transcripts
    # The plan
    # First loop through all the transcripts and extend the terminal exon by the extension length
    # Then cluster on genomic overlap
    # Go through the clusters, foreach transcript figure out if there is an extension possible
    # An extension is possible if:
    # 1) There is an overlap between the two transcript on the terminal intron
    # 2) There is an overlap between two transcript on the terminal exons
    #    For this one it doesn't matter if one of the exons overlaps an intron in the order (or both)
    # Join the transcript to:
    # 1) The candidate with the most exons
    # 2) The candidate that will extend the transcript the most
    # Tricks:
    # 1) Put all the exons on the forward strand temporarily, that way only a simple sort is needed instead of strand logic
    # 2) Resort to handle 5' and 3' in the same way if possible
=cut

sub join_transcripts {
    my ($transcripts_to_join) = @_;
    my $extension_length = 250;
    adjust_transcript_ends($transcripts_to_join,$extension_length);
    # Now after extending make them into genes to cluster on genomic overlap
    my $initial_genes = create_single_transcript_genes($transcripts_to_join);
    my $joined_transcripts = process_clusters_for_joining($initial_genes);
    $joined_transcripts =  restore_strand($joined_transcripts);
    return($joined_transcripts);
}

=head2 adjust_transcript_ends

Extends exon start and end by 250 bp.

=cut

sub adjust_transcript_ends {
    my ($transcripts_to_extend, $extension_length) = @_;
    # You should put a cloning step in here if it's better to keep the unmodified transcript
    foreach my $transcript (@$transcripts_to_extend) {
        my $exons = $transcript->get_all_Exons();
        my ($start_exon, $end_exon) = ($exons->[0], $exons->[-1]);
        my $new_start = $transcript->strand() == 1
            ? $start_exon->start() - ($extension_length - 1)
            : $end_exon->start() - ($extension_length - 1);
        my $new_end = $transcript->strand() == 1
            ? $end_exon->end() + ($extension_length - 1)
            : $start_exon->end() + ($extension_length - 1);
        $new_start = 1 if $new_start < 1;
        $new_end = $transcript->slice->length() if $new_end > $transcript->slice->length();
        $start_exon->start($new_start) if $transcript->strand() == 1;
        $end_exon->end($new_end) if $transcript->strand() == 1;
        $start_exon->start($new_start) if $transcript->strand() != 1;
        $end_exon->end($new_end) if $transcript->strand() != 1;        
    }
    return $transcripts_to_extend;
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

=head2 process_clusters_for_joining

Create gene clusters and for each cluster join transcripts.

=cut

sub process_clusters_for_joining {
    my ($genes) = @_;
    my $all_joined_transcripts = [];
    my $biotypes_hash;
    #my $biotypes_hash = get_all_biotypes([@$genes]);
    #my $biotypes_array = [keys(%$biotypes_hash)];
    my $unique_biotypes = get_unique_biotypes($genes);
    #$types_hash->{genes} = [@$biotypes_array];
    my $biotypes_hash = { genes => $unique_biotypes };
    $LoggerConfig::logger->info("Clustering genes from input_dbs...");
    #say "Clustering genes from input_dbs...";
    my ($clusters, $unclustered) = cluster_Genes($genes,$biotypes_hash);
    foreach my $cluster (@$clusters) {
        my $genes = $cluster->get_Genes();
        my $extracted_transcripts = get_transcripts_from_genes($genes);
        my $joined_transcripts = calculate_joins($extracted_transcripts);
        
        push @$all_joined_transcripts, @{@$joined_transcripts};
    }
    return($all_joined_transcripts);
}

=head2 get_unique_biotypes

not sure is needed as the input set has always one biotype; 
on the other hand it makes the sub usable by different kind of inputs

=cut

sub get_unique_biotypes {
    my ($genes) = @_;
    my %biotypes_hash;
    foreach my $gene (@$genes) {
        $biotypes_hash{$gene->biotype} = 1;
    }
    my @unique_biotypes = keys %biotypes_hash;
    return \@unique_biotypes;
}

=head2 get_transcripts_from_genes

Get the list of trnascript given a gene.

=cut

sub get_transcripts_from_genes {
    my ($genes) = @_;
    my @extracted_transcripts;
    for my $gene (@$genes) {
        push @extracted_transcripts, @{ $gene->get_all_Transcripts() };
    }
    return \@extracted_transcripts;
}

=head2 place_on_forward_strand

Convert all transcripts in forward strand keeping track of the original strand.

=cut

sub place_on_forward_strand {
    my ($transcripts_to_change) = @_;
    my $forward_transcripts = [];

    foreach my $transcript (@$transcripts_to_change) {
        my $original_strand = $transcript->strand;
        $transcript->{'original_strand'} = 1;

        if ($transcript->strand == 1) {
            push(@$forward_transcripts, $transcript);
            next;
        }
        my $forward_exons = [
            map {
                set_exon(
                    $_->start(),
                    $_->end(),
                    1,               # Strand set to 1 for forward strand
                    $_->slice(),
                )
            } @{ $transcript->get_all_Exons() }
        ];
        my $sorted_forward_exons = [sort { $a->start <=> $b->start } @{$forward_exons}];
        my $forward_transcript = set_transcript($sorted_forward_exons, 
                                    $transcript->stable_id(), 
                                    $transcript->biotype,#$forward_transcript->biotype($good_biotype);
                                    $transcript->slice()
                                    )
        $forward_transcript->strand(1);
        $forward_transcript->analysis($analysis);
        $forward_transcript->{'original_strand'} = $original_strand;

        push(@$forward_transcripts, $forward_transcript);
    }
    return $forward_transcripts;
}

=head2 calculate_joins

Compare two transcripts per time and keep at the end of the array the longest one or the one with most exons.
Revert back the original strand.

=cut

sub calculate_joins {
    my ($extracted_transcripts) = @_;
    my $joined_transcripts = [];
    my $exon_string_record = {};

    # Put all the transcripts on the forward strand
    $extracted_transcripts = place_on_forward_strand($extracted_transcripts);
    # 1: Find the join that will give the most length
    # 2: Find the join that will give the most introns
    # 3: Construct both of these. If they're identical remove one
    # 4: Pass any constructed transcripts onto the end of the array
    for (my $i = 0; $i < scalar(@$extracted_transcripts); $i++) {
        my $longest_seq_length = 0;
        my $most_exons_count   = 0;
        my $best_joined_transcripts = []; #hold the two best transcripts

        my $transcript_i = $extracted_transcripts->[$i];
        my $transcript_i_exons     = $transcript_i->get_all_Exons();
        my $transcript_i_exon_string = generate_exon_string($transcript_i_exons);

        # Skip if we've already processed a transcript with the same exon string
        next if $exon_string_record->{$transcript_i_exon_string};
        $exon_string_record->{$transcript_i_exon_string} = 1;

        my $single_exon = scalar(@$transcript_i_exons) == 1;

        for (my $j = 0; $j < scalar(@$extracted_transcripts); $j++) {
            next if $i == $j;  # No point in comparing to itself

            my $transcript_j = $extracted_transcripts->[$j];

            # Skip if there's no overlap between transcripts i and j
            next unless features_overlap($transcript_i, $transcript_j);

            my $joined_transcript = multi_exon_join($transcript_i, $transcript_j);
            # Check if the joined transcrpit has the longest seq or most exons
            if ($joined_transcript) {
                my $joined_transcript_exon_count = scalar(@{$joined_transcript->get_all_Exons()});
                my $joined_transcript_length    = $joined_transcript->length();

                # Check if the joined transcript has the longest sequence length and make it the longest
                if ($joined_transcript_length > $longest_seq_length) {
                    $longest_seq_length = $joined_transcript_length;
                    $best_joined_transcripts->[0] = $joined_transcript;
                }
                # Check if the joined transcript has more exons than the current best
                if ($joined_transcript_exon_count > $most_exons_count) {
                    $most_exons_count = $joined_transcript_exon_count;
                    $best_joined_transcripts->[1] = $joined_transcript;
                }
                # Check if the joined transcript has the most exons and is longer than the current one
                if ($joined_transcript_exon_count == $most_exons_count) {
                    if ($joined_transcript_length > $best_joined_transcripts->[1]->length()) {
                        ${$best_joined_transcripts}[1] = $joined_transcript;
                    }
                }
            }
        }

        # Identify and keep only distinct transcripts among the two best candidates at the end of the array
        # Make sure to assign an id, might also need to index the introns, and recond strand of original
        my $final_joined_transcripts = [];
        if ($best_joined_transcripts->[0]{
            ${$best_joined_transcripts}[0]->{'original_strand'} = $transcript_i->{'original_strand'};
            my $exon_string1 = generate_exon_string($best_joined_transcripts->[0]->get_all_Exons());
        }
        if $best_joined_transcripts->[1]) {
            ${$best_joined_transcripts}[1]->{'original_strand'} = $transcript_i->{'original_strand'};
            my $exon_string2 = generate_exon_string($best_joined_transcripts->[1]->get_all_Exons());
        }
        if ($exon_string1 && $exon_string1 eq $exon_string2) {
                push(@$final_joined_transcripts, $best_joined_transcripts->[0]);
            } elsif($exon_string1) {
                push(@$final_joined_transcripts,${$best_joined_transcripts}[0]);
            } elsif($exon_string2) {
                push(@$final_joined_transcripts,${$best_joined_transcripts}[1]);
            }

        # Index and add the final joined transcripts
        foreach my $final_joined_transcript (@$final_joined_transcripts) {
            my $final_transcript_exons = $final_joined_transcript->get_all_Exons();
            my $final_transcript_exon_string = generate_exon_string($final_transcript_exons);

            unless ($exon_string_record->{$final_transcript_exon_string}) {
                push(@$extracted_transcripts, $final_joined_transcript);
                push(@$joined_transcripts, $final_joined_transcript);
                $exon_string_record->{$final_transcript_exon_string} = 1;
            }
        }
    }
    return $joined_transcripts;
}

=head2 generate_exon_string

Encode exon start, end, strand in string format.

=cut

sub generate_exon_string {
    my ($exons) = @_;

    my $exon_string = "";
    foreach my $exon (@$exons) {
        $exon_string .= $exon->start.":".$exon->end.":".$exon->strand;
    }

    return($exon_string);
}

#OPTIONAL paths


=head2 multi_exon_join

    # Do this in two parts
    # First examine for matches on the terminal introns, if there are then:
    # 1: If t1 has more exons past the boundary, replace the end exon of t1 with the corresponding exon from t2
    #    and add extra exons past the boundary from t2
    # 2: If there are no more exons then pick the largest exon between t1 and t2
    # Consider just looking at the intron/exon boundary
    # Keep the transcript that adds the most exons and the one that adds the most length
    # If there are no intron matches then look for exon overlap with the terminal exon
    # Note this is probably not something that needs to be different for single/multi exon
    # As all comparisons are being done and judged, just look at it from the perspective
    # of what exon from exons2 is overlaps the 5'/3' boundary exon and is closest to or
    # crosses the boundary. When you figure out what exon that is for the 5' and 3' sides
    # you adjust the start/end to match
    # Just need the exons that are on or over the start/end to decide to join. Do not need
    # introns to make the decision

=cut

sub multi_exon_join {
    my ($transcript1,$transcript2) = @_;

    my $cloned_transcript1 = $transcript1;#clone_Transcript($transcript1);
    my $exons1 = $cloned_transcript1->get_all_Exons();
    my $exons2 = $transcript2->get_all_Exons();

    my $cloned_transcipt1_exon_string = generate_exon_string($exons1);
    my $transcript2_exon_string = generate_exon_string($exons2);

    if($cloned_transcipt1_exon_string eq $transcript2_exon_string) {
        return;
    }
    my $merged_exons = [];


    my $modified_five_prime = 0;
    my $modified_three_prime = 0;

    my $five_prime_exon2s = [];
    my $three_prime_exon2s = [];

    # Work out the internal exons for transcript1
    my $internal_exon1s = [];
    for(my $i=1; $i<scalar(@$exons1)-1; $i++) {
        push(@$internal_exon1s,${$exons1}[$i]);
    }

    my $terminal_exon1s = [${$exons1}[0]];
    if(scalar(@$exons1) > 1) {
        push(@$terminal_exon1s,${$exons1}[$#$exons1]);
    }

    for(my $i=0; $i<scalar(@$terminal_exon1s); $i++) {
        my $exon1 = ${$terminal_exon1s}[$i];
        for(my $j=0; $j<scalar(@$exons2); $j++) {
        my $exon2 = ${$exons2}[$j];
        if(features_overlap($exon1,$exon2)) {
            # Check the different conditions
            if($exon2->seq_region_start() > $exon1->seq_region_start() && $exon2->seq_region_end() < $exon1->seq_region_start()) {
            # This means the exon2 is completely in exon1 and can be ignored
            next;
            } elsif($exon2->seq_region_end() >= $exon1->seq_region_start() && $exon2->seq_region_start() <= $exon1->seq_region_start() && $i == 0) {
            # This means there is an overlap on the 5' boundary and exon1 should be updated to match exon2
            $exon1->start($exon2->start());
            $modified_five_prime = 1;
            } elsif($exon2->start() <= $exon1->end() && $exon2->end() >= $exon1->end() && $i == $#$terminal_exon1s) {
            # This means there is an overlap on the 3' boundary and exon1 should be updated to match exon2
            $exon1->end($exon2->end());
            $modified_three_prime = 1;
            }
            if   ($modified_three_prime || $modified_five_prime){
                push(@$merged_exons,$exon1);
            }
            
        } else {
            if($exon2->end() < $cloned_transcript1->start() && $i == 0) {
            #push(@{$five_prime_exon2s},$exon2);
            #$modified_five_prime = 1;
            push(@$merged_exons,$exon2);
            } elsif($exon2->start > $cloned_transcript1->end() && $i == $#$terminal_exon1s) {
            #push(@{$three_prime_exon2s},$exon2);
            #$modified_three_prime = 1;
            push(@$merged_exons,$exon2);
            }
        } # End else
        } # End for(my $j=0;
    } # End for(my $i=0;

    unless($modified_five_prime || $modified_three_prime) {
        return;
    }

=head2 TO DELETE IF THE NESTED LOGIC WORKS
  if($modified_five_prime && scalar(@{$five_prime_exon2s})) {
    push(@$merged_exons,@{$five_prime_exon2s});
  } else {
    push(@$merged_exons,${$terminal_exon1s}[0]);
  }
  if(scalar(@$internal_exon1s)) {
    push(@$merged_exons,@{$internal_exon1s});
  }
  if(scalar(@$terminal_exon1s) > 1) {
    push(@$merged_exons,${$terminal_exon1s}[$#$terminal_exon1s]);
  }
  if($modified_three_prime && scalar(@{$three_prime_exon2s})) {
    push(@$merged_exons,@{$three_prime_exon2s});
  }
=cut

    $merged_exons = [sort { $a->start <=> $b->start } @{$merged_exons}];

    # Since these exons might be modified in other merges, clone them before making the new transcript
    $merged_exons = clone_exon_array($merged_exons);

    my $merged_exon_string = generate_exon_string($merged_exons);

    if($merged_exon_string eq $cloned_transcipt1_exon_string || $merged_exon_string eq $transcript2_exon_string){
        return;
    }
    $LoggerConfig::logger->info("Creating joined transcript from single exon transcript. Joined transcript has ".scalar(@$merged_exons)." exons");
    #say "Creating joined transcript from single exon transcript. Joined transcript has ".scalar(@$merged_exons)." exons";
    my $joined_transcript = set_transcript($merged_exons, $transcript1->stable_id(), $transcript1->biotype, $transcript1->slice());
    $joined_transcript->analysis($analysis);
    copy_transcript_attribs($transcript1,$joined_transcript);
    return($joined_transcript);
}

=head2 clone_exon_array

Copy list of exons.

=cut

sub clone_exon_array {
    my ($exons) = @_;
    my $cloned_exons = [];
    foreach my $exon (@$exons) {
        push(@$cloned_exons,clone_Exon($exon));
    }
    return($cloned_exons);
}

=head2 restore_strand

Restore reverse strand in input transcript.

=cut

sub restore_strand {
    my ($transcripts_to_change) = @_;
    my $reverse_transcripts = [];

    foreach my $transcript (@$transcripts_to_change) {
        my $original_strand = $transcript->{'original_strand'};

        if ($original_strand == 1) {
            push(@$reverse_transcripts, $transcript);
            next;
        }

        my $reverse_exons = [map {
            set_exon($_->start(),$_->end(),-1,$_->slice());
        } @{$transcript->get_all_Exons()}];

        my $sorted_reverse_exons = [sort { $b->end <=> $a->end } @$reverse_exons];

        my $reverse_transcript = set_transcript( $sorted_reverse_exons, $transcript->stable_id(), $good_biotype, $transcript->slice())
        $reverse_transcript->strand(-1);
        $reverse_transcript->analysis($analysis);

        copy_transcript_attribs($transcript, $reverse_transcript);
        push(@$reverse_transcripts, $reverse_transcript);
    }

    return $reverse_transcripts;
}

=head2 copy_transcript_attribs

Copy attribs from transcript1 to transcript2

=cut

sub copy_transcript_attribs {
    my ($transcript1,$transcript2) = @_;
    $transcript2->{'original_strand'} = $transcript1->{'original_strand'};
}

=head2 clean_initial_transcripts

    # Removes initial transcripts if the cds  is short
    # It also removes short introns and recalculates the cds to see if removing them increases the cds length
    # then add the modified transcript back to the pile.     

=cut

sub clean_initial_transcripts {
    my ($transcripts) = @_;
    my $cleaned_transcripts = [];

    for my $index (0 .. $#{$transcripts}) {
        my $transcript = $transcripts->[$index];

        my $cds_seq = $transcript->translateable_seq();
        my $cds_length = $cds_seq ? length($cds_seq) : 0;

        my $introns = $transcript->get_all_Introns();
        if(scalar(@$introns) >= 1) {

            my $modified_transcript = check_introns($transcript, $introns, $min_intron_length);
            if($modified_transcript) {
                compute_translation($modified_transcript);
                my $modified_cds_seq = $modified_transcript->translateable_seq();
                if(length($modified_cds_seq) >= $cds_length) {
                push(@$transcripts,$modified_transcript);
                }
            }
        }
        my $valid_cds = scalar(@$introns) >= 1
        ? validate_cds_length($cds_seq, $min_multi_exon_cds_length)
        : validate_cds_length($cds_seq, $min_single_exon_cds_length);
        push(@$cleaned_transcripts, $transcript) if $valid_cds;
    }

    return $cleaned_transcripts;
}


=head2 validate_cds_length

    Checks the cds passes the length cutoff
    another check could be included
    # Checks the cds has a proper start/stop and passes the length cutoff
    #  unless($cds_seq =~ /^ATG/ && $cds_seq =~ /(TAA|TAG|TGA)$/ && length($cds_seq) >= $min_cds_length) {

=cut

sub validate_cds_length {
    my ($cds_seq, $min_cds_length) = @_;
    return length($cds_seq) >= $min_cds_length;
}

=head2 check_introns

    Checks the introns pass the length cutoff, if it is shorter the exon are merged and a new transcript is created

=cut

sub check_introns {
    my ($transcript, $introns, $min_intron_length) = @_;

    my $exons = $transcript->get_all_Exons();
    my $was_modified = 0;
    my $modified_exons = [];
    my $previous_exon = $exons->[0];
    push(@$modified_exons, clone_Exon($previous_exon));

    for my $intron_index (0 .. $#{$introns}) {
        my $intron = $introns->[$intron_index];
        my $next_exon = $exons->[$intron_index + 1];

        if ($intron->length() < $min_intron_length) {
            my $merged_exon = merge_exons($previous_exon, $next_exon);
            $was_modified = 1;
            $modified_exons->[-1] = $merged_exon;
            $previous_exon = $merged_exon;
        } else {
            push(@$modified_exons, clone_Exon($next_exon));
            $previous_exon = $next_exon;
        }
    }
    if($was_modified) {
        if($$modified_exons[0]->strand == 1) {
        $modified_exons = [sort { $a->start <=> $b->start } @{$modified_exons}];
        } else {
        $modified_exons = [sort { $b->end <=> $a->end } @{$modified_exons}];
        }

        my $modified_transcript = set_transcript( $modified_exons, $transcript->stable_id(), $good_biotype, $transcript->slice())
        $modified_transcript->analysis($analysis);
        compute_translation($modified_transcript);
        return($modified_transcript);
    }
}

=head2 merge_exons

Given two exons a new one is created

=cut

sub merge_exons {
    my ($exon1,$exon2) = @_;

    my @coords = ($exon1->seq_region_start(),$exon1->seq_region_end(),$exon2->seq_region_start(),$exon2->seq_region_end());
    my $start = min(@coords);
    my $end = max(@coords);
    my $merged_exon = set_exon($start, $end, $exon1->strand(), $exon1->slice());

    return($merged_exon);
}

=head2 process_transcript

    # Subroutine to process UTRs for potential ORFs
    # At this point the transcripts array has all the transcripts and each has a translation
    # Now we want to reprocess them to find additional ORFs
    # Iteratively, UTRs are made into new transcripts and added back to the pile
    # This reprocessing of UTRs continues until and of the following stop conditions are met:
    # 1: The UTR length is < 350bp
    # 2: The longest ORF from the computed translation is < 100aa
    # 3: The ORF does not start with a met or end with a stop

=cut    

sub process_transcript {
    my ($transcript, $transcripts, $processed_transcripts, $processed_exon_strings) = @_;
    #my $consider_cds_only = 0;
    if ($add_utr) {        
        process_with_utrs($transcript, $transcripts, $processed_transcripts, $processed_exon_strings);
    } else {
        process_cds_only($transcript, $processed_transcripts, $processed_exon_strings);
    }
}

=head2 process_cds_only

    # This is for pulling out cds seqs without UTR (for example in single cell stuff where it's mostly tightly packed single exon genes)
    # So the transcripts coming out of this are all then without UTR

=cut     

sub process_cds_only {
    my ($transcript, $processed_transcripts, $processed_exon_strings) = @_;
    my $cds_exons = $transcript->get_all_translateable_Exons();
    my $cds_exon_string = generate_exon_string($cds_exons);
    return if $processed_exon_strings->{$cds_exon_string};
    $processed_exon_strings->{$cds_exon_string} = 1;
    my $cds_transcript = set_transcript($cds_exons, 
                                    $transcript->stable_id(), 
                                    $good_biotype;
                                    $transcript->slice()
                                    )
    $cds_transcript->analysis($analysis);
    $cds_transcript->{'original_strand'} = $original_strand;
    compute_translation($cds_transcript);
    push(@$processed_transcripts, $cds_transcript);
    return $processed_transcripts;
}

=head2 process_with_utrs

    # This checks if the exon string has been seen and if not it's added to the processed_transcripts array and marked
    # If it has been seen already, then this just returns

=cut

sub process_with_utrs {
    my ($transcript, $transcripts, $processed_transcripts, $processed_exon_strings) = @_;
    # This checks if the exon string has been seen and if not it's added to the processed_transcripts array and marked
    # If it has been seen already, then this just returns
    my $exons = $transcript->get_all_Exons();
    my $exon_string = generate_exon_string($exons);
    return if $processed_exon_strings->{$exon_string};
    $processed_exon_strings->{$exon_string} = 1;
    my $utr_5p_features = $transcript->get_all_five_prime_UTRs();
    my $utr_3p_features = $transcript->get_all_three_prime_UTRs();
    my $utr_5p_exons = create_exons_from_utr($utr_5p_features);
    my $utr_3p_exons = create_exons_from_utr($utr_3p_features);
    my $min_5p_exon_count = 2;
    my $min_3p_exon_count = 2;
    my $min_5p_length = 500;
    my $min_3p_length = 500;
    my $cds_buffer_length = 100;
    my $utr_5p_length = cumulative_feature_length($utr_5p_exons, $cds_buffer_length);
    my $utr_3p_length = cumulative_feature_length($utr_3p_exons, $cds_buffer_length);
    # If these are exon strings that have not been seen before, we want to build a transcript out of them
    # and then assuming it passed out critera it will get added back to the transcripts array,
    # which means they'll come back into this method later and be further broke down if possible
    process_utr_transcript($utr_5p_exons, $transcript, $transcripts, $processed_exon_strings, $min_5p_exon_count, $min_5p_length);
    process_utr_transcript($utr_3p_exons, $transcript, $transcripts, $processed_exon_strings, $min_3p_exon_count, $min_3p_length);
    trim_transcript_utrs($transcript, $utr_3p_exons, $utr_5p_exons);

    push(@$processed_transcripts, $transcript);
    return $processed_transcripts;
}

=head2 create_exons_from_utr

    Get list of exons given utr feature 

=cut

sub create_exons_from_utr {
    my ($utr_features) = @_;
    my $exons = [];
    foreach my $utr (@$utr_features) {
        my $exon = set_exon($utr->start(),$utr->end(),$utr->strand(),$utr->slice());
        push(@$exons,$exon);
    }
    return($exons);
}

=head2 cumulative_feature_length

    Get cumulative length of CDS seq
    # If a buffer has been specificed take that into account, for example it would be unusual for one 
    # another gene to end closer than 100bp of the next so remove the buffer when considering
    # the availble length for a potential new CDS to fit into

=cut

sub cumulative_feature_length {
    my ($features, $buffer_length) = @_;
    # Set default buffer length to 0 if not provided
    $buffer_length //= 0;

    my $cumulative_length = 0;
    foreach my $feature (@$features) {
        $cumulative_length += $feature->length();
    }

    $cumulative_length -= $buffer_length;
    return $cumulative_length;
}

=head2 process_utr_transcript

    Buil a new transcript if it passes basic criteria

=cut

sub process_utr_transcript {
    my ($utr_exons, $transcript, $transcripts, $processed_exon_strings, $min_exon_count, $min_length) = @_;

    unless ($processed_exon_strings->{generate_exon_string($utr_exons)} || (scalar(@$utr_exons) < $min_exon_count) || (cumulative_feature_length($utr_exons) < $min_length)) {
        generate_new_transcript($utr_exons, $transcript, $transcripts);
    }
    return $transcripts;
}

=head2 generate_new_transcript

    Buil a new transcript ensure a minimum length for UTR transcripts

=cut

sub generate_new_transcript {
    my ($utr_exons, $transcript, $transcripts) = @_;

    my $min_cds_length = scalar(@$utr_exons) == 1 ? $min_single_exon_cds_length : $min_multi_exon_cds_length;

    my $utr_transcript = set_transcript(($utr_exons, $transcript->stable_id, $good_biotype, $transcript->slice()))
    $utr_transcript->analysis($transcript->analysis());
    copy_transcript_attribs($transcript, $utr_transcript);
    # Ensure a minimum length for UTR transcripts
    return unless $utr_transcript->length() >= ($min_single_exon_cds_length + 50);

    compute_translation($utr_transcript);
    my $cds_seq = $utr_transcript->translateable_seq();

    my $valid_cds = scalar(@$utr_exons) == 1
        ? validate_cds_length($cds_seq, $min_single_exon_cds_length)
        : validate_cds_length($cds_seq, $min_multi_exon_cds_length);

    return unless $valid_cds;
    push(@$transcripts, $utr_transcript);
    #return $transcripts;
}

=head2 trim_transcript_utrs

    Now for the original transcript, trim the 3' UTR as needed
    It is very unusual to have real introns in the 3' UTR. Usually an intron in the 3' UTR is incorrect, 
    because of the reconstructor or because of thing like transposons
    First cut off any additional introns (these will be searched for ORFs part of the recursive nature of this process)
    Then trim the remaining UTR based on looking for the PAS signal

=cut

sub trim_transcript_utrs {

    my ($transcript, $utr_3p_exons, $utr_5p_exons) = @_;

    if (scalar(@$utr_3p_exons)) {
        $transcript = remove_three_prime_exons($transcript, $utr_3p_exons);
        $transcript = trim_3prime_utr($transcript);
    }

    if (scalar(@$utr_5p_exons)) {
        $transcript = remove_five_prime_exons($transcript);
    }
    push(@$processed_transcripts,$transcript);
    #return $transcript;
}

=head2 remove_five_prime_exons

Trim 5 prime UTR removing short CDS

=cut

sub remove_three_prime_exons {
    my ($transcript, $utr_3p_features) = @_;
    my $exons         = $transcript->get_all_Exons();
    my $kept_exons    = [];
    my $skip          = 1;
    my $closest_utr_exon;
    for (my $i = @$exons - 1; $i >= 0; $i--) {
        my $exon = $exons->[$i];
        if ($exon->is_coding($transcript) and $skip) {
            $closest_utr_exon = $exons->[$i + 1] if $i < @$exons - 1;
            $skip              = 0;
        }
        push @$kept_exons, $exon unless $skip;
    }
    # Treat the first pure UTR exon (if it exists) as a special case. If it's close to the end of the CDS
    # and it's not abnormally short, then allow it
    my $max_first_intron_dist     = 100;
    my $min_first_intron_length   = 250;
    my $min_closest_exon_length   = 100;
    my $max_utr_length            = 1000;
    if (
        $closest_utr_exon
        and $utr_3p_features->[0]->length() < $max_utr_length
        and $closest_utr_exon->length() >= $min_closest_exon_length
    ) {
        my $intron_size = $closest_utr_exon->strand == 1
            ? $closest_utr_exon->seq_region_start() - $kept_exons->[0]->seq_region_end() + 1
            : $kept_exons->[0]->seq_region_start() - $closest_utr_exon->seq_region_end() + 1;
        push @$kept_exons, $closest_utr_exon if $intron_size >= $min_first_intron_length;
    }
    my @unsorted_exons = @$kept_exons;
    my @sorted_exons = sort {
    (${$kept_exons}[0]->strand == 1) ? $a->start <=> $b->start : $b->start <=> $a->start
    } @{$unsorted_exons};
    $kept_exons = \@sorted_exons;
    my $last_exon       = $kept_exons->[-1];
    my $coding_offset   = 0;
    # First figure out if there's a coding offset to take into account
    if ($last_exon->is_coding($transcript)) {
        #Returns the end position of the coding region of the exon in slice-relative coordinates on the forward strand.
        my $cds_end_coord = $last_exon->coding_region_end($transcript)
        $coding_offset = ($last_exon->strand() == 1)
        ? $cds_end_coord - $last_exon->seq_region_start()
        : $last_exon->seq_region_end() - $cds_end_coord;
    }
    my $last_exon_utr_length = $last_exon->length - $coding_offset;
    # Adjust the UTR of the final exon if its length exceeds the max
    if ($last_exon_utr_length > $max_utr_length) {
        my $diff = $last_exon_utr_length - $max_utr_length;
        if ($last_exon->strand() == 1 and $diff > 0) {
            $last_exon->end($last_exon->end() - $diff);
        } else {
            $last_exon->start($last_exon->start() + $diff);
        }
    }
    my $new_transcript = set_transcript($kept_exons, $transcript->stable_id, $transcript->biotype(), $transcript->slice())
    $new_transcript->analysis($transcript->analysis());
    copy_transcript_attribs($transcript, $new_transcript);
    return $new_transcript;
}

=head2 remove_five_prime_exons

Trim 5 prime UTR keeping coding and non coding exons until min length

=cut

sub remove_five_prime_exons {
    my ($transcript) = @_;
    compute_translation($transcript) unless $transcript->translation();
    my $exons          = $transcript->get_all_Exons();
    my $kept_exons     = [];
    my $min_utr_terminal_exon_length = 100;
    my $max_5p_utr_length  = 300;
    my $skip = 1;
    for my $i (0 .. $#$exons) {
        my $exon = $exons->[$i];
        $skip = 0, if $exon->is_coding($transcript);

        push(@$kept_exons, $exon) unless $skip || ($skip && $exon->length > $min_utr_terminal_exon_length);
    }
    my @unsorted_exons = @$kept_exons;
    my @sorted_exons = sort {
    (${$kept_exons}[0]->strand == 1) ? $a->start <=> $b->start : $b->start <=> $a->start
    } @{$unsorted_exons};
    $kept_exons = \@sorted_exons;
    # Now loop through the sorted exons and find the index if the first coding exon, 
    # then work back in the 5' direction
    my $cds_exon_index = 0;
    foreach my $exon (@$kept_exons) {
        last if($exon->is_coding($transcript));
        $cds_exon_index++;
    }
    my $cumulative_exon_length = 0;
    my $exons_5p_to_keep       = [];
    for my $i (reverse 0 .. $cds_exon_index) {
        my $exon = ${$kept_exons}[$i];
        my $coding_offset = $exon->is_coding($transcript) ? get_coding_offset($exon, $transcript) : 0;
        my $coding_offset = 0;
        if ($exon->is_coding($transcript)) {
            my $cds_start_coord = $exon->coding_region_start($transcript);
            $coding_offset = $exon->strand == 1 ? $exon->seq_region_end - $cds_start_coord : $cds_start_coord - $exon->seq_region_start;
        }
        my $exon_utr_length = $exon->length() - $coding_offset;
        # At this point we know the length of the current UTR and the cumulative length
        # If the length of the of the current utr exon added to the cumulative total > max_5p_utr_length 
        # (minus any coding offset)
        # then this exon needs to either be dropped or it's the last one (technically the first one).
        # First determine if once trimmed the length is >= min_remaining_length, if it is then trim and keep.
        # If it's not but it's a coding exon just keep as is. Otherwise just drop it.
        # Once any of these conditions are hit, stop and add all the exons after the coding exon index
        my $current_cumulative_length = $cumulative_exon_length + $exon_utr_length;
        my $min_remaining_length = 50;
        if ($current_cumulative_length > $max_5p_utr_len) {
            my $overlimit_length = $current_cumulative_length - $max_5p_utr_len;
            my $trimmed_length   = $exon_utr_length - $overlimit_length;
            if ($trimmed_length < $min_remaining_length && $exon->is_coding($transcript)) {
                push(@$exons_5p_to_keep, $exon);
                last;
            } elsif ($trimmed_length < $min_remaining_length) {
                last;
            } else {
                $exon->strand == 1 ? $exon->start($exon->start() + $overlimit_length) : $exon->end($exon->end() - $overlimit_length);
                push(@$exons_5p_to_keep, $exon);
                last;
            }
        }
        push(@$exons_5p_to_keep, $exon);
        $cumulative_exon_length = $current_cumulative_length;
    }
    my $final_exons = [@$exons_5p_to_keep, @$kept_exons[$coding_exon_index + 1 .. $#$kept_exons]];
    my $new_transcript = set_transcript($final_exons, $transcript->stable_id, $transcript->biotype(), $transcript->slice())
    $new_transcript->analysis($transcript->analysis());
    copy_transcript_attribs($transcript, $new_transcript);
    return $new_transcript;
}

=head2 trim_3prime_utr

Trim 3 prim UTR looking for the PAS signal

=cut

sub trim_3prime_utr {
    my ($transcript) = @_;

    my $pas_signal        = 'AATAAA';
    my $cleavage_signal   = 'CA';
    my $max_no_cleavage   = 1000;

    unless ($transcript->three_prime_utr) {
        warning("The trim_3prime_utr_short_read was called on a transcript with no 3 prime UTR. Nothing to trim");
        return $transcript;
    }

    my $exons                  = $transcript->get_all_Exons;
    my $final_exon             = $exons->[-1];
    my $translation_end_exon   = $transcript->translation->end_Exon;

    $coding_offset = $transcript->translation->end if $final_exon->start == $translation_end_exon->start:0;


    my $final_exon_seq = $final_exon->seq->seq;

    my ($pas_start, $pas_end) = (0, 0);
    my $found_pas = 0;

    while ($final_exon_seq =~ /$pas_signal/g && !$found_pas) {
        $pas_start = $-[0] + 1;
        $pas_end   = $+[0];

        next if $pas_start <= $coding_offset;
        $LoggerConfig::logger->info("Found PAS signal in final exon seq at the following coords: $pas_start..$pas_end");
        #say "Found PAS signal in final exon seq at the following coords: $pas_start..$pas_end";
        $found_pas = 1;
        last;
    }

    #my $cleavage_site = find_cleavage_site($final_exon_seq, $found_pas, $pas_end, $coding_offset);
    my $cleavage_site = 0;

    if ($found_pas) {
        # There are 15-30bp between the end of the pas signal and the cleavage site
        # as pas_end is already shifted to 1bp offset, just add 14
        my $post_pas_seq = substr($final_exon_seq, $pas_end + 14, 15);

        if ($post_pas_seq =~ /CA/) {
            $cleavage_site = $pas_end + 14 + $+[0];
        } else {
            $cleavage_site = $pas_end + 30;
        }
    } else {
        $cleavage_site = (length($final_exon_seq) - $coding_offset) > $max_no_cleavage
            ? $coding_offset + $max_no_cleavage
            : 0;
    }

    if ($cleavage_site >= length($final_exon_seq)) {
        $LoggerConfig::logger->info("Not cleaving as proposed cleavage site is at or over the end of the final exon");
        #say "Not cleaving as proposed cleavage site is at or over the end of the final exon";
        return $transcript;
    }
    if ($cleavage_site) {
        if ($final_exon->strand == 1 && $cleavage_site > $coding_offset) {
            $final_exon->end($final_exon->start + $cleavage_site - 1);
            $transcript->end($final_exon->end);
        } elsif ($final_exon->strand == -1 && $cleavage_site > $coding_offset) {
            $final_exon->start($final_exon->end - $cleavage_site + 1);
            $transcript->start($final_exon->start);
        }
    }
    return $transcript;
}

=head2 select_geneset

Filter geneset using parameters defined in stats

=cut

sub select_geneset {
    my ($genes) = @_;
    my @intron_sizes = ();
    my @scores = ();

    my $final_genes = [];
    my $initial_good_genes = [];
    my $initial_bad_genes = [];
    my $revised_good_genes = [];

    foreach my $gene (@$genes) {
    my $transcript = $gene->get_all_Transcripts->[0];
    my $exons = $transcript->get_all_Exons();
    my $cds_exons = $transcript->get_all_CDS();
    my $cds_exon_string = generate_exon_string($cds_exons);
    my $cds_exon_length = sum(map { $_->length } @$cds_exons) / scalar(@$cds_exons);
    $cds_exon_length //= 0;

    my $cds_seq = $transcript->translateable_seq;
    my $cds_length = length($cds_seq);
    my $genomic_span = $transcript->seq_region_end - $transcript->seq_region_start + 1;
    my $introns = $transcript->get_all_Introns();
    push(@intron_sizes, $_->length) foreach @$introns;

    my $score = 0;
    if (scalar(@$cds_exons) == 1 && $cds_length >= $stats->{'single_exon_cds_length'}) {
        $score++;
    } elsif (length($cds_seq) >= $stats->{'cds_length'}) {
        $score++;
    }
    push(@scores,$score);
    if($score < $min_score) {
        push(@$initial_bad_genes,$gene);
    } else {
        push(@$initial_good_genes,$gene);
    }
    
    }
    my $revised_good_genes = process_initial_good_genes($initial_good_genes, $initial_bad_genes);

    foreach my $gene (@$revised_good_genes) {
        my $transcript = ${$gene->get_all_Transcripts()}[0];
        $transcript->biotype($good_biotype);
        $gene->biotype($good_biotype);
        push(@$final_genes,$gene);
    }

    my $final_bad_genes = final_classification($revised_good_genes,$initial_bad_genes);
    foreach my $gene (@$final_bad_genes) {
        my $transcript = ${$gene->get_all_Transcripts()}[0];
        $transcript->biotype($bad_biotype);
        $gene->biotype($bad_biotype);
        push(@$final_genes,$gene);
    }

    return $final_genes;
}

=head2 process_initial_good_genes

Process initial geneset, calculates cluster and get the longest transcript or the one with more cds sequences

=cut

sub process_initial_good_genes {
    my ($good_genes, $bad_genes) = @_;

    my $revised_good_genes = [];
    my $unique_biotypes = get_unique_biotypes($genes);
    my $types_hash;
    $types_hash->{genes} = $biotypes_array;
    my ($clusters, $unclustered) = cluster_Genes_by_coding_exon_overlap($good_genes, $types_hash);

    foreach my $cluster (@$clusters) {
        my $intron_strings = {};
        foreach my $gene (@{$cluster->get_Genes()}) {
            my $transcript = ${$gene->get_all_Transcripts()}[0];
            my $cds_exons = $transcript->get_all_CDS();
            my $cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
            copy_transcript_attribs($transcript, $cds_transcript);
            my $introns = $cds_transcript->get_all_Introns();
            my $intron_string = generate_intron_string($introns);
            $intron_strings->{$intron_string} = 1;
        }
        my ($cluster_good_genes, $max_cluster_cds_exons, $max_cluster_orf_length) =
            process_cluster_genes($cluster, $intron_strings);
        update_revised_good_genes($revised_good_genes, $bad_genes, $cluster_good_genes,
            $max_cluster_cds_exons, $max_cluster_orf_length);
    }
    process_unclustered_genes($unclustered, $revised_good_genes);

    return $revised_good_genes;
}

=head2 generate_intron_string

Generate string with intron coordinates

=cut

sub generate_intron_string {
    # Assumes stranded data
    my ($intron_array) = @_;

    my $intron_string = "";
    my $count = 0;
    foreach my $intron (@{$intron_array}) {
        my $start = $intron->start();
        my $end = $intron->end();
        $intron_string .= $start."..".$end.":";
    }
    return($intron_string);
}

=head2 process_cluster_genes

Process gene cluster, move fragments into bad gene and calculates max cds exons, max orf length from the good genes 

=cut

sub process_cluster_genes {
    my ($cluster, $intron_strings) = @_;

    my $cluster_good_genes = [];
    my $max_cluster_cds_exons = 0;
    my $max_cluster_orf_length = 0;

    foreach my $gene (@{$cluster->get_Genes()}) {
        my $transcript = ${$gene->get_all_Transcripts()}[0];
        my $cds_exon_count = scalar(@{$transcript->get_all_CDS()});
        my $cds_exons = $transcript->get_all_CDS();
        # This is faster than calling to get all CDS introns from the API
        my $cds_transcript = Bio::EnsEMBL::Transcript->new(-exons => $cds_exons);
        copy_transcript_attribs($transcript,$cds_transcript);
        my $introns = $cds_transcript->get_all_Introns();
        my $intron_string = generate_intron_string($introns);

        
        my $cds_exon_count = scalar(@{$transcript->get_all_CDS()});
        my $cds_length = length($transcript->translateable_seq);

        my $is_substring =  is_substring($intron_string, $intron_strings);
        if ($is_substring) {
            push(@$bad_genes, $gene);
        } else {
            push(@$cluster_good_genes, $gene);
            $max_cluster_cds_exons = max($max_cluster_cds_exons, $cds_exon_count);
            $max_cluster_orf_length = max($max_cluster_orf_length, $cds_length);
        }
    }
    return ($cluster_good_genes, $max_cluster_cds_exons, $max_cluster_orf_length);
}

=head2 is_substring

Identifyy sequence fragments

=cut

sub is_substring {
    my ($intron_string, $unique_strings) = @_;

    foreach my $unique_string (keys %$unique_strings) {
        next if $unique_string eq $intron_string;
        return 1 if $unique_string =~ /$intron_string/;
    }
    return 0;
}

=head2 update_revised_good_genes

Filter geneset according to coding sequence length

=cut

sub update_revised_good_genes {
    my ($revised_good_genes, $bad_genes, $cluster_good_genes, $max_cluster_cds_exons, $max_cluster_orf_length) = @_;

    foreach my $gene (@$cluster_good_genes) {
        my $transcript = ${$gene->get_all_Transcripts()}[0];
        my $cds_exon_count = scalar(@{$transcript->get_all_CDS()});
        my $cds_length = length($transcript->translateable_seq);

        if ($cds_exon_count < ($max_cluster_cds_exons * 0.75) && $cds_length < ($max_cluster_orf_length * 0.75)) {
            push(@$bad_genes, $gene);
        } else {
            push(@$revised_good_genes, $gene);
        }
        return $bad_genes, $revised_good_genes;
    }

}

=head2 process_unclustered_genes

Add unclustered genes directly in the good batch

=cut

sub process_unclustered_genes {
    my ($unclustered, $revised_good_genes) = @_;

    foreach my $singleton (@$unclustered) {
        my $single_genes = $singleton->get_Genes();
        push @$revised_good_genes, @$single_genes;
    }
    return $revised_good_genes;
}

=head2 final_classification

cluster and filter the bad batch; add unclustered genes directly into the bad batch

=cut
sub final_classification {
    my ($good_genes, $bad_genes) = @_;
    my $final_bad_genes = [];
    my $bad_biotype = 'transcriptomic_check';
    $_->biotype($bad_biotype) for @$bad_genes;
    my $good_unique_biotypes = get_unique_biotypes($good_genes);
    my $bad_unique_biotypes = [$bad_biotype];
    my $types_hash = {genes => [@{$good_unique_biotypes}, $bad_biotype]};
    my $all_genes = [@{$good_genes}, @$bad_genes];

    $LoggerConfig::logger->info("Clustering genes from input_dbs...");
    #say "Clustering genes from input_dbs...";
    my ($clusters, $unclustered) = cluster_Genes_by_coding_exon_overlap($all_genes, $types_hash);

    my $final_bad_genes = [];
    for my $cluster (@$clusters) {
        my $good_cluster_genes = $cluster->get_Genes_by_Type($good_unique_biotypes);
        my $bad_cluster_genes  = $cluster->get_Genes_by_Type($bad_unique_biotypes);
        next unless scalar(@$bad_cluster_genes);
        if (!scalar(@$good_cluster_genes) and scalar(@$bad_cluster_genes)) {
            push @$final_bad_genes, @$bad_cluster_genes;
            next;
        }
        $_->biotype('transcriptomic_ignore') for @$bad_cluster_genes;
        push @$final_bad_genes, @$bad_cluster_genes;
    }
    for my $singleton (@$unclustered) {
        push @$final_bad_genes, @{ $singleton->get_Genes_by_Type($bad_unique_biotypes) };
    }

    return $final_bad_genes;
}

=head2 build_final_geneset

Run Genebuild module to filter the final geneset

=cut

sub build_final_geneset {
    my ($genes) = @_;

    my $final_genes = [];
    my $good_genes_by_slice = sort_features_by_slice($genes);

    foreach my $slice_name (keys(%$good_genes_by_slice)) {
    my $genes = $good_genes_by_slice->{$slice_name};
    my $slice = ${$genes}[0]->slice();
    my $biotype = ${$genes}[0]->biotype();
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
                        -query => $slice,
                        -analysis => $analysis,
                        -genes => $genes,
                        -output_biotype => $biotype,
                        -max_transcripts_per_cluster => 100,
                        -min_short_intron_len => 7,
                        -max_short_intron_len => 15,
                        -blessed_biotypes => {},
                        -skip_readthrough_check => 1,
                        -max_exon_length => 50000,
                        -coding_only => 1,
                    );
    $runnable->run;
    push(@$final_genes,@{$runnable->output});
    }
    return $final_genes;
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

=head2 prep_genes_for_writing

        # This array will store exactly what to write. By default, we would write the gene as constructed to this point.
        # However, if the run is just to process a particular type of transcripts into a cleaned set, we may want to
        # pull those back out into single transcript genes for further processing later. In that case, the array will
        # be replaced with single transcript genes

=cut

sub prep_genes_for_writing {
    my ($final_genes, $write_single_transcript_genes, $analysis) = @_;

    my $genes_to_write = [];

    foreach my $gene (@$final_genes) {


        my $transcripts = $gene->get_all_Transcripts();

        if ($write_single_transcript_genes && scalar(@$transcripts) > 1) {
            # If we want to write single transcript genes and the gene has multiple transcripts, split them into individual genes
            foreach my $transcript (@$transcripts) {
                my $new_gene = Bio::EnsEMBL::Gene->new(
                    -START  => $transcript->start(),
                    -END    => $transcript->end(),
                    -STRAND => $transcript->strand,
                    -SLICE  => $transcript->slice()
                );
                $new_gene->add_Transcript($transcript);
                $new_gene->analysis($analysis);
                push(@$genes_to_write, $new_gene);
            }
        } else {
            # Otherwise, just push the original gene to the array
            push(@$genes_to_write, $gene);
        }
    }

    return $genes_to_write;
}

=head2 write_to_gtf_file

    Write geneset into gtf, transcript sequences in cdna file and protein sequences in protein files.

=cut

sub write_to_gtf_file {
    my ($genes_to_write, $output_gtf_file, $analysis_name, $final_biotype) = @_;

    # Prepare additional output files for transcript and protein 
    my $output_transcript_seq_file = "$output_gtf_file.cdna";
    my $output_transcript_prot_file = "$output_gtf_file.prot";

    # Open output files
    open(my $out_gtf, ">", $output_gtf_file);
    open(my $out_cdna, ">", $output_transcript_seq_file);
    open(my $out_prot, ">", $output_transcript_prot_file);

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

            # Build and print GTF record
            my $gtf_record = build_gtf_record($transcript, $exons, $analysis_name, $gene_id, $transcript_id, $translation, $final_biotype);
            say $out_gtf join("\n", @$gtf_record);

            # Print transcript sequence
            say $out_cdna ">${transcript_id}\n" . $transcript->seq->seq();

            # Print protein sequence if available
            say $out_prot ">${transcript_id}\n" . ($translation ? $translation->seq() : "");

            $transcript_count++;
        }

        $gene_count++;
    }

    # Close output files
    close $out_gtf;
    close $out_cdna;
    close $out_prot;
}

=head2 build_gtf_record

    Define transcript and exon records for the GTF file

=cut

sub build_gtf_record {
    my ($transcript, $exons, $analysis_name, $gene_id, $transcript_id, $translation, $final_biotype) = @_;

    my $record = [];
    my $strand = ($transcript->strand() == -1) ? "-" : "+";

    # Build attributes for the transcript
    my $transcript_attribs = qq{gene_id "${gene_id}"; transcript_id "${transcript_id}"; biotype "${final_biotype}";};
    if ($translation) {
        # Encode translation information
        my $start_exon = $translation->start_Exon();
        my $end_exon = $translation->end_Exon();
        my $translation_coords = qq{ translation_coords "${start_exon->start()}:${start_exon->end()}:${translation->start()}:${end_exon->start()}:${end_exon->end()}:${translation->end()}";};
        $transcript_attribs .= $translation_coords;
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
