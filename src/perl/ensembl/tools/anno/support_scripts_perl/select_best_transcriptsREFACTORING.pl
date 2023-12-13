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

my ($host, $port, $user, $pass, $dbname, $dna_host, $dna_port, $dna_user, $dna_dbname);
my $analysis_name   = "anno";
my $module_name     = "anno";
my $region_details  = '';
my $specify_strand;
my $good_biotype    = 'transcriptomic';
my $bad_biotype     = 'transcriptomic_flagged';
my $all_cds_exons   = 0;
my $join_transcripts = 0;
my $clean_transcripts = 0;
my $cds_search       = 0;
my $write_protein_seqs    = 1;
my $write_transcript_seqs = 0;
my $write_single_transcript_genes = 1;
my $set_canonical    = 1;
my $skip_db_write    = 0;
my $output_path;
my $genome_file      = '';
my $input_gtf_file   = '';
my $output_gtf_file  = '';
my $final_biotype    = 'not_set';

GetOptions(
    'gtf_file=s'             => \$input_gtf_file,
    'region_details=s'       => \$region_details,
    'specify_strand=s'       => \$specify_strand,
    'cds_search!'            => \$cds_search,
    'write_protein_seqs!'    => \$write_protein_seqs,
    'write_transcript_seqs!' => \$write_transcript_seqs,
    'skip_db_write!'         => \$skip_db_write,
    'output_path=s'          => \$output_path,
    'input_gtf_file=s'       => \$input_gtf_file,
    'output_gtf_file=s'      => \$output_gtf_file,
    'genome_file=s'          => \$genome_file,
    'final_biotype=s'        => \$final_biotype,
    'join_transcripts!'      => \$join_transcripts,
    'clean_transcripts!'     => \$clean_transcripts,
    'all_cds_exons!'         => \$all_cds_exons
);

say "$input_gtf_file\n$output_gtf_file\n$region_details\n";

# Validate GTF file
validate_gtf_file($input_gtf_file);
# Setup fasta
setup_fasta(-FASTA => $genome_file);
# Parse region details
my ($region_name, $region_start, $region_end) = parse_region_details($region_details);
my $analysis = create_analysis_object($analysis_name, $module_name);
my $slice = set_slice($region_name, $region_start, $region_end, 1);

# Process GTF file and get transcripts in the slice
my $transcripts = process_gtf_file($input_gtf_file, $region_name, $slice, $specify_strand, $slice);

#my $transcripts_by_slice = sort_transcripts_by_slice($transcripts);
#foreach my $slice_name (keys %$transcripts_by_slice) {
#say "Processing region: $slice_name";

#my $initial_transcripts = $transcripts_by_slice->{$slice_name};
my $sorted_transcripts = remove_overlapping_exons($transcripts);

#if the condition is true builds a new array with the modified transcripts (ALWAYS DISABLED NOW)
my $cloned_transcripts = $join_transcripts ? join_transcripts([map { clone_Transcript($_->biotype('orig')) } @$sorted_transcripts]) : [];
# if the condition is true combine cloned transcripts with the original ones
my $joined_transcripts = $join_transcripts ? [@$cloned_transcripts, @$sorted_transcripts] : $sorted_transcripts;

say "Transcript count after joining: " . scalar(@$joined_transcripts);

# Continue with your existing code...
    

# More of your code...



=head1 METHODS

=head2 validate_gtf_file

Validate GTF file.

=cut

sub validate_gtf_file {
    my ($gtf_file) = @_;
    die "Could not open the GTF file, path used: $gtf_file" unless -e $gtf_file;
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


=head2 process_gtf_file

Read GTF file.

=cut
sub process_gtf_file {
    my ($input_gtf_file, $region_name, $region_start, $region_end, $specify_strand, $slice) = @_;
    my %strand_conversion = ('+' => '1', '-' => '-1', '.' => undef, '?' => undef);
    say "Reading in GTF";
    open(my $gtf_fh, '<', $input_gtf_file) or die "Could not open GTF file: $input_gtf_file $!";
    my $exons = [];
    my $transcripts = [];
    my $biotype = 'not_set';
    my $current_transcript_id = "";
    my $current_gtf_region_name = "";
    my $gtf_region_name;
    my $new_record = 0;

    while (<$gtf_fh>) {
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
        $biotype = $attributes->{'biotype'};
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
        my $transcript = Bio::EnsEMBL::Transcript->new( -EXONS => $exons );
        push( @{$transcripts}, $transcript );      
    }
    return $transcripts;
    say "Finished reading GTF";
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
  say "Clustering genes from input_dbs...";
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

            if ($joined_transcript) {
                my $joined_transcript_exon_count = scalar(@{$joined_transcript->get_all_Exons()});
                my $joined_transcript_length    = $joined_transcript->length();

                # Check if the joined transcript has the longest sequence length
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

        # Identify and keep only distinct transcripts among the two best candidates
        my $final_joined_transcripts = [];
        if ($best_joined_transcripts->[0] && $best_joined_transcripts->[1]) {
            ${$best_joined_transcripts}[0]->{'original_strand'} = $transcript_i->{'original_strand'};
            ${$best_joined_transcripts}[1]->{'original_strand'} = $transcript_i->{'original_strand'};
            my $exon_string1 = generate_exon_string($best_joined_transcripts->[0]->get_all_Exons());
            my $exon_string2 = generate_exon_string($best_joined_transcripts->[1]->get_all_Exons());

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
  # 1:

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

  say "Creating joined transcript from single exon transcript. Joined transcript has ".scalar(@$merged_exons)." exons";
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

Restore strand in input transcript.

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
            Bio::EnsEMBL::Exon->new(
                -START     => $_->start(),
                -END       => $_->end(),
                -STRAND    => -1,
                -SLICE     => $_->slice(),
                -PHASE     => -1,
                -END_PHASE => -1
            );
        } @{$transcript->get_all_Exons()}];

        my $sorted_reverse_exons = [sort { $b->end <=> $a->end } @$reverse_exons];

        my $reverse_transcript = Bio::EnsEMBL::Transcript->new(
            -EXONS    => $sorted_reverse_exons,
            -STRAND   => -1,
            -STABLE_ID => $transcript->stable_id(),
            -BIOTYPE  => $good_biotype,
            -SLICE    => $transcript->slice(),
            -ANALYSIS => $analysis
        );

        copy_transcript_attribs($transcript, $reverse_transcript);
        push(@$reverse_transcripts, $reverse_transcript);
    }

    return $reverse_transcripts;
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