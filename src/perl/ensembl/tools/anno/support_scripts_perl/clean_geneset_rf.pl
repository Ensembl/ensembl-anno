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
use Utils qw(check_gtf_file setup_fasta set_slice create_analysis_object set_transcript set_exon set_translation set_canonical_transcript parse_region_details write_to_gtf_file
  build_gtf_record set_attributes build_translation create_single_transcript_genes sort_features_by_slice
  run_gene_builder_for_slice);

my ($host, $port, $user, $pass, $dbname, $dna_host, $dna_port, $dna_user, $dna_dbname);

my $genome_file;
my $gtf_file;
my $slice_name;
my $compute_translations;
my $output_gtf_file;
my $set_canonical = 1;
my $analysis_name = "ensembl";
my $module_name = "Anno";


GetOptions(
            'host|dbhost|h=s'       => \$host,
            'port|dbport|P=s'       => \$port,
            'user|dbuser|u=s'       => \$user,
            'pass|dbpass|p=s'       => \$pass,
            'dbname|db|D=s'     => \$dbname,
            'dna_host=s'   => \$dna_host,
            'dna_port=s'   => \$dna_port,
            'dna_user=s'   => \$dna_user,
            'dna_dbname=s' => \$dna_dbname,
            'gtf_file=s' => \$gtf_file,
            'genome_file=s' => \$genome_file,
            'slice_name=s' => \$slice_name,
            'compute_translations!' => \$compute_translations,
            'output_gtf_file=s' => \$output_gtf_file);

check_gtf_file(input_gtf_file)
if ($genome_file && -e $genome_file) {
  setup_fasta(-FASTA => $genome_file);
}

my $slice_hash = {};
my $slices = fetch_slices_from_fasta($genome_file);
foreach my $slice (@{$slices}) {
  my $seq_region_name = $slice->seq_region_name;
  next if ($slice_name && $seq_region_name ne $slice_name);
  $slice_hash->{$seq_region_name} = $slice;
}
my $analysis = create_analysis_object($analysis_name, $module_name);
# Process GTF file and get transcripts in the slice
my $exons = process_gtf_file($input_gtf_file, $slice_hash, $analysis);



if(scalar(@$exons)) {
  if($$exons[0]->strand() == 1) {
    $exons = [sort { $a->start <=> $b->start } @{$exons}];
  } else {
    $exons = [sort { $b->start <=> $a->start } @{$exons}];
  }
  my $transcript = set_transcript($exons, $current_transcript_id, $current_biotype, $$exons[0]->slice())
  $transcript->analysis($analysis);

  if($current_translation_coords) {
    build_translation($transcript,$current_translation_coords);
  } elsif($compute_translations) {
    compute_translation($transcript);
  }

  push(@$transcripts,$transcript);
}

my $final_genes = [];
my $single_transcript_genes = create_single_transcript_genes($transcripts);
my $single_transcript_genes_by_slice = sort_features_by_slice($single_transcript_genes);
foreach my $slice_name (keys(%$single_transcript_genes_by_slice)) {
  my $slice_genes = $single_transcript_genes_by_slice->{$slice_name};
  my $slice = $slice_hash->{$slice_name};
  
  # Run GeneBuilder
  my $collapsed_slice_genes = run_gene_builder_for_slice($slice, $slice_genes, $analysis, 'gbuild');
  my $genes_with_unique_ids = set_ids($collapsed_slice_genes);
  my $genes_with_canonicals = set_canonical_transcripts($genes_with_unique_ids);
  my $genes_with_biotypes = set_gene_biotypes($genes_with_canonicals);
  my $post_lncrna_filtered_genes = remove_overlapping_lncrnas($genes_with_biotypes);
  my $collapsed_lncrna_genes = collapse_lncrna_genes($post_lncrna_filtered_genes,$analysis);
  my $genes_without_readthroughs = remove_potential_readthroughs($collapsed_lncrna_genes,$analysis);
  push(@$final_genes,@$genes_without_readthroughs);
}

$final_genes = set_ids($final_genes);
write_to_gtf_file($final_genes,$output_gtf_file,$analysis_name);

exit;


=head2 process_gtf_file

Read GTF file.

=cut
sub process_gtf_file {
    my ($input_gtf_file, $slice_hash, $analysis) = @_;
    my %strand_conversion = ('+' => '1', '-' => '-1', '.' => undef, '?' => undef);
    my $exons = [];
    my $transcripts = [];
    my $current_biotype = "";
    my $current_gene_id = "";
    my $current_transcript_id = "";
    my $current_translation_coords = "";
    my $current_gtf_region_name = "";
    my $gtf_region_name;
    my $new_record = 0;
    my $current_translation_coords = "";
    my $next_translation_coords = "";
    my $next_biotype = '';
    my $exon;
    my $biotype = 'not_set';
    my $genes_by_biotype = {};


    #say "Reading in GTF";
    $LoggerConfig::logger->info("Reading in GTF");
    open(IN,$input_gtf_file) or die "Could not open GTF file: $input_gtf_file $!";
    while (<IN>) {
        my $current_line = $_;
        my @columns = split("\t", $current_line);
        next unless scalar(@columns) == 9;
        $gtf_region_name = $columns[0];
        unless($slice_hash->{$gtf_region_name});
        my $slice = $slice_hash->{$gtf_region_name};
        #skip if the seq region name is different from the input one
        next if $gtf_region_name ne $region_name;
        my $type = $columns[2];
        my $start  = $columns[3];
        my $end    = $columns[4];
        my $strand = $strand_conversion{$columns[6]};  
        unless($strand) {
          throw("Issue with parsing strand")
        }
        my $attributes = set_attributes( $columns[8] );
        my $gene_id       = $attributes->{'gene_id'};
        my $transcript_id = $attributes->{'transcript_id'};
        
        if($type eq 'transcript' && !$current_transcript_id) {
          # This if for the very first transcript
          $current_transcript_id = $transcript_id;
          $current_gene_id = $gene_id;
          $biotype = $attributes->{'biotype'};
          $translation_coords = $attributes->{'translation_coords'};
          if($translation_coords) {
            unless($biotype) {
              $biotype = 'protein_coding';
            }
            $current_translation_coords = $translation_coords;
          }
          unless($biotype) {
            $biotype = 'not_set';
          }
          $current_biotype = $biotype;
          next;
        } elsif ($type eq 'transcript') {
          # This is for moving onto a new transcript, process the old one
          if($$exons[0]->strand() == 1) {
            $exons = [sort { $a->start <=> $b->start } @{$exons}];
          } else {
            $exons = [sort { $b->start <=> $a->start } @{$exons}];
          }
          my $transcript = set_transcript($exons, $current_transcript_id, $current_biotype, $$exons[0]->slice())
          $transcript->analysis($analysis);

          if($current_translation_coords) {
            build_translation($transcript,$current_translation_coords);
          } else {
            compute_translation($transcript);
          }
          push(@$transcripts,$transcript);
          
          # Now set the current variables to the ones from the new transcript
          $current_transcript_id = $transcript_id;
          $current_gene_id = $gene_id;
          my $biotype = $attributes->{'biotype'};
          my $translation_coords = $attributes->{'translation_coords'};
          if($translation_coords) {
            unless($biotype) {
              $biotype = 'protein_coding';
            }
            $current_translation_coords = $translation_coords;
          } else {
            $current_translation_coords = "";
          }

          unless($biotype) {
            $biotype = 'not_set';
          }
          $current_biotype = $biotype;
          $exons = [];
        } elsif($gtf_type eq 'exon') {
        # Build an exon and add to the current array of exons
        my $exon = set_exon($start, $end, $strand, $slice);
        push(@$exons,$exon);
        } else {
          throw("Found an unexpected type in the GTF file, expected transcript or exon, found: ".$gtf_type);
        }
    }
      return $exons;

      $LoggerConfig::logger->info("Finished reading GTF");
  }



# sub write_to_gtf_file {
#   my ($genes_to_write,$output_gtf_file,$analysis_name) = @_;

#   open(OUT1,">".$output_gtf_file);
#   my $gene_count = 1;
#   my $transcript_count = 1;
#   foreach my $gene (@$genes_to_write) {
#     my $gene_id = "gene_".$gene_count;
#     my $transcripts = $gene->get_all_Transcripts();
#     foreach my $transcript (@$transcripts) {
#       my $transcript_id = "transcript_".$transcript_count;
#       my $exons = $transcript->get_all_Exons();
#       my $translation = $transcript->translation();
#       my $record = build_gtf_record($transcript,$exons,$analysis_name,$gene_id,$transcript_id,$translation);
#       foreach my $line (@$record) {
#         say OUT1 $line;
#       }
#       $transcript_count++;
#     }
#     $gene_count++;
#   }
#   close OUT1;
# }




#sub set_attributes {
#  my ($attribute_string) = @_;
#  my $attribute_pairs = {};

#  my @attribute_array = split(";",$attribute_string);
#  foreach my $attribute (@attribute_array) {
#    my @pairs = split(" ",$attribute);
#    if(scalar(@pairs) == 2) {
#      $pairs[1] =~ s/\"//g;
#      $attribute_pairs->{$pairs[0]} = $pairs[1];
#    }
#  }

#  return($attribute_pairs);
#}

=head2 set_ids

Set unique gene/transcript stable ids. This is used initially on the slices to help process genes
and then re-applied on the final gene set to make completely unique ids

=cut

sub set_ids {
  my ($genes) = @_;

  my $gene_id_index = 1;
  my $transcript_id_index = 1;
  foreach my $gene (@$genes) {
    $gene->stable_id("gene_".$gene_id_index);
    $gene_id_index++;
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      $transcript->stable_id("transcript_".$transcript_id_index);
      $transcript_id_index++;
    }
  }
  return($genes);
}

=head2 set_gene_biotypes

Set the biotypes of the genes and transcripts based on the presence of a translation

=cut

sub set_gene_biotypes {
  my ($genes) = @_;

  # Note: This runs off the assumption that the input is just lncRNAs and protein coding genes
  foreach my $gene (@$genes) {
    $gene->biotype('lncRNA');
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      if($transcript->translation()) {
        $gene->biotype('protein_coding');
        $transcript->biotype('protein_coding');
      } else {
        $transcript->biotype('lncRNA');
      }
    }
  }
  return($genes);
}


=head2 collapse_lncrna_genes

Collapses overlapping lncRNA genes into single representative genes by 
running a gene-building process and returns a combined list of 
non-overlapping genes.

=cut

sub collapse_lncrna_genes {
    my ($genes, $analysis) = @_;
    my $collapsed_genes = []; 
    my $lncrna_genes = []; 

    # Separate non-lncRNA genes and lncRNA genes
    foreach my $gene (@$genes) {
        if ($gene->biotype() eq 'lncRNA') {
            push @$lncrna_genes, $gene;
        } else {
            push @$collapsed_genes, $gene;
        }
    }

    # Return non-lncRNA genes if no lncRNA genes are present
    return $collapsed_genes unless @$lncrna_genes;

    # Get slice and biotype information from the first lncRNA gene
    my $slice = $lncrna_genes->[0]->slice();
    my $biotype = $lncrna_genes->[0]->biotype();

    # Run the gene builder to collapse lncRNA genes
    my $collapsed_lncrna_genes = run_gene_builder_for_slice($slice, $lncrna_genes, $analysis, $biotype);

    # Add the collapsed lncRNA genes to the final list
    push @$collapsed_genes, @$collapsed_lncrna_genes;

    return $collapsed_genes;
}


=head2 remove_potential_readthroughs

Filters out potential readthrough transcripts in protein-coding genes 
by checking for introns that overlap coding exons beyond a specified 
threshold. If readthroughs are detected, a gene-building process is 
invoked to create clean single-transcript genes. Non-protein-coding 
genes are returned without modification.


=cut

sub remove_potential_readthroughs {
    my ($genes, $analysis) = @_;

    my $max_coding_exon_skip = 3;  # Maximum allowed skipped coding exons
    my $filtered_genes = [];       # Array to store filtered genes

    foreach my $gene (@$genes) {
        # Skip filtering for non-protein-coding genes
        unless ($gene->biotype() eq 'protein_coding') {
            push @$filtered_genes, $gene;
            next;
        }

        my $canonical_transcript = $gene->canonical_transcript();
        my $filtered_transcripts = [];
        my $transcripts = $gene->get_all_Transcripts();
        my $coding_exons = get_unique_coding_exons($transcripts);

        # Filter transcripts based on intron overlap with coding exons
        foreach my $transcript (@$transcripts) {
            my $remove_transcript = 0;
            my $introns = $transcript->get_all_Introns();

            # Check intron overlap with coding exons
            foreach my $intron (@$introns) {
                my $skipped_coding_exon_count = 0;
                foreach my $coding_exon (@$coding_exons) {
                    if (features_overlap($intron, $coding_exon)) {
                        $skipped_coding_exon_count++;
                        if ($skipped_coding_exon_count > $max_coding_exon_skip) {
                            $remove_transcript = 1;
                            last;
                        }
                    }
                }
                last if $skipped_coding_exon_count > $max_coding_exon_skip;
            }

            # Keep the transcript if it does not exceed the exon skip threshold
            push @$filtered_transcripts, $transcript unless $remove_transcript;
        }

        # Ensure at least one transcript is kept by adding the canonical transcript if necessary
        push @$filtered_transcripts, $canonical_transcript if @$filtered_transcripts == 0;

        # If transcripts were removed, rebuild the gene using the gene builder
        if (scalar(@$filtered_transcripts) != scalar(@$transcripts)) {
            my $single_transcript_genes = create_single_transcript_genes($filtered_transcripts);
            my $slice = $single_transcript_genes->[0]->slice();
            my $biotype = $single_transcript_genes->[0]->biotype();
            my $analysis = $single_transcript_genes->[0]->analysis();
            my $collapsed_genes = run_gene_builder_for_slice($slice, $single_transcript_genes, $analysis, $biotype);
            push @$filtered_genes, @$collapsed_genes;
        } else {
            # Add the original gene if no transcripts were removed
            push @$filtered_genes, $gene;
        }
    }

    return $filtered_genes;
}


=head2 get_unique_coding_exons

Returns a list of unique coding exons from a list of transcripts.

=cut

sub get_unique_coding_exons {
  my ($transcripts) = @_;

  my $seen_coding_exons = {};
  my $coding_exons = [];
  foreach my $transcript (@$transcripts) {
    my $cds_exons = $transcript->get_all_CDS();
    foreach my $cds_exon (@$cds_exons) {
      my $cds_exon_string = $cds_exon->start.":".$cds_exon->end;
      unless($seen_coding_exons->{$cds_exon_string}) {
        push(@$coding_exons,$cds_exon);
        $seen_coding_exons->{$cds_exon_string} = 1;
      }
    }
  }
  return($coding_exons);
}

=head2 remove_overlapping_lncrnas

Remove overlapping lncRNAs. This is done by sorting the genes by strand and then removing
any lncRNAs that overlap with protein coding genes.

=cut

sub remove_overlapping_lncrnas {
  my ($genes) = @_;

  my ($forward_genes,$reverse_genes) = sort_genes_by_strand($genes);

  my $filtered_forward_genes = remove_overlapping_lncrnas_by_strand($forward_genes);
  my $filtered_reverse_genes = remove_overlapping_lncrnas_by_strand($reverse_genes);

  return([@$filtered_forward_genes,@$filtered_reverse_genes]);
}

=head2 remove_overlapping_lncrnas_by_strand

Remove overlapping lncRNAs by strand. This is done by iterating through the genes and removing.

=cut
sub remove_overlapping_lncrnas_by_strand {
  my ($stranded_genes) = @_;

  my $gene_ids_to_remove = {};
  my $filtered_genes = [];
  for(my $i=0; $i<scalar(@$stranded_genes)-1; $i++) {
    my $gene_left = ${$stranded_genes}[$i];
    next if($gene_ids_to_remove->{$gene_left->stable_id});

    my $biotype_left = $gene_left->biotype();
    for(my $j = $i+1; $j<scalar(@$stranded_genes); $j++) {
      my $gene_right = ${$stranded_genes}[$j];
      next if($gene_ids_to_remove->{$gene_right->stable_id});

      my $biotype_right = $gene_right->biotype();
      last unless(features_overlap($gene_left,$gene_right));

      if(features_overlap($gene_left,$gene_right)) {
        if($biotype_left eq $biotype_right) {
          next;
        } if($biotype_left eq 'lncRNA') {
          $gene_ids_to_remove->{$gene_left->stable_id()} = 1;
        } else {
          $gene_ids_to_remove->{$gene_right->stable_id()} = 1;
        }
      }
    } # for(my $j = $i+1;
  } # End for(my $i=0;

  foreach my $gene (@$stranded_genes) {
    unless($gene_ids_to_remove->{$gene->stable_id()}) {
      push(@$filtered_genes,$gene);
    }
  }
  return($filtered_genes);
}


=head2 sort_genes_by_strand

Sort genes by strand and then sort the exons within the genes by start position

=cut

sub sort_genes_by_strand {
  my ($genes) = @_;

  # Separate forward strand genes and sort them by start position in ascending order
  my @sorted_forward_genes = sort { $a->start() <=> $b->start() } grep { $_->strand == 1 } @$genes;

  # Separate reverse strand genes and sort them by start position in ascending order
  my @sorted_reverse_genes = sort { $a->start() <=> $b->start() } grep { $_->strand == -1 } @$genes;

  # Return references to the sorted arrays for forward and reverse strand genes
  return (\@sorted_forward_genes, \@sorted_reverse_genes);
}

=head2 fetch_slices_from_fasta

Get slices from a fasta file.

=cut
sub fetch_slices_from_fasta {
  my ($genome_file) = @_;

  my @slice_info = ();
  my $slices = [];
  my ($current_slice_name, $current_seq) = (undef, '');
  open(IN,$genome_file) or die "Could not open GTF file: $genome_file $!";
    while (<IN>) {
        my $line = $_;
        chomp($line);
        if($line =~ /\>(.+)/ && defined $current_slice_name) {
          my $new_slice_name = $1;
          my $slice_length = length($current_seq);
          push(@slice_info,[$current_slice_name,$slice_length]);
          $current_seq = "";
          $current_slice_name = $new_slice_name;
          
        } elsif($line =~ /\>(.+)/) {
          $current_slice_name = $1;
        } else {
          $current_seq .= $line;
        }
  }
  close IN;

  if($current_seq) {
    my $slice_length = length($current_seq);
    push(@slice_info,[$current_slice_name,$slice_length]);
  }

  foreach my $slice_details (@slice_info) {
    my $slice = set_slice {
    my ($region_name, $region_length) = @$slice_details;
    my $slice = set_slice($region_name, 1, $region_end, 1);
    push(@$slices,$slice);
  }
  return($slices);
  }
}


=head2 set_canonical_transcript

Set the canonical transcript for each gene based on specific criteria:
- Prefer transcripts with a translation.
- Among coding transcripts, choose the one with the longest translation.
- If translations are of equal length, choose the transcript with the longest total length.
- If no coding transcripts are available, choose the longest non-coding transcript.

=cut

sub set_canonical_transcripts {
    my ($genes) = @_;

    foreach my $gene (@$genes) {
        my $transcripts = $gene->get_all_Transcripts();
        my $current_canonical = pop(@$transcripts);  # Initialize with the last transcript

        foreach my $transcript (@{$gene->get_all_Transcripts()}) {
            my $translation = $transcript->translation();

            if ($translation) {
                # Prefer coding transcripts over non-coding
                if (!$current_canonical->translation()) {
                    $current_canonical = $transcript;
                }
                # Among coding transcripts, prefer the one with the longest translation
                elsif ($current_canonical->translation->length() < $translation->length()) {
                    $current_canonical = $transcript;
                }
                # If translation lengths are equal, prefer the one with the longest total length
                elsif ($current_canonical->translation->length() == $translation->length() 
                        && $current_canonical->length() < $transcript->length()) {
                    $current_canonical = $transcript;
                }
            } 
            # Handle non-coding transcripts: prefer the longest one
            elsif (!$current_canonical->translation() 
                    && $current_canonical->length() < $transcript->length()) {
                $current_canonical = $transcript;
            }
        }

        # Set the chosen transcript as canonical for the current gene
        $gene->canonical_transcript($current_canonical);
    }

    return $genes;
}
