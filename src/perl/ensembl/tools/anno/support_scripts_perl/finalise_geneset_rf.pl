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
use List::Util qw( min max );

use File::Spec::Functions;
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
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases clone_Transcript exon_overlap coding_exon_overlap);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(contains_internal_stops compute_translation);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;

use FindBin;
use lib "$FindBin::Bin/../conf";

use LoggerConfig;
use Utils qw(check_gtf_file setup_fasta set_slice create_analysis_object set_transcript set_exon set_translation set_canonical_transcript parse_region_details write_to_gtf_file
  build_gtf_record set_attributes build_translation run_gene_builder_for_slice );

my ($host, $port, $user, $pass, $dbname, $dna_host, $dna_port, $dna_user, $dna_dbname);
my $specify_strand;


my $analysis_name = "ensembl";
my $module_name = "Anno";

my $region_details = '';

my $write_single_transcript_genes = 1;
my $set_canonical = 1;
my $output_path;
my $genome_file;
my $input_gtf_file;
my $output_gtf_file;
my $biotypes_hash = ['transcriptomic','busco','protein'];
my $good_biotype;
my $bad_biotype;
my $genes_by_biotype = {};
$genes_by_biotype->{'transcriptomic'} = [];
$genes_by_biotype->{'busco'} = [];
$genes_by_biotype->{'protein'} = [];
my $current_gene_id = "";
my $biotype = 'not_set';
my $genes = [];

GetOptions( 'gtf_file=s'       => \$input_gtf_file,
            'region_details=s' => \$region_details,
            'specify_strand=s'  => \$specify_strand,
            'output_path=s' => \$output_path,
            'input_gtf_file=s' => \$input_gtf_file,
            'output_gtf_file=s' => \$output_gtf_file,
            'genome_file=s' => \$genome_file);


check_gtf_file(input_gtf_file)
if ($genome_file && -e $genome_file) {
  setup_fasta(-FASTA => $genome_file);
}


my $analysis = create_analysis_object($analysis_name, $module_name);
# Parse region details
my ($region_name, $region_start, $region_end) = parse_region_details($region_details);
my $slice = set_slice($region_name, $region_start, $region_end, 1);
my $transcripts,$genes_by_biotype = process_gtf_file($input_gtf_file, $region_name, $slice, $specify_strand, $analysis);

my $output_genes = process_genes($genes_by_biotype);
say "Found ".scalar(@$output_genes)." genes for output";

if($set_canonical) {
  foreach my $gene (@$output_genes) {
    set_canonical_transcript($gene);
  }
}

write_to_gtf_file($output_genes,$output_gtf_file,$analysis_name,1);
exit;

=head2 process_gtf_file

Read GTF file.

=cut

sub process_genes {
    my ($genes_by_biotype) = @_;
    my $final_genes = [];
    my $layered_genes = [];

    # Define biotypes and extract genes
    my %biotypes = (
        transcriptomic => 'transcriptomic',
        busco          => 'busco',
        protein        => 'protein',
    );

    my %genes_by_biotype = map {
        $_ => $genes_by_biotype->{$biotypes{$_}}
    } keys %biotypes;

    # Define $types_hash for clustering
    my $types_hash = {
        transcriptomic => ['transcriptomic'],
        busco          => ['busco'],
        protein        => ['protein'],
    };

    # Include all genes into the final set initially
    foreach my $biotype (keys %genes_by_biotype) {
        push(@$final_genes, @{$genes_by_biotype{$biotype}});
    }

    $LoggerConfig::logger->info("Clustering genes from input_dbs...");

    # Cluster genes for each biotype
    my %clusters_by_biotype;
    foreach my $biotype (keys %genes_by_biotype) {
        my ($clusters, $unclustered) = cluster_Genes_by_coding_exon_overlap(
            $genes_by_biotype{$biotype},
            { $biotype => [$biotypes{$biotype}] } # $types_hash{$biotype}
        );
        $clusters_by_biotype{$biotype} = [@$clusters, @$unclustered];
        $LoggerConfig::logger->info("Found " . scalar(@{$clusters_by_biotype{$biotype}}) . " $biotype clusters");

    }

    # Process transcriptomic clusters
    foreach my $transcriptomic_cluster (@{$clusters_by_biotype{transcriptomic}}) {
        process_cluster_overlap(
            $transcriptomic_cluster,
            $clusters_by_biotype{busco},
            $biotypes{transcriptomic},
            $biotypes{busco},
            $final_genes,
            0.9,
            $types_hash
        );


        process_cluster_overlap(
            $transcriptomic_cluster,
            $clusters_by_biotype{protein},
            $biotypes{transcriptomic},
            $biotypes{protein},
            $final_genes,
            0.8,
            $types_hash
        );
    }
    # At this point all the transcriptomic clusters have been processed and BUSCO/protein canonicals added as appropriate
    # Now clusters that don't have overlapping transcriptomic CDS data need to be processed. Cluster everything again and
    # ignore clusters with transcriptomic genes

    # Process busco clusters
    foreach my $busco_cluster (@{$clusters_by_biotype{busco}}) {
        process_cluster_overlap(
            $busco_cluster,
            $clusters_by_biotype{protein},
            $biotypes{busco},
            $biotypes{protein},
            $final_genes,
            0.8,
            $types_hash
        );
    }

    # Now both the transcriptomic clusters and BUSCO clusters have been augmented by canonicals that represent
    # overlapping CDS that was not captured in the original clustered genes. Once this is done, layering can be
    # performed to select the final gene set. Layering will prioritise:
    # - Augmented transcriptomic set
    # - Augmented BUSCO set
    # - Protein set

    # Perform final clustering and layering
    my ($clusters, $unclustered) = cluster_Genes_by_coding_exon_overlap($final_genes, $types_hash);

    foreach my $cluster (@$clusters) {
        my $tc_genes = $cluster->get_Genes_by_Set($biotypes{transcriptomic});
        my $bc_genes = $cluster->get_Genes_by_Set($biotypes{busco});
        my $pc_genes = $cluster->get_Genes_by_Set($biotypes{protein});

        if (scalar(@$tc_genes)) {
            push(@$layered_genes, @$tc_genes);
        } elsif (scalar(@$bc_genes)) {
            push(@$layered_genes, @$bc_genes);
        } else {
            push(@$layered_genes, @$pc_genes);
        }
    }

    # Add unclustered genes
    foreach my $singleton (@$unclustered) {
        push(@$layered_genes, @{$singleton->get_Genes()});
    }
    $LoggerConfig::logger->info("Found " . scalar(@layered_genes) . " genes post-layering");

    # Run GeneBuilder
    my $slice_genes = run_gene_builder_for_slice($slice, $layered_genes, $analysis, 'gbuild');
    return $slice_genes;
}

=head2 process_cluster_overlap

First get all the genes from the secondary cluster that have a coding exon overlap with the canonical 
from the primary cluster.
Then select the canonical from the overlapping secondary genes and calculate the overlap
If there is less than X percent overlap on the secondary canonical seq, add it.

=cut
sub process_cluster_overlap {
    my ($primary_cluster, $secondary_clusters, $primary_biotype, $secondary_biotype, $final_genes, $overlap_threshold,$types_hash) = @_;

    my $primary_genes = $primary_cluster->get_Genes();
    my $primary_canonical_gene = select_canonical_gene($primary_genes);
    my $primary_canonical_transcript = ${$primary_canonical_gene->get_all_Transcripts}[0];

    foreach my $secondary_cluster (@$secondary_clusters) {
        next unless features_overlap($primary_cluster, $secondary_cluster);

        my $secondary_genes = $secondary_cluster->get_Genes();
        my ($clusters, $unclustered) = cluster_Genes_by_coding_exon_overlap(
            [$primary_canonical_gene, @$secondary_genes],
            $types_hash
        );

        foreach my $canonical_cluster (@$clusters) {
            my $overlapping_secondary_genes = $canonical_cluster->get_Genes_by_Set($secondary_biotype);
            my $overlapping_primary_genes = $canonical_cluster->get_Genes_by_Set($primary_biotype);

            my $secondary_canonical_gene = select_canonical_gene($overlapping_secondary_genes);
            my $secondary_canonical_transcript = ${$secondary_canonical_gene->get_all_Transcripts}[0];

            if (scalar(@$overlapping_primary_genes) && $primary_canonical_transcript->translation()) {
                my $coding_overlap = coding_exon_overlap(
                    $primary_canonical_transcript,
                    $secondary_canonical_transcript
                );
                my $secondary_canonical_length = length($secondary_canonical_transcript->translateable_seq());

                if ($coding_overlap / $secondary_canonical_length <= $overlap_threshold) {
                    my $clone_gene = clone_Gene($secondary_canonical_gene);
                    $clone_gene->biotype($primary_biotype);
                    push(@$final_genes, $clone_gene);
                }
            } else {
                my $clone_gene = clone_Gene($secondary_canonical_gene);
                $clone_gene->biotype($primary_biotype);
                push(@$final_genes, $clone_gene);
            }
        }
    }
    return $final_genes;
}



=head2 features_overlap

Return 1 if there featureA overlaps feature B
Return 0 otherwise

=cut

sub features_overlap {

  my ($featureA,$featureB) = @_;

  if (($featureA->start() <= $featureB->end()) and ($featureA->end() >= $featureB->start())) {
    return 1;
  }
  return 0;
}


=head2 build_single_transcript_gene

Create one gene per transcript.

=cut
sub build_single_transcript_gene {
  my ($transcript) = @_;

  my $gene = Bio::EnsEMBL::Gene->new(
              -START     => $transcript->start(),
              -END       => $transcript->end(),
              -STRAND    => $transcript->strand,
              -SLICE     => $transcript->slice());
  $gene->add_Transcript($transcript);
  $gene->biotype($transcript->biotype());
  return($gene);
}


=head2 select_canonical_gene

Select canonical gene from a set of single transcript genes.

=cut

sub select_canonical_gene {
  my ($genes) = @_;

  my $current_canonical_transcript;
  my $current_canonical_gene;
  # This method is to select a canonical between a set of single transcript genes
  foreach my $gene (@$genes) {
    my $transcripts = $gene->get_all_Transcripts();
    if(scalar(@$transcripts) > 1) {
      throw("Found a gene with more than one transcript where processing the clusters");
    }

    my $transcript = ${$transcripts}[0];
    unless($current_canonical_transcript) {
      $current_canonical_transcript = $transcript;
      $current_canonical_gene = $gene;
      next;
    }

    my $translation = $transcript->translation();
    if($translation) {
      if(!($current_canonical_transcript->translation())) {
        $current_canonical_transcript = $transcript;
        $current_canonical_gene = $gene;
      } elsif($current_canonical_transcript->translation->length() < $translation->length()) {
        $current_canonical_transcript = $transcript;
        $current_canonical_gene = $gene;
      } elsif($current_canonical_transcript->translation->length() == $translation->length() && $current_canonical_transcript->length() < $transcript->length()) {
        $current_canonical_transcript = $transcript;
        $current_canonical_gene = $gene;
      }
    } elsif(!($current_canonical_transcript->translation()) && $current_canonical_transcript->length() < $transcript->length()) {
      $current_canonical_transcript = $transcript;
      $current_canonical_gene = $gene;
    }
  } # End foreach my $gene

  return($current_canonical_gene);
}


=head2 process_gtf_file

Read GTF file.

=cut
sub process_gtf_file {
    my ($input_gtf_file, $region_name, $slice, $specify_strand, $analysis) = @_;
    my %strand_conversion = ('+' => '1', '-' => '-1', '.' => undef, '?' => undef);
    my $exons = [];
    my $transcripts = [];
    
    my $current_transcript_id = "";
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
        next if $current_line =~ /^\#/;
        my @columns = split("\t", $current_line);
        $gtf_region_name = $columns[0];
        #skip if the seq region name is different from the input one
        next if $gtf_region_name ne $region_name;
        my $type = $columns[2];
        if($type ne 'exon' && $type ne 'transcript') {
          next;
        }
        my $start  = $columns[3];
        my $end    = $columns[4];
        my $strand = $strand_conversion{$columns[6]} // 1;  
        if ( $specify_strand && $specify_strand != $strand ) {
        next;
        }
        my $attributes = set_attributes( $columns[8] );
        my $gene_id       = $attributes->{'gene_id'};
        my $transcript_id = $attributes->{'transcript_id'};
        
        if($type eq 'transcript' && $current_transcript_id) {
          $next_biotype = $attributes->{'biotype'};
          $next_translation_coords = $attributes->{'translation_coords'};
          $new_record = 1;
          next;
        } elsif ($type eq 'transcript') {
          $biotype = $attributes->{'biotype'};
          $current_translation_coords = $attributes->{'translation_coords'};
          next;
        }
        $exon = set_exon($start, $end, $strand, $slice);
        # This is weak if the transcript id is not unique
        if ($new_record) {
        $new_record = 0;
        $exons = [ sort { $a->start <=> $b->start } @{$exons} ] if $exons->[0]->strand == 1;
        $exons = [ sort { $b->start <=> $a->start } @{$exons} ] if $exons->[0]->strand != 1;
        my $transcript = set_transcript($exons, $current_transcript_id, $biotype, $slice)
        if($current_translation_coords) {
          build_translation($transcript,$current_translation_coords);
        } else {
          compute_translation($transcript);
        }
        push(@$transcripts,$transcript);
        
        my $gene = build_single_transcript_gene($transcript);
        push(@{$genes_by_biotype->{$biotype}},$gene);
        $exons = [];
        push(@$exons,$exon);
        $current_translation_coords = $next_translation_coords;
        $biotype = $next_biotype;
      } else {
        push(@$exons,$exon);
      }
      $current_transcript_id = $transcript_id;
      $current_gtf_region_name = $gtf_region_name;
    }
    # Process remaining exons
        if ( scalar(@$exons) ) {
        $exons = [ sort { $a->start <=> $b->start } @{$exons} ] if $exons->[0]->strand == 1;
        $exons = [ sort { $b->start <=> $a->start } @{$exons} ] if $exons->[0]->strand != 1;
        my $transcript = set_transcript($exons, $current_transcript_id, $biotype, $slice)
        if($current_translation_coords) {
          build_translation($transcript,$current_translation_coords);
        } else {
          compute_translation($transcript);
        }
        push( @{$transcripts}, $transcript );    
        my $gene = build_single_transcript_gene($transcript);
        push(@{$genes_by_biotype->{$biotype}},$gene);  
    }
    return $transcripts,$genes_by_biotype;

    $LoggerConfig::logger->info("Finished reading GTF");
}
