import argparse
import pandas as pd


def read_in_miniprot_gff(miniprot_gff_path):

    # Read in gtf as a pandas dataframe and set column names:
    miniprot_gtf = pd.read_csv(miniprot_gff_path, sep='\t', comment='#', header=None)
    miniprot_gtf.columns = ['chr', 'source', 'feature_type', 'start', 'end', 'score', \
                        'strand', 'frame', 'attributes']

    # Miniprot records as CDS and mRNA record for each alignment. Just keep the mRNA record as it
    # is easier to handle:
    miniprot_gtf = miniprot_gtf[miniprot_gtf['feature_type']=='mRNA']

    # Extract rank, target name and % identity from the attributes field:
    miniprot_gtf['rank'] = [[x for x in y.split(';') if 'Rank' in x][0] for y in miniprot_gtf['attributes']]
    miniprot_gtf['target'] = [[x.split('|')[1] for x in y.split(';') if 'Target' in x] for y in miniprot_gtf['attributes']]
    miniprot_gtf['target'] = miniprot_gtf['target'].explode()
    miniprot_gtf['identity'] = [[x.split('=')[1] for x in y.split(';') if 'Identity' in x] for y in miniprot_gtf['attributes']]
    miniprot_gtf['identity'] = pd.to_numeric(miniprot_gtf['identity'].explode())
    miniprot_gtf = miniprot_gtf.reset_index(drop=True)
    print(miniprot_gtf)
    return miniprot_gtf


def pre_filter_alignments(alignment_df):

    # Only keep Rank=1 alignments:
    alignment_df = alignment_df[alignment_df['rank']=='Rank=1']
    alignment_df = alignment_df.reset_index(drop=True)
    return alignment_df


def does_alignment_overlap_cluster(row, cluster_info):

    # Back out straight away if the chromosomes don't match
    if row['chr'] != cluster_info['chr']:
        return False

    # Extract coordinates to establish if there is an overlap:
    alignment_start = row['start']
    alignment_end = row['end']
    cluster_start = cluster_info['start']
    cluster_end = cluster_info['end']

    # If there is an overlap, this statement will be true:
    if (alignment_start <= cluster_end) and (cluster_start <= alignment_end):
        return True
    return False


def update_cluster_info(row, cluster_info):

    # Expand cluster interval if required:
    if row['start'] < cluster_info['start']:
        cluster_info['start'] = row['start']
    if row['end'] > cluster_info['end']:
        cluster_info['end'] = row['end']

    # Add the alignment record to the cluster specific dataframe:
    cluster_info['df'] = pd.concat([cluster_info['df'], row.to_frame().T]).reset_index(drop=True)

    return cluster_info


def merge_adjacent_clusters(cluster_info_dict, adjacent_cluster_info_dict):

    # cluster_info_dict is guaranteed to have a larger start coordinate than adjacent cluster_info_dict, so
    # no need to change the start coordinate.
    # Update the end coordinate if required:
    if cluster_info_dict['end'] > adjacent_cluster_info_dict['end']:
        adjacent_cluster_info_dict['end'] = cluster_info_dict['end']

    # Combine the two cluster dataframes (containing miniprot alignment records)
    adjacent_cluster_info_dict['df'] = pd.concat(
        [adjacent_cluster_info_dict['df'], cluster_info_dict['df']]
        ).reset_index(drop=True)
    return adjacent_cluster_info_dict


def merge_overlapping_clusters(cluster_dict):

    # Sort cluster dictionaries by chromosome and start coordinate
    sorted_clusters = sorted(cluster_dict.items(), key=lambda tup: (tup[1]["chr"], tup[1]["start"]))

    merged_clusters = []

    for cluster in sorted_clusters:

        cluster_info_dict = cluster[1]
        
        # If this is the first cluster, add it to the merged cluster list
        if len(merged_clusters)==0:
            merged_clusters.append(cluster_info_dict)
            continue       

        # Pull the last cluster out of the list of merged_clusters.
        # This will be the cluster we will see if we can merge with.
        # This logic works because we have sorted the cluster list by chr and start.
        adjacent_cluster_info_dict = merged_clusters[-1]

        # If the chromosomes don't match, we can't merge, so we just add the cluster as is to merged_clusters
        if adjacent_cluster_info_dict['chr'] != cluster_info_dict['chr']:
            merged_clusters.append(cluster_info_dict)
            continue

        # At this point, we know the chromosomes match, so we want to see if the intervals overlap.
        # Because we have sorted the cluster list by start coordinate, we know the start coordinate of
        # the current cluster comes after the start coordinate of the last merged cluster. Therefore,
        # the only way there can be an overlap with the last cluster in merged_clusters is if the start 
        # coordinate of the current cluster comes before the end coordinate of the last cluster in merged 
        # clusters:
        if cluster_info_dict['start'] <= adjacent_cluster_info_dict['end']:
            updated_cluster_info_dict = merge_adjacent_clusters(cluster_info_dict, adjacent_cluster_info_dict)
            merged_clusters[-1] = updated_cluster_info_dict
            continue

        # If we are here, there was no overlap, so we just add the cluster as is to merged_clusters
        merged_clusters.append(cluster_info_dict)

    return merged_clusters


def make_clusters(alignment_df):

    cluster_dict = {}
    cluster_id = 1

    for index, row in alignment_df.iterrows():
        for (cluster_id, cluster_info) in cluster_dict.items():

            # If the alignment overlaps a cluster, update the cluster dictionary to reflect this
            # and break out of the for loop so that no alignment is recorded in more than one cluster:
            if does_alignment_overlap_cluster(row, cluster_info):
                cluster_dict[cluster_id] = update_cluster_info(row, cluster_info)
                cluster_id += 1
                break
        else:
            # If there are no overlaps, add the alignment as a new cluster to the dictionary:
            cluster_dict[cluster_id] = {
                'chr': row['chr'],
                'start': row['start'],
                'end': row['end'],
                'df': pd.concat([row])
            }
            cluster_id += 1

    # The for loops above could generate a list of clusters where some clusters overlap
    # merge_overlapping_clusters consolidates overlapping clusters into one:
    merged_clusters = merge_overlapping_clusters(cluster_dict)
    return merged_clusters


def select_top_alignment_per_cluster(cluster_list):
    top_alignment_dfs = []

    for cluster in cluster_list:
        alignment_df = cluster['df']
        alignment_df = alignment_df.loc[alignment_df['identity'].idxmax()]
        top_alignment_dfs.append(alignment_df)

    top_alignment_df = pd.concat(top_alignment_dfs)
    return top_alignment_df 


def write_alignments_to_gtf(top_alignment_df, collapsed_gtf_path):
    gtf = top_alignment_df[['chr', 'source', 'feature_type', 'start', 'end', 'score', \
                            'strand', 'frame', 'attributes']]
    gtf.to_csv(collapsed_gtf_path, index=False, header=False, sep='\t')


def collapse_alignments(miniprot_gff_path, collapsed_gtf_path):

    # Read in gtf
    miniprot = read_in_miniprot_gff(miniprot_gff_path)

    # Filter out miniprot alignments that don't meet criteria
    miniprot = pre_filter_alignments(miniprot)

    # Cluster overlapping alignments
    cluster_list = make_clusters(miniprot)

    # For each cluster, pick the best alignment based on perc_identity
    top_alignment_df = select_top_alignment_per_cluster(cluster_list)

    # Write gtf to file
    write_alignments_to_gtf(top_alignment_df, collapsed_gtf_path)


def parse_args():
    parser = argparse.ArgumentParser(description="Arguments for script to cluster and collapse miniprot IGMT alignments")
    parser.add_argument("--miniprot_gff", nargs="+", help="Path to miniprot gtf containing IGMT alignments")
    parser.add_argument("--collapsed_gtf", help="Path to output gtf")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    collapse_alignments(args.miniprot_gtf, args.collapsed_gtf)