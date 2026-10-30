import pytest 
import pandas as pd
from bin.src.python.ensembl.tools.anno.nextflow_utils.ig_collapse import \
    read_in_miniprot_gff, \
    pre_filter_alignments, \
    does_alignment_overlap_cluster, \
    update_cluster_info, \
    merge_adjacent_clusters, \
    merge_overlapping_clusters, \
    make_clusters, \
    select_top_alignment_per_cluster, \
    write_alignments_to_gtf, \
    collapse_alignments


class TestIGCollapse:

    def test_read_in_gff(self):
        expected_df = pd.DataFrame({
            'chr': ['primary_assembly:T2T-mGirTip1v1.2.pri:12:1:150165456:1',
                    'primary_assembly:T2T-mGirTip1v1.2.pri:12:1:150165456:1',
                    'primary_assembly:T2T-mGirTip1v1.2.pri:12:1:150165456:1'],
            'source': ['miniprot', 'miniprot', 'miniprot'],
            'feature_type': ['mRNA', 'mRNA', 'mRNA'],
            'start': [370708, 370708, 370708],
            'end': [370993, 370964, 370993],
            'score': [183, 167, 177],
            'strand': ['+', '+', '+'],
            'frame': ['.','.', '.'],
            'attributes': ['ID=MP000001;Rank=1;Identity=0.5625;Positive=0.6667;Frameshift=2;StopCodon=1;Target=IG_V_gene|M99641_IGHV|IMGT 4 98',
                           'ID=MP000002;Rank=1;Identity=0.5349;Positive=0.6512;Frameshift=1;StopCodon=1;Target=IG_V_gene|X60503_IGHV|IMGT 4 88',
                           'ID=MP000003;Rank=1;Identity=0.5521;Positive=0.6562;Frameshift=2;StopCodon=1;Target=IG_V_gene|HM855463_IGHV|IMGT 4 98'],
            'rank': ['Rank=1', 'Rank=1', 'Rank=1'],
            'target': ['M99641_IGHV','X60503_IGHV', 'HM855463_IGHV'],
            'identity': [0.5625, 0.5349, 0.5521]
        })

        pd.testing.assert_frame_equal(
            read_in_miniprot_gff('tests/python/data/mini_miniprot.gff'),
            expected_df
        )

    
    def test_filter_alignments(self):
        alignments = pd.DataFrame({
            'record': [1,2,3,4],
            'rank': ['Rank=1', 'Rank=2', 'Rank=1', 'Rank=2']
        })
        expected_out = pd.DataFrame({
            'record': [1,3],
            'rank': ['Rank=1', 'Rank=1']
        })
        pd.testing.assert_frame_equal(
            pre_filter_alignments(alignments),
            expected_out)

    

    def test_does_alignment_overlap_cluster(self):

        row = {
            'chr': 1,
            'start': 5,
            'end' : 10
        }
        long_row = {
            'chr': 1,
            'start': 1,
            'end': 10
        }
        short_row = {
            'chr': 1,
            'start': 5,
            'end': 6
        }
        wrong_chr = {
            'chr': 2,
            'start': 5,
            'end' : 10
        }
        no_overlap = {
            'chr': 1,
            'start': 100,
            'end': 700
        }

        alignment = {
            'chr': 1,
            'start': 4,
            'end': 7
        }

        assert does_alignment_overlap_cluster(row, alignment) == True
        assert does_alignment_overlap_cluster(alignment, row) == True
        assert does_alignment_overlap_cluster(long_row, alignment) == True
        assert does_alignment_overlap_cluster(short_row, alignment) == True
        assert does_alignment_overlap_cluster(wrong_chr, alignment) == False
        assert does_alignment_overlap_cluster(no_overlap, alignment) == False


    def test_update_cluster_info(self):

        test_df = pd.DataFrame({
            'start': [10],
            'end': [20],
            'other_info': ['hello']
        })
        cluster_info = {
            'start': 15,
            'end': 25,
            'df' : pd.DataFrame({
                'start': [15],
                'end': [25],
                'other_info': ['world']
            })
        }
        expected_out_df = pd.DataFrame({
            'start': [15, 10],
            'end': [25, 20],
            'other_info': ['world', 'hello']
        })

        for index, row in test_df.iterrows():
            output_cluster_info = update_cluster_info(row, cluster_info)

        assert output_cluster_info['start'] == 10
        assert output_cluster_info['end'] == 25
        pd.testing.assert_frame_equal(output_cluster_info['df'], expected_out_df, check_dtype=False)
        

    def test_merge_adjacent_clusters(self):
        adjacent_cluster_info_dict = {
            'start': 15,
            'end': 25,
            'df' : pd.DataFrame({
                'start': [15],
                'end': [25],
                'other_info': ['hello']
            })
        }
        cluster_info_dict = {
            'start': 17,
            'end': 30,
            'df' : pd.DataFrame({
                'start': [17],
                'end': [30],
                'other_info': ['world']
            })
        }
        expected_output_dict = {
            'start': 15,
            'end': 30,
            'df' : pd.DataFrame({
                'start': [15, 17],
                'end': [25, 30],
                'other_info': ['hello', 'world']
            })
        }
        output_dict = merge_adjacent_clusters(cluster_info_dict, adjacent_cluster_info_dict)
        assert output_dict['end'] == 30
        pd.testing.assert_frame_equal(output_dict['df'], expected_output_dict['df'], check_dtype=False)

        return None

    def test_merge_overlapping_clusters(self):
        input_clusters = {
            1 : {
                'chr': 2,
                'start': 500,
                'end': 600,
                'df': pd.DataFrame({
                    'chr': [2, 2],
                    'start': [500, 555],
                    'end' : [550, 600]
                })
            },
            2: {
                'chr': 1,
                'start': 10,
                'end': 15,
                'df': pd.DataFrame({
                    'chr': [1, 1],
                    'start': [10, 11],
                    'end' : [15, 15]
                })
            },
            3: {
                'chr': 1,
                'start': 5,
                'end': 12,
                'df': pd.DataFrame({
                    'chr': [1, 1],
                    'start': [5, 7],
                    'end' : [10, 12]
                })
            }
        }
        expected_merged_output = [
            {
                'chr': 1,
                'start': 5,
                'end': 15,
                'df': pd.DataFrame({
                    'chr': [1, 1, 1, 1],
                    'start': [5, 7, 10, 11],
                    'end' : [10, 12, 15, 15]
                })
            },
            {
                'chr': 2,
                'start': 500,
                'end': 600,
                'df': pd.DataFrame({
                    'chr': [2, 2],
                    'start': [500, 555],
                    'end' : [550, 600]
                })
            }]

        merged_output = merge_overlapping_clusters(input_clusters)

        assert len(merged_output) == len(expected_merged_output)

        assert merged_output[0]['chr'] == expected_merged_output[0]['chr']
        assert merged_output[0]['start'] == expected_merged_output[0]['start']
        assert merged_output[0]['end'] == expected_merged_output[0]['end']
        pd.testing.assert_frame_equal(merged_output[0]['df'], merged_output[0]['df'], check_dtype=False)

        assert merged_output[1]['chr'] == expected_merged_output[1]['chr']
        assert merged_output[1]['start'] == expected_merged_output[1]['start']
        assert merged_output[1]['end'] == expected_merged_output[1]['end']
        pd.testing.assert_frame_equal(merged_output[1]['df'], merged_output[1]['df'], check_dtype=False)


    def test_make_clusters(self):
        return None

    def test_select_top_alignment_per_cluster(self):
        return None

    def test_write_alignments_to_gtf(self):
        return None

    def test_collapse_alignments(self):
        return None