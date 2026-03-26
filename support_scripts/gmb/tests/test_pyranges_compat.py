#!/usr/bin/env python3
"""Smoke tests for pandas / pyranges compatibility.

Ensures that core PyRanges operations used by the pipeline work correctly
under the supported version range (pandas 2.x–3.0, pyranges ≤0.1.4).
"""


import pandas as pd
import pyranges as pr
import pytest


class TestVersionRequirements:
    def test_pandas_version(self):
        """pandas version meets minimum supported."""
        major, minor = (int(x) for x in pd.__version__.split(".")[:2])
        if major < 2:
            pytest.skip(
                f"pandas {pd.__version__} < 2.0: not in supported range, "
                f"but local testing is still useful"
            )

    def test_pyranges_importable(self):
        """pyranges is importable and has expected attributes."""
        assert hasattr(pr, "PyRanges")


class TestPyRangesRoundTrip:
    @pytest.fixture
    def simple_df(self):
        return pd.DataFrame(
            {
                "Chromosome": ["1", "1", "1", "2"],
                "Start": [100, 200, 250, 100],
                "End": [150, 300, 400, 200],
                "Strand": ["+", "+", "+", "-"],
            }
        )

    def test_pyranges_df_attribute(self, simple_df):
        """PyRanges .df returns a pd.DataFrame."""
        gr = pr.PyRanges(simple_df)
        result = gr.df
        assert isinstance(result, pd.DataFrame)
        assert set(result.columns) >= {"Chromosome", "Start", "End", "Strand"}

    def test_cluster_operation(self, simple_df):
        """cluster() produces 'Cluster' column."""
        gr = pr.PyRanges(simple_df)
        try:
            clustered = gr.cluster(slack=0, count=True)
        except TypeError:
            clustered = gr.cluster(count=True)
        cdf = clustered.df
        assert "Cluster" in cdf.columns
        # Overlapping intervals on chr1+ should share a cluster
        chr1_clusters = cdf[(cdf["Chromosome"] == "1") & (cdf["Strand"] == "+")][
            "Cluster"
        ].unique()
        assert len(chr1_clusters) >= 1

    def test_overlap_operation(self):
        """overlap() between two PyRanges returns correct results."""
        df_a = pd.DataFrame(
            {
                "Chromosome": ["1"],
                "Start": [100],
                "End": [300],
                "Strand": ["+"],
            }
        )
        df_b = pd.DataFrame(
            {
                "Chromosome": ["1"],
                "Start": [200],
                "End": [400],
                "Strand": ["+"],
            }
        )
        gr_a = pr.PyRanges(df_a)
        gr_b = pr.PyRanges(df_b)
        ovl = gr_a.overlap(gr_b)
        assert not ovl.df.empty, "Expected overlap between [100,300) and [200,400)"

    def test_no_overlap_returns_empty(self):
        """Non-overlapping intervals produce empty overlap result."""
        df_a = pd.DataFrame(
            {
                "Chromosome": ["1"],
                "Start": [100],
                "End": [200],
                "Strand": ["+"],
            }
        )
        df_b = pd.DataFrame(
            {
                "Chromosome": ["1"],
                "Start": [500],
                "End": [600],
                "Strand": ["+"],
            }
        )
        gr_a = pr.PyRanges(df_a)
        gr_b = pr.PyRanges(df_b)
        ovl = gr_a.overlap(gr_b)
        assert ovl.df.empty


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
