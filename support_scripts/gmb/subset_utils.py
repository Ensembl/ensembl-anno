#!/usr/bin/env python3
"""Subset utils — legacy wrapper. All logic now in :mod:`gmb.pipeline.subset_utils`."""
import os, sys  # noqa: E401
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)
from gmb.pipeline.subset_utils import (  # noqa: E402, F401
    Region,
    add_subset_args,
    build_mapping,
    load_assembly_mapping,
    load_regions_file,
    load_seqname_map,
    parse_region,
    remap_df_seqnames,
    resolve_subset_regions,
    sample_loci,
    subset_df_by_regions,
    write_subset_manifest,
    _build_loci_from_exons,
)
