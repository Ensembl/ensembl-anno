#!/usr/bin/env python3
"""Configuration system — legacy wrapper.

All logic now lives in :mod:`gmb.pipeline.config`.
"""
import os
import sys

_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)

# Re-export everything for backwards compatibility
from gmb.pipeline.config import (  # noqa: E402, F401
    DedupConfig,
    ExportConfig,
    HelixerFilterConfig,
    OrfConfig,
    PipelineConfig,
    ProteinFilterConfig,
    ProteinValidationConfig,
    QcConfig,
    ReportingConfig,
    ScoringConfig,
    ScoringWeights,
    TranscriptomicFilterConfig,
    TranscriptSplittingConfig,
    UtrConfig,
    ValidationConfig,
    load_config,
)
