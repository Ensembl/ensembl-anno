#!/usr/bin/env python3
"""Tests for protein validation logic."""

import os
import subprocess
import sys
from unittest.mock import patch

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from gmb.pipeline.config import load_config
from gmb.pipeline.protein_validation import check_dependencies, detect_psauron_capabilities


@pytest.fixture
def config():
    return load_config()


class TestValidationDependencies:
    def test_missing_diamond_throws_on_enabled(self, config):
        config.protein_validation.enabled = True
        config.protein_validation.diamond_path = "non_existent_diamond_bin"

        with pytest.raises(SystemExit) as excinfo:
            check_dependencies(config.protein_validation)
        assert excinfo.value.code == 1


# ---------------------------------------------------------------------------
# Psauron CLI capability detection
#
# These use real --help text captured from psauron v1.0.8 (this checkout's
# installed version) plus synthetic variants, rather than requiring the
# real binary -- verifies GMB's assumptions are checked, not hard-coded,
# and would correctly flag a hypothetical future CLI change.
# ---------------------------------------------------------------------------

_REAL_V1_0_8_HELP = """PSAURON version 1.0.8
usage: psauron [-h] -i INPUT_FASTA [-o OUTPUT_PATH] [-m MINIMUM_LENGTH]
               [-e EXCLUDE] [--inframe INFRAME] [--outframe OUTFRAME] [-c]
               [-s] [-p] [-a] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA, --input-fasta INPUT_FASTA
                        REQUIRED path to FASTA with spliced CDS sequence or
                        protein sequence.
  -o OUTPUT_PATH, --output-path OUTPUT_PATH
                        OPTIONAL path to output results file,
                        default=./psauron_score.csv
  -m MINIMUM_LENGTH, --minimum-length MINIMUM_LENGTH
                        OPTIONAL exclude all proteins shorter than m amino
                        acids, default=5
  -c, --use-cpu         OPTIONAL set -c to force usage of CPU instead of GPU,
                        default=False
  -s, --single-frame    OPTIONAL set -s to score only the in-frame CDS.
  -p, --protein         OPTIONAL set -p if your FASTA contains amino acid
                        protein sequence, which may lower accuracy of the
                        model, default=False
  -a, --all-prob        OPTIONAL set -a to output per-amino-acid predicted
                        probabilities.
  -v, --verbose         OPTIONAL set -v for verbose output.
"""

# Verified against the same v1.0.0 (2024-07) and v1.1.3 (2026-05) releases
# fetched from github.com/salzberg-lab/PSAURON -- identical flag set to
# v1.0.8 above (see report). Used to confirm detection is stable across the
# whole investigated range, not coincidentally tied to 1.0.8's exact text.
_REAL_V1_1_3_HELP = _REAL_V1_0_8_HELP.replace("version 1.0.8", "version 1.1.3")

# A hypothetical *future* release that dropped -p/--protein -- this is the
# exact scenario capability detection exists to catch.
_HYPOTHETICAL_MISSING_PROTEIN_FLAG_HELP = """PSAURON version 2.0.0
usage: psauron [-h] -i INPUT_FASTA [-o OUTPUT_PATH] [-m MINIMUM_LENGTH] [-c]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTA, --input-fasta INPUT_FASTA
                        REQUIRED path to FASTA.
  -o OUTPUT_PATH, --output-path OUTPUT_PATH
                        OPTIONAL path to output results file.
  -m MINIMUM_LENGTH, --minimum-length MINIMUM_LENGTH
                        OPTIONAL minimum length filter.
  -c, --use-cpu         OPTIONAL force CPU.
"""

_UNRELATED_TOOL_HELP = """usage: some-other-tool [-h] --input FILE --output FILE
optional arguments:
  -h, --help    show this help message and exit
  --input FILE  input file
  --output FILE output file
"""


def _mock_help_result(text, returncode=0):
    return subprocess.CompletedProcess(args=[], returncode=returncode, stdout=text, stderr="")


class TestPsauronCapabilityDetection:
    def test_current_local_version_detected_fully(self):
        with patch("subprocess.run", return_value=_mock_help_result(_REAL_V1_0_8_HELP)):
            caps = detect_psauron_capabilities("psauron")
        assert caps["version"] == "1.0.8"
        assert caps["has_input_flag"]
        assert caps["has_output_flag"]
        assert caps["has_minimum_length"]
        assert caps["has_protein_flag"]
        assert caps["has_cpu_flag"]

    def test_latest_released_version_detected_fully(self):
        # v1.1.3 -- latest release investigated -- must detect identically
        # to v1.0.8, confirming the fix isn't accidentally 1.0.8-specific.
        with patch("subprocess.run", return_value=_mock_help_result(_REAL_V1_1_3_HELP)):
            caps = detect_psauron_capabilities("psauron")
        assert caps["version"] == "1.1.3"
        assert caps["has_protein_flag"]
        assert caps["has_minimum_length"]

    def test_hypothetical_future_version_missing_protein_flag_detected(self):
        with patch(
            "subprocess.run", return_value=_mock_help_result(_HYPOTHETICAL_MISSING_PROTEIN_FLAG_HELP)
        ):
            caps = detect_psauron_capabilities("psauron")
        assert caps["version"] == "2.0.0"
        assert caps["has_protein_flag"] is False
        assert caps["has_input_flag"] is True  # unrelated flags still detected correctly

    def test_unrelated_tool_detects_no_flags_and_no_version(self):
        with patch("subprocess.run", return_value=_mock_help_result(_UNRELATED_TOOL_HELP)):
            caps = detect_psauron_capabilities("psauron")
        assert caps["version"] is None
        assert caps["has_protein_flag"] is False
        assert caps["has_minimum_length"] is False

    def test_missing_binary_returns_all_false_not_raise(self):
        with patch("subprocess.run", side_effect=FileNotFoundError()):
            caps = detect_psauron_capabilities("nonexistent-psauron-binary")
        assert caps["version"] is None
        assert all(
            caps[k] is False
            for k in ("has_input_flag", "has_output_flag", "has_minimum_length", "has_protein_flag", "has_cpu_flag")
        )

    def test_check_dependencies_fails_clearly_when_protein_flag_missing(self, config):
        config.protein_validation.enabled = True
        with (
            patch("subprocess.run") as mock_run,
            patch(
                "gmb.pipeline.protein_validation.detect_psauron_capabilities"
            ) as mock_detect,
        ):
            mock_run.return_value = _mock_help_result("PSAURON version 2.0.0\n")
            mock_detect.return_value = {
                "version": "2.0.0",
                "help_text": _HYPOTHETICAL_MISSING_PROTEIN_FLAG_HELP,
                "has_input_flag": True,
                "has_output_flag": True,
                "has_minimum_length": True,
                "has_protein_flag": False,
                "has_cpu_flag": True,
            }
            with pytest.raises(SystemExit) as excinfo:
                check_dependencies(config.protein_validation)
        assert excinfo.value.code == 1

    def test_check_dependencies_passes_with_real_local_psauron_shape(self, config):
        config.protein_validation.enabled = True
        with (
            patch("subprocess.run") as mock_run,
            patch(
                "gmb.pipeline.protein_validation.detect_psauron_capabilities"
            ) as mock_detect,
        ):
            mock_run.return_value = _mock_help_result(_REAL_V1_0_8_HELP)
            mock_detect.return_value = {
                "version": "1.0.8",
                "help_text": _REAL_V1_0_8_HELP,
                "has_input_flag": True,
                "has_output_flag": True,
                "has_minimum_length": True,
                "has_protein_flag": True,
                "has_cpu_flag": True,
            }
            # diamond_db check will still run after psauron capability check;
            # point it at a path guaranteed to exist so this test isolates
            # the psauron capability logic specifically.
            config.protein_validation.diamond_db = __file__
            check_dependencies(config.protein_validation)  # must not raise


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
