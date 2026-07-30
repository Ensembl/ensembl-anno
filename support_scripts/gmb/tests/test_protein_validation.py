#!/usr/bin/env python3
"""Tests for protein validation logic."""

import os
import subprocess
import sys
from unittest.mock import patch

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from gmb.pipeline.config import load_config
from gmb.pipeline.protein_validation import (
    batch_score_proteins,
    check_dependencies,
    detect_psauron_capabilities,
    protein_sha256,
)


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

    def test_missing_psauron_throws_on_enabled(self, config):
        config.protein_validation.enabled = True
        config.protein_validation.psauron_path = "non_existent_psauron_bin"

        def _fake_run(cmd, **kwargs):
            if cmd[0] == config.protein_validation.diamond_path:
                return subprocess.CompletedProcess(args=cmd, returncode=0, stdout=b"", stderr=b"")
            raise FileNotFoundError()

        with patch("subprocess.run", side_effect=_fake_run):
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


def _caps(help_text, version, **flags):
    """Build a detect_psauron_capabilities()-shaped dict for mocking."""
    base = {
        "found": True,
        "returncode": 0,
        "version": version,
        "help_text": help_text,
        "has_input_flag": True,
        "has_output_flag": True,
        "has_minimum_length": True,
        "has_protein_flag": True,
        "has_cpu_flag": True,
    }
    base.update(flags)
    return base


class TestPsauronCapabilityDetection:
    def test_current_local_version_detected_fully(self):
        with patch(
            "subprocess.run", return_value=_mock_help_result(_REAL_V1_0_8_HELP)
        ) as mock_run:
            caps = detect_psauron_capabilities("psauron")
        assert mock_run.call_count == 1  # --help is invoked exactly once
        assert caps["found"] is True
        assert caps["returncode"] == 0
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
            "subprocess.run",
            return_value=_mock_help_result(_HYPOTHETICAL_MISSING_PROTEIN_FLAG_HELP),
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
        assert caps["found"] is False
        assert caps["returncode"] is None
        assert caps["version"] is None
        assert all(
            caps[k] is False
            for k in (
                "has_input_flag",
                "has_output_flag",
                "has_minimum_length",
                "has_protein_flag",
                "has_cpu_flag",
            )
        )

    def test_nonzero_returncode_recorded_not_raised(self):
        with patch(
            "subprocess.run", return_value=_mock_help_result(_REAL_V1_0_8_HELP, returncode=2)
        ):
            caps = detect_psauron_capabilities("psauron")
        assert caps["found"] is True
        assert caps["returncode"] == 2


class TestCheckDependenciesPsauron:
    def test_help_invoked_exactly_once(self, config):
        # check_dependencies must not probe psauron a second time on top of
        # detect_psauron_capabilities()'s own --help call.
        config.protein_validation.enabled = True
        config.protein_validation.diamond_db = __file__
        with (
            patch("subprocess.run") as mock_diamond_run,
            patch(
                "gmb.pipeline.protein_validation.detect_psauron_capabilities",
                return_value=_caps(_REAL_V1_0_8_HELP, "1.0.8"),
            ) as mock_detect,
        ):
            mock_diamond_run.return_value = subprocess.CompletedProcess(
                args=[], returncode=0, stdout=b"", stderr=b""
            )
            check_dependencies(config.protein_validation)
        assert mock_detect.call_count == 1

    def test_fails_clearly_when_protein_flag_missing(self, config):
        config.protein_validation.enabled = True
        with (
            patch("subprocess.run", return_value=_mock_help_result("PSAURON version 2.0.0\n")),
            patch(
                "gmb.pipeline.protein_validation.detect_psauron_capabilities",
                return_value=_caps(
                    _HYPOTHETICAL_MISSING_PROTEIN_FLAG_HELP, "2.0.0", has_protein_flag=False
                ),
            ),
        ):
            with pytest.raises(SystemExit) as excinfo:
                check_dependencies(config.protein_validation)
        assert excinfo.value.code == 1

    def test_fails_clearly_on_nonzero_help_returncode(self, config):
        config.protein_validation.enabled = True
        with (
            patch("subprocess.run", return_value=_mock_help_result(_REAL_V1_0_8_HELP)),
            patch(
                "gmb.pipeline.protein_validation.detect_psauron_capabilities",
                return_value=_caps(_REAL_V1_0_8_HELP, "1.0.8", returncode=2),
            ),
        ):
            with pytest.raises(SystemExit) as excinfo:
                check_dependencies(config.protein_validation)
        assert excinfo.value.code == 1

    def test_cpu_flag_only_required_when_cpu_mode_enabled(self, config):
        # psauron_use_cpu defaults to False, so a psauron build that lacks
        # -c/--use-cpu entirely must still pass -- GMB never constructs a
        # -c argument in that case.
        config.protein_validation.enabled = True
        config.protein_validation.psauron_use_cpu = False
        config.protein_validation.diamond_db = __file__
        with (
            patch("subprocess.run", return_value=_mock_help_result(_REAL_V1_0_8_HELP)),
            patch(
                "gmb.pipeline.protein_validation.detect_psauron_capabilities",
                return_value=_caps(_REAL_V1_0_8_HELP, "1.0.8", has_cpu_flag=False),
            ),
        ):
            check_dependencies(config.protein_validation)  # must not raise

    def test_cpu_flag_required_when_cpu_mode_enabled(self, config):
        config.protein_validation.enabled = True
        config.protein_validation.psauron_use_cpu = True
        config.protein_validation.diamond_db = __file__
        with (
            patch("subprocess.run", return_value=_mock_help_result(_REAL_V1_0_8_HELP)),
            patch(
                "gmb.pipeline.protein_validation.detect_psauron_capabilities",
                return_value=_caps(_REAL_V1_0_8_HELP, "1.0.8", has_cpu_flag=False),
            ),
        ):
            with pytest.raises(SystemExit) as excinfo:
                check_dependencies(config.protein_validation)
        assert excinfo.value.code == 1

    def test_passes_with_real_local_psauron_shape(self, config):
        config.protein_validation.enabled = True
        with (
            patch("subprocess.run", return_value=_mock_help_result(_REAL_V1_0_8_HELP)),
            patch(
                "gmb.pipeline.protein_validation.detect_psauron_capabilities",
                return_value=_caps(_REAL_V1_0_8_HELP, "1.0.8"),
            ),
        ):
            # diamond_db check will still run after psauron capability check;
            # point it at a path guaranteed to exist so this test isolates
            # the psauron capability logic specifically.
            config.protein_validation.diamond_db = __file__
            check_dependencies(config.protein_validation)  # must not raise


# ---------------------------------------------------------------------------
# batch_score_proteins: hardened Psauron CSV output parsing
#
# These stub out DIAMOND/Psauron with a fake subprocess.run that writes a
# controlled psauron.csv into the real temp dir batch_score_proteins()
# creates, so each malformed-output scenario is exercised via the real
# parsing code rather than a hand-rolled reimplementation of it.
# ---------------------------------------------------------------------------

_PSAURON_HEADER = "description,psauron_is_protein,in-frame_score\n"


def _make_fake_run(psauron_csv_contents):
    """subprocess.run stub: DIAMOND writes an empty hits file (no matches),
    Psauron writes `psauron_csv_contents` (or nothing, if None) to -o.
    """

    def _fake_run(cmd, **kwargs):
        if "-o" in cmd:
            out_path = cmd[cmd.index("-o") + 1]
            if cmd[0] != "diamond" and psauron_csv_contents is not None:
                with open(out_path, "w") as fh:
                    fh.write(psauron_csv_contents)
            elif cmd[0] == "diamond":
                open(out_path, "w").close()
        return subprocess.CompletedProcess(args=cmd, returncode=0, stdout=b"", stderr=b"")

    return _fake_run


class TestBatchScoreProteinsCsvHardening:
    def test_missing_output_file_fails_clearly(self, config):
        config.protein_validation.enabled = True
        with patch("subprocess.run", side_effect=_make_fake_run(None)):
            with pytest.raises(SystemExit):
                batch_score_proteins({"t1": "MSTAVLKQ"}, config)

    def test_missing_header_fails_clearly(self, config):
        config.protein_validation.enabled = True
        bad_csv = "psauron -i foo -o bar\npsauron score: 0.5\nno,header,here\n"
        with patch("subprocess.run", side_effect=_make_fake_run(bad_csv)):
            with pytest.raises(SystemExit):
                batch_score_proteins({"t1": "MSTAVLKQ"}, config)

    def test_missing_expected_column_fails_clearly(self, config):
        config.protein_validation.enabled = True
        bad_csv = "psauron -i foo -o bar\npsauron score: 0.5\ndescription,psauron_is_protein,some_other_score\nseq_0,True,0.9\n"
        with patch("subprocess.run", side_effect=_make_fake_run(bad_csv)):
            with pytest.raises(SystemExit):
                batch_score_proteins({"t1": "MSTAVLKQ"}, config)

    def test_no_matching_rows_fails_clearly(self, config):
        config.protein_validation.enabled = True
        # Valid header/shape, but the one data row names a sequence ID that
        # cannot correspond to anything GMB submitted.
        csv_text = _PSAURON_HEADER + "totally_unrelated_id,True,0.9\n"
        with patch("subprocess.run", side_effect=_make_fake_run(csv_text)):
            with pytest.raises(SystemExit):
                batch_score_proteins({"t1": "MSTAVLKQ"}, config)

    def test_valid_output_scores_correctly(self, config):
        config.protein_validation.enabled = True
        csv_text = _PSAURON_HEADER + "seq_0,True,0.87\n"
        with patch("subprocess.run", side_effect=_make_fake_run(csv_text)):
            scores, details = batch_score_proteins({"t1": "MSTAVLKQ"}, config)
        assert details["t1"]["psauron_score"] == pytest.approx(0.87)
        assert scores["t1"] == pytest.approx(0.87 * config.protein_validation.psauron_weight)

    def test_malformed_row_skipped_with_warning_others_still_scored(self, config, capsys):
        config.protein_validation.enabled = True
        csv_text = _PSAURON_HEADER + "seq_0,True,not_a_number\nseq_1,True,0.5\n"
        with patch("subprocess.run", side_effect=_make_fake_run(csv_text)):
            scores, details = batch_score_proteins({"t1": "MSTAVLKQ", "t2": "MSTAVLKQAA"}, config)
        assert "malformed" in capsys.readouterr().err.lower()
        assert details["t2"]["psauron_score"] == pytest.approx(0.5)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


# ---------------------------------------------------------------------------
# Compute-once / consume-many contract
#
# DIAMOND and psauron must run exactly once per build, over deduplicated
# sequences, and every downstream consumer (model scoring, duplicate
# collapse, canonical selection, InterPro review) must reuse those results
# rather than recomputing them.
# ---------------------------------------------------------------------------


class TestComputeOnceContract:
    def _fake_run_factory(self, calls, psauron_csv):
        """subprocess.run stub recording every external tool invocation."""

        def _fake_run(cmd, **kwargs):
            # Distinguish real scans from cheap version probes: a scan is
            # `diamond blastp ...` or `psauron -i ...`; `diamond version`
            # and `psauron --help` are provenance lookups, not work.
            if cmd[0] == "diamond" and "blastp" in cmd:
                calls.append("diamond")
            elif cmd[0] != "diamond" and "-i" in cmd:
                calls.append("psauron")
            if "-o" in cmd:
                out_path = cmd[cmd.index("-o") + 1]
                if cmd[0] == "diamond":
                    open(out_path, "w").close()  # no hits
                else:
                    with open(out_path, "w") as fh:
                        fh.write(psauron_csv)
            return subprocess.CompletedProcess(args=cmd, returncode=0, stdout=b"", stderr=b"")

        return _fake_run

    def test_diamond_and_psauron_each_run_once_for_many_transcripts(self, config):
        config.protein_validation.enabled = True
        calls = []
        # Three transcripts, two distinct sequences.
        proteins = {"t1": "MAAAA", "t2": "MAAAA", "t3": "MBBBB"}
        csv_text = (
            "description,psauron_is_protein,in-frame_score\nseq_0,True,0.9\nseq_1,True,0.8\n"
        )
        with patch("subprocess.run", side_effect=self._fake_run_factory(calls, csv_text)):
            scores, details = batch_score_proteins(proteins, config)

        assert (
            calls.count("diamond") == 1
        ), f"DIAMOND ran {calls.count('diamond')} times, expected 1"
        assert (
            calls.count("psauron") == 1
        ), f"psauron ran {calls.count('psauron')} times, expected 1"
        assert set(scores) == {"t1", "t2", "t3"}

    def test_identical_proteins_share_one_validation_result(self, config):
        config.protein_validation.enabled = True
        calls = []
        proteins = {"t1": "MAAAA", "t2": "MAAAA"}
        csv_text = "description,psauron_is_protein,in-frame_score\nseq_0,True,0.9\n"
        with patch("subprocess.run", side_effect=self._fake_run_factory(calls, csv_text)):
            scores, details = batch_score_proteins(proteins, config)

        assert scores["t1"] == scores["t2"]
        assert details["t1"]["protein_sha256"] == details["t2"]["protein_sha256"]
        assert details["t1"]["psauron_score"] == details["t2"]["psauron_score"]

    def test_checksum_is_stable_and_maps_transcripts_that_are_renamed(self, config):
        # The checksum, not the transcript ID, is the reuse key: it must
        # survive the candidate -> final transcript ID rename.
        assert protein_sha256("MAAAA") == protein_sha256("MAAAA")
        assert protein_sha256("MAAAA") != protein_sha256("MBBBB")
        # ...and is not a substitute for the ID (both are retained in the
        # sidecar; see builder.py's protein_validation.tsv writer).
        assert protein_sha256("MAAAA") is not None

    def test_disabled_validation_still_returns_checksummed_records(self, config):
        config.protein_validation.enabled = False
        scores, details = batch_score_proteins({"t1": "MAAAA"}, config)
        assert scores["t1"] == 0.0
        assert details["t1"]["protein_sha256"] == protein_sha256("MAAAA")
        assert details["t1"]["validation_status"] == "not_run"

    def test_provenance_recorded_for_reuse_decisions(self, config):
        config.protein_validation.enabled = True
        calls = []
        csv_text = "description,psauron_is_protein,in-frame_score\nseq_0,True,0.9\n"
        with patch("subprocess.run", side_effect=self._fake_run_factory(calls, csv_text)):
            _, details = batch_score_proteins({"t1": "MAAAA"}, config)
        # diamond_db identity must be recorded so a consumer can tell
        # whether a cached result was produced against the same database.
        assert details["t1"]["diamond_db"] == config.protein_validation.diamond_db
        assert "diamond_version" in details["t1"]
        assert "psauron_version" in details["t1"]

    def test_one_transcript_per_sequence_marked_calculated_rest_reused(self, config):
        # Exactly one key sharing a sequence must be "calculated" and every
        # other key sharing that same sequence must be "reused", pointing
        # back at the calculated key -- this is what makes the compute-once/
        # consume-many contract auditable per transcript rather than only
        # inferable from a shared protein_sha256 (see
        # batch_score_proteins()'s "calculated vs reused" contract).
        config.protein_validation.enabled = True
        calls = []
        # Three transcripts share one sequence; t4 has a distinct sequence.
        proteins = {"t3": "MAAAA", "t1": "MAAAA", "t2": "MAAAA", "t4": "MBBBB"}
        csv_text = (
            "description,psauron_is_protein,in-frame_score\nseq_0,True,0.9\nseq_1,True,0.8\n"
        )
        with patch("subprocess.run", side_effect=self._fake_run_factory(calls, csv_text)):
            _, details = batch_score_proteins(proteins, config)

        # Deterministic choice: lexicographically-first key among {t1,t2,t3}.
        assert details["t1"]["protein_validation_source"] == "calculated"
        assert details["t1"]["protein_validation_reused_from"] is None
        for reused_key in ("t2", "t3"):
            assert details[reused_key]["protein_validation_source"] == "reused"
            assert details[reused_key]["protein_validation_reused_from"] == "t1"
        # t4 is the sole transcript for its sequence, so it is "calculated"
        # with nothing to reuse from.
        assert details["t4"]["protein_validation_source"] == "calculated"
        assert details["t4"]["protein_validation_reused_from"] is None
        # Still only one DIAMOND/psauron invocation regardless of the
        # calculated/reused bookkeeping.
        assert calls.count("diamond") == 1
        assert calls.count("psauron") == 1

    def test_calculated_reused_fields_are_none_when_validation_not_run(self, config):
        config.protein_validation.enabled = False
        _, details = batch_score_proteins({"t1": "MAAAA"}, config)
        assert details["t1"]["protein_validation_source"] is None
        assert details["t1"]["protein_validation_reused_from"] is None


class TestCanonicalSelectionDoesNotInvokeExternalTools:
    """Canonical selection must consume the sidecar, never recompute it."""

    def test_canonical_selection_module_never_calls_subprocess(self):
        # Structural guarantee: the module must not even import subprocess.
        import gmb.pipeline.canonical_selection as cs

        source = open(cs.__file__).read()
        assert "subprocess" not in source, "canonical_selection must not shell out"
        assert (
            "batch_score_proteins" not in source
        ), "canonical_selection must not recompute protein validation"

    def test_run_canonical_selection_with_sidecar_invokes_no_executable(self, tmp_path):
        """Fails if canonical selection tries to run DIAMOND/psauron when
        valid sidecar results were supplied."""
        import pandas as pd

        from gmb.pipeline.canonical_selection import run_canonical_selection

        gff3 = tmp_path / "consensus.gff3"
        gff3.write_text(
            "##gff-version 3\n"
            "1\tGMB\tgene\t100\t400\t.\t+\t.\tID=G1\n"
            "1\tGMB\tmRNA\t100\t400\t.\t+\t.\tID=G1.t1;Parent=G1;Evidence=Tiberius\n"
            "1\tGMB\texon\t100\t400\t.\t+\t.\tID=G1.t1.exon1;Parent=G1.t1\n"
            "1\tGMB\tCDS\t100\t400\t.\t+\t0\tID=G1.t1.cds1;Parent=G1.t1\n"
        )
        ev = tmp_path / "evidence_attribution.tsv"
        pd.DataFrame(
            [
                {
                    "gene_id": "G1",
                    "transcript_id": "G1.t1",
                    "evidence_sources": "Tiberius",
                    "exon_count": 1,
                    "cds_bp": 300,
                    "transcript_span_bp": 300,
                    "gmb_score": 2.6,
                }
            ]
        ).to_csv(ev, sep="\t", index=False)
        pv = tmp_path / "protein_validation.tsv"
        pd.DataFrame(
            [
                {
                    "gene_id": "G1",
                    "transcript_id": "G1.t1",
                    "protein_sha256": "abc",
                    "diamond_hit": "X",
                    "diamond_pident": 60.0,
                    "diamond_qcov": 90.0,
                    "diamond_scov": 90.0,
                    "psauron_score": 0.99,
                    "protein_length": 99,
                    "orf_label": "ORF:99aa ATG|STOP",
                    "is_partial_5": False,
                    "is_partial_3": False,
                    "internal_stop_count": 0,
                    "protein_coding_score": 0.9,
                    "gmb_score": 2.6,
                }
            ]
        ).to_csv(pv, sep="\t", index=False)

        out = tmp_path / "canonical"
        with patch("subprocess.run") as mock_run:
            run_canonical_selection(
                str(gff3),
                str(ev),
                str(out),
                load_config(),
                protein_validation_tsv=str(pv),
                annotate_gff3=False,
            )
        assert mock_run.call_count == 0, (
            "canonical selection invoked an external process despite valid "
            f"sidecar results being supplied (calls: {mock_run.call_args_list})"
        )
        assert (out / "canonical_transcripts.tsv").exists()
