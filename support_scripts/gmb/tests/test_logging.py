#!/usr/bin/env python3
"""Tests for GMB logging setup."""

from gmb.utils.logging import resolve_log_file, restore_stdio, setup_logging


def test_resolve_log_file_defaults_to_output_dir(tmp_path):
    assert resolve_log_file(str(tmp_path), None) == str(tmp_path / "gmb.log")


def test_setup_logging_captures_stdout(tmp_path):
    log_file = tmp_path / "gmb.log"
    setup_logging(log_file=str(log_file), capture_stdio=True)
    try:
        print("captured progress line")
    finally:
        restore_stdio()

    assert "captured progress line" in log_file.read_text()

