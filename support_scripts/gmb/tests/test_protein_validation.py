#!/usr/bin/env python3
"""Tests for protein validation logic."""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from config import load_config
from protein_validation import check_dependencies


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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
