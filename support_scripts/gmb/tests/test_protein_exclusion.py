import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))

import subprocess
import tempfile

from test_integration import _create_test_data


def test_protein_exclusion():
    with tempfile.TemporaryDirectory() as tmpdir:
        files = _create_test_data(tmpdir)
        output_dir = os.path.join(tmpdir, "output")
        script = os.path.join(os.path.dirname(os.path.dirname(__file__)), "gene_model_builder.py")
        cmd = [
            sys.executable,
            script,
            "--scallop",
            files["scallop"],
            "--stringtie",
            files["stringtie"],
            "--helixer",
            files["helixer"],
            "--orthodb",
            files["orthodb"],
            "--uniprot",
            files["uniprot"],
            "--genome",
            files["genome"],
            "--output-dir",
            output_dir,
            "--gene-prefix",
            "TEST",
        ]

        subprocess.run(cmd, capture_output=True, text=True)
        gff3_path = os.path.join(output_dir, "consensus.gff3")

        with open(gff3_path) as f:
            content = f.read()

        # No OrthoDB_ or UniProt_ prefixed mRNA IDs should exist
        assert "ID=OrthoDB_" not in content
        assert "ID=UniProt_" not in content
