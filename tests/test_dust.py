"""Tests for Dust repeat annotation orchestration."""

import os
import shutil
import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

from src.python.ensembl.tools.anno.repeat_annotation import dust


class FakePool:
    """Minimal multiprocessing pool test double."""

    def __init__(self):
        self.tasks = []
        self.closed = False
        self.joined = False

    def apply_async(self, func, args=()):
        """Record scheduled work without starting a worker process."""
        self.tasks.append((func, args))

    def close(self):
        """Record pool closure."""
        self.closed = True

    def join(self):
        """Record pool join."""
        self.joined = True


class TestDust(unittest.TestCase):
    """Tests for Dust runner setup."""

    def setUp(self):
        """Create a temporary output directory."""
        self.cwd = Path.cwd()
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        """Restore the working directory and remove temporary files."""
        os.chdir(self.cwd)
        shutil.rmtree(self.temp_dir)

    def test_run_dust(self):
        """Run Dust orchestration without requiring the external dustmasker executable."""
        genome_file = "path/to/your/genome.fasta"
        output_dir = self.temp_dir
        dust_bin = "dustmasker"
        num_threads = 2
        pool = FakePool()

        def write_annotation(dust_dir, *_args):
            """Write the merged annotation file that slice_output_to_gtf would create."""
            (dust_dir / "annotation.gtf").write_text(
                'seq1\tDust\trepeat\t1\t10\t.\t+\t.\trepeat_id "1";\n',
                encoding="utf8",
            )

        with ExitStack() as stack:
            mock_check_exe = stack.enter_context(patch.object(dust, "check_exe"))
            stack.enter_context(patch.object(dust, "get_seq_region_length", return_value={"seq1": 5000}))
            stack.enter_context(patch.object(dust, "get_slice_id", return_value=[["seq1", "1", "5000"]]))
            mock_pool = stack.enter_context(patch.object(dust.multiprocessing, "Pool", return_value=pool))
            stack.enter_context(patch.object(dust, "slice_output_to_gtf", side_effect=write_annotation))

            dust.run_dust(genome_file, output_dir, dust_bin, num_threads)

        expected_gtf_file = output_dir / "dust_output" / "annotation.gtf"
        self.assertTrue(expected_gtf_file.exists())
        mock_check_exe.assert_called_once_with(dust_bin)
        mock_pool.assert_called_once_with(num_threads)
        self.assertEqual(pool.tasks[0][0].__name__, "_multiprocess_dust")
        self.assertTrue(pool.closed)
        self.assertTrue(pool.joined)

if __name__ == "__main__":
    unittest.main()
