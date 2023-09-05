import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch
from src.python.ensembl.tools.anno.repeat_annotation.dust import run_dust  # Replace 'your_module' with the actual module path

class TestDust(unittest.TestCase):

    def setUp(self):
        # Create a temporary directory for testing
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        # Remove the temporary directory and its contents
        shutil.rmtree(self.temp_dir)

    def test_run_dust(self):
        # Define test input values
        genome_file = "path/to/your/genome.fasta"
        output_dir = self.temp_dir
        dust_bin = "dust"
        num_threads = 2  # Set the number of threads for testing

        # Mock any external dependencies or subprocess calls
        # For example, if the `dustmasker` command is used, you can mock it:
        with patch("subprocess.run") as mock_run:
            mock_run.return_value.returncode = 0  # Mock a successful subprocess call
            run_dust(genome_file, output_dir, dust_bin, num_threads)

        # Perform assertions to check if the function behaves as expected
        # You can assert conditions based on the function's behavior
        # For example, check if the expected GTF file was created:
        expected_gtf_file = output_dir / "dust_output" / "annotation.gtf"
        self.assertTrue(expected_gtf_file.exists())

        # Add more assertions as needed based on your function's behavior

if __name__ == "__main__":
    unittest.main()
