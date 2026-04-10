"""
Unit and integration tests for the assembly module (Flye wrapper).
"""

from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import tempfile
import shutil
import pytest

from nanopore_pipeline.assembly.flye_wrapper import FlyeAssembler, AssemblyResult
from nanopore_pipeline.utils.fasta_parser import SequenceRecord


class TestFlyeAssembler:
    """Unit tests for FlyeAssembler class."""

    def test_init_default_params(self):
        """Test assembler initialization with default parameters."""
        assembler = FlyeAssembler()
        assert assembler.genome_size == "5m"
        assert assembler.min_overlap == 1000
        assert assembler.iterations == 1
        assert assembler.threads >= 1
        assert assembler.meta_mode is True
        assert assembler.min_read_length == 1000
        assert assembler.timeout == 3600

    def test_init_custom_params(self):
        """Test assembler initialization with custom parameters."""
        assembler = FlyeAssembler(
            genome_size="3m",
            min_overlap=500,
            iterations=2,
            threads=8,
            meta_mode=False,
            min_read_length=2000,
            timeout=7200,
        )
        assert assembler.genome_size == "3m"
        assert assembler.min_overlap == 500
        assert assembler.iterations == 2
        assert assembler.threads == 8
        assert assembler.meta_mode is False
        assert assembler.min_read_length == 2000
        assert assembler.timeout == 7200

    def test_check_flye_installed(self):
        """Test Flye installation check."""
        assembler = FlyeAssembler()

        # Mock shutil.which to simulate Flye not installed
        with patch("shutil.which", return_value=None):
            assert assembler._check_flye_installed() is False

        # Mock shutil.which to simulate Flye installed
        with patch("shutil.which", return_value="/usr/bin/flye"):
            assert assembler._check_flye_installed() is True

    def test_filter_short_reads(self, tmp_path):
        """Test filtering of short reads."""
        assembler = FlyeAssembler(min_read_length=1000)

        # Create test input FASTA with mixed lengths
        input_fasta = tmp_path / "input.fasta"
        with open(input_fasta, "w", encoding="utf-8") as f:
            f.write(">read1 short\n")
            f.write("A" * 500 + "\n")  # 500bp - should be filtered
            f.write(">read2 long\n")
            f.write("G" * 1500 + "\n")  # 1500bp - should pass
            f.write(">read3 exact\n")
            f.write("C" * 1000 + "\n")  # 1000bp - should pass (equal to threshold)
            f.write(">read4 short\n")
            f.write("T" * 100 + "\n")  # 100bp - should be filtered

        output_fasta = tmp_path / "output.fasta"

        # Run filter
        retained = assembler._filter_short_reads(input_fasta, output_fasta)

        # Check retained count
        assert retained == 2

        # Check output content
        with open(output_fasta, "r", encoding="utf-8") as f:
            content = f.read()
            assert "read2" in content
            assert "read3" in content
            assert "read1" not in content
            assert "read4" not in content

    def test_assemble_missing_input(self, tmp_path):
        """Test assembly with missing input file."""
        assembler = FlyeAssembler()
        missing_file = tmp_path / "nonexistent.fasta"
        output_dir = tmp_path / "output"

        result = assembler.assemble(
            input_fasta=missing_file,
            output_dir=output_dir,
            sample_name="test",
        )

        assert result.success is False
        assert "not found" in result.error_message.lower()

    def test_assemble_flye_not_installed(self, tmp_path):
        """Test assembly when Flye is not installed."""
        assembler = FlyeAssembler()

        # Create dummy input file
        input_fasta = tmp_path / "input.fasta"
        with open(input_fasta, "w", encoding="utf-8") as f:
            f.write(">read1\n")
            f.write("A" * 2000 + "\n")

        output_dir = tmp_path / "output"

        # Mock Flye not installed
        with patch.object(assembler, "_check_flye_installed", return_value=False):
            result = assembler.assemble(
                input_fasta=input_fasta,
                output_dir=output_dir,
                sample_name="test",
            )

        assert result.success is False
        assert "flye not found" in result.error_message.lower()

    def test_assemble_no_reads_after_filtering(self, tmp_path):
        """Test assembly when all reads are filtered out."""
        assembler = FlyeAssembler(min_read_length=1000)

        # Create input with only short reads
        input_fasta = tmp_path / "input.fasta"
        with open(input_fasta, "w", encoding="utf-8") as f:
            f.write(">read1\n")
            f.write("A" * 500 + "\n")
            f.write(">read2\n")
            f.write("G" * 300 + "\n")

        output_dir = tmp_path / "output"

        # Mock Flye installed
        with patch.object(assembler, "_check_flye_installed", return_value=True):
            result = assembler.assemble(
                input_fasta=input_fasta,
                output_dir=output_dir,
                sample_name="test",
            )

        assert result.success is False
        assert "no reads" in result.error_message.lower()

    def test_assemble_timeout(self, tmp_path):
        """Test assembly timeout handling."""
        assembler = FlyeAssembler(timeout=1)  # 1 second timeout

        # Create valid input
        input_fasta = tmp_path / "input.fasta"
        with open(input_fasta, "w", encoding="utf-8") as f:
            f.write(">read1\n")
            f.write("A" * 2000 + "\n")

        output_dir = tmp_path / "output"

        # Mock Flye installed and subprocess to timeout
        with patch.object(assembler, "_check_flye_installed", return_value=True):
            with patch("subprocess.run") as mock_run:
                import subprocess
                mock_run.side_effect = subprocess.TimeoutExpired(cmd="flye", timeout=1)

                result = assembler.assemble(
                    input_fasta=input_fasta,
                    output_dir=output_dir,
                    sample_name="test",
                )

        assert result.success is False
        assert "timed out" in result.error_message.lower()

    def test_assemble_nonzero_exit(self, tmp_path):
        """Test assembly when Flye exits with error."""
        assembler = FlyeAssembler()

        # Create valid input
        input_fasta = tmp_path / "input.fasta"
        with open(input_fasta, "w", encoding="utf-8") as f:
            f.write(">read1\n")
            f.write("A" * 2000 + "\n")

        output_dir = tmp_path / "output"

        # Mock Flye installed and subprocess returning error
        with patch.object(assembler, "_check_flye_installed", return_value=True):
            with patch("subprocess.run") as mock_run:
                mock_result = MagicMock()
                mock_result.returncode = 1
                mock_run.return_value = mock_result

                result = assembler.assemble(
                    input_fasta=input_fasta,
                    output_dir=output_dir,
                    sample_name="test",
                )

        assert result.success is False
        assert "exited with code 1" in result.error_message

    def test_cleanup_intermediate_files(self, tmp_path):
        """Test cleanup of intermediate assembly files."""
        assembler = FlyeAssembler()

        # Create mock assembly directory structure
        assembly_dir = tmp_path / "assembly"
        assembly_dir.mkdir()

        # Create intermediate directories
        intermediate_dirs = [
            "00-assembly",
            "10-consensus",
            "20-repeat",
            "30-contigger",
            "40-polishing",
        ]

        for dirname in intermediate_dirs:
            dirpath = assembly_dir / dirname
            dirpath.mkdir()
            # Add a dummy file
            (dirpath / "dummy.txt").write_text("test")

        # Create files to keep
        (assembly_dir / "assembly.fasta").write_text(">contig1\nATGC\n")
        (assembly_dir / "flye_run.log").write_text("log content")

        # Run cleanup
        assembler.cleanup_intermediate_files(assembly_dir)

        # Check intermediate dirs removed
        for dirname in intermediate_dirs:
            assert not (assembly_dir / dirname).exists()

        # Check kept files still exist
        assert (assembly_dir / "assembly.fasta").exists()
        assert (assembly_dir / "flye_run.log").exists()


@pytest.mark.integration
class TestFlyeAssemblerIntegration:
    """Integration tests requiring actual Flye installation."""

    @pytest.mark.skipif(
        shutil.which("flye") is None,
        reason="Flye not installed"
    )
    def test_assemble_small_dataset(self, tmp_path):
        """
        Integration test with real Flye on a tiny dataset.
        Note: This test requires Flye to be installed.
        """
        assembler = FlyeAssembler()

        # Create a small test dataset (3 reads, minimal overlap)
        input_fasta = tmp_path / "test_reads.fasta"
        with open(input_fasta, "w", encoding="utf-8") as f:
            # Read 1: 2000bp
            f.write(">read1\n")
            f.write("ATGC" * 500 + "\n")
            # Read 2: 2500bp (overlaps with read1)
            f.write(">read2\n")
            f.write("GCTA" * 625 + "\n")
            # Read 3: 3000bp
            f.write(">read3\n")
            f.write("CGAT" * 750 + "\n")

        output_dir = tmp_path / "assembly_output"

        result = assembler.assemble(
            input_fasta=input_fasta,
            output_dir=output_dir,
            sample_name="test_integration",
        )

        # Note: Flye may fail on such a tiny dataset, which is expected
        # This test mainly checks that the wrapper handles it gracefully
        if result.success:
            # If it succeeds (unlikely), validate output
            assert result.contig_path.exists()
            assert result.contigs > 0
            assert result.n50 > 0
        else:
            # Failure is acceptable for tiny dataset
            assert result.error_message != ""
            assert result.assembly_dir == output_dir
