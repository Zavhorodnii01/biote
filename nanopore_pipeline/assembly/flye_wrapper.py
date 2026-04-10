"""
Flye Assembler Wrapper
~~~~~~~~~~~~~~~~~~~~~~
Wraps Flye de novo assembler for Nanopore long-read assembly.
Reduces overlapping read redundancy before alignment to prevent hit inflation.

Flye parameters tuned for bacterial metagenomics with high error rate.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from config import settings
from nanopore_pipeline.utils.fasta_parser import compute_fasta_stats, parse_fasta, SequenceRecord
from nanopore_pipeline.utils.logging_config import setup_logging

logger = setup_logging("assembly.flye")


@dataclass
class AssemblyResult:
    """Result container for assembly operation."""
    success: bool
    contig_path: Optional[Path] = None
    contigs: int = 0
    n50: int = 0
    total_bases: int = 0
    mean_contig_length: float = 0.0
    assembly_dir: Optional[Path] = None
    error_message: str = ""
    flye_log: Optional[Path] = None


class FlyeAssembler:
    """
    Wrapper for Flye de novo assembler.

    Flye is optimized for long, noisy reads (Nanopore, PacBio).
    Using --meta mode for metagenomics (multiple strains in clinical samples).
    """

    def __init__(
        self,
        genome_size: str = settings.FLYE_GENOME_SIZE,
        min_overlap: int = settings.FLYE_MIN_OVERLAP,
        iterations: int = settings.FLYE_ITERATIONS,
        threads: int = settings.FLYE_THREADS,
        meta_mode: bool = settings.FLYE_META_MODE,
        min_read_length: int = settings.FLYE_MIN_READ_LENGTH,
        timeout: int = settings.FLYE_TIMEOUT,
    ):
        self.genome_size = genome_size
        self.min_overlap = min_overlap
        self.iterations = iterations
        self.threads = threads
        self.meta_mode = meta_mode
        self.min_read_length = min_read_length
        self.timeout = timeout
        logger.info(
            "FlyeAssembler initialized: genome_size=%s, min_overlap=%d, "
            "iterations=%d, threads=%d, meta=%s",
            genome_size, min_overlap, iterations, threads, meta_mode,
        )

    def _check_flye_installed(self) -> bool:
        """Check if Flye is available in PATH."""
        return shutil.which("flye") is not None

    def _filter_short_reads(
        self,
        input_fasta: Path,
        output_fasta: Path,
    ) -> int:
        """
        Filter out reads shorter than min_read_length.

        Returns:
            Number of reads retained after filtering.
        """
        logger.info("Filtering reads < %d bp from %s", self.min_read_length, input_fasta)
        retained = 0
        with open(output_fasta, "w", encoding="utf-8") as out_fh:
            for record in parse_fasta(input_fasta):
                if record.length >= self.min_read_length:
                    out_fh.write(f">{record.id} {record.description}\n")
                    out_fh.write(f"{record.sequence}\n")
                    retained += 1

        logger.info("Retained %d reads after filtering", retained)
        return retained

    def assemble(
        self,
        input_fasta: Path,
        output_dir: Path,
        sample_name: str,
    ) -> AssemblyResult:
        """
        Run Flye assembly on Nanopore reads.

        Args:
            input_fasta: Path to FASTA file with raw Nanopore reads
            output_dir: Directory for assembly output
            sample_name: Sample identifier (for logging)

        Returns:
            AssemblyResult with assembly metadata or error details
        """
        input_fasta = Path(input_fasta)
        output_dir = Path(output_dir)

        # Check prerequisites
        if not input_fasta.exists():
            error_msg = f"Input FASTA not found: {input_fasta}"
            logger.error(error_msg)
            return AssemblyResult(success=False, error_message=error_msg)

        if not self._check_flye_installed():
            error_msg = (
                "Flye not found in PATH. Install with: conda install -c bioconda flye"
            )
            logger.error(error_msg)
            return AssemblyResult(success=False, error_message=error_msg)

        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        # Filter short reads
        filtered_fasta = output_dir / f"{sample_name}_filtered.fasta"
        retained = self._filter_short_reads(input_fasta, filtered_fasta)

        if retained == 0:
            error_msg = f"No reads >= {self.min_read_length}bp remaining after filtering"
            logger.error(error_msg)
            return AssemblyResult(success=False, error_message=error_msg)

        # Build Flye command
        cmd = [
            "flye",
            "--nano-raw", str(filtered_fasta),
            "--out-dir", str(output_dir),
            "--genome-size", self.genome_size,
            "--min-overlap", str(self.min_overlap),
            "--iterations", str(self.iterations),
            "--threads", str(self.threads),
        ]

        if self.meta_mode:
            cmd.append("--meta")

        logger.info("Running Flye assembly for sample '%s'", sample_name)
        logger.info("Command: %s", " ".join(cmd))

        # Run Flye
        log_file = output_dir / "flye_run.log"
        try:
            with open(log_file, "w", encoding="utf-8") as log_fh:
                result = subprocess.run(
                    cmd,
                    stdout=log_fh,
                    stderr=subprocess.STDOUT,
                    timeout=self.timeout,
                    check=False,
                    text=True,
                )

            if result.returncode != 0:
                error_msg = f"Flye exited with code {result.returncode}. Check log: {log_file}"
                logger.error(error_msg)

                # Read last 20 lines of log for error context
                with open(log_file, "r", encoding="utf-8") as log_fh:
                    lines = log_fh.readlines()
                    error_context = "".join(lines[-20:])
                    logger.error("Flye error context:\n%s", error_context)

                return AssemblyResult(
                    success=False,
                    error_message=error_msg,
                    assembly_dir=output_dir,
                    flye_log=log_file,
                )

        except subprocess.TimeoutExpired:
            error_msg = f"Flye assembly timed out after {self.timeout}s"
            logger.error(error_msg)
            return AssemblyResult(
                success=False,
                error_message=error_msg,
                assembly_dir=output_dir,
                flye_log=log_file,
            )
        except Exception as e:
            error_msg = f"Flye assembly failed: {str(e)}"
            logger.exception(error_msg)
            return AssemblyResult(
                success=False,
                error_message=error_msg,
                assembly_dir=output_dir,
                flye_log=log_file,
            )

        # Check for output contigs
        contig_path = output_dir / "assembly.fasta"
        if not contig_path.exists():
            error_msg = f"Flye completed but no assembly.fasta produced in {output_dir}"
            logger.error(error_msg)
            return AssemblyResult(
                success=False,
                error_message=error_msg,
                assembly_dir=output_dir,
                flye_log=log_file,
            )

        # Compute assembly statistics
        try:
            stats = compute_fasta_stats(contig_path)

            if stats.total_reads == 0:
                error_msg = "Assembly produced empty FASTA file"
                logger.error(error_msg)
                return AssemblyResult(
                    success=False,
                    error_message=error_msg,
                    contig_path=contig_path,
                    assembly_dir=output_dir,
                    flye_log=log_file,
                )

            logger.info(
                "Assembly successful: %d contigs, N50=%d, total_bases=%d",
                stats.total_reads, stats.n50, stats.total_bases,
            )

            return AssemblyResult(
                success=True,
                contig_path=contig_path,
                contigs=stats.total_reads,
                n50=stats.n50,
                total_bases=stats.total_bases,
                mean_contig_length=stats.mean_read_length,
                assembly_dir=output_dir,
                flye_log=log_file,
            )

        except Exception as e:
            error_msg = f"Failed to compute assembly statistics: {str(e)}"
            logger.exception(error_msg)
            return AssemblyResult(
                success=False,
                error_message=error_msg,
                contig_path=contig_path,
                assembly_dir=output_dir,
                flye_log=log_file,
            )

    def cleanup_intermediate_files(self, assembly_dir: Path) -> None:
        """
        Remove large intermediate files to save disk space.

        Keeps:
            - assembly.fasta (final contigs)
            - flye_run.log (for debugging)

        Removes:
            - 00-assembly/ (raw assembly graphs)
            - 10-consensus/ (intermediate consensus)
            - 20-repeat/ (repeat detection)
            - 30-contigger/ (contig building)
            - 40-polishing/ (polishing data)
        """
        assembly_dir = Path(assembly_dir)
        if not assembly_dir.exists():
            logger.warning("Assembly directory not found: %s", assembly_dir)
            return

        # Directories to remove
        cleanup_patterns = [
            "00-assembly",
            "10-consensus",
            "20-repeat",
            "30-contigger",
            "40-polishing",
        ]

        for pattern in cleanup_patterns:
            target = assembly_dir / pattern
            if target.exists() and target.is_dir():
                try:
                    shutil.rmtree(target)
                    logger.info("Cleaned up intermediate directory: %s", target)
                except Exception as e:
                    logger.warning("Failed to remove %s: %s", target, e)

        logger.info("Cleanup complete for %s", assembly_dir)
