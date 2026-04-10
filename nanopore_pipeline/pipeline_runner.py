"""
End-to-end Pipeline Runner
~~~~~~~~~~~~~~~~~~~~~~~~~~
Orchestrates the full flow: raw FASTA in -> report out.

Steps:
  1. Compute read statistics from FASTA
  2. Register sample in database
  3. Align reads against VFDB (virulence factors)
  4. Align reads against CARD (antibiotic resistance)
  5. Classify hits into pathogenicity categories
  6. Generate visualisation reports + JSON export
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from config import settings
from nanopore_pipeline.db.manager import DatabaseManager
from nanopore_pipeline.alignment.wrapper import AlignmentWrapper
from nanopore_pipeline.classifier.pathogenicity import PathogenicityClassifier, CategoryProfile
from nanopore_pipeline.reporter.visualisation import PathogenicityReporter
from nanopore_pipeline.utils.fasta_parser import compute_fasta_stats, FastaStats
from nanopore_pipeline.utils.logging_config import setup_logging

logger = setup_logging("pipeline_runner")


@dataclass
class PipelineResult:
    sample_name: str
    sample_id: int
    fasta_stats: FastaStats
    vf_hits: int
    amr_hits: int
    profiles: list[CategoryProfile]
    report_path: Optional[Path]
    json_path: Optional[Path]


class PipelineRunner:

    def __init__(self, db_url: str = settings.DATABASE_URL):
        self.db = DatabaseManager(db_url)
        self.aligner = AlignmentWrapper(db_url)
        self.classifier = PathogenicityClassifier(db_url)
        self.reporter = PathogenicityReporter()
        logger.info("PipelineRunner initialised")

    def run(
        self,
        sample_name: str,
        fasta_path: Path,
        description: str = "",
        vfdb_blast_db: Optional[Path] = None,
        card_blast_db: Optional[Path] = None,
        alignment_tool: str = "blastn",
        skip_vf: bool = False,
        skip_amr: bool = False,
        word_size: int = 7,
    ) -> PipelineResult:
        fasta_path = Path(fasta_path)
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        # ── Step 1: Compute FASTA read statistics ─────────────────────────
        logger.info("Step 1/6: Computing read statistics for %s", fasta_path)
        stats = compute_fasta_stats(fasta_path)
        logger.info(
            "  Reads=%d  Bases=%d  MeanLen=%.0f  N50=%d  GC=%.1f%%",
            stats.total_reads, stats.total_bases,
            stats.mean_read_length, stats.n50, stats.mean_gc_content,
        )

        # ── Step 2: Register sample ──────────────────────────────────────
        logger.info("Step 2/6: Registering sample '%s'", sample_name)
        sample = self.db.register_sample(
            name=sample_name,
            file_path=str(fasta_path),
            description=description,
            total_reads=stats.total_reads,
            total_bases=stats.total_bases,
            mean_read_length=stats.mean_read_length,
            n50=stats.n50,
            mean_gc_content=stats.mean_gc_content,
        )

        # ── Step 3: Align against VFDB ───────────────────────────────────
        vf_hits = 0
        if not skip_vf:
            logger.info("Step 3/6: Aligning against VFDB")
            vf_db = vfdb_blast_db or (settings.VFDB_DIR / "VFDB_setA_nt.fas.blastdb")
            vf_hits, vf_hit_rows = self.aligner.align_sample(
                sample_path=fasta_path,
                blast_db_path=vf_db,
                sample_id=sample.id,
                hit_type="VF",
                tool=alignment_tool,
                word_size=word_size,
            )
            logger.info("  VF hits stored: %d", vf_hits)
            if vf_hit_rows:
                self.reporter.plot_alignment_quality(
                    hits_data=[{"pident": h.pident, "bitscore": h.bitscore} for h in vf_hit_rows],
                    sample_name=sample_name,
                    hit_type="VF",
                )
        else:
            logger.info("Step 3/6: Skipping VF alignment")

        # ── Step 4: Align against CARD ───────────────────────────────────
        amr_hits = 0
        if not skip_amr:
            logger.info("Step 4/6: Aligning against CARD")
            amr_db = card_blast_db or (
                settings.CARD_DIR / "nucleotide_fasta_protein_homolog_model.nucleotide_fasta_protein_homolog_model.fasta.blastdb"
            )
            amr_hits, amr_hit_rows = self.aligner.align_sample(
                sample_path=fasta_path,
                blast_db_path=amr_db,
                sample_id=sample.id,
                hit_type="AMR",
                tool=alignment_tool,
                word_size=word_size,
            )
            logger.info("  AMR hits stored: %d", amr_hits)
            if amr_hit_rows:
                self.reporter.plot_alignment_quality(
                    hits_data=[{"pident": h.pident, "bitscore": h.bitscore} for h in amr_hit_rows],
                    sample_name=sample_name,
                    hit_type="AMR",
                )
        else:
            logger.info("Step 4/6: Skipping AMR alignment")

        # ── Step 5: Classify ─────────────────────────────────────────────
        logger.info("Step 5/6: Classifying pathogenicity profile")
        profiles = self.classifier.classify_sample(sample.id)
        logger.info("  Categories found: %d", len(profiles))
        for p in profiles:
            logger.info(
                "    %-25s  genes=%d  hits=%d  TPM=%.1f  risk=%s",
                p.category, p.unique_genes, p.hit_count, p.tpm, p.risk_level,
            )

        # ── Step 6: Generate reports ─────────────────────────────────────
        logger.info("Step 6/6: Generating reports")
        report_path = None
        json_path = None
        if profiles:
            self.reporter.plot_category_bar(profiles, sample_name)
            self.reporter.plot_risk_dashboard(profiles, sample_name)
            json_path = self.reporter.export_json(profiles, sample_name)
            report_path = settings.RESULTS_DIR / f"{sample_name}_risk_dashboard.html"
            logger.info("  Reports saved to %s", settings.RESULTS_DIR)
        else:
            logger.warning("  No profiles to report (no alignment hits?)")

        logger.info("Pipeline complete for '%s'", sample_name)

        return PipelineResult(
            sample_name=sample_name,
            sample_id=sample.id,
            fasta_stats=stats,
            vf_hits=vf_hits,
            amr_hits=amr_hits,
            profiles=profiles,
            report_path=report_path,
            json_path=json_path,
        )

    def run_comparison(
        self,
        sample_name_a: str,
        sample_name_b: str,
    ) -> dict:
        """
        Compare two already-processed samples and generate comparison reports.

        Both samples must have been processed with `run()` first.
        """
        sa = self.db.get_sample(sample_name_a)
        sb = self.db.get_sample(sample_name_b)
        if not sa:
            raise ValueError(f"Sample '{sample_name_a}' not found in database")
        if not sb:
            raise ValueError(f"Sample '{sample_name_b}' not found in database")

        comparison = self.classifier.compare_samples(sa.id, sb.id)

        self.reporter.plot_comparison_heatmap(comparison, sample_name_a, sample_name_b)
        self.reporter.plot_fold_change(comparison, sample_name_a, sample_name_b)

        logger.info(
            "Comparison complete: %s vs %s (%d categories)",
            sample_name_a, sample_name_b, len(comparison),
        )
        return comparison
