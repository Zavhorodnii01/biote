"""
Ground Truth Evaluator
======================
Compares pipeline-detected genes against a known GFF annotation (ground truth)
to compute per-category recall and identify WHERE in the pipeline genes were missed.

Used for Test 4 (MRSA RefSeq reference) where we know exactly what genes are
present in the genome.

Failure point taxonomy:
  - "not_in_db"   : gene is annotated but not in VFDB/CARD (database gap)
  - "not_detected": gene is in the DB but pipeline produced zero hits for this category
  - "low_identity": BLAST ran but hits were filtered by identity/coverage threshold
  - "found"       : gene was correctly detected
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from nanopore_pipeline.utils.logging_config import setup_logging

logger = setup_logging("classifier.ground_truth_evaluator")


@dataclass
class CategoryEvaluation:
    category: str
    expected_genes: list[str]       # genes found in GFF annotation
    found_genes: list[str]          # genes detected by pipeline
    missed_genes: list[str]         # expected but not detected
    found_count: int = 0
    total_count: int = 0
    recall: float = 0.0             # found / total (0–1)
    failure_reason: str = ""        # dominant failure mode for this category


@dataclass
class EvaluationReport:
    sample_name: str
    categories: list[CategoryEvaluation] = field(default_factory=list)

    @property
    def overall_recall(self) -> float:
        total = sum(c.total_count for c in self.categories)
        found = sum(c.found_count for c in self.categories)
        return found / total if total > 0 else 0.0

    def summary_table(self) -> str:
        """Return a formatted console table."""
        lines = [
            f"\n{'Category':<30} {'Expected':>8} {'Found':>6} {'Recall':>8}  Failure",
            "-" * 68,
        ]
        for c in sorted(self.categories, key=lambda x: -x.recall):
            recall_pct = f"{c.recall * 100:.0f}%"
            lines.append(
                f"{c.category:<30} {c.total_count:>8} {c.found_count:>6} {recall_pct:>8}  {c.failure_reason}"
            )
        lines.append("-" * 68)
        overall = f"{self.overall_recall * 100:.0f}%"
        lines.append(f"{'OVERALL':<30} {sum(c.total_count for c in self.categories):>8} "
                     f"{sum(c.found_count for c in self.categories):>6} {overall:>8}")
        return "\n".join(lines)


class GroundTruthEvaluator:
    """
    Compares pipeline ClassificationResult profiles against ground truth
    built from a GFF3 annotation file.
    """

    def evaluate(
        self,
        ground_truth: dict[str, list[str]],
        pipeline_profiles: list,          # list[CategoryProfile] from classifier
        sample_name: str,
    ) -> EvaluationReport:
        """
        Args:
            ground_truth:      {category: [gene_names]} from GFFGroundTruth.build()
            pipeline_profiles: CategoryProfile objects from PipelineResult.profiles
            sample_name:       Used for labelling the report

        Returns:
            EvaluationReport with per-category recall and failure reasons
        """
        # Index profiles by category for fast lookup
        detected: dict[str, list[str]] = {
            p.category: p.gene_names for p in pipeline_profiles
        }

        report = EvaluationReport(sample_name=sample_name)

        for category, expected_genes in ground_truth.items():
            found_in_pipeline = detected.get(category, [])

            # Count how many expected genes appear in pipeline output
            # Use substring matching (gene names may have minor differences)
            found = []
            missed = []
            for expected_gene in expected_genes:
                eg = expected_gene.lower()
                hit = any(
                    eg in pg.lower() or pg.lower() in eg
                    for pg in found_in_pipeline
                )
                if hit:
                    found.append(expected_gene)
                else:
                    missed.append(expected_gene)

            found_count = len(found)
            total_count = len(expected_genes)
            recall = found_count / total_count if total_count > 0 else 0.0

            # Determine dominant failure reason
            failure_reason = self._classify_failure(
                category, found_count, total_count, found_in_pipeline
            )

            report.categories.append(CategoryEvaluation(
                category=category,
                expected_genes=expected_genes,
                found_genes=found,
                missed_genes=missed,
                found_count=found_count,
                total_count=total_count,
                recall=recall,
                failure_reason=failure_reason,
            ))

        logger.info(
            "Evaluation complete for '%s': overall recall=%.1f%%",
            sample_name, report.overall_recall * 100,
        )
        return report

    @staticmethod
    def _classify_failure(
        category: str,
        found: int,
        total: int,
        pipeline_genes: list[str],
    ) -> str:
        """Assign a human-readable failure reason for the category."""
        if total == 0:
            return "no_annotation"
        if found == total:
            return "all_found"
        if found == 0 and not pipeline_genes:
            return "not_detected"          # pipeline found nothing in this category
        if found == 0 and pipeline_genes:
            return "gene_name_mismatch"    # pipeline found something but names don't match
        if found < total:
            return "partial_detection"     # found some but not all
        return "unknown"
