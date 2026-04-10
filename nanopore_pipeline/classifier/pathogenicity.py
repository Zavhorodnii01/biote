"""
Module 3 -- Pathogenicity Classifier
Groups alignment hits into functional categories (Adherence, Motility,
Exotoxins, Endotoxins, Secretion Systems, AMR classes) and computes
TPM-normalised abundance.

Risk scoring:
    Critical  -  Exotoxin + Secretion System + AMR present
    High      -  Two of the above
    Medium    -  One major category, moderate hit counts
    Low       -  Minimal or non-specific hits
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional

from sqlalchemy.orm import Session

from config import settings
from nanopore_pipeline.models.database import (
    AlignmentHit, ClassificationResult, Sample,
    VirulenceFactorEntry, AMREntry,
    get_session_factory,
)
from nanopore_pipeline.utils.logging_config import setup_logging

logger = setup_logging("classifier.pathogenicity")


@dataclass
class CategoryProfile:
    category: str
    gene_names: list[str] = field(default_factory=list)
    hit_count: int = 0
    unique_genes: int = 0
    mean_identity: float = 0.0
    mean_coverage: float = 0.0
    tpm: float = 0.0
    risk_level: str = "Low"


class PathogenicityClassifier:
    """Classifies alignment hits into pathogenicity categories with TPM normalisation."""

    CRITICAL_CATEGORIES = {"Exotoxin", "Secretion System", "Beta-Lactam Resistance"}
    HIGH_CATEGORIES = {
        "Adherence", "Invasion", "Endotoxin", "Iron Uptake",
        "Aminoglycoside Resistance", "Fluoroquinolone Resistance",
    }

    def __init__(self, db_url: str = settings.DATABASE_URL):
        self._SessionFactory = get_session_factory(db_url)
        logger.info("PathogenicityClassifier initialised")

    def _session(self) -> Session:
        return self._SessionFactory()

    # ── TPM normalisation ─────────────────────────────────────────────────

    @staticmethod
    def compute_tpm(hit_lengths: list[tuple[int, int]]) -> dict[str, float]:
        """
        Compute TPM given list of (gene_index, alignment_length / gene_length) pairs
        grouped by category.

        hit_lengths: list of (alignment_length, gene_length) per hit
        Returns per-hit RPK and overall TPM scalar.
        """
        rpk_values = []
        for aln_len, gene_len in hit_lengths:
            if gene_len > 0:
                rpk = aln_len / (gene_len / 1000)
            else:
                rpk = 0
            rpk_values.append(rpk)

        total_rpk = sum(rpk_values)
        scaling = total_rpk / settings.TPM_SCALING_FACTOR if total_rpk > 0 else 1
        tpm_values = [rpk / scaling if scaling > 0 else 0 for rpk in rpk_values]
        return {
            "tpm_sum": sum(tpm_values),
            "tpm_mean": sum(tpm_values) / len(tpm_values) if tpm_values else 0,
        }

    # ── Classify a single sample ──────────────────────────────────────────

    def classify_sample(self, sample_id: int) -> list[CategoryProfile]:
        session = self._session()
        try:
            hits = session.query(AlignmentHit).filter_by(sample_id=sample_id).all()
            if not hits:
                logger.warning("No alignment hits for sample_id=%d", sample_id)
                return []

            # Group VF hits by category
            vf_groups: dict[str, list[AlignmentHit]] = defaultdict(list)
            amr_groups: dict[str, list[AlignmentHit]] = defaultdict(list)

            for hit in hits:
                if hit.hit_type == "VF" and hit.vf_entry:
                    cat = hit.vf_entry.category or "Other"
                    vf_groups[cat].append(hit)
                elif hit.hit_type == "AMR" and hit.amr_entry:
                    cat = hit.amr_entry.drug_class or "Other"
                    amr_groups[cat].append(hit)

            profiles: list[CategoryProfile] = []

            # Process VF categories
            for cat, cat_hits in vf_groups.items():
                profiles.append(self._build_profile(cat, cat_hits, "VF"))

            # Process AMR categories
            for cat, cat_hits in amr_groups.items():
                profiles.append(self._build_profile(cat, cat_hits, "AMR"))

            # Assign risk levels based on the overall profile
            self._assign_risk_levels(profiles)

            # Persist to DB
            self._store_classifications(session, sample_id, profiles)

            logger.info("Classified sample_id=%d into %d categories", sample_id, len(profiles))
            return profiles

        except Exception:
            session.rollback()
            logger.exception("Error classifying sample %d", sample_id)
            raise
        finally:
            session.close()

    def _build_profile(
        self,
        category: str,
        hits: list[AlignmentHit],
        hit_type: str,
    ) -> CategoryProfile:
        gene_names = list({
            (h.vf_entry.gene_name if hit_type == "VF" and h.vf_entry else
             h.amr_entry.gene_name if hit_type == "AMR" and h.amr_entry else
             h.subject_id)
            for h in hits
        })

        identities = [h.percent_identity for h in hits if h.percent_identity]
        coverages = [h.query_coverage for h in hits if h.query_coverage]

        hit_lengths = [
            (h.alignment_length or 0, h.subject_length or 1) for h in hits
        ]
        tpm_data = self.compute_tpm(hit_lengths)

        return CategoryProfile(
            category=category,
            gene_names=gene_names,
            hit_count=len(hits),
            unique_genes=len(gene_names),
            mean_identity=sum(identities) / len(identities) if identities else 0,
            mean_coverage=sum(coverages) / len(coverages) if coverages else 0,
            tpm=tpm_data["tpm_sum"],
        )

    def _assign_risk_levels(self, profiles: list[CategoryProfile]) -> None:
        category_names = {p.category for p in profiles}
        crit_present = category_names & self.CRITICAL_CATEGORIES
        high_present = category_names & self.HIGH_CATEGORIES

        for p in profiles:
            if p.category in self.CRITICAL_CATEGORIES:
                if len(crit_present) >= 2:
                    p.risk_level = "Critical"
                else:
                    p.risk_level = "High"
            elif p.category in self.HIGH_CATEGORIES:
                if len(crit_present) >= 1:
                    p.risk_level = "High"
                else:
                    p.risk_level = "Medium"
            else:
                if p.hit_count > 10:
                    p.risk_level = "Medium"
                else:
                    p.risk_level = "Low"

    def _store_classifications(
        self,
        session: Session,
        sample_id: int,
        profiles: list[CategoryProfile],
    ) -> None:
        # Clear previous results for this sample
        session.query(ClassificationResult).filter_by(sample_id=sample_id).delete()

        for p in profiles:
            cr = ClassificationResult(
                sample_id=sample_id,
                category=p.category,
                gene_count=p.unique_genes,
                total_hits=p.hit_count,
                mean_identity=p.mean_identity,
                mean_coverage=p.mean_coverage,
                tpm=p.tpm,
                risk_level=p.risk_level,
            )
            session.add(cr)
        session.commit()

    # ── Compare two samples ───────────────────────────────────────────────

    def compare_samples(
        self,
        sample_id_a: int,
        sample_id_b: int,
    ) -> dict:
        profiles_a = self.classify_sample(sample_id_a)
        profiles_b = self.classify_sample(sample_id_b)

        map_a = {p.category: p for p in profiles_a}
        map_b = {p.category: p for p in profiles_b}

        all_cats = sorted(set(map_a.keys()) | set(map_b.keys()))
        comparison = {}
        for cat in all_cats:
            pa = map_a.get(cat)
            pb = map_b.get(cat)
            tpm_a = pa.tpm if pa else 0
            tpm_b = pb.tpm if pb else 0
            fold_change = (tpm_a / tpm_b) if tpm_b > 0 else float("inf") if tpm_a > 0 else 0

            comparison[cat] = {
                "sample_a_hits": pa.hit_count if pa else 0,
                "sample_b_hits": pb.hit_count if pb else 0,
                "sample_a_genes": pa.unique_genes if pa else 0,
                "sample_b_genes": pb.unique_genes if pb else 0,
                "sample_a_tpm": tpm_a,
                "sample_b_tpm": tpm_b,
                "fold_change_a_vs_b": round(fold_change, 2),
                "sample_a_risk": pa.risk_level if pa else "None",
                "sample_b_risk": pb.risk_level if pb else "None",
            }

        logger.info(
            "Comparison complete: %d categories across samples %d vs %d",
            len(comparison), sample_id_a, sample_id_b,
        )
        return comparison
