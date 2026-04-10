"""
Tests for the Pathogenicity Classifier.
"""

import pytest
from pathlib import Path

from nanopore_pipeline.models.database import (
    init_db, get_session_factory, Base,
    Sample, VirulenceFactorEntry, AMREntry, AlignmentHit,
)
from nanopore_pipeline.classifier.pathogenicity import PathogenicityClassifier


@pytest.fixture
def db_url(tmp_path: Path) -> str:
    return f"sqlite:///{tmp_path / 'test_classifier.db'}"


@pytest.fixture
def populated_db(db_url: str):
    """Create a DB with sample, VF entries, AMR entries, and alignment hits."""
    init_db(db_url)
    SessionFactory = get_session_factory(db_url)
    session = SessionFactory()

    # Sample
    sample = Sample(sample_name="test_patient", file_path="/tmp/reads.fastq",
                    total_reads=10000)
    session.add(sample)
    session.flush()

    # VF entries
    vf_adherence = VirulenceFactorEntry(
        vfdb_id="VF001", gene_name="fimA", product="fimbrial adhesin",
        organism="E. coli", category="Adherence", subcategory="adhesin",
        sequence="ATCG", seq_length=4,
    )
    vf_toxin = VirulenceFactorEntry(
        vfdb_id="VF002", gene_name="stx1", product="shiga toxin",
        organism="E. coli", category="Exotoxin", subcategory="toxin",
        sequence="GCTA", seq_length=4,
    )
    vf_secretion = VirulenceFactorEntry(
        vfdb_id="VF003", gene_name="escN", product="type iii secretion",
        organism="E. coli", category="Secretion System", subcategory="t3ss",
        sequence="TTTT", seq_length=4,
    )
    session.add_all([vf_adherence, vf_toxin, vf_secretion])
    session.flush()

    # AMR entry
    amr = AMREntry(
        card_id="ARO001", gene_name="blaTEM-1", aro_accession="ARO:3000015",
        drug_class="Beta-Lactam Resistance", resistance_mechanism="hydrolysis",
        organism="E. coli", sequence="AAAA", seq_length=4,
    )
    session.add(amr)
    session.flush()

    # Alignment hits
    hits = [
        AlignmentHit(sample_id=sample.id, query_id="read_1", hit_type="VF",
                     reference_id=vf_adherence.id, subject_id="VF001",
                     percent_identity=95.0, alignment_length=100,
                     evalue=1e-50, bitscore=200, query_length=150, subject_length=120,
                     query_coverage=66.7),
        AlignmentHit(sample_id=sample.id, query_id="read_2", hit_type="VF",
                     reference_id=vf_toxin.id, subject_id="VF002",
                     percent_identity=92.0, alignment_length=80,
                     evalue=1e-40, bitscore=160, query_length=120, subject_length=100,
                     query_coverage=66.7),
        AlignmentHit(sample_id=sample.id, query_id="read_3", hit_type="VF",
                     reference_id=vf_secretion.id, subject_id="VF003",
                     percent_identity=88.0, alignment_length=60,
                     evalue=1e-30, bitscore=120, query_length=100, subject_length=80,
                     query_coverage=60.0),
        AlignmentHit(sample_id=sample.id, query_id="read_4", hit_type="AMR",
                     amr_reference_id=amr.id, subject_id="ARO001",
                     percent_identity=99.0, alignment_length=200,
                     evalue=1e-80, bitscore=400, query_length=250, subject_length=220,
                     query_coverage=80.0),
    ]
    session.add_all(hits)
    session.commit()
    session.close()

    return db_url, sample.id


class TestPathogenicityClassifier:
    def test_classify_returns_profiles(self, populated_db):
        db_url, sample_id = populated_db
        classifier = PathogenicityClassifier(db_url)
        profiles = classifier.classify_sample(sample_id)
        assert len(profiles) > 0

    def test_categories_present(self, populated_db):
        db_url, sample_id = populated_db
        classifier = PathogenicityClassifier(db_url)
        profiles = classifier.classify_sample(sample_id)
        cats = {p.category for p in profiles}
        assert "Adherence" in cats
        assert "Exotoxin" in cats
        assert "Secretion System" in cats
        assert "Beta-Lactam Resistance" in cats

    def test_risk_levels_assigned(self, populated_db):
        db_url, sample_id = populated_db
        classifier = PathogenicityClassifier(db_url)
        profiles = classifier.classify_sample(sample_id)
        for p in profiles:
            assert p.risk_level in ("Critical", "High", "Medium", "Low")

    def test_critical_risk_when_multiple_critical_cats(self, populated_db):
        db_url, sample_id = populated_db
        classifier = PathogenicityClassifier(db_url)
        profiles = classifier.classify_sample(sample_id)
        # Exotoxin, Secretion System, and Beta-Lactam are all CRITICAL_CATEGORIES
        critical_profiles = [p for p in profiles if p.risk_level == "Critical"]
        assert len(critical_profiles) >= 2

    def test_tpm_computed(self, populated_db):
        db_url, sample_id = populated_db
        classifier = PathogenicityClassifier(db_url)
        profiles = classifier.classify_sample(sample_id)
        for p in profiles:
            assert p.tpm >= 0

    def test_classification_stored_in_db(self, populated_db):
        db_url, sample_id = populated_db
        classifier = PathogenicityClassifier(db_url)
        classifier.classify_sample(sample_id)

        from nanopore_pipeline.models.database import ClassificationResult
        SessionFactory = get_session_factory(db_url)
        session = SessionFactory()
        results = session.query(ClassificationResult).filter_by(
            sample_id=sample_id
        ).all()
        assert len(results) > 0
        session.close()


class TestComparesamples:
    def test_compare_returns_dict(self, populated_db):
        db_url, sample_id = populated_db
        classifier = PathogenicityClassifier(db_url)
        # Compare sample with itself (trivial but tests the code path)
        comparison = classifier.compare_samples(sample_id, sample_id)
        assert isinstance(comparison, dict)
        assert len(comparison) > 0

    def test_fold_change_self_is_one(self, populated_db):
        db_url, sample_id = populated_db
        classifier = PathogenicityClassifier(db_url)
        comparison = classifier.compare_samples(sample_id, sample_id)
        for cat, data in comparison.items():
            assert data["fold_change_a_vs_b"] == pytest.approx(1.0)


class TestTPMComputation:
    def test_tpm_with_equal_lengths(self):
        classifier = PathogenicityClassifier.__new__(PathogenicityClassifier)
        result = classifier.compute_tpm([(100, 1000), (100, 1000)])
        assert result["tpm_sum"] > 0

    def test_tpm_empty_list(self):
        classifier = PathogenicityClassifier.__new__(PathogenicityClassifier)
        result = classifier.compute_tpm([])
        assert result["tpm_sum"] == 0
        assert result["tpm_mean"] == 0

    def test_tpm_zero_gene_length(self):
        classifier = PathogenicityClassifier.__new__(PathogenicityClassifier)
        result = classifier.compute_tpm([(100, 0)])
        assert result["tpm_sum"] == 0
