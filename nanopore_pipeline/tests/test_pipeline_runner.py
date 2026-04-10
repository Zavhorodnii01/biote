"""
Tests for the end-to-end PipelineRunner (using skip flags to avoid needing BLAST).
"""

from pathlib import Path

import pytest

from nanopore_pipeline.models.database import (
    init_db, get_session_factory,
    Sample, VirulenceFactorEntry, AMREntry, AlignmentHit,
)
from nanopore_pipeline.pipeline_runner import PipelineRunner


@pytest.fixture
def db_url(tmp_path: Path) -> str:
    return f"sqlite:///{tmp_path / 'test_runner.db'}"


@pytest.fixture
def sample_fasta(tmp_path: Path) -> Path:
    """A small FASTA file simulating Nanopore reads."""
    content = (
        ">read_001 length=500\n"
        + "ATCGATCG" * 62 + "ATCG\n"   # 500 bp
        + ">read_002 length=300\n"
        + "GCTAGCTA" * 37 + "GCTA\n"   # 300 bp
        + ">read_003 length=1000\n"
        + "ATGCATGC" * 125 + "\n"      # 1000 bp
    )
    fasta_file = tmp_path / "sample_reads.nucleotide_fasta_protein_homolog_model.fasta"
    fasta_file.write_text(content)
    return fasta_file


@pytest.fixture
def seeded_runner(db_url: str) -> PipelineRunner:
    """Runner with pre-seeded VF/AMR entries and hits (simulating BLAST output)."""
    init_db(db_url)
    SessionFactory = get_session_factory(db_url)
    session = SessionFactory()

    # Pre-seed VF entries
    vf = VirulenceFactorEntry(
        vfdb_id="VF_TEST", gene_name="testGene", product="test adhesin",
        organism="E. coli", category="Adherence", subcategory="adhesin",
        sequence="ATCG", seq_length=4,
    )
    session.add(vf)
    session.flush()

    # Pre-seed a sample
    sample = Sample(
        sample_name="pre_seeded", file_path="/tmp/fake.nucleotide_fasta_protein_homolog_model.fasta",
        total_reads=100, total_bases=50000,
    )
    session.add(sample)
    session.flush()

    # Pre-seed hits for that sample
    hit = AlignmentHit(
        sample_id=sample.id, query_id="read_1", hit_type="VF",
        reference_id=vf.id, subject_id="VF_TEST",
        percent_identity=95.0, alignment_length=100,
        evalue=1e-50, bitscore=200, query_length=150, subject_length=120,
        query_coverage=66.7,
    )
    session.add(hit)
    session.commit()
    session.close()

    return PipelineRunner(db_url)


class TestPipelineRunnerStats:
    def test_run_computes_stats(self, db_url: str, sample_fasta: Path):
        """Pipeline computes FASTA stats even when alignment is skipped."""
        runner = PipelineRunner(db_url)
        result = runner.run(
            sample_name="stats_test",
            fasta_path=sample_fasta,
            skip_vf=True,
            skip_amr=True,
        )
        assert result.fasta_stats.total_reads == 3
        assert result.fasta_stats.total_bases > 0
        assert result.fasta_stats.n50 > 0
        assert result.fasta_stats.mean_gc_content > 0

    def test_run_registers_sample(self, db_url: str, sample_fasta: Path):
        runner = PipelineRunner(db_url)
        result = runner.run(
            sample_name="register_test",
            fasta_path=sample_fasta,
            skip_vf=True,
            skip_amr=True,
        )
        assert result.sample_id > 0
        sample = runner.db.get_sample("register_test")
        assert sample is not None
        assert sample.total_reads == 3

    def test_run_file_not_found(self, db_url: str):
        runner = PipelineRunner(db_url)
        with pytest.raises(FileNotFoundError):
            runner.run("bad", Path("/nonexistent/file.nucleotide_fasta_protein_homolog_model.fasta"), skip_vf=True, skip_amr=True)


class TestPipelineRunnerClassification:
    def test_classify_seeded_data(self, seeded_runner: PipelineRunner):
        """Classifier works on pre-seeded alignment hits."""
        profiles = seeded_runner.classifier.classify_sample(1)
        assert len(profiles) > 0
        cats = {p.category for p in profiles}
        assert "Adherence" in cats

    def test_comparison_seeded(self, seeded_runner: PipelineRunner):
        """Compare pre-seeded sample with itself."""
        comparison = seeded_runner.run_comparison("pre_seeded", "pre_seeded")
        assert isinstance(comparison, dict)
        assert len(comparison) > 0
