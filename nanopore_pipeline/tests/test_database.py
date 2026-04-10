"""
Tests for database models and DatabaseManager.
"""

import tempfile
from pathlib import Path

import pytest

from nanopore_pipeline.models.database import (
    Base, Sample, VirulenceFactorEntry, AMREntry,
    init_db, get_session_factory,
)
from nanopore_pipeline.db.manager import DatabaseManager


@pytest.fixture
def db_url(tmp_path: Path) -> str:
    return f"sqlite:///{tmp_path / 'test.db'}"


@pytest.fixture
def db_manager(db_url: str) -> DatabaseManager:
    return DatabaseManager(db_url)


@pytest.fixture
def session_factory(db_url: str):
    init_db(db_url)
    return get_session_factory(db_url)


class TestDatabaseInit:
    def test_tables_created(self, db_url: str):
        engine = init_db(db_url)
        table_names = Base.metadata.tables.keys()
        assert "samples" in table_names
        assert "virulence_factors" in table_names
        assert "amr_entries" in table_names
        assert "alignment_hits" in table_names
        assert "classification_results" in table_names


class TestDatabaseManager:
    def test_register_sample(self, db_manager: DatabaseManager):
        sample = db_manager.register_sample(
            name="test_sample",
            file_path="/tmp/test.fastq",
            description="Test sample",
            total_reads=1000,
        )
        assert sample.id is not None
        assert sample.sample_name == "test_sample"

    def test_duplicate_sample_returns_existing(self, db_manager: DatabaseManager):
        s1 = db_manager.register_sample(name="dup", file_path="/tmp/a.fastq")
        s2 = db_manager.register_sample(name="dup", file_path="/tmp/a.fastq")
        assert s1.id == s2.id

    def test_get_sample(self, db_manager: DatabaseManager):
        db_manager.register_sample(name="lookup_test", file_path="/tmp/x.fq")
        found = db_manager.get_sample("lookup_test")
        assert found is not None
        assert found.sample_name == "lookup_test"

    def test_get_sample_missing(self, db_manager: DatabaseManager):
        assert db_manager.get_sample("nonexistent") is None

    def test_list_samples(self, db_manager: DatabaseManager):
        db_manager.register_sample(name="s1", file_path="/tmp/1.fq")
        db_manager.register_sample(name="s2", file_path="/tmp/2.fq")
        samples = db_manager.list_samples()
        assert len(samples) == 2


class TestVFDBLoading:
    def test_load_vfdb_from_file(self, db_manager: DatabaseManager, tmp_path: Path):
        fasta = tmp_path / "vfdb_test.nucleotide_fasta_protein_homolog_model.fasta"
        fasta.write_text(
            ">VFG000001 (fimA) fimbrial adhesin [Escherichia coli]\n"
            "ATCGATCGATCGATCGATCG\n"
            ">VFG000002 (hlyA) alpha-hemolysin exotoxin [E. coli]\n"
            "GCTAGCTAGCTAGCTAGCTA\n"
        )
        count = db_manager.load_vfdb(fasta)
        assert count == 2

    def test_categorisation_adherence(self, db_manager: DatabaseManager, tmp_path: Path):
        fasta = tmp_path / "vfdb_adh.nucleotide_fasta_protein_homolog_model.fasta"
        fasta.write_text(
            ">VFG_ADH (fimA) fimbrial adhesin [E. coli]\n"
            "ATCGATCG\n"
        )
        db_manager.load_vfdb(fasta)
        session = db_manager._session()
        entry = session.query(VirulenceFactorEntry).filter_by(vfdb_id="VFG_ADH").first()
        assert entry.category == "Adherence"
        session.close()

    def test_categorisation_exotoxin(self, db_manager: DatabaseManager, tmp_path: Path):
        fasta = tmp_path / "vfdb_tox.nucleotide_fasta_protein_homolog_model.fasta"
        fasta.write_text(
            ">VFG_TOX (stx) shiga toxin [E. coli]\n"
            "GCTAGCTA\n"
        )
        db_manager.load_vfdb(fasta)
        session = db_manager._session()
        entry = session.query(VirulenceFactorEntry).filter_by(vfdb_id="VFG_TOX").first()
        assert entry.category == "Exotoxin"
        session.close()

    def test_no_duplicates_on_reload(self, db_manager: DatabaseManager, tmp_path: Path):
        fasta = tmp_path / "vfdb_dup.nucleotide_fasta_protein_homolog_model.fasta"
        fasta.write_text(">VFG_DUP (gene) product [Org]\nATCG\n")
        db_manager.load_vfdb(fasta)
        db_manager.load_vfdb(fasta)
        session = db_manager._session()
        count = session.query(VirulenceFactorEntry).filter_by(vfdb_id="VFG_DUP").count()
        assert count == 1
        session.close()
