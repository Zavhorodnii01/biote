"""
Tests for the FASTA parser, read filtering, and stats computation.
"""

from pathlib import Path

import pytest

from nanopore_pipeline.utils.fasta_parser import (
    SequenceRecord,
    FastaStats,
    parse_fasta,
    filter_reads,
    compute_fasta_stats,
)


@pytest.fixture
def sample_fasta(tmp_path: Path) -> Path:
    content = (
        ">seq1 Test gene [E. coli]\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        ">seq2 Another gene [S. aureus]\n"
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n"
    )
    fasta_file = tmp_path / "test.nucleotide_fasta_protein_homolog_model.fasta"
    fasta_file.write_text(content)
    return fasta_file


@pytest.fixture
def multi_read_fasta(tmp_path: Path) -> Path:
    """FASTA with reads of varying lengths for stats testing."""
    reads = [
        (">read1", "A" * 100),
        (">read2", "GCGC" * 50),   # 200 bp, 100% GC
        (">read3", "ATCG" * 75),   # 300 bp, 50% GC
        (">short", "AT"),           # 2 bp, too short to pass filters
    ]
    content = "\n".join(f"{hdr}\n{seq}" for hdr, seq in reads) + "\n"
    fasta_file = tmp_path / "multi.nucleotide_fasta_protein_homolog_model.fasta"
    fasta_file.write_text(content)
    return fasta_file


class TestSequenceRecord:
    def test_length(self):
        rec = SequenceRecord(id="s1", description="test", sequence="ATCG")
        assert rec.length == 4

    def test_gc_content_50(self):
        rec = SequenceRecord(id="s1", description="test", sequence="ATGC")
        assert rec.gc_content == pytest.approx(50.0)

    def test_gc_content_100(self):
        rec = SequenceRecord(id="s1", description="test", sequence="GCGC")
        assert rec.gc_content == pytest.approx(100.0)

    def test_gc_content_zero(self):
        rec = SequenceRecord(id="s1", description="test", sequence="AATT")
        assert rec.gc_content == pytest.approx(0.0)

    def test_gc_content_empty(self):
        rec = SequenceRecord(id="s1", description="test", sequence="")
        assert rec.gc_content == 0.0


class TestParseFasta:
    def test_parses_two_records(self, sample_fasta: Path):
        records = list(parse_fasta(sample_fasta))
        assert len(records) == 2

    def test_first_record_id(self, sample_fasta: Path):
        records = list(parse_fasta(sample_fasta))
        assert records[0].id == "seq1"

    def test_first_record_description(self, sample_fasta: Path):
        records = list(parse_fasta(sample_fasta))
        assert "E. coli" in records[0].description

    def test_sequence_no_newlines(self, sample_fasta: Path):
        records = list(parse_fasta(sample_fasta))
        assert "\n" not in records[0].sequence

    def test_multiline_concatenated(self, sample_fasta: Path):
        records = list(parse_fasta(sample_fasta))
        assert records[0].length == 80  # 44 + 36


class TestFilterReads:
    def test_filters_short_reads(self, multi_read_fasta: Path):
        # min_length=50 should keep read1 (100), read2 (200), read3 (300) and drop short (2)
        records = list(filter_reads(multi_read_fasta, min_length=50))
        assert len(records) == 3

    def test_strict_filter(self, multi_read_fasta: Path):
        # min_length=250 should only keep read3 (300)
        records = list(filter_reads(multi_read_fasta, min_length=250))
        assert len(records) == 1
        assert records[0].id == "read3"

    def test_all_pass_relaxed(self, multi_read_fasta: Path):
        records = list(filter_reads(multi_read_fasta, min_length=1))
        assert len(records) == 4


class TestComputeFastaStats:
    def test_total_reads(self, multi_read_fasta: Path):
        stats = compute_fasta_stats(multi_read_fasta)
        assert stats.total_reads == 4

    def test_total_bases(self, multi_read_fasta: Path):
        stats = compute_fasta_stats(multi_read_fasta)
        assert stats.total_bases == 100 + 200 + 300 + 2  # 602

    def test_mean_read_length(self, multi_read_fasta: Path):
        stats = compute_fasta_stats(multi_read_fasta)
        assert stats.mean_read_length == pytest.approx(602 / 4)

    def test_min_max(self, multi_read_fasta: Path):
        stats = compute_fasta_stats(multi_read_fasta)
        assert stats.min_read_length == 2
        assert stats.max_read_length == 300

    def test_n50(self, multi_read_fasta: Path):
        stats = compute_fasta_stats(multi_read_fasta)
        # Sorted desc: 300, 200, 100, 2.  Total=602, half=301.
        # cumsum: 300 (< 301), 300+200=500 (>= 301) -> N50 = 200
        assert stats.n50 == 200

    def test_gc_content(self, multi_read_fasta: Path):
        stats = compute_fasta_stats(multi_read_fasta)
        # read1: 0% GC, read2: 100% GC, read3: 50% GC, short: 0% GC
        assert stats.mean_gc_content == pytest.approx((0 + 100 + 50 + 0) / 4)

    def test_empty_file(self, tmp_path: Path):
        empty = tmp_path / "empty.nucleotide_fasta_protein_homolog_model.fasta"
        empty.write_text("")
        stats = compute_fasta_stats(empty)
        assert stats.total_reads == 0
        assert stats.total_bases == 0
