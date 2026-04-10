"""
Tests for the Alignment Wrapper (parsing and filtering only -- no external tools).
"""

import tempfile
from pathlib import Path

import pytest

from nanopore_pipeline.alignment.wrapper import AlignmentWrapper, BlastHitRow


@pytest.fixture
def blast6_file(tmp_path: Path) -> Path:
    """BLAST tabular output (outfmt 6 with qlen/slen)."""
    rows = [
        "read1\tVF001\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200\t150\t120",
        "read2\tVF002\t60.0\t30\t12\t0\t1\t30\t1\t30\t1e-5\t50\t100\t80",
        "read3\tVF003\t85.0\t80\t12\t0\t1\t80\t1\t80\t1e-30\t150\t90\t100",
    ]
    out = tmp_path / "test.blast6"
    out.write_text("\n".join(rows) + "\n")
    return out


class TestParseBlast6:
    def test_parse_count(self, blast6_file: Path):
        wrapper = AlignmentWrapper.__new__(AlignmentWrapper)
        hits = wrapper.parse_blast6(blast6_file)
        assert len(hits) == 3

    def test_parse_fields(self, blast6_file: Path):
        wrapper = AlignmentWrapper.__new__(AlignmentWrapper)
        hits = wrapper.parse_blast6(blast6_file)
        h = hits[0]
        assert h.qseqid == "read1"
        assert h.sseqid == "VF001"
        assert h.pident == pytest.approx(95.0)
        assert h.length == 100
        assert h.evalue == pytest.approx(1e-50)
        assert h.qlen == 150
        assert h.slen == 120


class TestFilterHits:
    def test_identity_filter(self, blast6_file: Path):
        wrapper = AlignmentWrapper.__new__(AlignmentWrapper)
        hits = wrapper.parse_blast6(blast6_file)
        filtered = wrapper.filter_hits(hits, min_identity=70.0, min_coverage=0)
        # read1 (95%) and read3 (85%) pass; read2 (60%) fails
        assert len(filtered) == 2

    def test_coverage_filter(self, blast6_file: Path):
        wrapper = AlignmentWrapper.__new__(AlignmentWrapper)
        hits = wrapper.parse_blast6(blast6_file)
        filtered = wrapper.filter_hits(hits, min_identity=0, min_coverage=60)
        # read1: 100/150=66.7% pass; read2: 30/100=30% fail; read3: 80/90=88.9% pass
        assert len(filtered) == 2
        ids = {h.qseqid for h in filtered}
        assert "read1" in ids
        assert "read3" in ids

    def test_combined_filter(self, blast6_file: Path):
        wrapper = AlignmentWrapper.__new__(AlignmentWrapper)
        hits = wrapper.parse_blast6(blast6_file)
        filtered = wrapper.filter_hits(hits, min_identity=90.0, min_coverage=60)
        # Only read1 passes both
        assert len(filtered) == 1
        assert filtered[0].qseqid == "read1"

    def test_empty_input(self):
        wrapper = AlignmentWrapper.__new__(AlignmentWrapper)
        filtered = wrapper.filter_hits([], min_identity=70, min_coverage=50)
        assert filtered == []
