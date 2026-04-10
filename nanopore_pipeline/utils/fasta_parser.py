"""
FASTA parser for Nanopore reads and reference sequences.
Handles gzipped files transparently.

Primary input format: FASTA (Nanopore basecalled output).
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Generator

from config.settings import NANOPORE_MIN_READ_LENGTH


@dataclass(slots=True)
class SequenceRecord:
    id: str
    description: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def gc_content(self) -> float:
        if not self.sequence:
            return 0.0
        gc = sum(1 for c in self.sequence.upper() if c in ("G", "C"))
        return gc / len(self.sequence) * 100


@dataclass
class FastaStats:
    """Summary statistics computed from a FASTA file."""
    total_reads: int = 0
    total_bases: int = 0
    mean_read_length: float = 0.0
    min_read_length: int = 0
    max_read_length: int = 0
    mean_gc_content: float = 0.0
    n50: int = 0


def _open_file(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_fasta(path: Path) -> Generator[SequenceRecord, None, None]:
    """Parse a FASTA file (plain or gzipped) into SequenceRecord objects."""
    with _open_file(path) as fh:
        header = ""
        seq_parts: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    full = "".join(seq_parts)
                    parts = header.split(None, 1)
                    yield SequenceRecord(
                        id=parts[0],
                        description=parts[1] if len(parts) > 1 else "",
                        sequence=full,
                    )
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header:
            full = "".join(seq_parts)
            parts = header.split(None, 1)
            yield SequenceRecord(
                id=parts[0],
                description=parts[1] if len(parts) > 1 else "",
                sequence=full,
            )


def filter_reads(
    path: Path,
    min_length: int = NANOPORE_MIN_READ_LENGTH,
) -> Generator[SequenceRecord, None, None]:
    """Filter FASTA reads by minimum length."""
    for record in parse_fasta(path):
        if record.length >= min_length:
            yield record


def compute_fasta_stats(path: Path) -> FastaStats:
    """Compute summary statistics for a FASTA file without holding all reads in memory."""
    lengths: list[int] = []
    total_gc = 0.0

    for record in parse_fasta(path):
        lengths.append(record.length)
        total_gc += record.gc_content

    if not lengths:
        return FastaStats()

    lengths_sorted = sorted(lengths, reverse=True)
    total_bases = sum(lengths)
    half = total_bases / 2
    cumsum = 0
    n50 = 0
    for l in lengths_sorted:
        cumsum += l
        if cumsum >= half:
            n50 = l
            break

    return FastaStats(
        total_reads=len(lengths),
        total_bases=total_bases,
        mean_read_length=total_bases / len(lengths),
        min_read_length=min(lengths),
        max_read_length=max(lengths),
        mean_gc_content=total_gc / len(lengths),
        n50=n50,
    )
