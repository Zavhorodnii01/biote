"""
Module 2 -- Alignment Wrapper
Python interface to run BLAST+ / DIAMOND on Nanopore reads against the
local VFDB and CARD databases.  Handles:
  - makeblastdb / diamond makedb
  - blastn / diamond blastx
  - Output parsing into AlignmentHit ORM objects
Alignment parameters are tuned for Nanopore error profiles (~10% error).
"""

from __future__ import annotations

import csv
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from sqlalchemy.orm import Session

from config import settings
from nanopore_pipeline.models.database import (
    AlignmentHit, Sample, VirulenceFactorEntry, AMREntry,
    get_session_factory,
)
from nanopore_pipeline.utils.logging_config import setup_logging

logger = setup_logging("alignment.wrapper")


@dataclass
class BlastHitRow:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int = 0
    slen: int = 0


class AlignmentWrapper:
    """Runs BLAST+/DIAMOND against VFDB or CARD and stores hits in DB."""

    BLAST_FIELDS = (
        "qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen"
    )

    def __init__(self, db_url: str = settings.DATABASE_URL):
        self._SessionFactory = get_session_factory(db_url)
        logger.info("AlignmentWrapper initialised")

    def _session(self) -> Session:
        return self._SessionFactory()

    # ── BLAST+ database creation ──────────────────────────────────────────

    @staticmethod
    def make_blast_db(fasta_path: Path, db_type: str = "nucl", title: str = "", parse_seqids: bool = True) -> Path:
        out_db = fasta_path.with_suffix(".blastdb")
        cmd = [
            "makeblastdb",
            "-in", str(fasta_path),
            "-dbtype", db_type,
            "-out", str(out_db),
            "-title", title or fasta_path.stem,
        ]
        if parse_seqids:
            cmd.append("-parse_seqids")
        logger.info("Creating BLAST DB: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.debug("makeblastdb stdout: %s", result.stdout)
        return out_db

    # ── DIAMOND database creation ─────────────────────────────────────────

    @staticmethod
    def make_diamond_db(fasta_path: Path) -> Path:
        out_db = fasta_path.with_suffix(".dmnd")
        cmd = [
            "diamond", "makedb",
            "--in", str(fasta_path),
            "--db", str(out_db),
        ]
        logger.info("Creating DIAMOND DB: %s", " ".join(cmd))
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return out_db

    # ── Run BLAST ─────────────────────────────────────────────────────────

    def run_blastn(
        self,
        query_path: Path,
        db_path: Path,
        evalue: float = settings.BLAST_EVALUE,
        num_threads: int = settings.BLAST_NUM_THREADS,
        output_path: Optional[Path] = None,
        word_size: int = settings.BLAST_WORD_SIZE,
    ) -> Path:
        if output_path is None:
            output_path = Path(tempfile.mktemp(suffix=".blast6"))

        cmd = [
            "blastn",
            "-query", str(query_path),
            "-db", str(db_path),
            "-out", str(output_path),
            "-outfmt", f"6 {self.BLAST_FIELDS}",
            "-evalue", str(evalue),
            "-num_threads", str(num_threads),
            "-word_size", str(word_size),
            "-reward", str(settings.BLAST_REWARD),
            "-penalty", str(settings.BLAST_PENALTY),
            "-gapopen", str(settings.BLAST_GAPOPEN),
            "-gapextend", str(settings.BLAST_GAPEXTEND),
            "-max_target_seqs", "5",
            "-dust", "no",
        ]
        logger.info("Running blastn: query=%s db=%s", query_path, db_path)
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info("BLAST complete -> %s", output_path)
        return output_path

    # ── Run DIAMOND ───────────────────────────────────────────────────────

    def run_diamond_blastx(
        self,
        query_path: Path,
        db_path: Path,
        evalue: float = settings.DIAMOND_EVALUE,
        output_path: Optional[Path] = None,
    ) -> Path:
        if output_path is None:
            output_path = Path(tempfile.mktemp(suffix=".diamond6"))

        cmd = [
            "diamond", "blastx",
            "--query", str(query_path),
            "--db", str(db_path),
            "--out", str(output_path),
            "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
            "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
            "evalue", "bitscore", "qlen", "slen",
            "--evalue", str(evalue),
            "--id", str(settings.DIAMOND_MIN_IDENTITY),
            "--query-cover", str(settings.DIAMOND_MIN_COVERAGE),
            "--block-size", str(settings.DIAMOND_BLOCK_SIZE),
            "--index-chunks", str(settings.DIAMOND_INDEX_CHUNKS),
            "--max-target-seqs", "5",
            "--long-reads",  # Nanopore mode
        ]
        logger.info("Running diamond blastx: query=%s db=%s", query_path, db_path)
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info("DIAMOND complete -> %s", output_path)
        return output_path

    # ── Parse BLAST6 output ───────────────────────────────────────────────

    @staticmethod
    def parse_blast6(result_path: Path) -> list[BlastHitRow]:
        hits: list[BlastHitRow] = []
        with open(result_path, "r") as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                if len(row) < 12:
                    continue
                hits.append(BlastHitRow(
                    qseqid=row[0],
                    sseqid=row[1],
                    pident=float(row[2]),
                    length=int(row[3]),
                    mismatch=int(row[4]),
                    gapopen=int(row[5]),
                    qstart=int(row[6]),
                    qend=int(row[7]),
                    sstart=int(row[8]),
                    send=int(row[9]),
                    evalue=float(row[10]),
                    bitscore=float(row[11]),
                    qlen=int(row[12]) if len(row) > 12 else 0,
                    slen=int(row[13]) if len(row) > 13 else 0,
                ))
        logger.info("Parsed %d hits from %s", len(hits), result_path)
        return hits

    # ── Filter hits ───────────────────────────────────────────────────────

    @staticmethod
    def filter_hits(
        hits: list[BlastHitRow],
        min_identity: float = settings.BLAST_MIN_IDENTITY,
        min_coverage: float = settings.BLAST_MIN_COVERAGE,
    ) -> list[BlastHitRow]:
        filtered = []
        for h in hits:
            if h.pident < min_identity:
                continue
            if h.slen > 0:
                coverage = (h.length / h.slen) * 100
                if coverage < min_coverage:
                    continue
            filtered.append(h)
        logger.info("Filtered to %d / %d hits (identity>=%.1f%%, coverage>=%.1f%%)",
                     len(filtered), len(hits), min_identity, min_coverage)
        return filtered

    # ── Best-hit selection ────────────────────────────────────────────────

    @staticmethod
    def select_best_hits(hits: list[BlastHitRow]) -> list[BlastHitRow]:
        """
        For each query sequence, keep only the single best hit (highest bitscore).

        When a query gene matches multiple entries in VFDB/CARD, only the most
        confident match is retained. This prevents double-counting and makes the
        per-category statistics cleaner for the thesis analysis.

        Supervisor note (meeting): "you will have the same gene, maybe, from the
        query aligned against seven different hits from the FireLens database.
        How would you know which one is good?"
        """
        best: dict[str, BlastHitRow] = {}
        for h in hits:
            existing = best.get(h.qseqid)
            if existing is None or h.bitscore > existing.bitscore:
                best[h.qseqid] = h
        result = list(best.values())
        logger.info(
            "Best-hit selection: %d hits -> %d unique query sequences retained",
            len(hits), len(result),
        )
        return result

    # ── Store hits in DB ──────────────────────────────────────────────────

    def store_hits(
        self,
        hits: list[BlastHitRow],
        sample_id: int,
        hit_type: str,  # "VF" or "AMR"
    ) -> int:
        session = self._session()
        count = 0
        try:
            for h in hits:
                qcov = (h.length / h.qlen * 100) if h.qlen else None
                scov = (h.length / h.slen * 100) if h.slen else None

                ref_id = None
                amr_ref_id = None
                if hit_type == "VF":
                    vf = session.query(VirulenceFactorEntry).filter_by(
                        vfdb_id=h.sseqid
                    ).first()
                    ref_id = vf.id if vf else None
                elif hit_type == "AMR":
                    amr = session.query(AMREntry).filter_by(
                        card_id=h.sseqid
                    ).first()
                    amr_ref_id = amr.id if amr else None

                db_hit = AlignmentHit(
                    sample_id=sample_id,
                    query_id=h.qseqid,
                    hit_type=hit_type,
                    reference_id=ref_id,
                    amr_reference_id=amr_ref_id,
                    subject_id=h.sseqid,
                    percent_identity=h.pident,
                    alignment_length=h.length,
                    mismatches=h.mismatch,
                    gap_opens=h.gapopen,
                    query_start=h.qstart,
                    query_end=h.qend,
                    subject_start=h.sstart,
                    subject_end=h.send,
                    evalue=h.evalue,
                    bitscore=h.bitscore,
                    query_length=h.qlen,
                    subject_length=h.slen,
                    query_coverage=qcov,
                    subject_coverage=scov,
                )
                session.add(db_hit)
                count += 1

                if count % 1000 == 0:
                    session.commit()

            session.commit()
            logger.info("Stored %d %s hits for sample_id=%d", count, hit_type, sample_id)
        except Exception:
            session.rollback()
            logger.exception("Error storing hits")
            raise
        finally:
            session.close()
        return count

    # ── High-level: align sample against a DB ─────────────────────────────

    def align_sample(
        self,
        sample_path: Path,
        blast_db_path: Path,
        sample_id: int,
        hit_type: str = "VF",
        tool: str = "blastn",
        word_size: int = settings.BLAST_WORD_SIZE,
    ) -> tuple[int, list[BlastHitRow]]:
        """
        Run alignment, filter hits, apply best-hit selection, and store in DB.

        Returns:
            (hits_stored, best_hits) - count of stored hits and the filtered
            BlastHitRow list (used for quality-distribution plots).
        """
        if tool == "blastn":
            result_file = self.run_blastn(sample_path, blast_db_path, word_size=word_size)
        elif tool == "diamond":
            result_file = self.run_diamond_blastx(sample_path, blast_db_path)
        else:
            raise ValueError(f"Unknown alignment tool: {tool}")

        raw_hits = self.parse_blast6(result_file)
        filtered = self.filter_hits(raw_hits)
        best_hits = self.select_best_hits(filtered)
        count = self.store_hits(best_hits, sample_id, hit_type)
        return count, best_hits
