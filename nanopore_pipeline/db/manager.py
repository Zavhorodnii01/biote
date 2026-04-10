"""
Module 1 -- Database Manager
Downloads, parses, and indexes VFDB / CARD reference data into the local
relational database.  Also provides convenience methods for sample CRUD.
"""

from __future__ import annotations

import gzip
import io
import re
import shutil
import urllib.request
from pathlib import Path
from typing import Optional

from sqlalchemy.orm import Session

from config import settings
from nanopore_pipeline.models.database import (
    AMREntry, VirulenceFactorEntry, Sample,
    init_db, get_session_factory,
)
from nanopore_pipeline.utils.fasta_parser import parse_fasta, SequenceRecord
from nanopore_pipeline.utils.logging_config import setup_logging

logger = setup_logging("db.manager")


class DatabaseManager:
    """Handles all database lifecycle: init, reference loading, sample CRUD."""

    def __init__(self, db_url: str = settings.DATABASE_URL):
        self.db_url = db_url
        init_db(db_url)
        self._SessionFactory = get_session_factory(db_url)
        logger.info("DatabaseManager initialised (url=%s)", db_url)

    # ── helpers ───────────────────────────────────────────────────────────

    def _session(self) -> Session:
        return self._SessionFactory()

    @staticmethod
    def _download(url: str, dest: Path) -> Path:
        logger.info("Downloading %s -> %s", url, dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(url, str(dest))
        logger.info("Download complete: %s", dest)
        return dest

    @staticmethod
    def _decompress_gz(gz_path: Path) -> Path:
        out = gz_path.with_suffix("")
        with gzip.open(gz_path, "rb") as f_in, open(out, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        logger.info("Decompressed %s -> %s", gz_path, out)
        return out

    # ── VFDB ──────────────────────────────────────────────────────────────

    def _categorise_vf(self, description: str) -> tuple[str, str]:
        desc_lower = description.lower()
        for category, keywords in settings.VF_CATEGORIES.items():
            for kw in keywords:
                if kw in desc_lower:
                    return category, kw
        return "Other", ""

    def download_vfdb(self) -> Path:
        dest = settings.VFDB_DIR / "VFDB_setA_nt.fas.gz"
        self._download(settings.VFDB_FASTA_URL, dest)
        return self._decompress_gz(dest)

    def load_vfdb(self, fasta_path: Optional[Path] = None) -> int:
        if fasta_path is None:
            fasta_path = settings.VFDB_DIR / "VFDB_setA_nt.fas"
            if not fasta_path.exists():
                fasta_path = self.download_vfdb()

        session = self._session()
        count = 0
        try:
            for record in parse_fasta(fasta_path):
                existing = session.query(VirulenceFactorEntry).filter_by(
                    vfdb_id=record.id
                ).first()
                if existing:
                    continue

                category, subcat = self._categorise_vf(record.description)

                # Try to extract gene name from description pattern: (gene_name)
                gene_match = re.search(r"\(([^)]+)\)", record.description)
                gene_name = gene_match.group(1) if gene_match else record.id

                # Try to extract organism: [...organism...]
                org_match = re.search(r"\[([^\]]+)\]", record.description)
                organism = org_match.group(1) if org_match else "Unknown"

                entry = VirulenceFactorEntry(
                    vfdb_id=record.id,
                    gene_name=gene_name,
                    product=record.description,
                    organism=organism,
                    category=category,
                    subcategory=subcat,
                    sequence=record.sequence,
                    seq_length=record.length,
                )
                session.add(entry)
                count += 1

                if count % 500 == 0:
                    session.commit()
                    logger.info("Loaded %d VFDB entries...", count)

            session.commit()
            logger.info("VFDB loading complete: %d new entries", count)
        except Exception:
            session.rollback()
            logger.exception("Error loading VFDB data")
            raise
        finally:
            session.close()
        return count

    # ── CARD ──────────────────────────────────────────────────────────────

    def _categorise_amr(self, description: str) -> tuple[str, str]:
        desc_lower = description.lower()
        for category, keywords in settings.AMR_CATEGORIES.items():
            for kw in keywords:
                if kw in desc_lower:
                    return category, kw
        return "Other", ""

    def download_card(self) -> Path:
        dest = settings.CARD_DIR / "card-data.tar.bz2"
        self._download(settings.CARD_DATA_URL, dest)
        return dest

    def load_card(self, fasta_path: Optional[Path] = None) -> int:
        if fasta_path is None:
            fasta_path = settings.CARD_DIR / "nucleotide_fasta_protein_homolog_model.nucleotide_fasta_protein_homolog_model.fasta"
            if not fasta_path.exists():
                logger.warning(
                    "CARD FASTA not found at %s. "
                    "Download and extract CARD data first, or provide a path.",
                    fasta_path,
                )
                return 0

        session = self._session()
        count = 0
        try:
            for record in parse_fasta(fasta_path):
                existing = session.query(AMREntry).filter_by(
                    card_id=record.id
                ).first()
                if existing:
                    continue

                aro_match = re.search(r"ARO:(\d+)", record.description)
                aro = f"ARO:{aro_match.group(1)}" if aro_match else ""

                drug_class, _ = self._categorise_amr(record.description)

                gene_match = re.search(r"\|([^|]+)\|", record.description)
                gene_name = gene_match.group(1).strip() if gene_match else record.id

                org_match = re.search(r"\[([^\]]+)\]", record.description)
                organism = org_match.group(1) if org_match else "Unknown"

                entry = AMREntry(
                    card_id=record.id,
                    gene_name=gene_name,
                    aro_accession=aro,
                    drug_class=drug_class,
                    resistance_mechanism="",
                    organism=organism,
                    sequence=record.sequence,
                    seq_length=record.length,
                )
                session.add(entry)
                count += 1

                if count % 500 == 0:
                    session.commit()
                    logger.info("Loaded %d CARD entries...", count)

            session.commit()
            logger.info("CARD loading complete: %d new entries", count)
        except Exception:
            session.rollback()
            logger.exception("Error loading CARD data")
            raise
        finally:
            session.close()
        return count

    # ── Samples ───────────────────────────────────────────────────────────

    def register_sample(
        self,
        name: str,
        file_path: str,
        description: str = "",
        total_reads: int = 0,
        total_bases: int = 0,
        mean_read_length: float = 0.0,
        n50: int = 0,
        mean_gc_content: float = 0.0,
    ) -> Sample:
        session = self._session()
        try:
            sample = session.query(Sample).filter_by(sample_name=name).first()
            if sample:
                logger.info("Sample '%s' already exists (id=%d)", name, sample.id)
                return sample

            sample = Sample(
                sample_name=name,
                file_path=file_path,
                description=description,
                total_reads=total_reads,
                total_bases=total_bases,
                mean_read_length=mean_read_length,
                n50=n50,
                mean_gc_content=mean_gc_content,
            )
            session.add(sample)
            session.commit()
            logger.info("Registered sample '%s' (id=%d)", name, sample.id)
            return sample
        except Exception:
            session.rollback()
            raise
        finally:
            session.close()

    def get_sample(self, name: str) -> Optional[Sample]:
        session = self._session()
        try:
            return session.query(Sample).filter_by(sample_name=name).first()
        finally:
            session.close()

    def list_samples(self) -> list[Sample]:
        session = self._session()
        try:
            return session.query(Sample).all()
        finally:
            session.close()
