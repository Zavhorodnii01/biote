"""
SQLAlchemy ORM models for the Nanopore Pathogenicity Pipeline.
Stores VFDB/CARD reference entries, sample metadata, alignment hits,
and classification results.
"""

from __future__ import annotations

from datetime import datetime, timezone

from sqlalchemy import (
    Column, Integer, String, Float, Text, DateTime, Boolean,
    ForeignKey, Index, create_engine, Enum,
)
from sqlalchemy.orm import declarative_base, relationship, sessionmaker

from config.settings import DATABASE_URL

Base = declarative_base()


# ── Reference tables ──────────────────────────────────────────────────────────

class VirulenceFactorEntry(Base):
    __tablename__ = "virulence_factors"

    id = Column(Integer, primary_key=True, autoincrement=True)
    vfdb_id = Column(String(64), unique=True, nullable=False, index=True)
    gene_name = Column(String(128), nullable=False)
    product = Column(Text)
    organism = Column(String(256))
    category = Column(String(64))           # Adherence, Motility, Exotoxin, etc.
    subcategory = Column(String(128))
    sequence = Column(Text, nullable=False)
    seq_length = Column(Integer, nullable=False)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))

    hits = relationship("AlignmentHit", back_populates="vf_entry",
                        foreign_keys="AlignmentHit.reference_id")

    __table_args__ = (
        Index("ix_vf_category", "category"),
        Index("ix_vf_organism", "organism"),
    )


class AMREntry(Base):
    __tablename__ = "amr_entries"

    id = Column(Integer, primary_key=True, autoincrement=True)
    card_id = Column(String(64), unique=True, nullable=False, index=True)
    gene_name = Column(String(128), nullable=False)
    aro_accession = Column(String(64))
    drug_class = Column(String(128))
    resistance_mechanism = Column(String(256))
    organism = Column(String(256))
    sequence = Column(Text, nullable=False)
    seq_length = Column(Integer, nullable=False)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))

    hits = relationship("AlignmentHit", back_populates="amr_entry",
                        foreign_keys="AlignmentHit.amr_reference_id")

    __table_args__ = (
        Index("ix_amr_drug_class", "drug_class"),
        Index("ix_amr_organism", "organism"),
    )


# ── Sample & results tables ──────────────────────────────────────────────────

class Sample(Base):
    __tablename__ = "samples"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sample_name = Column(String(256), unique=True, nullable=False, index=True)
    description = Column(Text)
    file_path = Column(String(1024))
    total_reads = Column(Integer, default=0)
    total_bases = Column(Integer, default=0)
    mean_read_length = Column(Float)
    n50 = Column(Integer)
    mean_gc_content = Column(Float)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))

    # Assembly metadata
    assembly_performed = Column(Boolean, default=False)
    assembly_path = Column(String(1024), nullable=True)
    assembly_contigs = Column(Integer, default=0)
    assembly_n50 = Column(Integer, nullable=True)
    assembly_total_bases = Column(Integer, nullable=True)
    assembly_mean_contig_length = Column(Float, nullable=True)
    assembly_status = Column(String(64), nullable=True)  # "success", "failed", "skipped"

    hits = relationship("AlignmentHit", back_populates="sample")
    classifications = relationship("ClassificationResult", back_populates="sample")


class AlignmentHit(Base):
    __tablename__ = "alignment_hits"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sample_id = Column(Integer, ForeignKey("samples.id"), nullable=False, index=True)
    query_id = Column(String(256), nullable=False)
    hit_type = Column(Enum("VF", "AMR", name="hit_type_enum"), nullable=False)

    # VF reference
    reference_id = Column(Integer, ForeignKey("virulence_factors.id"), nullable=True)
    # AMR reference
    amr_reference_id = Column(Integer, ForeignKey("amr_entries.id"), nullable=True)

    subject_id = Column(String(256), nullable=False)
    percent_identity = Column(Float, nullable=False)
    alignment_length = Column(Integer, nullable=False)
    mismatches = Column(Integer)
    gap_opens = Column(Integer)
    query_start = Column(Integer)
    query_end = Column(Integer)
    subject_start = Column(Integer)
    subject_end = Column(Integer)
    evalue = Column(Float, nullable=False)
    bitscore = Column(Float, nullable=False)
    query_length = Column(Integer)
    subject_length = Column(Integer)
    query_coverage = Column(Float)
    subject_coverage = Column(Float)

    sample = relationship("Sample", back_populates="hits")
    vf_entry = relationship("VirulenceFactorEntry", back_populates="hits",
                            foreign_keys=[reference_id])
    amr_entry = relationship("AMREntry", back_populates="hits",
                             foreign_keys=[amr_reference_id])

    __table_args__ = (
        Index("ix_hit_sample_type", "sample_id", "hit_type"),
    )


class ClassificationResult(Base):
    __tablename__ = "classification_results"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sample_id = Column(Integer, ForeignKey("samples.id"), nullable=False, index=True)
    category = Column(String(64), nullable=False)
    gene_count = Column(Integer, default=0)
    total_hits = Column(Integer, default=0)
    mean_identity = Column(Float)
    mean_coverage = Column(Float)
    tpm = Column(Float)
    risk_level = Column(String(32))         # Low, Medium, High, Critical
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))

    sample = relationship("Sample", back_populates="classifications")

    __table_args__ = (
        Index("ix_class_sample_cat", "sample_id", "category"),
    )


# ── Engine / session factory ──────────────────