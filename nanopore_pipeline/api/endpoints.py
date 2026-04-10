"""
FastAPI endpoints for the Nanopore Pathogenicity Diagnostic Pipeline.
Provides REST API for:
  - Full pipeline execution (FASTA in -> report out)
  - Sample management (upload, list, detail)
  - Running alignments
  - Retrieving classification results
  - Comparative analysis between samples
  - Generating visualisation reports
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, UploadFile, File, HTTPException, Query
from fastapi.responses import HTMLResponse
from pydantic import BaseModel

from config import settings
from nanopore_pipeline.db.manager import DatabaseManager
from nanopore_pipeline.alignment.wrapper import AlignmentWrapper
from nanopore_pipeline.classifier.pathogenicity import PathogenicityClassifier
from nanopore_pipeline.reporter.visualisation import PathogenicityReporter
from nanopore_pipeline.pipeline_runner import PipelineRunner
from nanopore_pipeline.utils.logging_config import setup_logging

logger = setup_logging("api.endpoints")

app = FastAPI(
    title="Nanopore Pathogenicity Diagnostic Pipeline",
    description=(
        "REST API for classifying Nanopore sequencing results (FASTA) against "
        "VFDB (virulence factors) and CARD (antibiotic resistance) databases."
    ),
    version="0.2.0",
)

# ── Lazy-initialised singletons ───────────────────────────────────────────────
_db_manager: Optional[DatabaseManager] = None
_aligner: Optional[AlignmentWrapper] = None
_classifier: Optional[PathogenicityClassifier] = None
_reporter: Optional[PathogenicityReporter] = None
_runner: Optional[PipelineRunner] = None


def _get_db() -> DatabaseManager:
    global _db_manager
    if _db_manager is None:
        _db_manager = DatabaseManager()
    return _db_manager


def _get_aligner() -> AlignmentWrapper:
    global _aligner
    if _aligner is None:
        _aligner = AlignmentWrapper()
    return _aligner


def _get_classifier() -> PathogenicityClassifier:
    global _classifier
    if _classifier is None:
        _classifier = PathogenicityClassifier()
    return _classifier


def _get_reporter() -> PathogenicityReporter:
    global _reporter
    if _reporter is None:
        _reporter = PathogenicityReporter()
    return _reporter


def _get_runner() -> PipelineRunner:
    global _runner
    if _runner is None:
        _runner = PipelineRunner()
    return _runner


# ── Pydantic schemas ──────────────────────────────────────────────────────────

class SampleOut(BaseModel):
    id: int
    sample_name: str
    description: str | None
    total_reads: int
    total_bases: int
    mean_read_length: float | None

    class Config:
        from_attributes = True


class ClassificationOut(BaseModel):
    category: str
    gene_count: int
    total_hits: int
    mean_identity: float | None
    mean_coverage: float | None
    tpm: float | None
    risk_level: str | None


class ComparisonOut(BaseModel):
    category: str
    sample_a_hits: int
    sample_b_hits: int
    sample_a_tpm: float
    sample_b_tpm: float
    fold_change: float
    sample_a_risk: str
    sample_b_risk: str


class PipelineRunOut(BaseModel):
    sample_name: str
    sample_id: int
    total_reads: int
    total_bases: int
    mean_read_length: float
    n50: int
    mean_gc_content: float
    vf_hits: int
    amr_hits: int
    categories_found: int
    profiles: list[ClassificationOut]
    report_path: str | None
    json_path: str | None
    assembly_performed: bool = False
    assembly_contigs: int | None = None
    assembly_n50: int | None = None
    assembly_status: str | None = None


# ── Health check ──────────────────────────────────────────────────────────────

@app.get("/", tags=["health"])
async def root():
    return {"status": "ok", "service": "Nanopore Pathogenicity Pipeline", "version": "0.2.0"}


# ── Full pipeline (single endpoint) ──────────────────────────────────────────

@app.post("/run", tags=["pipeline"], response_model=PipelineRunOut)
async def run_pipeline(
    name: str,
    file: UploadFile = File(...),
    description: str = "",
    tool: str = Query("blastn", enum=["blastn", "diamond"]),
    skip_vf: bool = False,
    skip_amr: bool = False,
    perform_assembly: bool = False,
):
    """
    Full pipeline: upload a FASTA file and get back pathogenicity classification.

    1. Computes read statistics (count, N50, GC content)
    1.5. (Optional) Performs Flye assembly to reduce read overlap inflation
    2. Registers sample metadata in database
    3. Aligns reads (or contigs) against VFDB (virulence factors)
    4. Aligns reads (or contigs) against CARD (antibiotic resistance)
    5. Classifies hits into pathogenicity categories with TPM
    6. Generates interactive HTML reports and JSON export
    """
    # Save uploaded FASTA
    dest = settings.SAMPLES_DIR / file.filename
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as out:
        shutil.copyfileobj(file.file, out)

    runner = _get_runner()
    result = runner.run(
        sample_name=name,
        fasta_path=dest,
        description=description,
        alignment_tool=tool,
        skip_vf=skip_vf,
        skip_amr=skip_amr,
        perform_assembly=perform_assembly,
    )

    return PipelineRunOut(
        sample_name=result.sample_name,
        sample_id=result.sample_id,
        total_reads=result.fasta_stats.total_reads,
        total_bases=result.fasta_stats.total_bases,
        mean_read_length=result.fasta_stats.mean_read_length,
        n50=result.fasta_stats.n50,
        mean_gc_content=result.fasta_stats.mean_gc_content,
        vf_hits=result.vf_hits,
        amr_hits=result.amr_hits,
        categories_found=len(result.profiles),
        profiles=[
            ClassificationOut(
                category=p.category,
                gene_count=p.unique_genes,
                total_hits=p.hit_count,
                mean_identity=p.mean_identity,
                mean_coverage=p.mean_coverage,
                tpm=p.tpm,
                risk_level=p.risk_level,
            )
            for p in result.profiles
        ],
        report_path=str(result.report_path) if result.report_path else None,
        json_path=str(result.json_path) if result.json_path else None,
        assembly_performed=perform_assembly,
        assembly_contigs=result.assembly_result.contigs if result.assembly_result and result.assembly_result.success else None,
        assembly_n50=result.assembly_result.n50 if result.assembly_result and result.assembly_result.success else None,
        assembly_status="success" if result.assembly_result and result.assembly_result.success else ("failed" if result.assembly_result else "skipped"),
    )


# ── Reference data management ────────────────────────────────────────────────

@app.post("/db/init", tags=["database"])
async def init_database():
    db = _get_db()
    return {"message": "Database initialised", "url": db.db_url}


@app.post("/db/load-vfdb", tags=["database"])
async def load_vfdb(fasta_path: Optional[str] = None):
    db = _get_db()
    path = Path(fasta_path) if fasta_path else None
    count = db.load_vfdb(path)
    return {"message": f"Loaded {count} VFDB entries"}


@app.post("/db/load-card", tags=["database"])
async def load_card(fasta_path: Optional[str] = None):
    db = _get_db()
    path = Path(fasta_path) if fasta_path else None
    count = db.load_card(path)
    return {"message": f"Loaded {count} CARD entries"}


# ── Sample management ────────────────────────────────────────────────────────

@app.post("/samples/upload", tags=["samples"])
async def upload_sample(
    name: str,
    description: str = "",
    file: UploadFile = File(...),
):
    """Upload a FASTA file and register it as a sample (without running the pipeline)."""
    dest = settings.SAMPLES_DIR / file.filename
    dest.parent.mkdir(parents=True, exist_ok=True)

    with open(dest, "wb") as out:
        shutil.copyfileobj(file.file, out)

    db = _get_db()
    sample = db.register_sample(
        name=name,
        file_path=str(dest),
        description=description,
    )
    logger.info("Uploaded sample '%s' -> %s", name, dest)
    return {"sample_id": sample.id, "file_path": str(dest)}


@app.get("/samples", tags=["samples"], response_model=list[SampleOut])
async def list_samples():
    db = _get_db()
    return db.list_samples()


@app.get("/samples/{sample_name}", tags=["samples"])
async def get_sample(sample_name: str):
    db = _get_db()
    sample = db.get_sample(sample_name)
    if not sample:
        raise HTTPException(404, f"Sample '{sample_name}' not found")
    return SampleOut.model_validate(sample)


# ── Alignment (step-by-step) ─────────────────────────────────────────────────

@app.post("/align", tags=["alignment"])
async def align_sample(
    sample_name: str,
    db_type: str = Query("VF", enum=["VF", "AMR"]),
    tool: str = Query("blastn", enum=["blastn", "diamond"]),
    blast_db_path: str = "",
):
    db = _get_db()
    sample = db.get_sample(sample_name)
    if not sample:
        raise HTTPException(404, f"Sample '{sample_name}' not found")

    if not blast_db_path:
        if db_type == "VF":
            blast_db_path = str(settings.VFDB_DIR / "VFDB_setA_nt.fas.blastdb")
        else:
            blast_db_path = str(settings.CARD_DIR / "nucleotide_fasta_protein_homolog_model.nucleotide_fasta_protein_homolog_model.fasta.blastdb")

    aligner = _get_aligner()
    hit_count, _ = aligner.align_sample(
        sample_path=Path(sample.file_path),
        blast_db_path=Path(blast_db_path),
        sample_id=sample.id,
        hit_type=db_type,
        tool=tool,
    )
    return {"sample": sample_name, "hit_type": db_type, "hits_stored": hit_count}


# ── Classification (step-by-step) ────────────────────────────────────────────

@app.post("/classify/{sample_name}", tags=["classification"])
async def classify_sample(sample_name: str):
    db = _get_db()
    sample = db.get_sample(sample_name)
    if not sample:
        raise HTTPException(404, f"Sample '{sample_name}' not found")

    classifier = _get_classifier()
    profiles = classifier.classify_sample(sample.id)
    return [
        ClassificationOut(
            category=p.category,
            gene_count=p.unique_genes,
            total_hits=p.hit_count,
            mean_identity=p.mean_identity,
            mean_coverage=p.mean_coverage,
            tpm=p.tpm,
            risk_level=p.risk_level,
        )
        for p in profiles
    ]


# ── Comparison ────────────────────────────────────────────────────────────────

@app.get("/compare", tags=["comparison"])
async def compare_samples(sample_a: str, sample_b: str):
    db = _get_db()
    sa = db.get_sample(sample_a)
    sb = db.get_sample(sample_b)
    if not sa:
        raise HTTPException(404, f"Sample '{sample_a}' not found")
    if not sb:
        raise HTTPException(404, f"Sample '{sample_b}' not found")

    classifier = _get_classifier()
    comparison = classifier.compare_samples(sa.id, sb.id)

    return [
        ComparisonOut(
            category=cat,
            sample_a_hits=data["sample_a_hits"],
            sample_b_hits=data["sample_b_hits"],
            sample_a_tpm=data["sample_a_tpm"],
            sample_b_tpm=data["sample_b_tpm"],
            fold_change=data["fold_change_a_vs_b"],
            sample_a_risk=data["sample_a_risk"],
            sample_b_risk=data["sample_b_risk"],
        )
        for cat, data in comparison.items()
    ]


# ── Visualisation ─────────────────────────────────────────────────────────────

@app.get("/report/{sample_name}", tags=["reports"], response_class=HTMLResponse)
async def generate_report(sample_name: str):
    db = _get_db()
    sample = db.get_sample(sample_name)
    if not sample:
        raise HTTPException(404, f"Sample '{sample_name}' not found")

    classifier = _get_classifier()
    profiles = classifier.classify_sample(sample.id)

    reporter = _get_reporter()
    fig = reporter.plot_risk_dashboard(profiles, sample_name, save_html=True)
    reporter.export_json(profiles, sample_name)

    return fig.to_html(full_html=True)


@app.get("/report/compare/{sample_a}/{sample_b}", tags=["reports"], response_class=HTMLResponse)
async def generate_comparison_report(sample_a: str, sample_b: str):
    db = _get_db()
    sa = db.get_sample(sample_a)
    sb = db.get_sample(sample_b)
    if not sa:
        raise HTTPException(404, f"Sample '{sample_a}' not found")
    if not sb:
        raise HTTPException(404, f"Sample '{sample_b}' not found")

    classifier = _get_classifier()
    comparison = classifier.compare_samples(sa.id, sb.id)

    reporter = _get_reporter()
    fig = reporter.plot_comparison_heatmap(comparison, sample_a, sample_b)
    return fig.to_html(full_html=True)
