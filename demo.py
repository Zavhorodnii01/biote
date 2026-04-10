"""
Demo Runner -- Full Pipeline on E. coli O157:H7 Reference Genome
=================================================================
Runs the complete pipeline on the real E. coli O157:H7 reference genome.
This is Test 1 (perfect reference) of the supervisor's testing framework.
BLAST+ must be installed and on PATH before running this.

Steps:
  1. Verify reference genome exists
  2. Download + load VFDB into SQLite, build BLAST index
  3. Load CARD into SQLite, build BLAST index
  4. Run PipelineRunner.run() -- real blastn alignment + classification + reports

Usage:
    python demo.py
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import settings
from nanopore_pipeline.db.manager import DatabaseManager
from nanopore_pipeline.alignment.wrapper import AlignmentWrapper
from nanopore_pipeline.pipeline_runner import PipelineRunner

SAMPLE_NAME = "ecoli_O157H7_reference"
DEMO_FASTA  = settings.SAMPLES_DIR / "ecoli_O157H7_reference.fasta"


# ── Step 2+3: setup reference databases ──────────────────────────────────────

def setup_databases() -> tuple[Path, Path]:
    """
    Returns (vfdb_blast_db_path, card_blast_db_path).
    Downloads VFDB if needed, loads both into SQLite, builds BLAST indices.
    """
    db = DatabaseManager(settings.DATABASE_URL)
    aligner = AlignmentWrapper(settings.DATABASE_URL)

    # VFDB
    vfdb_fasta = settings.VFDB_DIR / "VFDB_setA_nt.fas"
    vfdb_blast = vfdb_fasta.with_suffix(".blastdb")
    if not vfdb_fasta.exists():
        print("  Downloading VFDB (this may take a moment)...")
        vfdb_fasta = db.download_vfdb()
    else:
        print(f"  VFDB FASTA already exists: {vfdb_fasta}")

    print("  Loading VFDB into database...")
    n_vf = db.load_vfdb(vfdb_fasta)
    print(f"  VFDB entries loaded: {n_vf}")

    if not Path(str(vfdb_blast) + ".nin").exists():
        print("  Building VFDB BLAST index...")
        aligner.make_blast_db(vfdb_fasta, db_type="nucl", title="VFDB")
    else:
        print("  VFDB BLAST index already exists")

    # CARD
    card_fasta = settings.CARD_DIR / "nucleotide_fasta_protein_homolog_model.fasta"
    card_blast = card_fasta.with_suffix(".blastdb")
    if not card_fasta.exists():
        raise FileNotFoundError(
            f"CARD FASTA not found at {card_fasta}\n"
            "Download from https://card.mcmaster.ca/download -> 'Data' -> extract the .tar.bz2"
        )

    print("  Loading CARD into database...")
    n_amr = db.load_card(card_fasta)
    print(f"  CARD entries loaded: {n_amr}")

    if not Path(str(card_blast) + ".nin").exists():
        print("  Building CARD BLAST index...")
        aligner.make_blast_db(card_fasta, db_type="nucl", title="CARD", parse_seqids=False)
    else:
        print("  CARD BLAST index already exists")

    return vfdb_blast, card_blast


# ── Print results ─────────────────────────────────────────────────────────────

def print_results(result) -> None:
    stats = result.fasta_stats
    print(f"\n  Sample : {result.sample_name}")
    print(f"  Reads  : {stats.total_reads}   Bases: {stats.total_bases:,}")
    print(f"  N50    : {stats.n50}   Mean len: {stats.mean_read_length:.0f}   GC: {stats.mean_gc_content:.1f}%")
    print(f"  VF hits: {result.vf_hits}   AMR hits: {result.amr_hits}")

    if not result.profiles:
        print("\n  No pathogenicity profiles detected.")
        return

    RISK_COLORS = {
        "Critical": "\033[91m", "High": "\033[93m",
        "Medium":   "\033[94m", "Low":  "\033[92m",
    }
    RESET = "\033[0m"
    header = f"\n  {'Category':<30} {'Genes':>5} {'Hits':>5} {'TPM':>9} Risk"
    print(header)
    print("  " + "-" * 55)
    for p in sorted(result.profiles,
                    key=lambda x: ["Low", "Medium", "High", "Critical"].index(x.risk_level),
                    reverse=True):
        color = RISK_COLORS.get(p.risk_level, "")
        print(f"  {p.category:<30} {p.unique_genes:>5} {p.hit_count:>5} "
              f"{p.tpm:>9.1f} {color}{p.risk_level}{RESET}")

    if result.report_path:
        print(f"\n  Reports -> {result.report_path.parent}/")
        print(f"    {result.sample_name}_profile.html")
        print(f"    {result.sample_name}_risk_dashboard.html")
    if result.json_path:
        print(f"    {result.json_path.name}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("\n=== Nanopore Pathogenicity Pipeline -- Demo (Test 1: Perfect Reference) ===\n")
    settings.RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/3] Checking input...")
    if not DEMO_FASTA.exists():
        print(f"ERROR: Reference genome not found at:\n  {DEMO_FASTA}")
        print("Download E. coli O157:H7 from NCBI and place it there.")
        sys.exit(1)
    size_mb = DEMO_FASTA.stat().st_size / 1_000_000
    print(f"  Found: {DEMO_FASTA.name}  ({size_mb:.1f} MB)")

    print("\n[2/3] Setting up reference databases...")
    vfdb_blast, card_blast = setup_databases()

    print("\n[3/3] Running pipeline...")
    runner = PipelineRunner()
    result = runner.run(
        sample_name=SAMPLE_NAME,
        fasta_path=DEMO_FASTA,
        description="E. coli O157:H7 reference genome -- Test 1 perfect reference",
        vfdb_blast_db=vfdb_blast,
        card_blast_db=card_blast,
        alignment_tool="blastn",
    )

    print("\n=== Results ===")
    print_results(result)
    print()


if __name__ == "__main__":
    main()