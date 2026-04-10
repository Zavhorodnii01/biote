"""
3-Way Performance Comparison
=============================
Runs the full pipeline on all three test inputs and prints a side-by-side
comparison of results -- the core of the supervisor's testing framework.

Test 1 -- Perfect reference  : ecoli_O157H7_reference.fasta
Test 2 -- Simulated Nanopore : ecoli_O157H7_simulated.fasta  (from simulate_reads.py)
Test 3 -- Real Nanopore SRA  : ecoli_O157H7_real_sra.fasta   (from SRA DRR198813)

Usage:
    python run_comparison.py

Skip a test by removing/renaming its FASTA -- the script skips missing files.
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
from nanopore_pipeline.pipeline_runner import PipelineRunner, PipelineResult

# ── Test inputs ───────────────────────────────────────────────────────────────
TESTS = [
    {
        "name":      "test1_perfect_reference",
        "label":     "Test 1 -- Perfect Reference",
        "fasta":     settings.SAMPLES_DIR / "ecoli_O157H7_reference.fasta",
        "word_size": 11,   # standard BLAST — perfect DNA, no Nanopore errors
    },
    {
        "name":      "test2_simulated_nanopore",
        "label":     "Test 2 -- Simulated Nanopore",
        "fasta":     settings.SAMPLES_DIR / "ecoli_O157H7_simulated_5k.fasta",
        "word_size": 7,    # Nanopore-optimised — ~10% error rate
    },
    {
        "name":      "test3_real_sra",
        "label":     "Test 3 -- Real Nanopore (SRA)",
        "fasta":     settings.SAMPLES_DIR / "ecoli_O157H7_real_sra_5k.fasta",
        "word_size": 7,    # Nanopore-optimised — real sequencing errors
    },
]


# ── Setup reference databases (once) ─────────────────────────────────────────

def setup_databases() -> tuple[Path, Path]:
    db = DatabaseManager(settings.DATABASE_URL)
    aligner = AlignmentWrapper(settings.DATABASE_URL)

    vfdb_fasta = settings.VFDB_DIR / "VFDB_setA_nt.fas"
    vfdb_blast = vfdb_fasta.with_suffix(".blastdb")

    if not vfdb_fasta.exists():
        print("  Downloading VFDB...")
        vfdb_fasta = db.download_vfdb()
    if not db._session().query(
        __import__("nanopore_pipeline.models.database", fromlist=["VirulenceFactorEntry"]).VirulenceFactorEntry
    ).first():
        print("  Loading VFDB into database...")
        db.load_vfdb(vfdb_fasta)
    if not Path(str(vfdb_blast) + ".nin").exists():
        print("  Building VFDB BLAST index...")
        aligner.make_blast_db(vfdb_fasta, db_type="nucl", title="VFDB")

    card_fasta = settings.CARD_DIR / "nucleotide_fasta_protein_homolog_model.fasta"
    card_blast = card_fasta.with_suffix(".blastdb")

    if not card_fasta.exists():
        raise FileNotFoundError(f"CARD FASTA not found at {card_fasta}")
    if not Path(str(card_blast) + ".nin").exists():
        print("  Building CARD BLAST index...")
        aligner.make_blast_db(card_fasta, db_type="nucl", title="CARD", parse_seqids=False)
    if not db._session().query(
        __import__("nanopore_pipeline.models.database", fromlist=["AMREntry"]).AMREntry
    ).first():
        print("  Loading CARD into database...")
        db.load_card(card_fasta)

    return vfdb_blast, card_blast


# ── Print comparison table ────────────────────────────────────────────────────

def print_comparison(results: list[tuple[str, PipelineResult]]) -> None:
    RISK_COLORS = {
        "Critical": "\033[91m", "High": "\033[93m",
        "Medium":   "\033[94m", "Low":  "\033[92m", "None": "",
    }
    RESET = "\033[0m"

    print("\n" + "=" * 80)
    print("PERFORMANCE COMPARISON -- E. coli O157:H7")
    print("=" * 80)

    # per-sample summary
    print(f"\n{'Metric':<25}", end="")
    for label, _ in results:
        print(f"  {label:<26}", end="")
    print()
    print("-" * (25 + 28 * len(results)))

    def row(metric, getter):
        print(f"{metric:<25}", end="")
        for _, r in results:
            print(f"  {getter(r):<26}", end="")
        print()

    row("Total reads",      lambda r: f"{r.fasta_stats.total_reads:,}")
    row("Total bases",      lambda r: f"{r.fasta_stats.total_bases:,}")
    row("Mean read length", lambda r: f"{r.fasta_stats.mean_read_length:.0f} bp")
    row("N50",              lambda r: f"{r.fasta_stats.n50:,} bp")
    row("GC content",       lambda r: f"{r.fasta_stats.mean_gc_content:.1f}%")
    row("VF hits",          lambda r: str(r.vf_hits))
    row("AMR hits",         lambda r: str(r.amr_hits))

    # collect all categories across all results
    all_cats = sorted({p.category for _, r in results for p in r.profiles})
    if not all_cats:
        print("\n  No pathogenicity profiles detected in any test.")
        return

    print(f"\n{'Category':<28}", end="")
    for label, _ in results:
        print(f"  {label[:12]:<12} {'TPM':>8} {'Risk':<8}", end="")
    print()
    print("-" * (28 + 30 * len(results)))

    for cat in all_cats:
        print(f"{cat:<28}", end="")
        for _, r in results:
            profile = next((p for p in r.profiles if p.category == cat), None)
            if profile:
                color = RISK_COLORS.get(profile.risk_level, "")
                print(f"  {profile.hit_count:<12} {profile.tpm:>8.1f} "
                      f"{color}{profile.risk_level:<8}{RESET}", end="")
            else:
                print(f"  {'0':<12} {'0.0':>8} {'None':<8}", end="")
        print()

    print("=" * 80)

    # report paths
    print("\nReports saved to data/results/:")
    for _, r in results:
        if r.report_path:
            print(f"  {r.sample_name}_risk_dashboard.html")
            print(f"  {r.sample_name}_profile.html")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("\n=== 3-Way Pipeline Performance Comparison ===\n")
    settings.RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("[Setup] Initialising reference databases...")
    vfdb_blast, card_blast = setup_databases()
    print("  Databases ready.\n")

    runner  = PipelineRunner()
    results = []

    for test in TESTS:
        fasta: Path = test["fasta"]
        if not fasta.exists():
            print(f"[SKIP] {test['label']} -- file not found: {fasta.name}")
            continue

        print(f"[Running] {test['label']}  ({fasta.name})")
        result = runner.run(
            sample_name=test["name"],
            fasta_path=fasta,
            description=test["label"],
            vfdb_blast_db=vfdb_blast,
            card_blast_db=card_blast,
            alignment_tool="blastn",
            word_size=test["word_size"],
        )
        results.append((test["label"].split("--")[1].strip(), result))
        print(f"  VF hits: {result.vf_hits}  AMR hits: {result.amr_hits}\n")

    if not results:
        print("No test inputs found. Check data/samples/ directory.")
        return

    print_comparison(results)


if __name__ == "__main__":
    main()