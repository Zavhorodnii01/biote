"""
6-Way Pipeline Performance Benchmark
======================================
Runs the full pipeline on all six test inputs and produces:
  - Side-by-side comparison table (console)
  - Ground truth recall for Test 4 (MRSA RefSeq with GFF annotation)
  - Per-category error analysis (where does the pipeline fail?)
  - Interactive benchmark heatmap (data/results/benchmark_heatmap.html)

Test 1 -- Perfect reference  : ecoli_O157H7_reference.fasta     (word_size=11)
Test 2 -- Simulated Nanopore : ecoli_O157H7_simulated_5k.fasta  (word_size=7)
Test 3 -- Real Nanopore SRA  : ecoli_O157H7_real_sra_5k.fasta   (SRA DRR198813)
Test 4 -- MRSA RefSeq genome : mrsa_reference.fasta + .gff      (ground truth)
Test 5 -- MRSA Nanopore SRA  : mrsa_nanopore_sra.fasta          (single strain)
Test 6 -- Metagenomic Nanopore: metagenomic_nanopore_sra.fasta  (multi-species)

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
from nanopore_pipeline.utils.gff_parser import GFFGroundTruth
from nanopore_pipeline.classifier.ground_truth_evaluator import GroundTruthEvaluator
from nanopore_pipeline.reporter.visualisation import PathogenicityReporter

# ── Test definitions ──────────────────────────────────────────────────────────
TESTS = [
    {
        "name":         "test1_perfect_reference",
        "label":        "T1 Perfect Ref",
        "description":  "Test 1 — Perfect Reference (E. coli O157:H7)",
        "fasta":        settings.SAMPLES_DIR / "ecoli_O157H7_reference.fasta",
        "word_size":    11,    # standard BLAST — perfect DNA, no Nanopore errors
        "gff":          None,
        "ground_truth": False,
    },
    {
        "name":         "test2_simulated_nanopore",
        "label":        "T2 Sim Nanopore",
        "description":  "Test 2 — Simulated Nanopore (E. coli O157:H7, ~10% error)",
        "fasta":        settings.SAMPLES_DIR / "ecoli_O157H7_simulated.fasta",
        "word_size":    7,
        "gff":          None,
        "ground_truth": False,
    },
    {
        "name":         "test3_real_sra",
        "label":        "T3 Real SRA",
        "description":  "Test 3 — Real Nanopore SRA (E. coli O157:H7, DRR198813)",
        "fasta":        settings.SAMPLES_DIR / "ecoli_O157H7_real_sra.fasta",
        "word_size":    7,
        "gff":          None,
        "ground_truth": False,
    },
    {
        "name":         "test4_mrsa_reference",
        "label":        "T4 MRSA RefSeq",
        "description":  "Test 4 — MRSA RefSeq genome (S. aureus MRSA252, NC_002952) + GFF ground truth",
        "fasta":        settings.SAMPLES_DIR / "mrsa_reference.fasta",
        "word_size":    11,    # already assembled reference — no Nanopore errors
        "gff":          settings.SAMPLES_DIR / "mrsa_reference.gff",
        "ground_truth": True,
    },
    {
        "name":         "test5_mrsa_nanopore",
        "label":        "T5 MRSA Nanopore",
        "description":  "Test 5 — Real Nanopore SRA, single-strain MRSA",
        "fasta":        settings.SAMPLES_DIR / "mrsa_nanopore_sra.fasta",
        "word_size":    7,
        "gff":          None,
        "ground_truth": False,
    },
    {
        "name":         "test6_metagenomic",
        "label":        "T6 Metagenomic",
        "description":  "Test 6 — Real Nanopore SRA, environmental/hospital metagenome",
        "fasta":        settings.SAMPLES_DIR / "metagenomic_nanopore_sra.fasta",
        "word_size":    7,
        "gff":          None,
        "ground_truth": False,
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


# ── Ground truth evaluation (Test 4) ─────────────────────────────────────────

def run_ground_truth_eval(test: dict, result: PipelineResult) -> dict | None:
    """Run GFF-based ground truth evaluation for tests that have a GFF file."""
    if not test.get("ground_truth") or not test.get("gff"):
        return None
    gff_path = Path(test["gff"])
    if not gff_path.exists():
        print(f"  [WARN] GFF not found, skipping ground truth: {gff_path.name}")
        return None

    print(f"  Running ground truth evaluation from {gff_path.name}...")
    gt_builder = GFFGroundTruth()
    ground_truth = gt_builder.build(gff_path, settings.VF_CATEGORIES, settings.AMR_CATEGORIES)

    evaluator = GroundTruthEvaluator()
    report = evaluator.evaluate(ground_truth, result.profiles, result.sample_name)

    print(report.summary_table())

    # Return recall dict {category: recall_float} for heatmap
    return {c.category: c.recall for c in report.categories}


# ── Console comparison table ──────────────────────────────────────────────────

def print_comparison(results: list[tuple[str, PipelineResult, dict | None]]) -> None:
    RISK_COLORS = {
        "Critical": "\033[91m", "High": "\033[93m",
        "Medium":   "\033[94m", "Low":  "\033[92m", "None": "",
    }
    RESET = "\033[0m"

    print("\n" + "=" * 90)
    print("6-WAY PIPELINE BENCHMARK")
    print("=" * 90)

    print(f"\n{'Metric':<22}", end="")
    for label, _, _ in results:
        print(f"  {label:<16}", end="")
    print()
    print("-" * (22 + 18 * len(results)))

    def row(metric, getter):
        print(f"{metric:<22}", end="")
        for _, r, _ in results:
            try:
                val = getter(r)
            except Exception:
                val = "N/A"
            print(f"  {str(val):<16}", end="")
        print()

    row("Total reads",      lambda r: f"{r.fasta_stats.total_reads:,}")
    row("Total bases (Mb)", lambda r: f"{r.fasta_stats.total_bases / 1e6:.1f}")
    row("N50",              lambda r: f"{r.fasta_stats.n50:,} bp")
    row("GC content",       lambda r: f"{r.fasta_stats.mean_gc_content:.1f}%")
    row("VF hits",          lambda r: str(r.vf_hits))
    row("AMR hits",         lambda r: str(r.amr_hits))

    # Category table
    all_cats = sorted({p.category for _, r, _ in results for p in r.profiles})
    if not all_cats:
        print("\n  No pathogenicity profiles detected in any test.")
        return

    print(f"\n{'Category':<26} {'Diff':>4}", end="")
    for label, _, _ in results:
        print(f"  {label[:8]:<8} {'Hits':>4} {'Risk':<6}", end="")
    print()
    print("-" * (32 + 20 * len(results)))

    for cat in all_cats:
        diff = settings.CATEGORY_DIFFICULTY.get(cat, "?")[0]   # E/M/H
        diff_colors = {"E": "\033[92m", "M": "\033[93m", "H": "\033[91m"}
        diff_color = diff_colors.get(diff, "")
        print(f"{cat:<26} {diff_color}{diff:>4}{RESET}", end="")
        for _, r, recall_dict in results:
            profile = next((p for p in r.profiles if p.category == cat), None)
            if profile:
                color = RISK_COLORS.get(profile.risk_level, "")
                # Show recall % if available
                if recall_dict and cat in recall_dict:
                    hits_str = f"{recall_dict[cat] * 100:.0f}%"
                else:
                    hits_str = str(profile.hit_count)
                print(f"  {hits_str:<8} {profile.hit_count:>4} "
                      f"{color}{profile.risk_level[:5]:<6}{RESET}", end="")
            else:
                print(f"  {'—':<8} {'0':>4} {'None':<6}", end="")
        print()

    print("=" * 90)
    print("\nNotes:")
    print("  Diff column: E=Easy  M=Medium  H=Hard (bioinformatic detectability)")
    print("  T4 hits column shows ground-truth recall % when GFF is available")


# ── Benchmark heatmap ─────────────────────────────────────────────────────────

def generate_heatmap(results: list[tuple[str, PipelineResult, dict | None]]) -> None:
    reporter = PathogenicityReporter()
    test_data = []
    for label, r, recall_dict in results:
        test_data.append({
            "label":      label,
            "profiles":   r.profiles,
            "recall":     recall_dict,
            "difficulty": settings.CATEGORY_DIFFICULTY,
        })
    fig = reporter.plot_benchmark_heatmap(test_data)
    out = settings.RESULTS_DIR / "benchmark_heatmap.html"
    print(f"\nBenchmark heatmap saved -> {out}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("\n=== 6-Way Pipeline Performance Benchmark ===\n")
    settings.RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("[Setup] Initialising reference databases...")
    vfdb_blast, card_blast = setup_databases()
    print("  Databases ready.\n")

    runner  = PipelineRunner()
    results = []   # (label, PipelineResult, recall_dict|None)

    for test in TESTS:
        fasta: Path = test["fasta"]
        if not fasta.exists():
            print(f"[SKIP] {test['label']} — file not found: {fasta.name}")
            continue

        print(f"[Running] {test['description']}")
        result = runner.run(
            sample_name=test["name"],
            fasta_path=fasta,
            description=test["description"],
            vfdb_blast_db=vfdb_blast,
            card_blast_db=card_blast,
            alignment_tool="blastn",
            word_size=test["word_size"],
        )
        print(f"  VF hits: {result.vf_hits}  AMR hits: {result.amr_hits}")

        recall_dict = run_ground_truth_eval(test, result)

        results.append((test["label"], result, recall_dict))
        print()

    if not results:
        print("No test inputs found. Check data/samples/ directory.")
        return

    print_comparison(results)
    generate_heatmap(results)

    print("\nReports saved to data/results/:")
    for label, r, _ in results:
        if r.report_path:
            print(f"  [{label}] {r.sample_name}_risk_dashboard.html")
    print("  benchmark_heatmap.html")


if __name__ == "__main__":
    main()
