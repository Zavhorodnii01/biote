"""
Nanopore Read Simulator
=======================
Generates simulated Nanopore reads from a reference genome using badread.
Produces Test 2 input for the performance testing framework.

Usage:
    pip install badread
    python simulate_reads.py

Input:
    data/samples/ecoli_O157H7_reference.fasta

Output:
    data/samples/ecoli_O157H7_simulated.fasta
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import settings

REFERENCE = settings.SAMPLES_DIR / "ecoli_O157H7_reference.fasta"
SIMULATED = settings.SAMPLES_DIR / "ecoli_O157H7_simulated.fasta"

# Nanopore R9.4 typical profile
QUANTITY    = "50x"          # 50x coverage of the genome
READ_LENGTH = "8000,6000"    # mean ~8 kb, stdev 6 kb
ERROR_MODEL = "nanopore2023" # built-in Nanopore error model
IDENTITY    = "85,95"        # ~5-15% error rate, realistic for Nanopore R9.4


def check_badread() -> None:
    try:
        subprocess.run(
            [sys.executable, "-m", "badread", "--version"],
            capture_output=True, check=True,
        )
        print("  badread: OK")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERROR: badread not found.\nInstall with:  pip install badread")
        sys.exit(1)


def simulate(reference: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "-m", "badread", "simulate",
        "--reference",   str(reference),
        "--quantity",    QUANTITY,
        "--length",      READ_LENGTH,
        "--error_model", ERROR_MODEL,
        "--identity",    IDENTITY,
        "--seed",        "42",
    ]
    print(f"  Running badread simulate ...")
    with open(output, "w") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"ERROR:\n{result.stderr}")
        sys.exit(1)

    n_reads = sum(1 for line in open(output) if line.startswith(">"))
    size_mb = output.stat().st_size / 1_000_000
    print(f"  Done: {n_reads:,} simulated reads, {size_mb:.1f} MB -> {output.name}")


def main():
    print("\n=== Nanopore Read Simulator (badread) ===\n")
    if not REFERENCE.exists():
        print(f"ERROR: Reference not found at:\n  {REFERENCE}")
        print("Download E. coli O157:H7 from NCBI (Genome sequences FASTA) and place it there.")
        sys.exit(1)

    ref_size = REFERENCE.stat().st_size / 1_000_000
    print(f"  Reference  : {REFERENCE.name}  ({ref_size:.1f} MB)")
    print(f"  Coverage   : {QUANTITY}")
    print(f"  Read length: {READ_LENGTH} bp (mean, stdev)")
    print(f"  Error model: {ERROR_MODEL}  identity range={IDENTITY}\n")

    check_badread()
    simulate(REFERENCE, SIMULATED)

    print(f"\nNext: run  python run_comparison.py  to compare all 3 test inputs")


if __name__ == "__main__":
    main()