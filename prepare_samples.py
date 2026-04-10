"""
Sample Preparation
==================
Prepares properly-sized FASTA subsets for Test 2 and Test 3.

Test 2 (simulated): filters to reads <= 30 kb, keeps up to 5000 reads
Test 3 (real SRA):  keeps first 5000 reads

Run this once before run_comparison.py.

Usage:
    python prepare_samples.py
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import settings

MAX_READS      = 5000
MAX_READ_LEN   = 30_000   # filter out badread giant reads


def subsample_fasta(
    input_path: Path,
    output_path: Path,
    max_reads: int = MAX_READS,
    max_len: int = 0,          # 0 = no length filter
) -> int:
    """Write up to max_reads records from input_path to output_path.
    If max_len > 0, skips reads longer than max_len bp.
    Returns number of reads written.
    """
    count = 0
    current_id = None
    current_seq: list[str] = []

    def write_record(fh, rec_id, seq_lines):
        seq = "".join(seq_lines)
        if max_len > 0 and len(seq) > max_len:
            return False
        fh.write(rec_id + "\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")
        return True

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(input_path) as fh_in, open(output_path, "w") as fh_out:
        for line in fh_in:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None and count < max_reads:
                    if write_record(fh_out, current_id, current_seq):
                        count += 1
                if count >= max_reads:
                    break
                current_id = line
                current_seq = []
            else:
                current_seq.append(line)

        # last record
        if current_id is not None and count < max_reads:
            if write_record(fh_out, current_id, current_seq):
                count += 1

    return count


def main():
    print("\n=== Sample Preparation ===\n")

    # Test 2 — simulated reads (filter giant reads, keep 5000)
    sim_in  = settings.SAMPLES_DIR / "ecoli_O157H7_simulated.fasta"
    sim_out = settings.SAMPLES_DIR / "ecoli_O157H7_simulated_5k.fasta"

    if not sim_in.exists():
        print(f"  [SKIP] Simulated FASTA not found: {sim_in.name}")
    else:
        print(f"  Preparing Test 2: {sim_in.name} -> {sim_out.name}")
        print(f"    (max {MAX_READS} reads, max length {MAX_READ_LEN:,} bp)")
        n = subsample_fasta(sim_in, sim_out, max_reads=MAX_READS, max_len=MAX_READ_LEN)
        size_mb = sim_out.stat().st_size / 1_000_000
        print(f"    Done: {n} reads, {size_mb:.1f} MB")

    # Test 3 — real SRA (keep first 5000 reads, no length filter)
    sra_in  = settings.SAMPLES_DIR / "ecoli_O157H7_real_sra.fasta"
    sra_out = settings.SAMPLES_DIR / "ecoli_O157H7_real_sra_5k.fasta"

    if not sra_in.exists():
        print(f"  [SKIP] SRA FASTA not found: {sra_in.name}")
    else:
        print(f"\n  Preparing Test 3: {sra_in.name} -> {sra_out.name}")
        print(f"    (first {MAX_READS} reads)")
        n = subsample_fasta(sra_in, sra_out, max_reads=MAX_READS, max_len=0)
        size_mb = sra_out.stat().st_size / 1_000_000
        print(f"    Done: {n} reads, {size_mb:.1f} MB")

    print("\nDone. Now run:  python run_comparison.py")


if __name__ == "__main__":
    main()
