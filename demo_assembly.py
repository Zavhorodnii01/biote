"""
Demo script showing how to use the assembly feature in the pipeline.

Usage:
    python demo_assembly.py <fasta_file> [--assembly]

Examples:
    # Run without assembly (default)
    python demo_assembly.py data/samples/ecoli_O157H7_reference.fasta

    # Run with assembly
    python demo_assembly.py data/samples/ecoli_O157H7_reference.fasta --assembly
"""

import argparse
import sys
from pathlib import Path

from nanopore_pipeline.pipeline_runner import PipelineRunner
from config import settings


def main():
    parser = argparse.ArgumentParser(
        description="Run nanopore pipeline with optional assembly step"
    )
    parser.add_argument(
        "fasta_file",
        type=Path,
        help="Path to input FASTA file",
    )
    parser.add_argument(
        "--assembly",
        action="store_true",
        help="Perform Flye assembly before alignment (reduces read overlap inflation)",
    )
    parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name (defaults to filename)",
    )

    args = parser.parse_args()

    # Validate input file
    if not args.fasta_file.exists():
        print(f"Error: FASTA file not found: {args.fasta_file}")
        sys.exit(1)

    # Generate sample name if not provided
    sample_name = args.sample_name or args.fasta_file.stem
    if args.assembly:
        sample_name += "_assembled"

    print("=" * 70)
    print(f"Nanopore Pathogenicity Pipeline")
    print("=" * 70)
    print(f"Sample name:     {sample_name}")
    print(f"Input FASTA:     {args.fasta_file}")
    print(f"Assembly:        {'Yes (Flye)' if args.assembly else 'No (raw reads)'}")
    print("=" * 70)
    print()

    # Initialize pipeline
    runner = PipelineRunner()

    # Run pipeline
    try:
        result = runner.run(
            sample_name=sample_name,
            fasta_path=args.fasta_file,
            perform_assembly=args.assembly,
        )

        # Print results
        print("\n" + "=" * 70)
        print("RESULTS")
        print("=" * 70)
        print(f"Sample ID:       {result.sample_id}")
        print()
        print("FASTA Statistics:")
        print(f"  Total reads:   {result.fasta_stats.total_reads:,}")
        print(f"  Total bases:   {result.fasta_stats.total_bases:,}")
        print(f"  Mean length:   {result.fasta_stats.mean_read_length:.0f} bp")
        print(f"  N50:           {result.fasta_stats.n50:,} bp")
        print(f"  GC content:    {result.fasta_stats.mean_gc_content:.1f}%")

        if result.assembly_result and result.assembly_result.success:
            print()
            print("Assembly Statistics:")
            print(f"  Contigs:       {result.assembly_result.contigs:,}")
            print(f"  N50:           {result.assembly_result.n50:,} bp")
            print(f"  Total bases:   {result.assembly_result.total_bases:,}")
            print(f"  Mean length:   {result.assembly_result.mean_contig_length:.0f} bp")
            print(f"  Output:        {result.assembly_result.contig_path}")
            print()
            print(f"  Read reduction: {result.fasta_stats.total_reads:,} → {result.assembly_result.contigs:,} "
                  f"({100 * (1 - result.assembly_result.contigs / result.fasta_stats.total_reads):.1f}% reduction)")
        elif result.assembly_result and not result.assembly_result.success:
            print()
            print(f"Assembly failed: {result.assembly_result.error_message}")
            print("Pipeline continued with raw reads.")

        print()
        print("Alignment Results:")
        print(f"  VF hits:       {result.vf_hits:,}")
        print(f"  AMR hits:      {result.amr_hits:,}")

        print()
        print(f"Categories found: {len(result.profiles)}")
        for profile in result.profiles[:10]:  # Show top 10
            print(f"  {profile.category:30s}  "
                  f"genes={profile.unique_genes:3d}  "
                  f"TPM={profile.tpm:8.1f}  "
                  f"risk={profile.risk_level}")

        if len(result.profiles) > 10:
            print(f"  ... and {len(result.profiles) - 10} more categories")

        if result.report_path:
            print()
            print(f"HTML report:     {result.report_path}")
        if result.json_path:
            print(f"JSON export:     {result.json_path}")

        print("=" * 70)

    except Exception as e:
        print(f"\nError running pipeline: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
