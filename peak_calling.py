#!/usr/bin/env python3
import os
import subprocess
import glob
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Consensus Peak Calling Script using MACS2"
    )
    parser.add_argument(
        "-i", "--input", required=True, 
        help="Input directory containing .fragments.tsv.gz files"
    )
    parser.add_argument(
        "-o", "--output", required=True, 
        help="Output directory for BED and MACS2 peaks"
    )
    parser.add_argument(
        "-g", "--genome", default="mm", 
        help="Genome size for MACS2 (default: mm)"
    )

    args = parser.parse_args()

    frag_dir = args.input
    out_dir = args.output
    genome = args.genome

    print("=== STARTING PEAK CALLING SCRIPT ===", flush=True)
    print(f"Creating output directory: {out_dir}", flush=True)
    os.makedirs(out_dir, exist_ok=True)

    # Find all fragment files
    print(f"Looking for fragment files in: {frag_dir}", flush=True)
    frag_files = glob.glob(os.path.join(frag_dir, "*.fragments.tsv.gz"))
    print(f"Found {len(frag_files)} fragment files", flush=True)

    for i, frag in enumerate(frag_files):
        sample = os.path.basename(frag).replace(".fragments.tsv.gz", "")
        print(f">>> [{i+1}/{len(frag_files)}] Processing {sample}", flush=True)

        bed_file = os.path.join(out_dir, f"{sample}.bed")

        # Convert fragments to BED WITHOUT SORTING (MACS2 can handle unsorted)
        print("  Step 1: Converting fragments to BED (no sort)...", flush=True)
        cmd_bed = f"zcat {frag} | awk 'BEGIN{{OFS=\"\\t\"}} {{print $1,$2,$3}}' > {bed_file}"
        print(f"  Running: {cmd_bed}", flush=True)

        try:
            subprocess.run(cmd_bed, shell=True, check=True, timeout=3600)  # 1 hour timeout
            file_size = os.path.getsize(bed_file) / (1024 * 1024)  # Size in MB
            print(f"  BED created: {file_size:.2f} MB", flush=True)
        except subprocess.TimeoutExpired:
            print("  TIMEOUT - skipping to next sample", flush=True)
            continue
        except Exception as e:
            print(f"  Conversion failed: {e}", flush=True)
            continue

        # Call MACS2
        print("  Step 2: Calling peaks with MACS2...", flush=True)
        cmd_macs = [
            "macs2", "callpeak",
            "-t", bed_file,
            "-f", "BED",
            "-n", sample,
            "--outdir", out_dir,
            "--nomodel",
            "--shift", "-100",
            "--extsize", "200",
            "-g", genome,
            "--keep-dup", "all"  # Important for ATAC-seq
        ]
        print(f"  Running: {' '.join(cmd_macs)}", flush=True)

        try:
            subprocess.run(cmd_macs, check=True, timeout=7200)  # 2 hour timeout
            print(f"  MACS2 completed for {sample}", flush=True)
        except subprocess.TimeoutExpired:
            print("  MACS2 TIMEOUT", flush=True)
        except Exception as e:
            print(f"  MACS2 failed: {e}", flush=True)

        print(f">>> Done {sample}", flush=True)
        print("-" * 50, flush=True)

    print("=== ALL SAMPLES PROCESSED ===", flush=True)


if __name__ == "__main__":
    main()

