#!/usr/bin/env python3
import os
import glob
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Merge MACS2 narrowPeak files and generate consensus peaks"
    )

    parser.add_argument(
        "-i", "--macs2_dir", required=True,
        help="Directory containing MACS2 narrowPeak files"
    )
    parser.add_argument(
        "-o", "--output_dir", required=True,
        help="Directory to store combined and consensus peak files"
    )
    parser.add_argument(
        "-c", "--combined_bed", default="combined_peaks.bed",
        help="Filename for combined narrowPeak file (default: combined_peaks.bed)"
    )
    parser.add_argument(
        "-p", "--consensus_peaks", default="consensus_peaks.bed",
        help="Filename for consensus peaks file (default: consensus_peaks.bed)"
    )

    args = parser.parse_args()

    macs2_dir = args.macs2_dir
    output_dir = args.output_dir
    combined_bed = os.path.join(output_dir, args.combined_bed)
    consensus_bed = os.path.join(output_dir, args.consensus_peaks)

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # -----------------------------
    # Merge all narrowPeak files
    # -----------------------------
    narrow_files = glob.glob(os.path.join(macs2_dir, "*.narrowPeak"))
    if not narrow_files:
        print(f"No narrowPeak files found in {macs2_dir}")
        return

    # Combine all narrowPeak files into one
    with open(combined_bed, "w") as out_f:
        for f in narrow_files:
            with open(f) as in_f:
                for line in in_f:
                    out_f.write(line)
    print(f"Combined narrowPeak file created: {combined_bed}")

    # Sort and merge peaks using bedtools
    try:
        subprocess.run(
            f"bedtools sort -i {combined_bed} | bedtools merge > {consensus_bed}",
            shell=True,
            check=True
        )
        print(f"Consensus peak set created: {consensus_bed}")
    except subprocess.CalledProcessError as e:
        print(f"Error during bedtools processing: {e}")

if __name__ == "__main__":
    main()

