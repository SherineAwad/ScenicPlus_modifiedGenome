#!/usr/bin/env python3
import os
import argparse
import pickle

def main():
    parser = argparse.ArgumentParser(description="Generate pycisTopic QC command file for multiple samples")
    parser.add_argument("--out_dir", required=True, help="Base output directory")
    parser.add_argument("--consensus_dir", required=True, help="Directory containing consensus peak calling files")
    parser.add_argument("--tss_bed", required=True, help="Path to TSS BED file (e.g., tss_mm10.bed)")
    parser.add_argument("--th1_fragments", required=True, help="Path to TH1 fragment file (TSV.gz)")
    parser.add_argument("--th2_fragments", required=True, help="Path to TH2 fragment file (TSV.gz)")
    parser.add_argument("--qc_commands_filename", default="pycistopic_qc_commands.txt",
                        help="Filename to save QC commands")
    args = parser.parse_args()

    # 1. Create fragments_dict
    fragments_dict = {
        "TH1": args.th1_fragments,
        "TH2": args.th2_fragments
    }

    # Save fragments_dict as pickle (optional, for later use)
    fragments_dict_path = os.path.join(args.out_dir, "fragments_dict.pkl")
    os.makedirs(args.out_dir, exist_ok=True)
    with open(fragments_dict_path, "wb") as f:
        pickle.dump(fragments_dict, f)
    print(f"Saved fragments_dict to {fragments_dict_path}")

    # 2. Define regions BED from consensus directory
    consensus_regions_dir = args.consensus_dir
    if not os.path.exists(consensus_regions_dir):
        raise FileNotFoundError(f"Consensus regions directory not found: {consensus_regions_dir}")
    
    # Assuming consensus BED file is named something like consensus_regions.bed inside the folder
    # If multiple files exist, pick the first BED file
    consensus_bed_files = [f for f in os.listdir(consensus_regions_dir) if f.endswith(".bed")]
    if len(consensus_bed_files) == 0:
        raise FileNotFoundError(f"No BED files found in {consensus_regions_dir}")
    
    regions_bed_filename = os.path.join(consensus_regions_dir, consensus_bed_files[0])
    print(f"Using consensus BED file: {regions_bed_filename}")

    # 3. TSS BED file
    tss_bed_filename = args.tss_bed
    if not os.path.exists(tss_bed_filename):
        raise FileNotFoundError(f"TSS BED file not found: {tss_bed_filename}")

    # Make sure QC output directory exists
    qc_output_dir = os.path.join(args.out_dir, "QC")
    os.makedirs(qc_output_dir, exist_ok=True)

    # 4. Create text file with pycisTopic QC command lines
    with open(args.qc_commands_filename, "w") as fh:
        for sample, fragment_filename in fragments_dict.items():
            sample_qc_dir = os.path.join(qc_output_dir, sample)
            os.makedirs(sample_qc_dir, exist_ok=True)
            command = (
                f"pycistopic qc "
                f"--fragments {fragment_filename} "
                f"--regions {regions_bed_filename} "
                f"--tss {tss_bed_filename} "
                f"--output {sample_qc_dir}"
            )
            print(command, file=fh)

    print(f"QC commands written to {args.qc_commands_filename}")
    print("You can now run them in the terminal, e.g.:")
    print(f"  bash {args.qc_commands_filename}")

if __name__ == "__main__":
    main()

