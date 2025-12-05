#!/usr/bin/env python3
import os
import argparse
import pickle
from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
import matplotlib.pyplot as plt

def main(fragments_dict_path, qc_output_dir, plots_output_dir, qc_results_pickle=None, barcode_plots_output_dir=None):
    """
    Generate QC plots for each sample and save as PNG.
    
    fragments_dict_path: pickle file containing sample_id -> fragment path
    qc_output_dir: directory where pycisTopic QC outputs are located
    plots_output_dir: directory to save sample-level QC plots (PNG)
    qc_results_pickle: optional pickle file with barcodes + thresholds for barcode plots
    barcode_plots_output_dir: optional directory to save barcode-level QC plots (PNG)
    """

    # Load fragments_dict
    if not os.path.exists(fragments_dict_path):
        raise FileNotFoundError(f"Fragments dictionary not found: {fragments_dict_path}")

    with open(fragments_dict_path, "rb") as f:
        fragments_dict = pickle.load(f)

    # Create plots output directory
    os.makedirs(plots_output_dir, exist_ok=True)
    
    # Create barcode plots output directory if needed
    if barcode_plots_output_dir:
        os.makedirs(barcode_plots_output_dir, exist_ok=True)

    # Load QC results for barcode plots if provided
    if qc_results_pickle:
        if not os.path.exists(qc_results_pickle):
            raise FileNotFoundError(f"QC results pickle not found: {qc_results_pickle}")
        with open(qc_results_pickle, "rb") as f:
            qc_results = pickle.load(f)
        sample_id_to_barcodes_passing_filters = qc_results["barcodes"]
        sample_id_to_thresholds = qc_results["thresholds"]

    # Generate plots for each sample
    for sample_id in fragments_dict:
        print(f"Processing sample: {sample_id}")
        
        # Generate sample-level QC plots
        print(f"  Generating sample-level QC plot...")
        fig = plot_sample_stats(
            sample_id=sample_id,
            pycistopic_qc_output_dir=qc_output_dir
        )
        plot_file = os.path.join(plots_output_dir, f"{sample_id}_qc.png")
        fig.savefig(plot_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved sample-level plot: {plot_file}")

        # Generate barcode-level QC plots if QC results are available
        if qc_results_pickle and barcode_plots_output_dir:
            print(f"  Generating barcode-level QC plot...")
            fig_barcode = plot_barcode_stats(
                sample_id=sample_id,
                pycistopic_qc_output_dir=qc_output_dir,
                bc_passing_filters=sample_id_to_barcodes_passing_filters[sample_id],
                detailed_title=False,
                **sample_id_to_thresholds[sample_id]
            )
            barcode_plot_file = os.path.join(barcode_plots_output_dir, f"{sample_id}_barcode_qc.png")
            fig_barcode.savefig(barcode_plot_file, dpi=300, bbox_inches="tight")
            plt.close(fig_barcode)
            print(f"  Saved barcode-level plot: {barcode_plot_file}")

    print("All plots generated successfully!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate unified pycisTopic QC plots")
    parser.add_argument("--fragments_dict", required=True,
                        help="Pickle file containing sample_id -> fragment file path")
    parser.add_argument("--qc_output_dir", required=True,
                        help="Directory containing pycisTopic QC output for all samples")
    parser.add_argument("--plots_output_dir", required=True,
                        help="Directory to save sample-level QC plots as PNG")
    parser.add_argument("--qc_results_pickle", required=False, default=None,
                        help="Optional: Pickle file with barcodes + thresholds for barcode plots")
    parser.add_argument("--barcode_plots_output_dir", required=False, default=None,
                        help="Optional: Directory to save barcode-level QC plots as PNG")
    args = parser.parse_args()

    main(
        args.fragments_dict,
        args.qc_output_dir,
        args.plots_output_dir,
        args.qc_results_pickle,
        args.barcode_plots_output_dir
    )
