#!/usr/bin/env python3
import os
import argparse
import pickle
from pycisTopic.qc import get_barcodes_passing_qc_for_sample

def main(fragments_dict_path, qc_output_dir, output_pickle,
         unique_fragments_threshold, tss_enrichment_threshold, frip_threshold):
    """
    Collect barcodes passing QC for each sample.
    
    fragments_dict_path: pickle file containing sample_id -> fragment path
    qc_output_dir: directory where pycisTopic QC outputs are located
    output_pickle: path to save results (barcodes + thresholds)
    unique_fragments_threshold: minimum unique fragments to pass QC
    tss_enrichment_threshold: minimum TSS enrichment to pass QC
    frip_threshold: minimum fraction of reads in peaks to pass QC
    """

    # Load fragments_dict
    if not os.path.exists(fragments_dict_path):
        raise FileNotFoundError(f"Fragments dictionary not found: {fragments_dict_path}")

    with open(fragments_dict_path, "rb") as f:
        fragments_dict = pickle.load(f)

    sample_id_to_barcodes_passing_filters = {}
    sample_id_to_thresholds = {}

    # Loop over each sample
    for sample_id in fragments_dict:
        print(f"Processing sample: {sample_id}")
        barcodes, thresholds = get_barcodes_passing_qc_for_sample(
            sample_id=sample_id,
            pycistopic_qc_output_dir=qc_output_dir,
            unique_fragments_threshold=unique_fragments_threshold,
            tss_enrichment_threshold=tss_enrichment_threshold,
            frip_threshold=frip_threshold,
            use_automatic_thresholds=True
        )
        sample_id_to_barcodes_passing_filters[sample_id] = barcodes
        sample_id_to_thresholds[sample_id] = thresholds

    # Save results
    results = {
        "barcodes": sample_id_to_barcodes_passing_filters,
        "thresholds": sample_id_to_thresholds
    }
    with open(output_pickle, "wb") as f:
        pickle.dump(results, f)

    print(f"QC barcodes and thresholds saved to {output_pickle}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect barcodes passing pycisTopic QC")
    parser.add_argument("--fragments_dict", required=True,
                        help="Pickle file containing sample_id -> fragment file path")
    parser.add_argument("--qc_output_dir", required=True,
                        help="Directory containing pycisTopic QC output for all samples")
    parser.add_argument("--output_pickle", required=True,
                        help="Path to save barcodes and thresholds pickle")
    parser.add_argument("--unique_fragments_threshold", type=int, default=500,
                        help="Minimum unique fragments (default: 500)")
    parser.add_argument("--tss_enrichment_threshold", type=float, default=None,
                        help="Minimum TSS enrichment (default: None, automatic)")
    parser.add_argument("--frip_threshold", type=float, default=0,
                        help="Minimum fraction of reads in peaks (default: 0)")
    args = parser.parse_args()

    main(args.fragments_dict, args.qc_output_dir, args.output_pickle,
         args.unique_fragments_threshold, args.tss_enrichment_threshold, args.frip_threshold)

