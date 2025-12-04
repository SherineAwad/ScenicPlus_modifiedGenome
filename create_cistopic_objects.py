#!/usr/bin/env python3
import os
import argparse
import pickle
import polars as pl
import scipy.sparse as sp
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments

def fix_transposed_matrices(cistopic_obj):
    """
    Fix transposed matrices in cisTopic object
    """
    print(f"  Checking matrix dimensions...")
    
    # Check and fix fragment_matrix
    if hasattr(cistopic_obj, 'fragment_matrix'):
        if cistopic_obj.fragment_matrix.shape[0] == len(cistopic_obj.region_names):
            print(f"    Fixing fragment_matrix: {cistopic_obj.fragment_matrix.shape} -> {cistopic_obj.fragment_matrix.shape[1], cistopic_obj.fragment_matrix.shape[0]}")
            cistopic_obj.fragment_matrix = cistopic_obj.fragment_matrix.T
    
    # Check and fix binary_matrix  
    if hasattr(cistopic_obj, 'binary_matrix'):
        if cistopic_obj.binary_matrix.shape[0] == len(cistopic_obj.region_names):
            print(f"    Fixing binary_matrix: {cistopic_obj.binary_matrix.shape} -> {cistopic_obj.binary_matrix.shape[1], cistopic_obj.binary_matrix.shape[0]}")
            cistopic_obj.binary_matrix = cistopic_obj.binary_matrix.T
    
    return cistopic_obj

def main(fragments_dict_path, qc_results_pickle, regions_bed, blacklist_bed,
         qc_output_dir, output_pickle, n_cpu):
    """
    Create CistopicObjects from fragments after QC filtering.

    fragments_dict_path: pickle file mapping sample_id -> fragment path
    qc_results_pickle: pickle file with barcodes + thresholds (from QC filtering)
    regions_bed: consensus regions BED file
    blacklist_bed: blacklist BED file (mm10 in this case)
    qc_output_dir: directory containing pycisTopic QC outputs (with per-cell stats parquet)
    output_pickle: output pickle to save list of CistopicObjects
    n_cpu: number of CPUs to use
    """

    # Load fragments_dict
    if not os.path.exists(fragments_dict_path):
        raise FileNotFoundError(f"Fragments dictionary not found: {fragments_dict_path}")
    with open(fragments_dict_path, "rb") as f:
        fragments_dict = pickle.load(f)

    # Load QC results
    if not os.path.exists(qc_results_pickle):
        raise FileNotFoundError(f"QC results pickle not found: {qc_results_pickle}")
    with open(qc_results_pickle, "rb") as f:
        qc_results = pickle.load(f)
    sample_id_to_barcodes_passing_filters = qc_results["barcodes"]

    # Create list of CistopicObjects
    cistopic_obj_list = []
    for sample_id in fragments_dict:
        print(f"Creating CistopicObject for sample: {sample_id}")

        # Load per-cell metrics parquet
        parquet_path = os.path.join(qc_output_dir, f"{sample_id}.fragments_stats_per_cb.parquet")
        if not os.path.exists(parquet_path):
            raise FileNotFoundError(f"QC parquet not found for sample {sample_id}: {parquet_path}")

        sample_metrics = (
            pl.read_parquet(parquet_path)
              .to_pandas()
              .set_index("CB")
              .loc[sample_id_to_barcodes_passing_filters[sample_id]]
        )

        # Create CistopicObject
        cistopic_obj = create_cistopic_object_from_fragments(
            path_to_fragments=fragments_dict[sample_id],
            path_to_regions=regions_bed,
            path_to_blacklist=blacklist_bed,
            metrics=sample_metrics,
            valid_bc=sample_id_to_barcodes_passing_filters[sample_id],
            n_cpu=n_cpu,
            project=sample_id,
            split_pattern="-"
        )

        # FIX: Transpose matrices if they're in wrong orientation
        print(f"  Checking and fixing matrix orientation...")
        cistopic_obj = fix_transposed_matrices(cistopic_obj)
        
        # Verify dimensions are correct
        if hasattr(cistopic_obj, 'fragment_matrix'):
            expected_shape = (len(cistopic_obj.cell_names), len(cistopic_obj.region_names))
            actual_shape = cistopic_obj.fragment_matrix.shape
            if actual_shape == expected_shape:
                print(f"    ✓ fragment_matrix shape correct: {actual_shape}")
            else:
                print(f"    ✗ fragment_matrix shape wrong: {actual_shape} vs expected {expected_shape}")

        # FIX BARCODES: Remove -TH1/-TH2 suffixes to match scRNA-seq format
        print(f"  Fixing barcodes for {sample_id}...")
        cistopic_obj.cell_data.index = ['-'.join(barcode.split('-')[0:2]) for barcode in cistopic_obj.cell_data.index]

        # ADD THIS LINE: Ensure barcodes are also in the barcode column for consistency
        cistopic_obj.cell_data['barcode'] = cistopic_obj.cell_data.index

        # NEW: Remove cells with nan celltype annotations if they exist
        if 'celltype' in cistopic_obj.cell_data.columns:
            before_count = len(cistopic_obj.cell_data)
            cistopic_obj.cell_data = cistopic_obj.cell_data[cistopic_obj.cell_data['celltype'].notna()]
            after_count = len(cistopic_obj.cell_data)
            if before_count != after_count:
                print(f"  Removed {before_count - after_count} cells with nan celltype annotations")

        cistopic_obj_list.append(cistopic_obj)

    # Save all CistopicObjects
    with open(output_pickle, "wb") as f:
        pickle.dump(cistopic_obj_list, f)

    print(f"Saved {len(cistopic_obj_list)} CistopicObjects to {output_pickle}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create CistopicObjects from fragments after QC")
    parser.add_argument("--fragments_dict", required=True,
                        help="Pickle file with sample_id -> fragment path mapping")
    parser.add_argument("--qc_results_pickle", required=True,
                        help="Pickle file with barcodes + thresholds from QC filtering")
    parser.add_argument("--regions_bed", required=True,
                        help="Consensus regions BED file")
    parser.add_argument("--blacklist_bed", required=True,
                        help="Blacklist BED file (mm10)")
    parser.add_argument("--qc_output_dir", required=True,
                        help="Directory containing pycisTopic QC outputs (with parquet files)")
    parser.add_argument("--output_pickle", required=True,
                        help="Path to save CistopicObjects as pickle")
    parser.add_argument("--n_cpu", type=int, default=1,
                        help="Number of CPUs to use for processing")
    args = parser.parse_args()

    main(
        args.fragments_dict,
        args.qc_results_pickle,
        args.regions_bed,
        args.blacklist_bed,
        args.qc_output_dir,
        args.output_pickle,
        args.n_cpu
    )
