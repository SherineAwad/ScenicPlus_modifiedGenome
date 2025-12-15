#!/usr/bin/env python3
import argparse
import os
import pickle
import numpy as np

from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)

def run_dar(cistopic_pickle, output_dir, var_column, scale_factor_impute=1e7, scale_factor_norm=1e4,
            n_cpu=5, temp_dir=None, adjpval_thr=0.1, log2fc_thr=np.log2(1.2)):
    os.makedirs(output_dir, exist_ok=True)

    # Load CistopicObject
    with open(cistopic_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # Validate and correct fragment_matrix orientation if needed
    if hasattr(cistopic_obj, "fragment_matrix"):
        n_cells = len(cistopic_obj.cell_names)
        n_regions = len(cistopic_obj.region_names)
        if cistopic_obj.fragment_matrix.shape == (n_cells, n_regions):
            print("[INFO] Matrix orientation OK (cells x regions).")
        elif cistopic_obj.fragment_matrix.shape == (n_regions, n_cells):
            print("[FIX] Transposing fragment_matrix (regions x cells -> cells x regions).")
            cistopic_obj.fragment_matrix = cistopic_obj.fragment_matrix.T
        else:
            raise ValueError(
                f"[ERROR] Unexpected fragment_matrix shape: {cistopic_obj.fragment_matrix.shape}, "
                f"expected ({n_cells}, {n_regions}) or ({n_regions}, {n_cells})."
            )
    else:
        raise AttributeError("[ERROR] CistopicObject is missing 'fragment_matrix' attribute.")

    print("[INFO] Running imputation of accessibility...")
    imputed_acc_obj = impute_accessibility(
        cistopic_obj,
        selected_cells=None,
        selected_regions=None,
        scale_factor=scale_factor_impute
    )

    print("[INFO] Normalizing imputed accessibility...")
    normalized_imputed_acc_obj = normalize_scores(
        imputed_acc_obj,
        scale_factor=scale_factor_norm
    )

    print("[INFO] Finding highly variable regions...")
    variable_regions = find_highly_variable_features(
        normalized_imputed_acc_obj,
        min_disp=0.01,
        min_mean=0.001,
        max_mean=3,
        max_disp=np.inf,
        n_bins=20,
        n_top_features=None,
        plot=False  # Disable plotting for CLI runs
    )
    print(f"[INFO] Number of highly variable regions: {len(variable_regions)}")

    print("[INFO] Finding differential accessibility regions (DARs)...")
    markers_dict = find_diff_features(
        cistopic_obj,
        imputed_acc_obj,
        variable=var_column,
        var_features=variable_regions,
        contrasts=None,
        adjpval_thr=adjpval_thr,
        log2fc_thr=log2fc_thr,
        n_cpu=n_cpu,
        _temp_dir=temp_dir if temp_dir else "/tmp",
        split_pattern='-'
    )

    # Save marker tables
    for celltype, df in markers_dict.items():
        out_file = os.path.join(output_dir, f"markers_{celltype}.tsv")
        df.to_csv(out_file, sep="\t")
        print(f"[INFO] Saved markers for {celltype} -> {out_file}")

    # Attach markers_dict safely
    setattr(cistopic_obj, "markers_dict", markers_dict)

    # Save updated object
    cistopic_obj_file = os.path.join(output_dir, "cistopic_obj_with_DARs.pkl")
    with open(cistopic_obj_file, "wb") as f:
        pickle.dump(cistopic_obj, f)
    print(f"[INFO] Saved CistopicObject with DARs -> {cistopic_obj_file}")

    # Print summary
    print("\n[SUMMARY] Number of DARs found per group:")
    for group, df in markers_dict.items():
        print(f"  {group}: {len(df)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Differential Accessibility Analysis (DAR) with pycisTopic")
    parser.add_argument("-i", "--input_pickle", required=True, help="CistopicObject pickle")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save DAR results and plots")
    parser.add_argument("-v", "--var_column", default="celltype", help="Metadata column in cell_data for group comparisons")
    parser.add_argument("--scale_impute", type=float, default=1e7, help="Scale factor for impute_accessibility")
    parser.add_argument("--scale_norm", type=float, default=1e4, help="Scale factor for normalize_scores")
    parser.add_argument("--n_cpu", type=int, default=5, help="Number of CPUs for find_diff_features")
    parser.add_argument("--temp_dir", default=None, help="Temporary directory for intermediate files")
    parser.add_argument("--adjpval_thr", type=float, default=0.1, help="Adjusted p-value threshold for DAR")
    parser.add_argument("--log2fc_thr", type=float, default=np.log2(1.2), help="Log2 fold-change threshold for DAR")
    args = parser.parse_args()

    run_dar(
        cistopic_pickle=args.input_pickle,
        output_dir=args.output_dir,
        var_column=args.var_column,
        scale_factor_impute=args.scale_impute,
        scale_factor_norm=args.scale_norm,
        n_cpu=args.n_cpu,
        temp_dir=args.temp_dir,
        adjpval_thr=args.adjpval_thr,
        log2fc_thr=args.log2fc_thr
    )

