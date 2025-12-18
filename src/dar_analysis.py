#!/usr/bin/env python3
import argparse
import os
import pickle
import numpy as np
import pandas as pd
from types import SimpleNamespace
import scipy.sparse as sp

from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)

# FIX: Define ModelObject at module level for pickle compatibility
class ModelObject:
    def __init__(self, model_dict):
        # Copy all dictionary items to object attributes
        for key, value in model_dict.items():
            setattr(self, key, value)

def run_dar(cistopic_pickle, output_dir, var_column, scale_factor_impute=1e7,
            scale_factor_norm=1e4, n_cpu=5, temp_dir=None,
            adjpval_thr=0.1, log2fc_thr=np.log2(1.2)):

    os.makedirs(output_dir, exist_ok=True)

    # Load the binarized Cistopic object
    with open(cistopic_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # FIX: Convert selected_model from dict to object for pycisTopic compatibility
    if hasattr(cistopic_obj, 'selected_model') and isinstance(cistopic_obj.selected_model, dict):
        print("[INFO] Converting selected_model dict to object for pycisTopic compatibility")
        
        # Create object version using the module-level ModelObject
        cistopic_obj.selected_model = ModelObject(cistopic_obj.selected_model)
        print("[INFO] selected_model converted to object")

    # CRITICAL FIX: Check and fix topic_region orientation
    # pycisTopic expects: topic_region rows = region names, columns = topic names
    print("[INFO] Checking topic_region orientation...")
    
    if hasattr(cistopic_obj, 'selected_model') and hasattr(cistopic_obj.selected_model, 'topic_region'):
        tr = cistopic_obj.selected_model.topic_region
        print(f"  topic_region shape: {tr.shape}")
        
        if isinstance(tr, pd.DataFrame):
            # Check if region names are in rows (index) or columns
            region_names = cistopic_obj.region_names
            
            # Sample check: are region names in columns?
            if any(region in tr.columns for region in region_names[:10]):
                print(f"[FIX] topic_region has region names as COLUMNS - transposing!")
                cistopic_obj.selected_model.topic_region = tr.T
                print(f"  After transpose shape: {cistopic_obj.selected_model.topic_region.shape}")

    # CRITICAL FIX: Transpose cell_topic to match pycisTopic's expectation
    # pycisTopic expects: topics as rows, cells as columns
    # Your data has: cells as rows, topics as columns
    print("[INFO] Checking cell_topic orientation...")
    
    if hasattr(cistopic_obj, 'cell_topic'):
        ct = cistopic_obj.cell_topic
        print(f"  Original cell_topic shape: {ct.shape}")
        
        # Check if we need to transpose
        # If first column name looks like a topic name (Topic1, Topic2, etc.)
        if isinstance(ct, pd.DataFrame) and len(ct.columns) > 0:
            first_col = str(ct.columns[0])
            if first_col.startswith('Topic') or first_col.startswith('topic'):
                print(f"[FIX] Transposing cell_topic (cells×topics -> topics×cells)")
                cistopic_obj.cell_topic = ct.T
                print(f"  Transposed cell_topic shape: {cistopic_obj.cell_topic.shape}")
    
    # Also fix selected_model.cell_topic
    if hasattr(cistopic_obj, 'selected_model') and hasattr(cistopic_obj.selected_model, 'cell_topic'):
        sm_ct = cistopic_obj.selected_model.cell_topic
        if isinstance(sm_ct, pd.DataFrame) and len(sm_ct.columns) > 0:
            first_col = str(sm_ct.columns[0])
            if first_col.startswith('Topic') or first_col.startswith('topic'):
                print(f"[FIX] Transposing selected_model.cell_topic")
                cistopic_obj.selected_model.cell_topic = sm_ct.T

    # Validate that it has the necessary attributes
    for attr in ["fragment_matrix", "cell_topic", "cell_names", "region_names"]:
        if not hasattr(cistopic_obj, attr):
            raise AttributeError(f"[ERROR] cistopic_obj missing required attribute: {attr}")

    # Check fragment_matrix orientation
    fm = cistopic_obj.fragment_matrix
    n_cells = len(cistopic_obj.cell_names)
    n_regions = len(cistopic_obj.region_names)
    if fm.shape == (n_cells, n_regions):
        print("[INFO] Matrix orientation OK (cells x regions).")
    elif fm.shape == (n_regions, n_cells):
        print("[FIX] Transposing fragment_matrix (regions x cells -> cells x regions).")
        cistopic_obj.fragment_matrix = fm.T
    else:
        raise ValueError(f"[ERROR] Unexpected fragment_matrix shape: {fm.shape}")

    # DEBUG: Check original fragment matrix statistics
    print("\n[DEBUG] Original fragment matrix stats:")
    if sp.issparse(fm):
        print(f"  Non-zero entries: {fm.nnz}")
        print(f"  Sparsity: {1 - (fm.nnz / (fm.shape[0] * fm.shape[1])):.4f}")
        print(f"  Min value: {fm.data.min() if len(fm.data) > 0 else 0}")
        print(f"  Max value: {fm.data.max() if len(fm.data) > 0 else 0}")
        print(f"  Mean value: {fm.data.mean() if len(fm.data) > 0 else 0}")

    # Run imputation
    print("[INFO] Running imputation of accessibility...")
    imputed_acc_obj = impute_accessibility(
        cistopic_obj,
        selected_cells=None,
        selected_regions=None,
        scale_factor=scale_factor_impute
    )

    # DEBUG: Check imputed data statistics
    print("\n[DEBUG] Imputed accessibility stats:")
    # FIX: Use .mtx attribute for CistopicImputedFeatures object
    imputed_matrix = imputed_acc_obj.mtx
    print(f"  Shape: {imputed_matrix.shape}")
    if sp.issparse(imputed_matrix):
        print(f"  Non-zero entries: {imputed_matrix.nnz}")
        print(f"  Sparsity: {1 - (imputed_matrix.nnz / (imputed_matrix.shape[0] * imputed_matrix.shape[1])):.4f}")
        if imputed_matrix.nnz > 0:
            print(f"  Min value: {imputed_matrix.data.min():.6f}")
            print(f"  Max value: {imputed_matrix.data.max():.6f}")
            print(f"  Mean value: {imputed_matrix.data.mean():.6f}")
            # Check for any non-zero values
            non_zero_count = (imputed_matrix.data > 0).sum()
            print(f"  Values > 0: {non_zero_count}/{imputed_matrix.nnz}")
    else:
        print(f"  Min value: {imputed_matrix.min():.6f}")
        print(f"  Max value: {imputed_matrix.max():.6f}")
        print(f"  Mean value: {imputed_matrix.mean():.6f}")

    # Normalize imputed accessibility
    print("[INFO] Normalizing imputed accessibility...")
    normalized_imputed_acc_obj = normalize_scores(
        imputed_acc_obj,
        scale_factor=scale_factor_norm
    )

    # DEBUG: Check normalized data statistics
    print("\n[DEBUG] Normalized imputed accessibility stats:")
    # FIX: Use .mtx attribute for normalized object too
    norm_matrix = normalized_imputed_acc_obj.mtx
    print(f"  Shape: {norm_matrix.shape}")
    if sp.issparse(norm_matrix):
        print(f"  Non-zero entries: {norm_matrix.nnz}")
        if norm_matrix.nnz > 0:
            print(f"  Min value: {norm_matrix.data.min():.6f}")
            print(f"  Max value: {norm_matrix.data.max():.6f}")
            print(f"  Mean value: {norm_matrix.data.mean():.6f}")
    else:
        print(f"  Min value: {norm_matrix.min():.6f}")
        print(f"  Max value: {norm_matrix.max():.6f}")
        print(f"  Mean value: {norm_matrix.mean():.6f}")

    # Find highly variable regions with relaxed parameters
    print("[INFO] Finding highly variable regions (with relaxed parameters)...")
    variable_regions = find_highly_variable_features(
        normalized_imputed_acc_obj,
        min_disp=0.001,      # RELAXED: was 0.01
        min_mean=0.0001,     # RELAXED: was 0.001  
        max_mean=10,         # RELAXED: was 3
        max_disp=np.inf,
        n_bins=20,
        n_top_features=5000,  # NEW: Force top 5000 regions
        plot=False
    )
    print(f"[INFO] Number of highly variable regions: {len(variable_regions)}")
    
    # If still no variable regions, investigate further
    if len(variable_regions) == 0:
        print("\n[ERROR] CRITICAL: Still no highly variable regions found!")
        print("[ERROR] Possible causes:")
        print("  1. Imputation failed (all values are 0 or NaN)")
        print("  2. Normalization failed (scale factor too large/small)")
        print("  3. Topic model doesn't capture meaningful variation")
        print("  4. Cell_topic matrix has problems")
        
        # Try manual calculation as last resort
        print("\n[WARNING] Attempting manual workaround...")
        # Calculate mean and variance manually
        if sp.issparse(norm_matrix):
            # For sparse matrix, calculate mean per region
            region_means = np.array(norm_matrix.mean(axis=0)).flatten()
            # Simple variance approximation for sparse data
            region_vars = np.array((norm_matrix.power(2).mean(axis=0) - np.square(region_means))).flatten()
        else:
            region_means = norm_matrix.mean(axis=0)
            region_vars = norm_matrix.var(axis=0)
        
        print(f"  Region means - min: {region_means.min():.6f}, max: {region_means.max():.6f}")
        print(f"  Region vars - min: {region_vars.min():.6f}, max: {region_vars.max():.6f}")
        
        # Select top 5000 regions by variance
        if len(region_vars) > 5000:
            top_indices = np.argsort(region_vars)[-5000:]
            variable_regions = [normalized_imputed_acc_obj.region_names[i] for i in top_indices]
            print(f"  Selected top 5000 variable regions manually")
        else:
            variable_regions = normalized_imputed_acc_obj.region_names
            print(f"  Using all regions as fallback")
        
        print(f"[INFO] Now using {len(variable_regions)} regions for DAR")

    # Find DARs
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

    # Attach markers_dict
    cistopic_obj.markers_dict = markers_dict

    # FIX: Convert ModelObject back to dict before saving for pickle compatibility
    if hasattr(cistopic_obj, 'selected_model') and isinstance(cistopic_obj.selected_model, ModelObject):
        print("[INFO] Converting selected_model back to dict for pickle compatibility")
        model_dict = {}
        for attr in dir(cistopic_obj.selected_model):
            if not attr.startswith('_'):
                value = getattr(cistopic_obj.selected_model, attr)
                model_dict[attr] = value
        cistopic_obj.selected_model = model_dict

    # Save updated object
    out_obj_file = os.path.join(output_dir, "cistopic_obj_with_DARs.pkl")
    with open(out_obj_file, "wb") as f:
        pickle.dump(cistopic_obj, f)
    print(f"[INFO] Saved CistopicObject with DARs -> {out_obj_file}")

    # Summary
    print("\n[SUMMARY] Number of DARs found per group:")
    for group, df in markers_dict.items():
        print(f"  {group}: {len(df)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Differential Accessibility Analysis (DAR) with pycisTopic")
    parser.add_argument("-i", "--input_pickle", required=True, help="Binarized CistopicObject pickle")
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
