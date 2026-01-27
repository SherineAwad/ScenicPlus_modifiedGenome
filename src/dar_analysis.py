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

    # CREATE OUTPUT DIRECTORY FIRST WITH PROPER ERROR HANDLING
    print(f"[INFO] Creating output directory: {output_dir}")
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"[INFO] Output directory ready: {output_dir}")
    except Exception as e:
        print(f"[ERROR] Failed to create output directory: {e}")
        raise

    # Load the binarized Cistopic object
    print(f"[INFO] Loading cisTopic object from {cistopic_pickle}")
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

    # Run imputation
    print("[INFO] Running imputation of accessibility...")
    imputed_acc_obj = impute_accessibility(
        cistopic_obj,
        selected_cells=None,
        selected_regions=None,
        scale_factor=scale_factor_impute
    )

    # Normalize imputed accessibility
    print("[INFO] Normalizing imputed accessibility...")
    normalized_imputed_acc_obj = normalize_scores(
        imputed_acc_obj,
        scale_factor=scale_factor_norm
    )

    # ATTACH IMPUTED OBJECTS FOR GENE ACTIVITY ANALYSIS
    cistopic_obj.imputed_acc_obj = imputed_acc_obj
    cistopic_obj.normalized_imputed_acc_obj = normalized_imputed_acc_obj

    # Find highly variable regions with relaxed parameters
    print("[INFO] Finding highly variable regions (with relaxed parameters)...")
    variable_regions = find_highly_variable_features(
        normalized_imputed_acc_obj,
        min_disp=0.001,
        min_mean=0.0001,
        max_mean=10,
        max_disp=np.inf,
        n_bins=20,
        n_top_features=5000,
        plot=False
    )
    print(f"[INFO] Number of highly variable regions: {len(variable_regions)}")

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
        try:
            df.to_csv(out_file, sep="\t")
            print(f"[INFO] Saved markers for {celltype} -> {out_file}")
        except Exception as e:
            print(f"[ERROR] Failed to save markers for {celltype}: {e}")

    # Attach markers_dict
    cistopic_obj.markers_dict = markers_dict

    # ADD PLOTTING HERE - FIXED VERSION
    print("\n" + "="*60)
    print("[INFO] Creating UMAP plots for top DARs...")
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        # Create plot directory
        plot_dir = os.path.join(output_dir, "umap_plots")
        os.makedirs(plot_dir, exist_ok=True)
        print(f"[INFO] Plot directory: {plot_dir}")
        
        # Check if UMAP needs to be computed
        if "cell" not in cistopic_obj.projections or "UMAP" not in cistopic_obj.projections["cell"]:
            print("[INFO] Computing UMAP for plotting...")
            
            try:
                # Try umap-learn with simplified parameters
                from umap import UMAP
                
                # Get cell_topic and ensure correct orientation
                cell_topic = cistopic_obj.cell_topic
                if cell_topic.shape[1] == len(cistopic_obj.cell_names):
                    cell_topic = cell_topic.T
                
                print(f"  Using cell_topic shape: {cell_topic.shape}")
                print("  Computing UMAP...")
                
                # Use Euclidean distance to avoid pynndescent issues
                umap_model = UMAP(
                    n_neighbors=15,
                    min_dist=0.1,
                    random_state=42,
                    n_components=2,
                    metric='euclidean',
                    n_jobs=1,
                    verbose=False
                )
                
                umap_embedding = umap_model.fit_transform(cell_topic)
                
            except Exception as e:
                print(f"  UMAP failed ({e}), using PCA instead...")
                from sklearn.decomposition import PCA
                
                cell_topic = cistopic_obj.cell_topic
                if cell_topic.shape[1] == len(cistopic_obj.cell_names):
                    cell_topic = cell_topic.T
                
                pca = PCA(n_components=2, random_state=42)
                umap_embedding = pca.fit_transform(cell_topic)
                print(f"  Using PCA (explained variance: {pca.explained_variance_ratio_})")
            
            # Add to projections
            if "cell" not in cistopic_obj.projections:
                cistopic_obj.projections["cell"] = {}
            
            cistopic_obj.projections["cell"]["UMAP"] = pd.DataFrame(
                umap_embedding,
                columns=["UMAP_1", "UMAP_2"],
                index=cistopic_obj.cell_names[:cell_topic.shape[0]]
            )
            print("[INFO] Embedding computed successfully")
        
        # Get UMAP coordinates
        umap_df = cistopic_obj.projections["cell"]["UMAP"]
        x = umap_df["UMAP_1"].values
        y = umap_df["UMAP_2"].values
        
        # Get region names from the cistopic object, not imputed_acc_obj
        region_names = cistopic_obj.region_names
        
        # Plot top DAR for each group
        plotted = 0
        for group, df in markers_dict.items():
            if df.shape[0] == 0:
                continue
            
            top_dar = df.index.tolist()[0]
            safe_dar = str(top_dar).replace(':', '_').replace('/', '_').replace('\\', '_')
            safe_group = group.replace(' ', '_')
            out_png = os.path.join(plot_dir, f"UMAP_{safe_group}_{safe_dar}.png")
            
            print(f"  Plotting top DAR for '{group}': {top_dar}")
            
            try:
                # Find region index in region_names
                if top_dar not in region_names:
                    print(f"    [ERROR] Region {top_dar} not found in region_names")
                    continue
                
                idx = list(region_names).index(top_dar)
                
                # Get accessibility values from imputed data
                # Try different ways to access the imputed matrix
                imputed_matrix = None
                
                # Method 1: Check if it has .mtx attribute
                if hasattr(imputed_acc_obj, 'mtx'):
                    imputed_matrix = imputed_acc_obj.mtx
                # Method 2: Check if it's the matrix itself
                elif hasattr(imputed_acc_obj, 'shape'):
                    imputed_matrix = imputed_acc_obj
                # Method 3: Try to access as DataFrame
                elif hasattr(imputed_acc_obj, 'columns') and hasattr(imputed_acc_obj, 'iloc'):
                    imputed_matrix = imputed_acc_obj
                else:
                    print(f"    [ERROR] Cannot access imputed data structure")
                    continue
                
                # Extract accessibility values
                if sp.issparse(imputed_matrix):
                    if imputed_matrix.shape[1] > idx:
                        acc_values = imputed_matrix[:, idx].toarray().flatten()
                    else:
                        print(f"    [ERROR] Index {idx} out of bounds for matrix shape {imputed_matrix.shape}")
                        continue
                elif hasattr(imputed_matrix, 'iloc'):  # DataFrame
                    if idx < len(imputed_matrix.columns):
                        acc_values = imputed_matrix.iloc[:, idx].values
                    else:
                        print(f"    [ERROR] Index {idx} out of bounds for DataFrame")
                        continue
                elif hasattr(imputed_matrix, 'shape'):  # numpy array
                    if imputed_matrix.shape[1] > idx:
                        acc_values = imputed_matrix[:, idx]
                    else:
                        print(f"    [ERROR] Index {idx} out of bounds for array")
                        continue
                else:
                    print(f"    [ERROR] Unknown imputed data type: {type(imputed_matrix)}")
                    continue
                
                # Create plot
                fig, ax = plt.subplots(figsize=(10, 8))
                scatter = ax.scatter(x, y, c=acc_values, cmap='viridis', s=5, alpha=0.7, edgecolors='none')
                plt.colorbar(scatter, ax=ax, label='Imputed Accessibility')
                ax.set_xlabel('UMAP 1')
                ax.set_ylabel('UMAP 2')
                ax.set_title(f"Top DAR for {group}\n{top_dar}", fontsize=14)
                plt.tight_layout()
                plt.savefig(out_png, dpi=300, bbox_inches='tight')
                plt.close(fig)
                
                plotted += 1
                print(f"    ✓ Saved to {out_png}")
                
            except Exception as e:
                print(f"    [ERROR] Failed to plot: {e}")
        
        print(f"[INFO] Successfully plotted {plotted} out of {len(markers_dict)} DARs")
        
    except ImportError as e:
        print(f"[WARNING] Could not import plotting dependencies: {e}")
        print("[INFO] Skipping plotting...")
    except Exception as e:
        print(f"[WARNING] Plotting failed: {e}")
        print("[INFO] Continuing without plots...")
    
    print("="*60)

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
