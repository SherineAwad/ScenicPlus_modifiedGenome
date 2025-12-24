#!/usr/bin/env python3
import argparse
import pickle
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot_top_dars_fixed(dar_pickle, cluster_pickle, output_dir):
    """Correct version - matrix is transposed!"""
    os.makedirs(output_dir, exist_ok=True)

    # Load objects
    with open(dar_pickle, "rb") as f:
        dar_obj = pickle.load(f)
    
    with open(cluster_pickle, "rb") as f:
        cluster_obj = pickle.load(f)

    # UMAP coordinates
    umap_df = cluster_obj.cell_data[['UMAP_1', 'UMAP_2']]
    x = umap_df['UMAP_1'].values
    y = umap_df['UMAP_2'].values
    cluster_cells = list(umap_df.index)

    # Get matrix - IT'S TRANSPOSED: regions × cells
    mtx = dar_obj.imputed_acc_obj.mtx  # Shape: (regions, cells)
    print(f"Matrix shape: {mtx.shape} (regions × cells)")
    
    # Create region to ROW index mapping (not column!)
    region_to_row = {}
    for row_idx, region in enumerate(dar_obj.region_names):
        region_to_row[region] = row_idx
    
    # Create cell to COLUMN index mapping
    cell_to_col = {}
    for col_idx, cell in enumerate(dar_obj.cell_names):
        cell_to_col[cell] = col_idx

    # Plot top DAR for each group
    for group, df in dar_obj.markers_dict.items():
        if len(df) == 0:
            continue
        
        # Get top DAR
        region = df.index[0]
        
        if region not in region_to_row:
            print(f"Region {region} not found in matrix rows")
            continue
        
        row_idx = region_to_row[region]
        
        # Get accessibility values for ALL cells (this row)
        row_values = mtx[row_idx, :]
        if hasattr(row_values, "toarray"):
            row_values = row_values.toarray().flatten()
        else:
            row_values = np.array(row_values).flatten()
        
        # Map to cluster cells
        colors = np.zeros(len(cluster_cells))
        for i, cell in enumerate(cluster_cells):
            if cell in cell_to_col:
                colors[i] = row_values[cell_to_col[cell]]
        
        # Debug: check distribution
        print(f"\n{group}: {region}")
        print(f"  Non-zero cells: {np.sum(colors > 0)}/{len(colors)}")
        print(f"  Mean: {np.mean(colors):.4f}, Max: {np.max(colors):.4f}")
        print(f"  95th percentile: {np.percentile(colors, 95):.4f}")
        
        # Plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Use 95th percentile as vmax to see variation
        vmax = np.percentile(colors, 95)
        
        sc = ax.scatter(x, y, c=colors, cmap='viridis', s=5, alpha=0.7,
                       vmin=0, vmax=vmax)
        
        plt.colorbar(sc, ax=ax, label="Accessibility")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(f"{group}: {region}")
        plt.tight_layout()
        
        # Save
        safe_region = region.replace(':', '_').replace('-', '_')
        out_png = os.path.join(output_dir, f"{group}_{safe_region}.png")
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  Saved: {out_png}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dar_pickle", required=True)
    parser.add_argument("-c", "--cluster_pickle", required=True)
    parser.add_argument("-o", "--output_dir", required=True)
    args = parser.parse_args()
    
    plot_top_dars_fixed(args.dar_pickle, args.cluster_pickle, args.output_dir)
