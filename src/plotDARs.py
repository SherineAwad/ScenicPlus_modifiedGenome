#!/usr/bin/env python3
import argparse
import pickle
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import PowerNorm

def plot_top_dars_on_cluster_umap(dar_pickle, cluster_pickle, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Load DAR object
    with open(dar_pickle, "rb") as f:
        dar_obj = pickle.load(f)

    # Load Cluster object
    with open(cluster_pickle, "rb") as f:
        cluster_obj = pickle.load(f)

    # UMAP coordinates
    if not {'UMAP_1', 'UMAP_2'}.issubset(cluster_obj.cell_data.columns):
        raise ValueError("[ERROR] No UMAP columns found in cluster_obj.cell_data")
    umap_df = cluster_obj.cell_data[['UMAP_1', 'UMAP_2']]
    x = umap_df['UMAP_1'].values
    y = umap_df['UMAP_2'].values
    cluster_cells = list(umap_df.index)

    # Top DARs
    if not hasattr(dar_obj, "markers_dict"):
        raise ValueError("[ERROR] DAR object has no markers_dict")
    top_dars = {group: df.index[0] for group, df in dar_obj.markers_dict.items() if len(df) > 0}

    # Imputed matrix
    if not hasattr(dar_obj, "imputed_acc_obj"):
        raise ValueError("[ERROR] DAR object has no imputed_acc_obj")
    imputed_obj = dar_obj.imputed_acc_obj

    if not hasattr(imputed_obj, "mtx"):
        raise ValueError("[ERROR] imputed_acc_obj has no mtx attribute")
    mtx = imputed_obj.mtx
    imputed_cells = list(dar_obj.cell_names)
    imputed_regions = list(dar_obj.region_names)

    for group, region in top_dars.items():
        if region not in imputed_regions:
            print(f"[WARNING] Region {region} not found, skipping {group}")
            continue
        region_idx = imputed_regions.index(region)

        try:
            cell_indices = [imputed_cells.index(c) for c in cluster_cells]
        except ValueError:
            print(f"[WARNING] Some cluster cells not found in DAR object, skipping {group}")
            continue

        values = mtx[cell_indices, region_idx]
        if hasattr(values, "toarray"):
            values = values.toarray().flatten()
        else:
            values = np.array(values).flatten()

        # ---- FIX: better color scaling ----
        # PowerNorm stretches low values while keeping high values distinguishable
        norm = PowerNorm(gamma=0.3, vmin=np.min(values), vmax=np.max(values))

        fig, ax = plt.subplots(figsize=(10, 8))
        sc = ax.scatter(x, y, c=values, cmap='viridis', s=5, alpha=0.7, norm=norm)
        plt.colorbar(sc, ax=ax, label="Imputed Accessibility")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(f"Top DAR: {group}\n{region}")
        plt.tight_layout()
        out_png = os.path.join(output_dir, f"UMAP_{group}_{region.replace(':','_')}.png")
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"[INFO] Plotted {group} -> {out_png}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot top DARs on cluster UMAP with visible color scale")
    parser.add_argument("-d", "--dar_pickle", required=True, help="DAR object pickle")
    parser.add_argument("-c", "--cluster_pickle", required=True, help="Cluster object pickle (with UMAP in cell_data)")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for plots")
    args = parser.parse_args()

    plot_top_dars_on_cluster_umap(args.dar_pickle, args.cluster_pickle, args.output_dir)

