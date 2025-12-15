#!/usr/bin/env python3
import argparse
import pickle
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np

def cluster_cistopic_manual(input_pickle, output_pickle, outdir, resolutions=[0.6, 1.2, 3], k=10):
    os.makedirs(outdir, exist_ok=True)

    # Load object
    with open(input_pickle, "rb") as f:
        obj = pickle.load(f)

    if 'celltype' not in obj.cell_data.columns:
        raise ValueError("Missing 'celltype' column in obj.cell_data — cannot plot or cluster.")

    print(f"Available cell types: {obj.cell_data['celltype'].unique()}")

    print("Performing manual clustering with Scanpy...")

    # Get topic matrix
    if hasattr(obj, 'selected_model') and hasattr(obj.selected_model, 'cell_topic'):
        cell_topic = obj.selected_model.cell_topic
        print(f"Using cell_topic from selected_model: {cell_topic.shape}")
    elif hasattr(obj, 'cell_topic'):
        cell_topic = obj.cell_topic
        print(f"Using cell_topic from object: {cell_topic.shape}")
    else:
        raise ValueError("No cell_topic found in object or selected_model.")

    # Orient matrix: cells × topics
    if cell_topic.shape[0] == obj.cell_data.shape[0]:
        topic_matrix = cell_topic
        print("Using cell_topic as cells × topics")
    else:
        topic_matrix = cell_topic.T
        print("Using cell_topic.T as cells × topics")

    # Create AnnData
    adata = sc.AnnData(X=topic_matrix)
    adata.obs_names = obj.cell_data.index.astype(str)
    adata.var_names = [f"Topic_{i}" for i in range(adata.n_vars)]

    # Add metadata
    for col in obj.cell_data.columns:
        adata.obs[col] = obj.cell_data[col].astype(str).values

    print(f"AnnData shape: {adata.shape}")
    print(f"AnnData obs_names match cell_data index: {adata.obs_names.equals(obj.cell_data.index)}")

    # Compute neighbors and clustering
    sc.pp.neighbors(adata, n_neighbors=k, use_rep='X', random_state=0)

    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f'pycisTopic_leiden_{res}')

    # UMAP
    sc.tl.umap(adata, random_state=0)

    # Transfer cluster + UMAP back
    for res in resolutions:
        cluster_col = f'pycisTopic_leiden_{res}'
        obj.cell_data[cluster_col] = adata.obs[cluster_col].values

    obj.cell_data['UMAP_1'] = adata.obsm['X_umap'][:, 0]
    obj.cell_data['UMAP_2'] = adata.obsm['X_umap'][:, 1]

    # -------- PLOTS --------

    # 1. Cell type UMAP
    plt.figure(figsize=(12, 10))
    unique_celltypes = obj.cell_data['celltype'].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_celltypes)))
    color_map = dict(zip(unique_celltypes, colors))

    centroids = {}
    for celltype in unique_celltypes:
        mask = obj.cell_data['celltype'] == celltype
        if mask.sum() > 0:
            x_data = obj.cell_data.loc[mask, 'UMAP_1']
            y_data = obj.cell_data.loc[mask, 'UMAP_2']
            plt.scatter(x_data, y_data, c=[color_map[celltype]],
                        label=celltype, s=10, alpha=0.7)
            centroids[celltype] = (x_data.mean(), y_data.mean())

    for celltype, (x, y) in centroids.items():
        plt.annotate(celltype, (x, y), xytext=(5, 5),
                     textcoords='offset points', fontsize=9,
                     bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('Cell Types - UMAP')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "celltype_umap.png"), bbox_inches='tight', dpi=300)
    plt.close()
    print("✓ Created celltype UMAP with labels")

    # 2. Annotated clusters
    for res in resolutions:
        cluster_col = f'pycisTopic_leiden_{res}'
        plt.figure(figsize=(12, 10))
        unique_clusters = sorted(obj.cell_data[cluster_col].unique())
        cluster_centroids = {}

        for cluster in unique_clusters:
            mask = obj.cell_data[cluster_col] == cluster
            celltype_counts = obj.cell_data.loc[mask, 'celltype'].value_counts()
            if len(celltype_counts) > 0:
                dominant_celltype = celltype_counts.index[0]
                dominant_fraction = celltype_counts.iloc[0] / mask.sum()
            else:
                dominant_celltype = 'Unknown'
                dominant_fraction = 0
            plt.scatter(obj.cell_data.loc[mask, 'UMAP_1'],
                        obj.cell_data.loc[mask, 'UMAP_2'],
                        c=[color_map.get(dominant_celltype, 'gray')],
                        s=10, alpha=0.7)
            x_centroid = obj.cell_data.loc[mask, 'UMAP_1'].mean()
            y_centroid = obj.cell_data.loc[mask, 'UMAP_2'].mean()
            cluster_centroids[cluster] = (x_centroid, y_centroid, dominant_celltype, dominant_fraction)

        for cluster, (x, y, celltype, fraction) in cluster_centroids.items():
            plt.annotate(f'C{cluster}: {celltype}\n({fraction:.1%})',
                         (x, y), xytext=(5, 5),
                         textcoords='offset points', fontsize=9, fontweight='bold',
                         bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
        plt.title(f'Clusters annotated with dominant cell types (resolution {res})')
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"annotated_clusters_res_{res}.png"),
                    bbox_inches='tight', dpi=300)
        plt.close()

        print(f"\n📊 Resolution {res} - Cluster composition:")
        for cluster in unique_clusters:
            mask = obj.cell_data[cluster_col] == cluster
            celltype_composition = obj.cell_data.loc[mask, 'celltype'].value_counts(normalize=True)
            dominant = celltype_composition.index[0] if len(celltype_composition) > 0 else 'Unknown'
            dominant_pct = celltype_composition.iloc[0] * 100 if len(celltype_composition) > 0 else 0
            print(f"  Cluster {cluster}: {mask.sum()} cells, dominant: {dominant} ({dominant_pct:.1f}%)")

        print(f"✓ Created annotated clusters for resolution {res}")

    # 3. Topic-celltype heatmap
    try:
        if cell_topic.shape[0] == obj.cell_data.shape[0]:
            cell_topic_df = pd.DataFrame(cell_topic, index=obj.cell_data.index)
        else:
            cell_topic_df = pd.DataFrame(cell_topic.T, index=obj.cell_data.index)

        cell_topic_df.columns = [f"Topic_{i}" for i in range(cell_topic_df.shape[1])]
        cell_topic_df['celltype'] = obj.cell_data['celltype'].values
        cell_topic_df = cell_topic_df.dropna(subset=['celltype'])
        topic_by_celltype = cell_topic_df.groupby('celltype').mean()

        plt.figure(figsize=(12, 8))
        sns.heatmap(topic_by_celltype, cmap='viridis', annot=False, linewidths=0.5)
        plt.title('Topic Enrichment by Cell Type')
        plt.xlabel('Topics')
        plt.ylabel('Cell Types')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "topic_celltype_heatmap.png"), dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Created topic-celltype heatmap")
    except Exception as e:
        print(f"⚠️ Topic heatmap failed: {e}")

    # Save updated object
    with open(output_pickle, "wb") as f:
        pickle.dump(obj, f)
    print(f"\n✓ Manual clustering completed and saved to: {output_pickle}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster CistopicObject with manual Scanpy method")
    parser.add_argument("-i", "--input_pickle", required=True)
    parser.add_argument("-o", "--output_pickle", required=True)
    parser.add_argument("-d", "--outdir", required=True)
    parser.add_argument("--resolutions", nargs="+", type=float, default=[0.6, 1.2, 3])
    parser.add_argument("--k", type=int, default=10)
    args = parser.parse_args()
    cluster_cistopic_manual(args.input_pickle, args.output_pickle, args.outdir, args.resolutions, args.k)

