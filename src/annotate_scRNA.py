import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import os

# Fix for importlib metadata
sys.modules['importlib.metadata'] = importlib_metadata

# -----------------------------
# Parse arguments
# -----------------------------
parser = argparse.ArgumentParser()
parser.add_argument('myObject', help="Input AnnData file (.h5ad)")
parser.add_argument('annotations', help="Text file with annotations")
args = parser.parse_args()

myObject = args.myObject
newObject = "annotated_" + myObject
annot_file = args.annotations

# Base name for figure outputs
base_name = os.path.splitext(os.path.basename(newObject))[0]

# -----------------------------
# Load AnnData
# -----------------------------
combined_adata = sc.read_h5ad(myObject)

cluster_to_celltype_dict = {}
with open(annot_file, "r") as f:
    for line in f:
        cluster, celltype = line.strip().split(',')
        cluster_to_celltype_dict[cluster] = celltype

cluster_to_celltype_dict = {str(key): value for key, value in cluster_to_celltype_dict.items()}

combined_adata = combined_adata[combined_adata.obs["leiden"] != "80"].copy()
combined_adata.obs["celltype"] = combined_adata.obs["leiden"].map(cluster_to_celltype_dict)

figure_name = "figures/"+base_name+"_annotationsON.png"
combined_adata.obs_names_make_unique()
fig = sc.pl.umap(combined_adata, color='celltype',legend_loc="on data", show=False,return_fig=True)
fig.savefig(figure_name, dpi=600,bbox_inches='tight')
plt.close(fig)

figure_name = "figures/"+base_name+"_annotations.png"
fig = sc.pl.umap(combined_adata, color='celltype', show=False,return_fig=True)
fig.savefig(figure_name, dpi=600,bbox_inches='tight')
plt.close(fig)

# -----------------------------
# Save clustered object
# -----------------------------
combined_adata.obs_names_make_unique()
combined_adata.write(newObject, compression="gzip")
print(f"✅ Saved clustered AnnData to {newObject}")

