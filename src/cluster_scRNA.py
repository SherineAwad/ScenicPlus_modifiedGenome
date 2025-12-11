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
parser.add_argument('markers', help="Text file with one marker gene per line")
args = parser.parse_args()

myObject = args.myObject
newObject = "clustered_" + myObject
markers = args.markers
print("Markers file:", markers)

# Base name for figure outputs
base_name = os.path.splitext(os.path.basename(newObject))[0]

# -----------------------------
# Load AnnData
# -----------------------------
combined_adata = sc.read(myObject)

# -----------------------------
# CRITICAL FIX: Preserve raw counts for SCENIC+
# -----------------------------
# Store the raw counts before normalization
combined_adata.raw = combined_adata
print("✅ Raw counts preserved in adata.raw for SCENIC+ compatibility")

# -----------------------------
# Rename samples before anything else
# -----------------------------
rename_dict = {"TH1": "Control", "TH2": "OE"}
if "sample" in combined_adata.obs.columns:
    combined_adata.obs["sample"] = combined_adata.obs["sample"].replace(rename_dict)
    print("✅ Renamed samples:", combined_adata.obs["sample"].unique())
else:
    print("⚠️ Column 'sample' not found. Available columns:", combined_adata.obs.columns.tolist())

# -----------------------------
# Preprocessing & clustering
# -----------------------------
sc.pp.normalize_total(combined_adata, target_sum=1e4)
sc.pp.log1p(combined_adata)
combined_adata.raw = combined_adata 
sc.pp.scale(combined_adata)

sc.tl.pca(combined_adata, svd_solver="arpack")
sc.pp.neighbors(combined_adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(combined_adata)

sc.tl.leiden(combined_adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)

# Save cluster plot
figure_name = f"_{base_name}_Clusters.png"
sc.pl.umap(combined_adata, color=["leiden"], legend_loc="on data", save=figure_name)
sc.pl.umap(combined_adata, color="sample", legend_loc="on data", save=f"_{base_name}_Samples.png")

# -----------------------------
# Marker gene plotting
# -----------------------------
with open(markers) as f:
    marker_genes = [line.strip() for line in f]
figurename = f"figures/{base_name}_Dotplot.png"
combined_adata.obs["leiden"] = combined_adata.obs["leiden"].astype("category")
marker_genes_present = [g for g in marker_genes if g in combined_adata.var_names]

print([g for g in combined_adata.var_names if "Neurog2" in g])


fig = sc.pl.dotplot(
    combined_adata,
    var_names=marker_genes_present,
    groupby="leiden",
    standard_scale="var",
    show=False,
    return_fig=True
)
# Save figure
fig.savefig(figurename, dpi=600, bbox_inches="tight")

# -----------------------------
# Individual gene UMAP plots
# -----------------------------
with open(markers) as f:
    marker_genes = [line.strip() for line in f]

for gene in marker_genes:
    if gene in combined_adata.var_names:
        sc.pl.umap(
            combined_adata,
            color=gene,
            title=gene,
            save=f"_{base_name}_{gene}.png"
        )
    else:
        print(f"⚠️ Skipping {gene}: not in adata.var_names")

# -----------------------------
# Save clustered object
# -----------------------------
combined_adata.obs_names_make_unique()
combined_adata.write(newObject, compression="gzip")
print(f"✅ Saved clustered AnnData to {newObject}")
print(f"✅ Raw counts preserved for SCENIC+ compatibility")
