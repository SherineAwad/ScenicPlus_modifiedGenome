import scanpy as sc
import pandas as pd
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
parser.add_argument('barcode_csv', help="Output CSV file for barcode annotations")

args = parser.parse_args()

myObject = args.myObject
barcode_mapping = args.barcode_csv

# FIXED: Use the input argument instead of hardcoded filename
adata = sc.read_h5ad(myObject)

# Make sure the annotations are in adata.obs
print(adata.obs.columns)

# Suppose your cell type annotations are stored in 'celltype'
# If instead they are in 'leiden' or 'clusters', adjust accordingly.
barcode_annotations = adata.obs[['celltype']].copy()
barcode_annotations.index.name = 'barcode'
barcode_annotations.reset_index(inplace=True)

# Save to CSV for SCENIC+
barcode_annotations.to_csv(barcode_mapping, index=False)
