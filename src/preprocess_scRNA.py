import scanpy as sc
import sys
import importlib_metadata
import argparse
import warnings
import logging

# Suppress annoying warnings while keeping functionality
warnings.filterwarnings('ignore', category=UserWarning, module='anndata')
logging.getLogger('scanpy').setLevel(logging.ERROR)

sys.modules['importlib.metadata'] = importlib_metadata

def read_samples(file_path):
    samples = {}
    with open(file_path, 'r') as file:
        for line in file:
            sample_id = line.strip()
            if sample_id:
                filename = f"{sample_id}_filtered_feature_bc_matrix.h5"
                samples[sample_id] = filename
    return samples

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('obj_name')
    parser.add_argument('inputfile')
    args = parser.parse_args()

    inputfile = args.inputfile
    obj_name = args.obj_name

    adatas = {}
    samples = read_samples(inputfile)
    print(samples)
    
    for sample_id, filename in samples.items():
        print(f"Reading {filename}...")
        
        # Read with suppressed warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata = sc.read_10x_h5(filename)
        
        # CRITICAL: Fix duplicate names immediately after reading
        adata.var_names_make_unique()  # Fixes duplicate gene names
        adata.obs_names_make_unique()  # Ensures unique cell barcodes
        
        adata.obs['sample'] = sample_id
        adatas[sample_id] = adata

    # Concatenate all samples
    combined_adata = sc.concat(adatas.values(), label='sample', keys=adatas.keys())
    
    # CRITICAL: Fix duplicate names after concatenation
    combined_adata.var_names_make_unique()
    combined_adata.obs_names_make_unique()

    # Calculate QC metrics
    combined_adata.var["mt"] = combined_adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
         combined_adata, qc_vars=["mt"], inplace=True, log1p=True
    )

    # Plot QC with suppressed warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.violin(
             combined_adata,
             ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
             jitter=0.4,
             multi_panel=True, 
             save="_QC.png"
        )

    # Filter cells with explicit copy to avoid view warnings
    # Less stringent filtering of cells
    combined_adata = combined_adata[
    (combined_adata.obs['n_genes_by_counts'] > 200) &    # keep cells with at least 200 genes
    (combined_adata.obs['n_genes_by_counts'] < 8000) &   # remove extreme high gene counts (doublets)
    (combined_adata.obs['total_counts'] > 500) &         # keep cells with at least 500 UMIs
    (combined_adata.obs['total_counts'] < 30000) &       # remove extreme high UMI counts
    (combined_adata.obs['pct_counts_mt'] < 20),          # remove cells with >20% mitochondrial counts
    :
    ].copy()


    # Basic filtering
    sc.pp.filter_cells(combined_adata, min_genes=100)
    sc.pp.filter_genes(combined_adata, min_cells=3)

    # Plot after QC with suppressed warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.violin(
            combined_adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True, 
            save="_AfterQC.png"
        )

    # Save the final object
    filename = obj_name + ".h5ad"
    combined_adata.write(filename)
    print(f"Successfully created {filename}")

if __name__ == "__main__":
    main()
