import os
import pandas as pd
import pycisTopic
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk

print("pycisTopic version:", pycisTopic.__version__)

# -----------------------------
# 1. Define fragments dictionary
# -----------------------------
fragments_dict = {
    "Control": "TH1_atac_fragments.tsv.gz",
    "OE": "TH2_atac_fragments.tsv.gz"
}

# -----------------------------
# 2. Load chromsizes FROM YOUR CUSTOM GENOME
# -----------------------------
# FIRST, create chrom.sizes file from your neurog2.fa:
# Run this command in terminal: 
#   samtools faidx neurog2.fa
#   cut -f1,2 neurog2.fa.fai > neurog2.chrom.sizes

chromsizes = pd.read_table(
    "neurog2.chrom.sizes",  # <-- CHANGE TO YOUR FILE
    header=None,
    names=["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
print("Chromsizes head:")
print(chromsizes.head())

# -----------------------------
# 3. Create output directories
# -----------------------------
out_dir = "scenicOuts"
consensus_dir = os.path.join(out_dir, "consensus_peak_calling")
bed_dir       = os.path.join(consensus_dir, "pseudobulk_bed_files")
bw_dir        = os.path.join(consensus_dir, "pseudobulk_bw_files")

os.makedirs(consensus_dir, exist_ok=True)
os.makedirs(bed_dir, exist_ok=True)
os.makedirs(bw_dir, exist_ok=True)

# -----------------------------
# 4. Create cell data from fragments files
# -----------------------------
cell_data_list = []

for sample, fragments_file in fragments_dict.items():
    print(f"Reading barcodes from {fragments_file}...")

    if not os.path.exists(fragments_file):
        print(f"  ERROR: File not found!")
        continue

    try:
        fragments_df = pd.read_csv(
            fragments_file,
            sep='\t',
            header=None,
            comment='#',
            names=['chrom', 'start', 'end', 'barcode', 'count'],
            usecols=[0, 1, 2, 3, 4]
        )

        print(f"  Data shape: {fragments_df.shape}")

        barcodes = fragments_df['barcode'].unique()
        print(f"  Found {len(barcodes)} unique barcodes")

        sample_df = pd.DataFrame({
            'barcode': barcodes.astype(str),
            'sample': str(sample),
            'celltype': str(sample)
        })

        cell_data_list.append(sample_df)

    except Exception as e:
        print(f"  Error processing {fragments_file}: {e}")

# Check if we have any data
if len(cell_data_list) == 0:
    print("\nERROR: No barcode data found!")
    exit(1)

# Combine all samples
cell_data = pd.concat(cell_data_list, ignore_index=True)
cell_data_simple = cell_data[['barcode', 'sample', 'celltype']].copy()

for col in ['barcode', 'sample', 'celltype']:
    cell_data_simple[col] = cell_data_simple[col].astype(str)

print(f"\nFinal cell data summary:")
print(f"Total barcodes: {len(cell_data_simple)}")
print("Sample distribution:")
print(cell_data_simple['sample'].value_counts())

# -----------------------------
# 5. Run pseudobulk export
# -----------------------------
print("\nStarting pseudobulk export...")

try:
    bw_paths, bed_paths = export_pseudobulk(
        input_data=cell_data_simple,
        variable="celltype",
        sample_id_col="sample",
        chromsizes=chromsizes,  # <-- NOW USING YOUR MODIFIED GENOME
        bed_path=bed_dir,
        bigwig_path=bw_dir,
        path_to_fragments=fragments_dict,
        n_cpu=4,
        normalize_bigwig=True,
        temp_dir="/tmp"
    )

    print("Pseudobulk export completed successfully!")

    # Save paths
    with open(os.path.join(consensus_dir, "bw_paths.tsv"), "w") as f:
        for v in bw_paths:
            f.write(f"{v}\t{bw_paths[v]}\n")

    with open(os.path.join(consensus_dir, "bed_paths.tsv"), "w") as f:
        for v in bed_paths:
            f.write(f"{v}\t{bed_paths[v]}\n")

    print(f"\nOutput directories:")
    print(f"  - BED files: {bed_dir}")
    print(f"  - BigWig files: {bw_dir}")

except Exception as e:
    print(f"Error during pseudobulk export: {e}")
