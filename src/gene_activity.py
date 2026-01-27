#!/usr/bin/env python3
import argparse
import os
import pickle
import numpy as np
import pandas as pd
import scipy.sparse as sp
import pyranges as pr
from pycisTopic.gene_activity import get_gene_activity

def run_gene_activity_fixed(
    cistopic_pickle,
    annotation_file,
    chromsizes_file,
    output_dir,
    chromsizes_sep="\t",
    annotation_sep="\t"
):
    """
    Compute gene activity scores from imputed accessibility data.
    FIXED VERSION: Handles both sparse and dense matrices.
    """

    os.makedirs(output_dir, exist_ok=True)

    # Load cistopic object
    print("[INFO] Loading cistopic object...")
    with open(cistopic_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    if not hasattr(cistopic_obj, 'imputed_acc_obj'):
        raise AttributeError("[ERROR] cistopic_obj missing imputed_acc_obj attribute.")

    imputed_acc_obj = cistopic_obj.imputed_acc_obj

    # Load chromosome sizes as PyRanges
    print("[INFO] Loading chromosome sizes...")
    chromsizes_df = pd.read_csv(
        chromsizes_file,
        sep=chromsizes_sep,
        header=None,
        names=["Chromosome", "Length"]
    )

    # Create PyRanges object with Chromosome, Start, End
    chromsizes_pr = pr.PyRanges(
        chromsizes_df.assign(
            Start=0,
            End=chromsizes_df["Length"]
        )[["Chromosome", "Start", "End"]]
    )

    # Load annotation
    print("[INFO] Loading annotation...")
    annotation_df = pd.read_table(annotation_file, sep=annotation_sep, header=None)

    if annotation_df.shape[1] < 6:
        raise ValueError(f"Annotation file must have at least 6 columns, got {annotation_df.shape[1]}")

    # Format annotation for PyRanges
    pr_annotation = pd.DataFrame({
        "Chromosome": annotation_df[0].astype(str),
        "Start": annotation_df[1],
        "End": annotation_df[2],
        "Gene": annotation_df[3],
        "Strand": annotation_df[5]
    })

    # Add required columns
    pr_annotation["Score"] = 0
    pr_annotation["Transcription_Start_Site"] = pr_annotation["Start"]
    pr_annotation = pr_annotation[["Chromosome", "Start", "End", "Gene", "Score", "Strand", "Transcription_Start_Site"]]
    pr_annotation = pr.PyRanges(pr_annotation)

    print(f"[INFO] Processing {len(pr_annotation)} genes")
    print(f"[INFO] Imputed accessibility shape: {imputed_acc_obj.mtx.shape}")

    # Compute gene activity WITHOUT Gini weights for speed
    print("[INFO] Computing gene activity scores (without Gini weights for speed)...")

    try:
        gene_act, weights = get_gene_activity(
            imputed_acc_obj,
            pr_annotation,
            chromsizes_pr,
            use_gene_boundaries=True,
            upstream=[1000, 100000],
            downstream=[1000, 100000],
            distance_weight=True,
            decay_rate=1,
            extend_gene_body_upstream=10000,
            extend_gene_body_downstream=500,
            gene_size_weight=False,
            gene_size_scale_factor='median',
            remove_promoters=False,
            average_scores=True,
            scale_factor=1,
            extend_tss=[10, 10],
            gini_weight=False,  # DISABLED for speed
            return_weights=True,
            project='Gene_activity'
        )

        print("[SUCCESS] Gene activity computed successfully!")

    except Exception as e:
        print(f"[ERROR] Failed to compute gene activity: {e}")
        raise

    # Save outputs
    print("[INFO] Saving outputs...")

    # Determine the correct attribute names
    if hasattr(gene_act, 'feature_names'):
        feature_names = gene_act.feature_names
    elif hasattr(gene_act, 'gene_names'):
        feature_names = gene_act.gene_names
    else:
        raise AttributeError("Could not find feature names attribute in gene_act object")
    
    cell_names = gene_act.cell_names
    
    print(f"[INFO] Features: {len(feature_names)}")
    print(f"[INFO] Cells: {len(cell_names)}")

    # Save gene activity matrix as TSV - optimized for large matrices
    print("[INFO] Converting matrix to DataFrame...")
    
    # Check if matrix is sparse or dense
    if sp.issparse(gene_act.mtx):
        print(f"[INFO] Sparse matrix with {gene_act.mtx.nnz} non-zero entries")
        # For sparse matrices
        gene_act_df = pd.DataFrame.sparse.from_spmatrix(
            gene_act.mtx,
            index=feature_names,
            columns=cell_names
        )
        is_sparse = True
    else:
        print("[INFO] Dense matrix")
        gene_act_df = pd.DataFrame(
            gene_act.mtx,
            index=feature_names,
            columns=cell_names
        )
        is_sparse = False

    tsv_path = os.path.join(output_dir, "gene_activity.tsv")
    
    # Save in chunks if matrix is very large
    if gene_act_df.shape[0] > 10000 or gene_act_df.shape[1] > 10000:
        print("[INFO] Large matrix detected, saving in chunks...")
        chunk_size = 1000
        with open(tsv_path, 'w') as f:
            # Write header
            f.write("Gene\t" + "\t".join(map(str, cell_names)) + "\n")
            
            # Write data in chunks
            total_chunks = (len(feature_names) - 1) // chunk_size + 1
            for i in range(0, len(feature_names), chunk_size):
                chunk_idx = i // chunk_size + 1
                print(f"[INFO] Saving chunk {chunk_idx}/{total_chunks}")
                chunk = gene_act_df.iloc[i:i+chunk_size]
                chunk.to_csv(f, sep='\t', header=False, mode='a')
    else:
        gene_act_df.to_csv(tsv_path, sep="\t")
    
    print(f"[INFO] Gene activity TSV saved to: {tsv_path}")

    # Save pickle
    pickle_path = os.path.join(output_dir, "gene_activity_obj.pkl")
    with open(pickle_path, "wb") as f:
        pickle.dump(gene_act, f)
    print(f"[INFO] Gene activity pickle saved to: {pickle_path}")

    # Save weights if they exist
    if weights is not None:
        weights_path = os.path.join(output_dir, "gene_activity_weights.pkl")
        with open(weights_path, "wb") as f:
            pickle.dump(weights, f)
        print(f"[INFO] Weights saved to: {weights_path}")

    print("\n[SUMMARY]")
    print(f"  Genes/Features: {gene_act.mtx.shape[0]}")
    print(f"  Cells: {gene_act.mtx.shape[1]}")
    
    # Calculate sparsity correctly for both sparse and dense matrices
    if is_sparse:
        total_elements = gene_act.mtx.shape[0] * gene_act.mtx.shape[1]
        non_zero = gene_act.mtx.nnz
        sparsity = 100 * (1 - non_zero / total_elements)
        print(f"  Non-zero entries: {non_zero}")
        print(f"  Sparsity: {sparsity:.2f}%")
    else:
        # For dense matrices, calculate non-zero elements
        non_zero = np.count_nonzero(gene_act.mtx)
        total_elements = gene_act.mtx.size
        sparsity = 100 * (1 - non_zero / total_elements)
        print(f"  Non-zero entries: {non_zero}")
        print(f"  Sparsity: {sparsity:.2f}%")
    
    print(f"\n[INFO] Outputs saved to: {output_dir}")

    return gene_act

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute gene activity scores (FIXED version)"
    )

    parser.add_argument(
        "-i", "--input_pickle",
        required=True,
        help="Path to cistopic_obj pickle file"
    )

    parser.add_argument(
        "-a", "--annotation",
        required=True,
        help="Path to gene annotation file"
    )

    parser.add_argument(
        "-c", "--chromsizes",
        required=True,
        help="Path to chromosome sizes file"
    )

    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="Directory to save outputs"
    )

    parser.add_argument(
        "--chromsizes_sep",
        default="\t",
        help="Separator for chromosome sizes file"
    )

    parser.add_argument(
        "--annotation_sep",
        default="\t",
        help="Separator for annotation file"
    )

    args = parser.parse_args()

    run_gene_activity_fixed(
        cistopic_pickle=args.input_pickle,
        annotation_file=args.annotation,
        chromsizes_file=args.chromsizes,
        output_dir=args.output_dir,
        chromsizes_sep=args.chromsizes_sep,
        annotation_sep=args.annotation_sep
    )
