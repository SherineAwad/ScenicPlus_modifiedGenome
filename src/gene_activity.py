#!/usr/bin/env python3
import argparse
import os
import pickle
import numpy as np
import pandas as pd
import scipy.sparse as sp
import pyranges as pr
from pycisTopic.gene_activity import get_gene_activity
from pycisTopic.utils import region_names_to_coordinates

def run_gene_activity_simple(
    cistopic_pickle,
    annotation_file,
    chromsizes_file,
    output_dir,
    chromsizes_sep="\t",
    annotation_sep="\t"
):
    """
    Compute gene activity scores from imputed accessibility data.
    Simple version following SCENIC+ tutorial.
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load cistopic object
    print("[INFO] Loading cistopic object...")
    with open(cistopic_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)
    
    if not hasattr(cistopic_obj, 'imputed_acc_obj'):
        raise AttributeError("[ERROR] cistopic_obj missing imputed_acc_obj attribute.")
    
    imputed_acc_obj = cistopic_obj.imputed_acc_obj
    
    # Load chromosome sizes
    print("[INFO] Loading chromosome sizes...")
    chromsizes_df = pd.read_csv(
        chromsizes_file,
        sep=chromsizes_sep,
        header=None,
        names=["Chromosome", "Length"]
    )
    
    # Convert chromsizes to dict for get_gene_activity
    chromsizes = dict(zip(chromsizes_df["Chromosome"], chromsizes_df["Length"]))
    
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
    
    # THE FIX: Disable gini_weight to make it run faster
    print("[INFO] Computing gene activity scores (without Gini weights for speed)...")
    
    try:
        # Use the tutorial approach but without Gini weights
        gene_act, weights = get_gene_activity(
            imputed_acc_obj,
            pr_annotation,
            chromsizes,
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
            gini_weight=False,  # <-- THE KEY CHANGE: Disable Gini weights
            return_weights=True,
            project='Gene_activity'
        )
        
        print("[SUCCESS] Gene activity computed successfully!")
        
    except Exception as e:
        print(f"[ERROR] get_gene_activity failed: {e}")
        print("[INFO] Trying an even simpler approach...")
        
        # Minimal approach if the above still fails
        gene_act, weights = get_gene_activity(
            imputed_acc_obj,
            pr_annotation,
            chromsizes,
            use_gene_boundaries=True,
            upstream=[1000, 100000],
            downstream=[1000, 100000],
            distance_weight=True,  # Keep distance weights
            decay_rate=1,
            extend_gene_body_upstream=10000,
            extend_gene_body_downstream=500,
            gene_size_weight=False,
            remove_promoters=False,
            average_scores=True,
            gini_weight=False,  # Definitely disable Gini
            return_weights=True,
            project='Gene_activity'
        )
    
    # Save outputs
    print("[INFO] Saving outputs...")
    
    # Save gene activity matrix as TSV
    gene_act_df = pd.DataFrame(
        gene_act.mtx.toarray() if sp.issparse(gene_act.mtx) else gene_act.mtx,
        index=gene_act.gene_names,
        columns=gene_act.cell_names
    )
    
    tsv_path = os.path.join(output_dir, "gene_activity.tsv")
    gene_act_df.to_csv(tsv_path, sep="\t")
    print(f"[INFO] Gene activity TSV saved to: {tsv_path}")
    
    # Save pickle
    pickle_path = os.path.join(output_dir, "gene_activity_obj.pkl")
    with open(pickle_path, "wb") as f:
        pickle.dump(gene_act, f)
    print(f"[INFO] Gene activity pickle saved to: {pickle_path}")
    
    # Save weights if you want them
    weights_path = os.path.join(output_dir, "gene_activity_weights.pkl")
    with open(weights_path, "wb") as f:
        pickle.dump(weights, f)
    
    print("\n[SUMMARY]")
    print(f"  Genes: {gene_act.mtx.shape[0]}")
    print(f"  Cells: {gene_act.mtx.shape[1]}")
    if sp.issparse(gene_act.mtx):
        print(f"  Non-zero entries: {gene_act.mtx.nnz}")
    print(f"\n[INFO] Outputs saved to: {output_dir}")
    
    return gene_act

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute gene activity scores (simple tutorial version)"
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
    
    run_gene_activity_simple(
        cistopic_pickle=args.input_pickle,
        annotation_file=args.annotation,
        chromsizes_file=args.chromsizes,
        output_dir=args.output_dir,
        chromsizes_sep=args.chromsizes_sep,
        annotation_sep=args.annotation_sep
    )
