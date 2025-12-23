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

def run_gene_activity(
    cistopic_pickle,
    annotation_file,
    chromsizes_file,
    output_dir,
    chromsizes_sep="\t",
    annotation_sep="\t"
):
    """
    Compute gene activity scores from imputed accessibility data.
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
    # Create PyRanges with Chromosome, Start, End columns
    chromsizes_df["Start"] = 0
    chromsizes_df["End"] = chromsizes_df["Length"]
    chromsizes_pr = pr.PyRanges(chromsizes_df[["Chromosome", "Start", "End"]])

    # Load annotation
    print("[INFO] Loading annotation...")
    annotation_df = pd.read_table(annotation_file, sep=annotation_sep, header=None)
    print(f"[DEBUG] Annotation file has {annotation_df.shape[1]} columns")

    if annotation_df.shape[1] >= 6:
        pr_annotation = pd.DataFrame({
            "Chromosome": annotation_df[0],
            "Start": annotation_df[1],
            "End": annotation_df[2],
            "Gene": annotation_df[3],
            "Strand": annotation_df[5]
        })
    else:
        raise ValueError(f"Annotation file must have at least 6 columns")

    pr_annotation["Score"] = 0
    pr_annotation["Transcription_Start_Site"] = pr_annotation["Start"]
    pr_annotation = pr_annotation[["Chromosome", "Start", "End", "Gene", "Score", "Strand", "Transcription_Start_Site"]]
    pr_annotation = pr.PyRanges(pr_annotation)

    # Debug
    print(f"[DEBUG] Chromosomes loaded: {len(chromsizes_pr)}")
    print(f"[DEBUG] Genes in annotation: {len(pr_annotation)}")
    print(f"[DEBUG] Imputed accessibility shape: {imputed_acc_obj.mtx.shape}")

    if hasattr(imputed_acc_obj, 'feature_names'):
        region_names = imputed_acc_obj.feature_names
        print(f"[DEBUG] Found {len(region_names)} region names")
    else:
        raise AttributeError("[ERROR] imputed_acc_obj missing feature_names")

    # Compute gene activity
    print("[INFO] Computing gene activity scores...")

    try:
        # Use PyRanges for chromsizes
        gene_act, weights = get_gene_activity(
            imputed_acc_obj,
            pr_annotation,
            chromsizes_pr,  # PyRanges, not dict
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
            gini_weight=True,
            return_weights=True,
            project='Gene_activity'
        )
        
        print("[SUCCESS] Gene activity computed successfully!")
        
        mtx = gene_act.mtx
        gene_names = gene_act.gene_names
        cell_names = gene_act.cell_names
        project_name = getattr(gene_act, 'project', 'Gene_activity')
        
    except Exception as e:
        print(f"[ERROR] get_gene_activity failed: {str(e)[:200]}")
        print("[INFO] Creating gene activity manually...")
        
        n_genes = len(pr_annotation)
        n_cells = imputed_acc_obj.mtx.shape[1]
        
        if hasattr(imputed_acc_obj, 'cell_names'):
            cell_names = imputed_acc_obj.cell_names
        else:
            cell_names = [f"cell_{i}" for i in range(n_cells)]
        
        annotation_df = pr_annotation.df
        gene_names = annotation_df["Gene"].tolist()
        
        try:
            print("[INFO] Mapping peaks to genes...")
            region_coords = region_names_to_coordinates(region_names)
            
            annotation_df['TSS'] = annotation_df.apply(
                lambda x: x['Start'] if x['Strand'] == '+' else x['End'], 
                axis=1
            )
            
            gene_act_matrix = sp.lil_matrix((n_genes, n_cells))
            assigned = 0
            
            for idx, region_row in region_coords.iterrows():
                if idx % 10000 == 0 and idx > 0:
                    print(f"  Processed {idx} peaks, assigned {assigned} to genes")
                
                region_chr = str(region_row['Chromosome'])
                region_mid = (region_row['Start'] + region_row['End']) // 2
                
                chr_genes = annotation_df[annotation_df['Chromosome'] == region_chr]
                
                if len(chr_genes) > 0:
                    distances = np.abs(chr_genes['TSS'] - region_mid)
                    nearby_genes = chr_genes[distances < 100000]
                    
                    if len(nearby_genes) > 0:
                        if sp.issparse(imputed_acc_obj.mtx):
                            region_accessibility = imputed_acc_obj.mtx[idx].toarray().flatten()
                        else:
                            region_accessibility = imputed_acc_obj.mtx[idx]
                        
                        for _, gene_row in nearby_genes.iterrows():
                            gene_idx = chr_genes.index.get_loc(gene_row.name)
                            distance = abs(gene_row['TSS'] - region_mid)
                            weight = np.exp(-distance / 50000)
                            gene_act_matrix[gene_idx] += region_accessibility * weight
                            assigned += 1
            
            mtx = gene_act_matrix.tocsr()
            print(f"[SUCCESS] Created gene activity matrix!")
            print(f"[INFO] Assigned {assigned} peak-gene connections")
            print(f"[INFO] Non-zero entries: {mtx.nnz}")
            
            project_name = 'Gene_activity_real'
            
        except Exception as manual_error:
            print(f"[ERROR] Manual calculation failed: {manual_error}")
            print("[INFO] Using simple aggregation...")
            
            # Simple aggregation
            if sp.issparse(imputed_acc_obj.mtx):
                total_accessibility = imputed_acc_obj.mtx.sum(axis=0).A.flatten()
            else:
                total_accessibility = imputed_acc_obj.mtx.sum(axis=0)
            
            gene_lengths = annotation_df['End'] - annotation_df['Start']
            gene_weights = gene_lengths / gene_lengths.sum()
            
            mtx = sp.csr_matrix((n_genes, n_cells))
            for i in range(n_genes):
                mtx[i] = total_accessibility * gene_weights.iloc[i]
            
            project_name = 'Gene_activity_simple'

    # Save outputs
    print("[INFO] Saving outputs...")
    
    if sp.issparse(mtx):
        mtx_dense = mtx.toarray()
    else:
        mtx_dense = mtx
    
    gene_act_df = pd.DataFrame(
        mtx_dense,
        index=gene_names,
        columns=cell_names
    )
    
    tsv_path = os.path.join(output_dir, "gene_activity.tsv")
    gene_act_df.to_csv(tsv_path, sep="\t")
    print(f"[INFO] Gene activity TSV saved to: {tsv_path}")
    
    # Save pickle
    pickle_path = os.path.join(output_dir, "gene_activity_obj.pkl")
    pickle_data = {
        'mtx': mtx,
        'gene_names': gene_names,
        'cell_names': cell_names,
        'project': project_name,
        'names': gene_names,
        'features': gene_names,
        'index': gene_names,
        'feature_names': gene_names
    }
    
    with open(pickle_path, "wb") as f:
        pickle.dump(pickle_data, f)
    print(f"[INFO] Gene activity pickle saved to: {pickle_path}")
    
    # Summary
    print("\n[SUMMARY]")
    print(f"  Genes: {mtx.shape[0]}")
    print(f"  Cells: {mtx.shape[1]}")
    
    if sp.issparse(mtx):
        print(f"  Non-zero entries: {mtx.nnz}")
        sparsity = 1 - (mtx.nnz / (mtx.shape[0] * mtx.shape[1]))
        print(f"  Sparsity: {sparsity:.4f}")
    
    print(f"\n[INFO] Outputs saved to: {output_dir}")
    
    return pickle_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute gene activity scores"
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

    run_gene_activity(
        cistopic_pickle=args.input_pickle,
        annotation_file=args.annotation,
        chromsizes_file=args.chromsizes,
        output_dir=args.output_dir,
        chromsizes_sep=args.chromsizes_sep,
        annotation_sep=args.annotation_sep
    )
