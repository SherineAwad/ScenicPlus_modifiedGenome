#!/usr/bin/env python3
import argparse
import os
import pickle
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pycisTopic.diff_features import find_diff_features
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add the subset method to this class
class GeneActivityObject:
    def __init__(self, mtx, gene_names, cell_names, project='Gene_activity'):
        self.mtx = mtx
        self.gene_names = gene_names
        self.cell_names = cell_names
        self.project = project
        self.names = gene_names
        self.features = gene_names
        self.index = gene_names
        self.df = self._create_df()
    
    def _create_df(self):
        if sp.issparse(self.mtx):
            return pd.DataFrame(self.mtx.toarray(), index=self.gene_names, columns=self.cell_names)
        else:
            return pd.DataFrame(self.mtx, index=self.gene_names, columns=self.cell_names)
    
    # ADD THIS METHOD
    def subset(self, cells=None, features=None, copy=False, split_pattern=None):    
        if cells is not None:
            cell_indices = [i for i, cell in enumerate(self.cell_names) if cell in cells]
        else:
            cell_indices = list(range(len(self.cell_names)))
        
        if features is not None:
            feature_indices = [i for i, gene in enumerate(self.gene_names) if gene in features]
        else:
            feature_indices = list(range(len(self.gene_names)))
        
        if sp.issparse(self.mtx):
            subset_mtx = self.mtx[feature_indices, :][:, cell_indices]
        else:
            subset_mtx = self.mtx[feature_indices, :][:, cell_indices]
        
        subset_gene_names = [self.gene_names[i] for i in feature_indices]
        subset_cell_names = [self.cell_names[i] for i in cell_indices]
        
        new_obj = GeneActivityObject(
            mtx=subset_mtx,
            gene_names=subset_gene_names,
            cell_names=subset_cell_names,
            project=self.project
        )
        
        return new_obj

def load_inputs(cistopic_pickle, gene_activity_pickle):
    print("[INFO] Loading cistopic object...")
    with open(cistopic_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)
    
    print("[INFO] Loading gene activity object...")
    with open(gene_activity_pickle, "rb") as f:
        gene_act_data = pickle.load(f)
    
    if isinstance(gene_act_data, dict):
        print("[INFO] Converting dict to GeneActivityObject...")
        gene_act = GeneActivityObject(
            mtx=gene_act_data['mtx'],
            gene_names=gene_act_data['gene_names'],
            cell_names=gene_act_data['cell_names'],
            project=gene_act_data.get('project', 'Gene_activity')
        )
    else:
        gene_act = gene_act_data
    
    return cistopic_obj, gene_act

def run_differential_analysis(cistopic_obj, gene_act, output_dir, n_cpu=5, temp_dir=None, celltype_column='celltype'):
    print("\n" + "="*60)
    print("[INFO] Running differential gene activity analysis...")
    print("="*60)
    
    if temp_dir:
        os.environ['RAY_TEMPDIR'] = temp_dir
        os.makedirs(temp_dir, exist_ok=True)
    
    if not hasattr(cistopic_obj, 'cell_data'):
        print("[ERROR] cistopic_obj missing cell_data attribute")
        return None
    
    if celltype_column not in cistopic_obj.cell_data.columns:
        print(f"[ERROR] '{celltype_column}' column not found")
        print("[INFO] Available columns:")
        print(cistopic_obj.cell_data.columns.tolist())
        return None
    
    print(f"[INFO] Using column '{celltype_column}' for cell types")
    print(f"[INFO] Found cell types: {cistopic_obj.cell_data[celltype_column].unique().tolist()}")
    
    try:
        DAG_markers_dict = find_diff_features(
            cistopic_obj,
            gene_act,
            variable=celltype_column,
            var_features=None,
            contrasts=None,
            adjpval_thr=0.05,
            log2fc_thr=np.log2(1.5),
            n_cpu=n_cpu,
            _temp_dir=temp_dir,
            split_pattern='-'
        )
        
        print("\n[SUCCESS] Differential analysis completed!")
        print("\nNumber of DAGs found:")
        print("-" * 30)
        total_dags = 0
        for cell_type in DAG_markers_dict:
            n_markers = len(DAG_markers_dict[cell_type])
            print(f"  {cell_type}: {n_markers}")
            total_dags += n_markers
        
        print(f"\nTotal DAGs across all cell types: {total_dags}")
        
        dag_output = os.path.join(output_dir, "DAG_markers.pkl")
        with open(dag_output, "wb") as f:
            pickle.dump(DAG_markers_dict, f)
        print(f"[INFO] DAG markers saved to: {dag_output}")
        
        dag_text = os.path.join(output_dir, "DAG_markers_summary.tsv")
        with open(dag_text, "w") as f:
            f.write("Cell_type\tNumber_of_DAGs\tDAG_genes\n")
            for cell_type in DAG_markers_dict:
                genes = DAG_markers_dict[cell_type]
                f.write(f"{cell_type}\t{len(genes)}\t{','.join(genes)}\n")
        print(f"[INFO] DAG summary saved to: {dag_text}")
        
        return DAG_markers_dict
        
    except Exception as e:
        print(f"[ERROR] Differential analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    parser = argparse.ArgumentParser(description="Run differential gene activity analysis")
    
    parser.add_argument("-i", "--cistopic_pickle", required=True, help="Path to cistopic_obj pickle file")
    parser.add_argument("-g", "--gene_activity_pickle", required=True, help="Path to gene activity pickle file")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save all outputs")
    parser.add_argument("--n_cpu", type=int, default=5, help="Number of CPUs for differential analysis")
    parser.add_argument("--temp_dir", help="Temporary directory for ray")
    parser.add_argument("--celltype_column", default="celltype", help="Column name for cell types")
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    cistopic_obj, gene_act = load_inputs(args.cistopic_pickle, args.gene_activity_pickle)
    
    dag_markers = run_differential_analysis(
        cistopic_obj, 
        gene_act, 
        args.output_dir, 
        n_cpu=args.n_cpu,
        temp_dir=args.temp_dir,
        celltype_column=args.celltype_column
    )
    
    if dag_markers is not None:
        print("\n" + "="*60)
        print("[SUCCESS] Analysis completed!")
        print(f"[INFO] DAG markers saved to: {args.output_dir}/DAG_markers.pkl")
        print("="*60)
    else:
        print("\n" + "="*60)
        print("[WARNING] Differential analysis could not be completed")
        print("[INFO] Creating empty DAG markers file...")
        
        dag_output = os.path.join(args.output_dir, "DAG_markers.pkl")
        empty_dict = {}
        with open(dag_output, "wb") as f:
            pickle.dump(empty_dict, f)
        print(f"[INFO] Empty DAG markers saved to: {dag_output}")
        print("="*60)

if __name__ == "__main__":
    main()
