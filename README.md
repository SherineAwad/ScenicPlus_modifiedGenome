# 🌟 Scenic+ GRN Overview

This is **Scenic+ GRN** based on the **modified genome `neurog2_S9A`**.  

For a full explanation of **Scenic+** using a custom genome, refer to our pipeline: [multiomicsGRN](https://github.com/SherineAwad/multiomicsGRN).


# 🔧 Modifying cisTarget for a Custom Genome (`neurog2_S9A`)

This workflow describes how to adapt **cisTarget** to a **modified genome** for use with **Scenic+**:

1. **Create a cisTarget database**: Download motif collections and prepare them for motif enrichment analysis.  
2. **Prepare the modified genome**: Generate a reference genome FASTA and annotation files, then define genomic windows for background sequences.  
3. **Generate padded sequences**: Create a padded FASTA to serve as input for cisTarget, ensuring motifs are evaluated in context.  
4. **Run motif scanning**: Use tools like Cluster-Buster to scan the genome with the motif collection.  
5. **Build the custom cisTarget database**: Combine the padded genome, motif collection, and scanning results to produce a database tailored to the modified genome.  

The resulting **custom cisTarget database** can then be used in **Scenic+** to perform **gene regulatory network inference** specifically for the modified genome `neurog2_S9A`.  

All steps are automated and reproducible using the [Makefile](https://github.com/SherineAwad/ScenicPlus_modifiedGenome/blob/master/Makefile).



# How cisTarget DB Creation Includes the New Gene 

### Step 1: Our FASTA contains the gene

The FASTA passed to cisTarget **contains our new gene sequence**, e.g.:

```fasta
>GENE_X_custom
ATGCGTACGATCG...
```

✔ At this moment, cisTarget **knows this sequence exists**.

---

### Step 2: FASTA entries become “regions”

cisTarget treats each FASTA record as **one independent region**:

```text
Region_ID = GENE_X_custom
Sequence  = ATGCGTACGATCG...
```

✔ The new gene is now a **region** in the database.

---

### Step 3: Motif scanning on each region

For every motif:

```text
for motif in motifs:
    for region in fasta_entries:
        score(motif, region.sequence)
```

✔ The new gene sequence is scanned
✔ Motif scores are computed
✔ No filtering or skipping occurs

#### ⚡⚡⚡ How cisTarget Scores a New Gene Region

1. **Region definition**: The new gene (or its 1kb windows) exists as a FASTA entry. cisTarget treats each entry as a “region” with a sequence.

2. **Motif scanning**: For each motif in the collection, cisTarget uses Cluster-Buster (CBUST) to slide the PWM along the region, calculate match scores at each position, and integrate them across the region.

3. **Score assignment**: CBUST produces a single numeric score per motif per region, representing how strongly the motif matches the sequence.

4. **Score storage**: Each motif score is stored in the cisTarget database as `region_id (the gene) → motif_id → score`. Every motif is scored independently, and no enrichment or interpretation occurs at this stage.


---

### Step 4: Scores are written to the database

For each motif, cisTarget stores:

```text
motif_id → list of (region_id, score)
```

Example for our gene:

```text
motif_X → (GENE_X_custom, 0.87)
motif_Y → (GENE_X_custom, 0.12)
```

✔ Motif scores for our gene exist on disk.

---

### Step 5: Ranking conversion 

Scores are sorted per motif:

```text
motif_X:
  region_1
  region_2
  ...
  GENE_X_custom   ← some rank
```

✔ The new gene has a rank for each motif
✔ No biological interpretation
✔ No enrichment logic


> **If a new gene sequence is present as a FASTA entry in the input FASTA, cisTarget DB creation WILL include it, scan it, compute motif scores, and store those scores in the database.**


### Critical concept (this is the key point)

**Motif scanning is associated with DNA sequences, NOT with genes.**

- cisTarget scans **raw DNA sequence**.
- Motifs are detected in **genomic regions (windows)**.
- Genes are **not used or referenced** during motif scanning.

If a DNA sequence exists in the FASTA and is scanned, motifs in that sequence **can be detected**, regardless of whether the sequence corresponds to a known, annotated, endogenous, or synthetic gene.

---

## What rebuilding cisTarget does NOT do

- ❌ It does **not** use gene names.
- ❌ It does **not** use the GTF annotation.
- ❌ It does **not** require genes to be annotated or “known”.

Gene annotations are only used **later**, during downstream steps, when genomic regions are assigned to genes.

----- 

# scRNA preprocessing 

##### Filtering 
```
combined_adata = combined_adata[
    (combined_adata.obs['n_genes_by_counts'] > 200) &    # keep cells with at least 200 genes
    (combined_adata.obs['n_genes_by_counts'] < 8000) &   # remove extreme high gene counts (doublets)
    (combined_adata.obs['total_counts'] > 500) &         # keep cells with at least 500 UMIs
    (combined_adata.obs['total_counts'] < 30000) &       # remove extreme high UMI counts
    (combined_adata.obs['pct_counts_mt'] < 20),          # remove cells with >20% mitochondrial counts
    :
    ].copy()

``` 

##### Before filtering 

![Before filtering](figures/violin_QC.png?v=2)

##### After filtering 
![After filtering](figures/violin_AfterQC.png?v=2)

<img src="figures/umap_clustered_Neurog2_Abca8a.png?v=2" alt="Abca8a" width="33%"><img src="figures/umap_clustered_Neurog2_Chat.png?v=2" alt="Chat" width="33%"><img src="figures/umap_clustered_Neurog2_Insm1.png?v=2" alt="Insm1" width="33%">
<img src="figures/umap_clustered_Neurog2_Nrl.png?v=2" alt="Nrl" width="33%"><img src="figures/umap_clustered_Neurog2_Notch1.png?v=2" alt="Notch1" width="33%"><img src="figures/umap_clustered_Neurog2_Acta2.png?v=2" alt="Acta2" width="33%">
<img src="figures/umap_clustered_Neurog2_Rpe65.png?v=2" alt="Rpe65" width="33%"><img src="figures/umap_clustered_Neurog2_Isl1.png?v=2" alt="Isl1" width="33%"><img src="figures/umap_clustered_Neurog2_Olig2.png?v=2" alt="Olig2" width="33%">
<img src="figures/umap_clustered_Neurog2_Sebox.png?v=2" alt="Sebox" width="33%"><img src="figures/umap_clustered_Neurog2_Apoe.png?v=2" alt="Apoe" width="33%"><img src="figures/umap_clustered_Neurog2_Csf1r.png?v=2" alt="Csf1r" width="33%">
<img src="figures/umap_clustered_Neurog2_Kcnj8.png?v=2" alt="Kcnj8" width="33%"><img src="figures/umap_clustered_Neurog2_Otx2.png?v=2" alt="Otx2" width="33%"><img src="figures/umap_clustered_Neurog2_Slc17a7.png?v=2" alt="Slc17a7" width="33%">
<img src="figures/umap_clustered_Neurog2_Aqp4.png?v=2" alt="Aqp4" width="33%"><img src="figures/umap_clustered_Neurog2_Elavl3.png?v=2" alt="Elavl3" width="33%"><img src="figures/umap_clustered_Neurog2_Lhx1.png?v=2" alt="Lhx1" width="33%">
<img src="figures/umap_clustered_Neurog2_Pax2.png?v=2" alt="Pax2" width="33%"><img src="figures/umap_clustered_Neurog2_Slc1a3.png?v=2" alt="Slc1a3" width="33%"><img src="figures/umap_clustered_Neurog2_Arr3.png?v=2" alt="Arr3" width="33%">
<img src="figures/umap_clustered_Neurog2_Elavl4.png?v=2" alt="Elavl4" width="33%"><img src="figures/umap_clustered_Neurog2_Lhx2.png?v=2" alt="Lhx2" width="33%"><img src="figures/umap_clustered_Neurog2_Pax6.png?v=2" alt="Pax6" width="33%">
<img src="figures/umap_clustered_Neurog2_Slc6a9.png?v=2" alt="Slc6a9" width="33%"><img src="figures/umap_clustered_Neurog2_Ascl1.png?v=2" alt="Ascl1" width="33%"><img src="figures/umap_clustered_Neurog2_Emx1.png?v=2" alt="Emx1" width="33%">
<img src="figures/umap_clustered_Neurog2_Lhx4.png?v=2" alt="Lhx4" width="33%"><img src="figures/umap_clustered_Neurog2_Pou4f2.png?v=2" alt="Pou4f2" width="33%"><img src="figures/umap_clustered_Neurog2_Sox11.png?v=2" alt="Sox11" width="33%">
<img src="figures/umap_clustered_Neurog2_Atoh7.png?v=2" alt="Atoh7" width="33%"><img src="figures/umap_clustered_Neurog2_Foxn4.png?v=2" alt="Foxn4" width="33%"><img src="figures/umap_clustered_Neurog2_Malat1.png?v=2" alt="Malat1" width="33%">
<img src="figures/umap_clustered_Neurog2_Prdm1.png?v=2" alt="Prdm1" width="33%"><img src="figures/umap_clustered_Neurog2_Prdx6.png?v=2" alt="Prdx6" width="33%"><img src="figures/umap_clustered_Neurog2_Sox9.png?v=2" alt="Sox9" width="33%">
<img src="figures/umap_clustered_Neurog2_Bsn.png?v=2" alt="Bsn" width="33%"><img src="figures/umap_clustered_Neurog2_Gad1.png?v=2" alt="Gad1" width="33%"><img src="figures/umap_clustered_Neurog2_mScarlet3.png?v=2" alt="mScarlet3" width="33%">
<img src="figures/umap_clustered_Neurog2_mt-Atp6.png?v=2" alt="mt-Atp6" width="33%"><img src="figures/umap_clustered_Neurog2_Rbfox3.png?v=2" alt="Rbfox3" width="33%"><img src="figures/umap_clustered_Neurog2_Tfap2a.png?v=2" alt="Tfap2a" width="33%">
<img src="figures/umap_clustered_Neurog2_Cabp5.png?v=2" alt="Cabp5" width="33%"><img src="figures/umap_clustered_Neurog2_Gfap.png?v=2" alt="Gfap" width="33%"><img src="figures/umap_clustered_Neurog2_Tie1.png?v=2" alt="Tie1" width="33%">
<img src="figures/umap_clustered_Neurog2_Calb1.png?v=2" alt="Calb1" width="33%"><img src="figures/umap_clustered_Neurog2_Glul.png?v=2" alt="Glul" width="33%"><img src="figures/umap_clustered_Neurog2_Neurog2_9SA.png?v=2" alt="Neurog2_9SA" width="33%">
<img src="figures/umap_clustered_Neurog2_Rho.png?v=2" alt="Rho" width="33%"><img src="figures/umap_clustered_Neurog2_Rlbp1.png?v=2" alt="Rlbp1" width="33%"><img src="figures/umap_clustered_Neurog2_Vim.png?v=2" alt="Vim" width="33%">
<img src="figures/umap_clustered_Neurog2_Calb2.png?v=2" alt="Calb2" width="33%"><img src="figures/umap_clustered_Neurog2_Hes1.png?v=2" alt="Hes1" width="33%"><img src="figures/umap_clustered_Neurog2_Neurog2.png?v=2" alt="Neurog2" width="33%">
<img src="figures/umap_clustered_Neurog2_Ccr2.png?v=2" alt="Ccr2" width="33%"><img src="figures/umap_clustered_Neurog2_Hes5.png?v=2" alt="Hes5" width="33%">




![Dotplot](figures/clustered_Neurog2_Dotplot.png?v=3)
![Clusters UMAP](figures/umap_clustered_Neurog2_Clusters.png?v=3)

#### scRNA annotations 
![annotations](figures/annotated_clustered_Neurog2_annotationsON.png?v=3)



# ATAC preprocessing 


![](scenicOuts/QC/TH1_barcode_qc.png?v=2)  
![](scenicOuts/QC/TH1_qc.png?v=2)  
![](scenicOuts/QC/TH2_barcode_qc.png?v=2)
![](scenicOuts/QC/TH2_qc.png?v=2) 





# Co-accessible regions using Mallet 

![](scenicOuts/umap_clusters/celltype_umap.png?v=3) 

![](scenicOuts/umap_clusters/annotated_clusters_res_0.5.png?v=3)

![](scenicOuts/umap_clusters/annotated_clusters_res_0.8.png?v=3)

![](scenicOuts/umap_clusters/annotated_clusters_res_1.0.png?v=3)

![](scenicOuts/umap_clusters/annotated_clusters_res_1.5.png?v=3)

![](scenicOuts/umap_clusters/annotated_clusters_res_2.0.png?v=3)

![](scenicOuts/umap_clusters/annotated_clusters_res_3.0.png?v=3)



# Differential Accessible Regions (DARs)

Differential Accessible Regions (DARs) are chromatin regions whose accessibility **differs between groups of cells**. These regions can include promoters, enhancers, or other regulatory elements.  

DARs are identified by comparing accessibility values across groups, such as cell types or conditions. Regions that are significantly more open in one group compared to another are called **DARs**.  

**Intuition:**  
- DARs highlight regulatory regions that change accessibility across conditions.  
- They are the **region-level differences** that can contribute to differential gene activity if they are linked to a gene.


### Summary: Number of DARs Found per Group

| Group        | Number of DARs |
|-------------|----------------|
| AC          | 40             |
| BC          | 40             |
| Cones       | 40             |
| MG          | 40             |
| MGPC        | 40             |
| Microglia   | 20             |
| Rod         | 20             |



![](scenicOuts/DARs/UMAP_AC_chr10_128821288-128822049.png?v=3)
![](scenicOuts/DARs/UMAP_BC_chr11_54122132-54123756.png?v=3)
![](scenicOuts/DARs/UMAP_Cones_chr10_128821288-128822049.png?v=3)
![](scenicOuts/DARs/UMAP_MG_chr10_21082176-21084166.png?v=3)
![](scenicOuts/DARs/UMAP_MGPC_chr11_54122132-54123756.png?v=3)
![](scenicOuts/DARs/UMAP_Microglia_chr10_44380210-44381883.png?v=3)
![](scenicOuts/DARs/UMAP_Rod_chr10_128821288-128822049.png?v=3)


# Gene activity 

Gene activity is computed by first identifying **co-accessible regions** — regions of chromatin that tend to be open together across cells. For a specific gene, SCENIC only considers the regions that are **physically close to or linked to that gene**, such as its promoter, nearby enhancers, or distal elements connected via co-accessibility.  

The gene activity score depends on how open **these linked regions** are in a given cell. Open regions that are linked to gene X increase the activity score of gene X, but do **not** affect gene Y, because each gene’s score is computed from its own set of regulatory regions.  

**Intuition:** a gene is more likely to be expressed if the chromatin around it and its linked regulatory regions is accessible, which is why gene activity scores are gene-specific.



# Differential gene activity 


Differential gene activity identifies genes whose activity differs **between groups of cells**. After computing gene activity scores for all cells, SCENIC compares these scores across a **metadata category**, such as cell type, condition, or treatment.  

For each gene, it tests whether the activity scores are significantly different between groups. This is still **based on chromatin accessibility**: genes are differentially active because the accessibility of their linked regulatory regions varies between groups. A gene whose regulatory regions are more open in one group than another will have a higher activity score in that group, making it a **differentially active gene**.  

**Intuition:** differential gene activity highlights genes whose regulatory potential changes across conditions, revealing which genes are likely “turned on” in one group of cells compared to another.

