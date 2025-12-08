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


# scRNA preprocessing 


<img src="figures/umap_clustered_Neurog2_Abca8a.png?v=1" alt="Abca8a" width="33%"><img src="figures/umap_clustered_Neurog2_Chat.png?v=1" alt="Chat" width="33%"><img src="figures/umap_clustered_Neurog2_Insm1.png?v=1" alt="Insm1" width="33%">
<img src="figures/umap_clustered_Neurog2_Nrl.png?v=1" alt="Nrl" width="33%"><img src="figures/umap_clustered_Neurog2_Samples.png?v=1" alt="Samples" width="33%"><img src="figures/umap_clustered_Neurog2_Acta2.png?v=1" alt="Acta2" width="33%">
<img src="figures/umap_clustered_Neurog2_Clusters.png?v=1" alt="Clusters" width="33%"><img src="figures/umap_clustered_Neurog2_Isl1.png?v=1" alt="Isl1" width="33%"><img src="figures/umap_clustered_Neurog2_Olig2.png?v=1" alt="Olig2" width="33%">
<img src="figures/umap_clustered_Neurog2_Sebox.png?v=1" alt="Sebox" width="33%"><img src="figures/umap_clustered_Neurog2_Apoe.png?v=1" alt="Apoe" width="33%"><img src="figures/umap_clustered_Neurog2_Csf1r.png?v=1" alt="Csf1r" width="33%">
<img src="figures/umap_clustered_Neurog2_Kcnj8.png?v=1" alt="Kcnj8" width="33%"><img src="figures/umap_clustered_Neurog2_Otx2.png?v=1" alt="Otx2" width="33%"><img src="figures/umap_clustered_Neurog2_Slc17a7.png?v=1" alt="Slc17a7" width="33%">
<img src="figures/umap_clustered_Neurog2_Aqp4.png?v=1" alt="Aqp4" width="33%"><img src="figures/umap_clustered_Neurog2_Elavl3.png?v=1" alt="Elavl3" width="33%"><img src="figures/umap_clustered_Neurog2_Lhx1.png?v=1" alt="Lhx1" width="33%">
<img src="figures/umap_clustered_Neurog2_Pax2.png?v=1" alt="Pax2" width="33%"><img src="figures/umap_clustered_Neurog2_Slc1a3.png?v=1" alt="Slc1a3" width="33%"><img src="figures/umap_clustered_Neurog2_Arr3.png?v=1" alt="Arr3" width="33%">
<img src="figures/umap_clustered_Neurog2_Elavl4.png?v=1" alt="Elavl4" width="33%"><img src="figures/umap_clustered_Neurog2_Lhx2.png?v=1" alt="Lhx2" width="33%"><img src="figures/umap_clustered_Neurog2_Pax6.png?v=1" alt="Pax6" width="33%">
<img src="figures/umap_clustered_Neurog2_Slc6a9.png?v=1" alt="Slc6a9" width="33%"><img src="figures/umap_clustered_Neurog2_Ascl1.png?v=1" alt="Ascl1" width="33%"><img src="figures/umap_clustered_Neurog2_Emx1.png?v=1" alt="Emx1" width="33%">
<img src="figures/umap_clustered_Neurog2_Lhx4.png?v=1" alt="Lhx4" width="33%"><img src="figures/umap_clustered_Neurog2_Pou4f2.png?v=1" alt="Pou4f2" width="33%"><img src="figures/umap_clustered_Neurog2_Sox11.png?v=1" alt="Sox11" width="33%">
<img src="figures/umap_clustered_Neurog2_Atoh7.png?v=1" alt="Atoh7" width="33%"><img src="figures/umap_clustered_Neurog2_Foxn4.png?v=1" alt="Foxn4" width="33%"><img src="figures/umap_clustered_Neurog2_Malat1.png?v=1" alt="Malat1" width="33%">
<img src="figures/umap_clustered_Neurog2_Prdm1.png?v=1" alt="Prdm1" width="33%"><img src="figures/umap_clustered_Neurog2_Prdx6.png?v=1" alt="Prdx6" width="33%"><img src="figures/umap_clustered_Neurog2_Sox9.png?v=1" alt="Sox9" width="33%">
<img src="figures/umap_clustered_Neurog2_Bsn.png?v=1" alt="Bsn" width="33%"><img src="figures/umap_clustered_Neurog2_Gad1.png?v=1" alt="Gad1" width="33%"><img src="figures/umap_clustered_Neurog2_mScarlet3.png?v=1" alt="mScarlet3" width="33%">
<img src="figures/umap_clustered_Neurog2_mt-Atp6.png?v=1" alt="mt-Atp6" width="33%"><img src="figures/umap_clustered_Neurog2_Rbfox3.png?v=1" alt="Rbfox3" width="33%"><img src="figures/umap_clustered_Neurog2_Tfap2a.png?v=1" alt="Tfap2a" width="33%">
<img src="figures/umap_clustered_Neurog2_Cabp5.png?v=1" alt="Cabp5" width="33%"><img src="figures/umap_clustered_Neurog2_Gfap.png?v=1" alt="Gfap" width="33%"><img src="figures/umap_clustered_Neurog2_Tie1.png?v=1" alt="Tie1" width="33%">
<img src="figures/umap_clustered_Neurog2_Calb1.png?v=1" alt="Calb1" width="33%"><img src="figures/umap_clustered_Neurog2_Glul.png?v=1" alt="Glul" width="33%"><img src="figures/umap_clustered_Neurog2_Neurog2_9SA.png?v=1" alt="Neurog2_9SA" width="33%">
<img src="figures/umap_clustered_Neurog2_Rho.png?v=1" alt="Rho" width="33%"><img src="figures/umap_clustered_Neurog2_Rlbp1.png?v=1" alt="Rlbp1" width="33%"><img src="figures/umap_clustered_Neurog2_Vim.png?v=1" alt="Vim" width="33%">
<img src="figures/umap_clustered_Neurog2_Calb2.png?v=1" alt="Calb2" width="33%"><img src="figures/umap_clustered_Neurog2_Hes1.png?v=1" alt="Hes1" width="33%"><img src="figures/umap_clustered_Neurog2_Neurog2.png?v=1" alt="Neurog2" width="33%">
<img src="figures/umap_clustered_Neurog2_Ccr2.png?v=1" alt="Ccr2" width="33%"><img src="figures/umap_clustered_Neurog2_Hes5.png?v=1" alt="Hes5" width="33%"><img src="figures/umap_clustered_Neurog2_Notch1.png?v=1" alt="Notch1" width="33%">
<img src="figures/umap_clustered_Neurog2_Rpe65.png?v=1" alt="Rpe65" width="33%">




![Dotplot](figures/clustered_Neurog2_Dotplot.png?v=1)
![Clusters UMAP](figures/umap_clustered_Neurog2_Clusters.png?v=)




