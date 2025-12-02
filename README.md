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




