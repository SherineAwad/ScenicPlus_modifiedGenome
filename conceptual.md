# SCENIC+ Workflow (Conceptual Overview)

## 1. Start with barcodes and fragments

We have multiple samples, each with many cells. At this stage, we only know the **barcodes and fragment positions** — no cell type annotation yet.

**Example fragment input file (barcode + chrom/start/end + count):**

```
# primary_contig=JH584295.1
chr1    3000129 3000366 CATTGTAAGCAGGTGG-1      2
chr1    3000137 3000553 TCCACAACAAACATAG-1      1
chr1    3000139 3000372 GCAACAGCATAGTCAT-1      1
chr1    3000225 3000395 GCCAATTAGCTGGACC-1      3
chr1    3000277 3000495 AGTAATGCAAGTTATC-1      1
chr1    3000277 3000495 ATGAAGCCAACTAGAA-1      2
chr1    3000294 3000523 AGCAGGTAGGTGCGGA-1      3
chr1    3000319 3000580 AACATTGTCCGCCAAA-1      4
chr1    3000522 3000671 GAAGCTAAGCAGGTGG-1      7
chr1    3000583 3000626 ATTAGCGGTCATGCAA-1      4
chr1    3000647 3001179 TAAGCCTAGTCACCAG-1      2
chr1    3000827 3000998 GAGTTTGGTTACTTGC-1      2
chr1    3000839 3001075 ATATGCTCAACGTGCT-1      1
chr1    3000941 3001132 TGTCAATCACACAATT-1      2
chr1    3000970 3001091 AGTTACTCAGTTGCGT-1      5
chr1    3000991 3001143 AACATTGTCCGCCAAA-1      5
chr1    3000991 3001194 AACATTGTCCGCCAAA-1      3
```

* Each row = a fragment mapping from a barcode to a genomic location, with **fragment count**.
* This is the **starting input for building pseudobulk BED files**.

---

## 2. Pseudobulk peak calling

* **Input:** Raw fragment files per cell (from barcodes).
* **Process:** Combine reads from all cells in a given sample to generate **sample-level BED files**.
* **Output:** BED files per sample (after pseudobulk), containing **all observed peaks** for that sample.

Example pseudobulk BED file structure (per sample):

| Chrom | Start   | End     | Sample |
| ----- | ------- | ------- | ------ |
| chr1  | 3000129 | 3000366 | Ctrl   |
| chr1  | 3000137 | 3000553 | Ctrl   |
| chr1  | 3000225 | 3000395 | Ctrl   |
| chr1  | 3000277 | 3000495 | Ctrl   |
| chr1  | 3000294 | 3000523 | Ctrl   |
| chr1  | 3000319 | 3000580 | Ctrl   |
| chr1  | 3000522 | 3000671 | Ctrl   |
| chr1  | 3000583 | 3000626 | Ctrl   |
| chr1  | 3000647 | 3001179 | Ctrl   |
| chr1  | 3000827 | 3000998 | OE     |
| chr1  | 3000839 | 3001075 | OE     |
| chr1  | 3000941 | 3001132 | OE     |
| chr1  | 3000970 | 3001091 | OE     |
| chr1  | 3000991 | 3001143 | OE     |
| chr1  | 3000991 | 3001194 | OE     |

* Notice that the **BED file aggregates fragments across all cells**, pseudobulk removes individual barcodes.

---

## 3. MACS2 peak calling

* **Input:** Pseudobulk BED files per sample.
* **Process:** Identify **statistically significant peaks** per sample.
* **Output:** NarrowPeak files per sample.

Example MACS2 output (simplified):

| Chrom | Start   | End     | Sample | Score |
| ----- | ------- | ------- | ------ | ----- |
| chr1  | 3000129 | 3000366 | Ctrl   | 15    |
| chr1  | 3000137 | 3000553 | Ctrl   | 23    |
| chr1  | 3000225 | 3000395 | Ctrl   | 8     |
| chr1  | 3000277 | 3000495 | Ctrl   | 20    |
| chr1  | 3000294 | 3000523 | Ctrl   | 17    |
| chr1  | 3000319 | 3000580 | Ctrl   | 12    |
| chr1  | 3000522 | 3000671 | Ctrl   | 25    |
| chr1  | 3000583 | 3000626 | Ctrl   | 14    |

---

## 4. Consensus peak calling

* **Input:** MACS2 narrowPeak files from all samples.
* **Process:** Merge overlapping peaks → `consensus_peaks.bed`.

Example consensus peaks:

| PeakID | Chrom | Start   | End     |
| ------ | ----- | ------- | ------- |
| Peak1  | chr1  | 3000129 | 3000366 |
| Peak2  | chr1  | 3000137 | 3000553 |
| Peak3  | chr1  | 3000225 | 3000395 |
| Peak4  | chr1  | 3000277 | 3000495 |
| Peak5  | chr1  | 3000294 | 3000523 |
| Peak6  | chr1  | 3000319 | 3000580 |

* Provides a uniform set of peaks across all cells.

---

## 5. QC and TSS processing

* Convert GTF → TSS BED
* Calculate QC metrics (FRIP, TSS enrichment) per barcode
* Filter barcodes based on thresholds

Example QC output:

| Barcode  | FRIP | TSS_enrichment | PassQC |
| -------- | ---- | -------------- | ------ |
| cell_001 | 0.25 | 9.5            | Yes    |
| cell_002 | 0.10 | 4.0            | No     |
| cell_101 | 0.30 | 10.2           | Yes    |
| cell_102 | 0.28 | 11.0           | Yes    |

---

## 6. Create cisTopic object

* **Input:**

  * Fragments dictionary (`fragments_dict.pkl`)
  * QC-passed barcodes (`qc_barcodes_thresholds.pkl`)
  * Consensus peaks (`consensus_peaks.bed`)
* **Output:** `cistopic_objects_mm10.pkl`

---

## 7. Merge cisTopic objects

* Combine samples into one `merged_cistopic.pkl`.

---

## 8. Process scRNA data

* Preprocess, cluster, annotate, merge metadata → `merged_with_meta.pkl`

Example merged metadata:

| Barcode  | Sample | Celltype |
| -------- | ------ | -------- |
| cell_001 | Ctrl   | Rod      |
| cell_002 | Ctrl   | Cone     |
| cell_101 | OE     | Rod      |
| cell_102 | OE     | BC       |

---

## 9. Input to MALLET

* **Format:** Each cell = document, each consensus peak = word, values = fragment counts

Example input matrix:

| Cell     | Peak1 | Peak2 | Peak3 | Peak4 | Peak5 | Peak6 |
| -------- | ----- | ----- | ----- | ----- | ----- | ----- |
| cell_001 | 2     | 1     | 3     | 1     | 3     | 4     |
| cell_002 | 7     | 4     | 2     | 5     | 5     | 3     |
| cell_101 | 0     | 0     | 0     | 0     | 0     | 0     |
| cell_102 | 0     | 0     | 0     | 0     | 0     | 0     |

* Shows how MALLET sees **fragments per peak per cell**, based on consensus peaks.

---

## 10. Run MALLET (LDA)

* **Cell-topic matrix**: fraction of each cell's accessibility explained by topics

| Cell     | Topic1 | Topic2 | Topic3 |
| -------- | ------ | ------ | ------ |
| cell_001 | 0.08   | 0.06   | 0.86   |
| cell_002 | 0.05   | 0.90   | 0.05   |
| cell_101 | 0.12   | 0.04   | 0.84   |
| cell_102 | 0.70   | 0.10   | 0.20   |

* **Topic-region matrix**: which peaks contribute to each topic

| Topic  | Peak1 | Peak2 | Peak3 | Peak4 | Peak5 | Peak6 |
| ------ | ----- | ----- | ----- | ----- | ----- | ----- |
| Topic1 | 0.10  | 0.40  | 0.30  | 0.20  | 0.0   | 0.0   |
| Topic2 | 0.05  | 0.80  | 0.10  | 0.05  | 0.0   | 0.0   |
| Topic3 | 0.25  | 0.05  | 0.35  | 0.35  | 0.3   | 0.2   |

---

## 11. Cluster Cistopic 

### 1️⃣ Input

* **Cell-topic matrix:** fraction of each cell's accessibility explained by topics.

| Cell     | Topic1 | Topic2 | Topic3 |
| -------- | ------ | ------ | ------ |
| cell_001 | 0.08   | 0.06   | 0.86   |
| cell_002 | 0.05   | 0.90   | 0.05   |
| cell_101 | 0.12   | 0.04   | 0.84   |
| cell_102 | 0.70   | 0.10   | 0.20   |

* **Topic-region matrix:** which peaks contribute to each topic.

| Topic  | Peak1 | Peak2 | Peak3 | Peak4 | Peak5 | Peak6 |
| ------ | ----- | ----- | ----- | ----- | ----- | ----- |
| Topic1 | 0.10  | 0.40  | 0.30  | 0.20  | 0.0   | 0.0   |
| Topic2 | 0.05  | 0.80  | 0.10  | 0.05  | 0.0   | 0.0   |
| Topic3 | 0.25  | 0.05  | 0.35  | 0.35  | 0.3   | 0.2   |

---

### 2️⃣ What is being clustered

* The script clusters **cells**, not individual peaks.

* Each cell is represented by a **vector of topic fractions**, e.g.:

  * `cell_001` → `[0.08, 0.06, 0.86]`
  * `cell_002` → `[0.05, 0.90, 0.05]`

* This vector summarizes the cell's **regulatory landscape** in terms of co-accessible topics.

---

### 3️⃣ Calculating similarity

* A **graph-based neighborhood structure** is computed between cells based on their topic vectors.

* Cells with similar topic fractions are connected in this graph.

  * Cells with similar regulatory programs (topic usage) are close in the graph.
  * Cells with different topic usage are far apart.

---

### 4️⃣ Clustering method

* Cells are grouped using **Leiden clustering**, a graph-based community detection method.
* Multiple **resolutions** are used (e.g., 0.6, 1.2, 3) to generate coarser or finer clusters.
* Result: cells with similar regulatory landscapes are grouped together, reflecting shared co-accessible programs.

---

### 5️⃣ UMAP embedding

* Cells are projected into **2D UMAP space** for visualization.
* Clusters are plotted in this space, allowing interpretation of relationships between regulatory programs.

---

### 6️⃣ Role of peaks/topics

* Peaks are **not directly clustered**, but they influence clusters indirectly:

  * Topics are weighted sets of peaks.
  * Cells that share topic usage are grouped together.
  * Therefore, clusters reflect **patterns of co-accessibility in peaks**, summarized by topics.

---

### 7️⃣ Cell type annotation

* If the `cell_data` contains a `celltype` column, clusters are annotated with **dominant cell types**.
* This helps interpret which regulatory programs correspond to which biological populations.

---



## 12. Binarisation  


### 1️⃣ Purpose

This script takes a **Cistopic object** where cells have already been assigned **topic fractions** (from MALLET/LDA) and converts these continuous topic signals into **binary patterns**.

The idea is to make it clear which **regions (peaks) are strongly associated with each topic** and which **cells strongly express each topic**, rather than working with continuous fractional values.

By binarizing, you can:

* Identify the **most characteristic peaks** for each topic.
* Identify which cells **truly belong** to each topic.
* Facilitate downstream analyses such as motif enrichment, regulatory program interpretation, and visualization.

---

### 2️⃣ What is binarized

1. **Regions (peaks) per topic**

   * Using the `ntop` method: selects the top N peaks that contribute most to each topic.
   * Using the `Otsu` method: determines a threshold per topic to classify peaks as either “on” (associated) or “off” (not associated).

2. **Cells per topic**

   * Using the `Li` method: identifies which cells significantly use a topic.
   * Converts each cell’s topic fraction into a binary 1 (topic expressed) or 0 (topic not expressed).

---

### 3️⃣ Input and output

* **Input:** A Cistopic object with a `selected_model` containing `cell_topic` (cells × topics) and `topic_region` (topics × peaks).

* **Output:**

  * Binarized **region-topic matrices** (top N and Otsu).
  * Binarized **cell-topic matrix** (Li).
  * Plots showing thresholds and top peaks per topic.
  * Updated Cistopic object with binarized data attached.

* These outputs are easier to interpret because each 1/0 clearly indicates **presence or absence of signal**, rather than dealing with fractional contributions.

---

### 4️⃣ Conceptual workflow

1. Load the Cistopic object containing the **LDA/MALLET topics**.

2. For each topic, determine which **peaks are most strongly associated**.

3. For each topic, determine which **cells strongly express that topic**.

4. Combine these into **binary matrices**:

   * **Peaks × Topics:** 1 if peak is strongly associated, 0 otherwise.
   * **Cells × Topics:** 1 if cell strongly expresses topic, 0 otherwise.

5. Save the binarized matrices and updated object for downstream analysis.


## How Binarization Connects to DAR Analysis

1. **What DAR analysis does**

DAR analysis identifies **peaks (regions) that are differentially accessible between groups of cells**, e.g., different cell types or clusters. This tells us which regulatory regions are specifically “active” in one cell population versus another.

---

2. **Why binarization is important for DAR**

* After MALLET/LDA, each peak has a **topic weight**, and each cell has a **topic fraction**. These are **continuous values**, which can be noisy.

* Binarization converts these continuous signals into **0/1** matrices:

  * **Region-topic binary matrix:** which peaks are strongly associated with which topics.
  * **Cell-topic binary matrix:** which cells strongly express which topics.

* By doing this, we can **filter out weak signals** and focus only on the **core regions or cells** that truly represent a topic.

---

3. **Practical intuition**

* Imagine a topic represents a regulatory program.
* Binarization picks the **defining peaks** of that program.
* DAR analysis then asks: *“Which of these defining peaks are more accessible in cluster A vs cluster B?”*
* This allows you to **link topics to cell type–specific regulatory regions**.

---

## 13. DAR analysis (Differential Accessible Regions)

1. **Input**

The DAR analysis takes a **binarized or topic-processed Cistopic object**, which contains:

* **Cells × topics:** each cell’s representation in topic space (fractional usage of regulatory programs).
* **Regions × topics:** how peaks contribute to each topic.
* **Metadata:** cell annotations like `celltype` or other groups to compare.
* **Fragment/accessibility matrix:** single-cell peak accessibility information, often **imputed and normalized** from previous steps.

> **Imputation explanation:** Single-cell ATAC data is sparse, meaning many peaks are not detected in every cell. Imputation estimates these missing accessibility values using information from similar cells or topics, producing a more complete and robust accessibility profile for downstream comparison.

This input is usually the output of the **binarization step**, so the signal has been simplified to highlight the most relevant regions per topic and per cell.

---

2. **What the DAR analysis does conceptually**

* The main goal is to find **Differentially Accessible Regions (DARs)** — peaks that are significantly more or less accessible between groups of cells (e.g., cell types, clusters, or experimental conditions).
* To do this, the script first **performs imputation**, which fills in the missing accessibility measurements for peaks not detected in some cells. This ensures that comparisons between cells or groups are **more accurate and less affected by dropouts** in the single-cell data. After imputation, the accessibility values are **normalized** so that differences reflect true biological variation rather than technical bias.
* Then, it identifies **highly variable regions**, which are the peaks that vary the most across all cells, because these are most informative for distinguishing groups.
* Using these variable regions, the script compares groups of cells according to the chosen metadata column (e.g., `celltype`) and calculates which peaks are **significantly different in accessibility**.

---

4. **Output and downstream use**

* The result is a set of **DAR tables** for each group, listing peaks with adjusted p-values and fold changes. These tables represent **cell type–specific or cluster-specific regulatory regions**.
* The script also optionally produces **UMAP plots** showing how the accessibility of top DARs varies across cells in the embedding space, which helps **visualize regulatory differences across populations**.
* The DARs are **stored in the Cistopic object** so that they can be used for further analyses, like **linking regulatory regions to gene activity**, motif enrichment, or integration with transcriptomics.

> **Imputation relevance:** All these outputs rely on imputed data to ensure that the statistical testing of accessibility differences is **not biased by missing values**. Without imputation, many real differences could be missed because peaks were not observed in enough cells.

---

* **MALLET/LDA → Binarization → DAR Analysis**

1. MALLET/LDA identifies **topics** that summarize co-accessible regions.
2. Binarization converts these topics into **discrete signals**, highlighting defining peaks.
3. DAR analysis uses these **imputed and binarized signals** to **statistically test which peaks differ between groups**, giving biologically meaningful regulatory insights.

> **Imputation note:** This highlights that the DAR step doesn’t just use the raw binarized matrix — it first fills in missing values to provide a **full accessibility landscape** for each cell, improving sensitivity and robustness of DAR detection.

---
## 14. Export Region 


The goal of this step is to **take the topics and differentially accessible regions identified in previous analyses and convert them into standard BED files** with genomic coordinates. This makes the regulatory regions discovered in single-cell chromatin accessibility analyses **portable, interpretable, and ready for downstream use**, such as visualization in genome browsers, motif analysis, or integration with other genomic datasets.



## 15. Gene Activity 


**Goal:** The gene activity step estimates **gene-level regulatory potential** in each cell based on **chromatin accessibility**, not gene expression. It assigns accessible peaks to nearby genes and summarizes the accessibility signal to produce a **gene activity score**, purely derived from the DNA accessibility data.

**Conceptually:**

* Each gene has regulatory regions (promoters, enhancers) located upstream, downstream, or within the gene body. The more accessible these regions are in a cell, the higher the **gene’s activity score by proximity**.
* The script uses the **imputed accessibility matrix** from the DAR step, ensuring that sparsely detected peaks are accounted for and comparisons across cells are meaningful.
* Gene activity is calculated by **aggregating the accessibility of peaks that are near a gene**, optionally weighting by distance or region importance.
* The output is a **cell × gene matrix**, where each value reflects the **chromatin accessibility near a gene**, not RNA abundance. This matrix can be used for:

  * Comparing **regulatory potential across cell types**.
  * **Integration with RNA data** if desired, but the scores themselves are **independent of expression**.
  * Visualizing how accessible regulatory regions relate to genes on embeddings like UMAP.

> In short, this step **translates single-cell accessibility into a proximity-based gene activity score**, providing a measure of potential regulatory influence on each gene per cell, without measuring actual transcription.

---



