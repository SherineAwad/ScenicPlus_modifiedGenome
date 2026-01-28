Perfect — I understand exactly what you want. I will **rewrite the previous Markdown explanation exactly as it was**, keeping **all the explanations and logic from the entire thread**, **no changes anywhere else**, but **replace the very first example with the fragment input you just provided**. I’ll make sure the tables and flow stay the same otherwise.

---

# Single-cell ATAC/SCENIC+ Workflow (Conceptual Overview)

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

Here’s your latest explanation written neatly in **Markdown** format:

---

## How the Cluster Cistopic Script Works

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

### 2️⃣ What is beingg clustered

* The script clusters **cells**, not individual peaks.

* Each cell is represented by a **vector of topic fractions**, e.g.:

  * `cell_001` → `[0.08, 0.06, 0.86]`
  * `cell_002` → `[0.05, 0.90, 0.05]`

* This vector summarizes the cell's **regulatory landscape** in terms of co-accessible topics.

---

### 3️⃣ Calculating similarity

* A **distance or similarity metric** (e.g., Euclidean distance or cosine similarity) is computed between cells based on their topic vectors.

**Example (Euclidean distance):**

* Distance between `cell_001` `[0.08, 0.06, 0.86]` and `cell_101` `[0.12, 0.04, 0.84]` → small → similar cells.
* Distance between `cell_001` `[0.08, 0.06, 0.86]` and `cell_102` `[0.70, 0.10, 0.20]` → large → very different cells.

---

### 4️⃣ Clustering method

* After computing distances, a clustering algorithm groups cells with similar topic distributions:

  * **Hierarchical clustering** → builds a tree and cuts it into clusters.
  * **K-means** → assigns cells to `k` clusters based on vectors.
  * **Graph-based clustering (Louvain/Leiden)** → builds a neighbor graph and detects communities.

* Result: cells with similar regulatory landscapes are grouped together.

---

### 5️⃣ Role of peaks/topics

* Peaks are **not directly clustered**, but they influence clusters indirectly:

  * Topics are weighted sets of peaks.
  * Cells that share topic usage are grouped together.
  * Therefore, clusters reflect **co-accessibility patterns in peaks**, summarized by topics.

---

### 🔑 TL;DR

1. Each cell → vector of topic fractions.
2. Compute similarity/distance between cells.
3. Apply clustering algorithm.
4. Cells with similar topic distributions → same cluster.
5. Peaks influence clusters indirectly via topics.

---





