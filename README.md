# Arabidopsis miRNA-seq Analysis – Multiple Independent Pipelines

This repository includes **three independent workflows** for analyzing small RNA-seq data in *Arabidopsis thaliana*, using:

- `miRDeep2` for known miRNA quantification
- `Bowtie2 + featureCounts` for read alignment and counting
- `DESeq2` and `edgeR` for differential expression
- Optional: `clusterProfiler` for GSEA-based biological interpretation

---

## 🧪 Pipeline A – miRDeep2 for Known miRNA Quantification

**Tools**: `cutadapt`, `mapper.pl`, `quantifier.pl`, `bowtie`, `miRDeep2`

### Steps:

1. **Adapter Trimming**  
   - Adapter: `AGATCGGAAGAGC`  
   - Output: `trimmed/*.fastq`

2. **Read Collapsing & Mapping** using `mapper.pl`  
   - Adapter for collapse: `TGGAATTCTCGGGTGCCAAGG`  
   - Output:
     - `collapsed/*.fa` (collapsed reads)
     - `collapsed/*.arf` (Bowtie mappings)

3. **miRNA Quantification** using `quantifier.pl`  
   - Reference: `mature_ath.fa`, `hairpin_ath.fa`  
   - Output: `quantifier_out/` (read count files per sample)

4. **Merge Counts** with R:
   - All `.csv` outputs from `quantifier.pl` are merged into `merged_counts.csv` using:
     ```r
     tools::file_path_sans_ext() + dplyr::merge()
     ```
---

## 🧮 Pipeline B – Bowtie2 + featureCounts

**Tools**: `fastp`, `FastQC`, `bowtie2`, `samtools`, `featureCounts`

### Steps:

1. **QC + Trimming**  
   - `fastp` generates trimmed FASTQ + quality reports

2. **Mapping**  
   - Aligns reads to reference genome using `bowtie2`  
   - Index: `references/mirna_index`  
   - Output: `alignment/*.sam`

3. **SAM → BAM Processing**
   - Sorting, indexing, and `flagstat` via `samtools`

4. **Feature Quantification**
   - `featureCounts` counts reads using a miRNA-specific `GTF`  
   - Output: `counts/miRNA_counts.txt`

---

## 📊 Pipeline C – Differential Expression & GSEA

Both `DESeq2` and `edgeR` were applied **separately** using `merged_counts.csv`.

### A. DESeq2

- **Script**: `01_deseq2_analysis.R`
- Filters: keep miRNAs with count > 5 in ≥2 samples  
- Output:
  - `deseq2_all_results.csv`
  - `deseq2_significant_results.csv` (padj < 0.05)

### B. edgeR

- **Script**: `03_edger_analysis.R`
- Uses `filterByExpr()` + TMM normalization
- Output:
  - `edgeR_results.csv`

---

## 🔍 Optional: GSEA (GO Enrichment)

- **Script**: `02_gsea_annotation.R`
- Uses ranked log2FC list from DESeq2
- Tool: `clusterProfiler` + `org.At.tair.db`
- Output:
  - `gsea_go_results.csv`
  - `dotplot()` visualization

---
## 🛠 Requirements

### Command-line:
- `cutadapt`, `fastp`, `FastQC`, `multiqc`
- `miRDeep2`, `mapper.pl`, `quantifier.pl`
- `bowtie`, `bowtie2`
- `samtools`, `featureCounts`

### R packages:
```r
BiocManager::install(c("DESeq2", "edgeR", "clusterProfiler", "org.At.tair.db", "AnnotationDbi"))
install.packages("tidyverse")


