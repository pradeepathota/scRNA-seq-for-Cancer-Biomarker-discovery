# 🧬 scRNA-seq Tumor vs Normal Analysis

This repository presents a complete pipeline for single-cell RNA sequencing (scRNA-seq) data analysis, comparing tumor and normal samples. It demonstrates a reproducible workflow from raw FASTQ data preprocessing to downstream bioinformatics and visualization steps.

---

## 🔍 Overview

- **Objective:** To identify differentially expressed genes and potential biomarkers that distinguish tumor cells from normal cells using single-cell RNA-seq data.
- **Tools Used:** FASTQC, fastp, STAR, featureCounts, Seurat, clusterProfiler, EnhancedVolcano, pROC.
- **Environment:** UNIX shell, SLURM cluster (for alignment), R (for downstream analysis).

---

## 📋 Workflow Summary

1. **Data Acquisition:** Download raw sequencing data from SRA.
2. **Quality Control:** Perform pre- and post-trimming quality checks using FastQC and fastp.
3. **Read Trimming:** Trim low-quality reads and adapters using fastp.
4. **Genome Indexing:** Prepare STAR genome index with the human genome.
5. **Alignment:** Align reads to the reference genome using STAR.
6. **Quantification:** Generate gene expression matrix using featureCounts.
7. **Seurat Analysis:**
   - Create Seurat object
   - Filter low-quality cells
   - Normalize and scale data
   - Perform PCA, clustering, UMAP, and t-SNE
   - Identify DEGs between tumor and normal samples
   - Visualize key marker genes
8. **Functional Analysis:**
   - GO and KEGG enrichment of DEGs
   - ROC analysis for biomarker validation

---

## 📁 Outputs

- **Quality reports**: FastQC/fastp HTML files
- **Trimmed FASTQ files**
- **Aligned BAM files**: Sorted and indexed
- **Gene count matrix**: Output from featureCounts
- **Seurat results**: Clustering, UMAP, DEGs
- **Plots**: Volcano, violin, UMAP, t-SNE, heatmaps, ROC
- **Enrichment results**: GO/KEGG dotplots and barplots
- **CSV tables**: DEGs with adjusted p-values and fold-changes

---

## 👩‍💻 Author

**venkata pradeep kumar athota**  
---

## 📚 References

- Seurat: https://satijalab.org/seurat/
- STAR: https://github.com/alexdobin/STAR
- featureCounts: http://bioinf.wehi.edu.au/featureCounts/
- fastp: https://github.com/OpenGene/fastp
- clusterProfiler: https://yulab-smu.top/biomedical-knowledge-mining-book/
- EnhancedVolcano: https://github.com/kevinblighe/EnhancedVolcano
- pROC: https://cran.r-project.org/web/packages/pROC/

---

## 📌 Note

This project showcases an end-to-end scRNA-seq pipeline built for research and learning purposes. Please adapt file paths and parameters as needed for your environment or dataset.
