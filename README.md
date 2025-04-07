# ScRNA-Seq Analysis Pipeline

This repository documents an **ScRNA-Seq analysis pipeline** including **FASTQ preprocessing, quality control, genome indexing, alignment**, and **downstream analysis** using **Seurat**, **DESeq2**, and **clusterProfiler**.

---

## Table of Contents
- [Overview](#overview)
- [Requirements](#requirements)
- [Data Preparation](#data-preparation)
- [Quality Control](#quality-control)
- [Trimming Reads](#trimming-reads)
- [Genome Indexing](#genome-indexing)
- [Read Alignment](#read-alignment)
- [Job Submission on SLURM](#job-submission-on-slurm)
- [Gene Quantification with featureCounts](#gene-quantification-with-featurecounts)
- [Downstream Analysis](#downstream-analysis)
  - [Creating Seurat Object](#creating-seurat-object)
  - [Quality Control & Filtering](#quality-control--filtering)
  - [Normalization & Feature Selection](#normalization--feature-selection)
  - [Dimensionality Reduction & Clustering](#dimensionality-reduction--clustering)
  - [Differential Expression Analysis](#differential-expression-analysis)
  - [Visualization & Biomarker Discovery](#visualization--biomarker-discovery)
  - [Enrichment Analysis](#enrichment-analysis)
  - [ROC Analysis](#roc-analysis)
- [References](#references)

---

## Overview

This pipeline processes **ScRNA-Seq** data through:
1. **Data acquisition**
2. **Read preprocessing and QC**
3. **Genome alignment**
4. **Feature quantification**
5. **Downstream analysis** using Seurat and other tools.

---

## Requirements

Install necessary tools via conda or load via module:

```bash
conda install -c bioconda fastqc fastp star bwa samtools subread
```
---

## Data Preparation
Download **ScRNA-Seq reads** from **SRA** and convert them to fastq files:

```bash
prefetch SRRxxxxxxx
fasterq-dump --split-files SRRxxxxxxx
```

Move files to the working directory:

```bash
mv SRRxxxxxxx_1.fastq SRRxxxxxxx_2.fastq raw_data/
```

---

## Quality Control
Perform **FastQC** on raw reads:

```bash
fastqc raw_data/*.fastq -o fastqc_reports/
```

---

## Trimming Reads
Trim reads using **fastp**:

```bash
fastp -i raw_data/SRRxxxxxxx_1.fastq -I raw_data/SRRxxxxxxx_2.fastq \
      -o trimmed_data/SRRxxxxxxx_1_trimmed.fastq \
      -O trimmed_data/SRRxxxxxxx_2_trimmed.fastq \
      -h fastp_reports/SRRxxxxxxx_fastp.html -j fastp_reports/SRRxxxxxxx_fastp.json
```

Verify trimmed reads:

```bash
fastqc trimmed_data/*.fastq -o trimmed_fastqc_reports/
```

---

## Genome Indexing
### **STAR Indexing**

```bash
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir star_genome_index/ \
     --genomeFastaFiles reference/Homo_sapiens.GRCh38.fa \
     --sjdbGTFfile reference/Homo_sapiens.GRCh38.109.gtf \
     --sjdbOverhang 100
```

### **BWA Indexing**

```bash
bwa index reference/Homo_sapiens.GRCh38.fa
```

---

## Read Alignment
### **STAR Alignment**

```bash
STAR --runThreadN 16 \
     --genomeDir star_genome_index/ \
     --readFilesIn trimmed_data/SRRxxxxxxx_1_trimmed.fastq trimmed_data/SRRxxxxxxx_2_trimmed.fastq \
     --outFileNamePrefix alignments/SRRxxxxxxx_ \
     --outSAMtype BAM SortedByCoordinate
```

Convert SAM to BAM:

```bash
samtools view -bS alignments/SRRxxxxxxx_aligned.sam | samtools sort -o alignments/SRRxxxxxxx_sorted.bam
```

Index BAM files:

```bash
samtools index alignments/SRRxxxxxxx_sorted.bam
```

---

## Job Submission on SLURM
### **STAR Alignment SLURM Job**

Create a SLURM job script:

```bash
cat <<EOT > star_alignment_job.sh
#!/bin/bash
#SBATCH -J STAR_align
#SBATCH -p general
#SBATCH -o star_align_output.log_%j.txt
#SBATCH -e star_align_error.log_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=08:00:00
#SBATCH --mem=100G
#SBATCH -A r00000 #project account code

module load star

STAR --runThreadN 16 \
     --genomeDir star_genome_index/ \
     --readFilesIn trimmed_data/SRRxxxxxxx_1_trimmed.fastq trimmed_data/SRRxxxxxxx_2_trimmed.fastq \
     --outFileNamePrefix alignments/SRRxxxxxxx_ \
     --outSAMtype BAM SortedByCoordinate
EOT
```

Submit the job:

```bash
sbatch star_alignment_job.sh
```

---

## Gene quantification with featureCounts:

```
featureCounts -T 8 \
  -a Homo_sapiens.GRCh38.109.gtf \
  -o featurecounts_output/gene_counts.txt \
  -g gene_id \
  -t exon \
  -s 2 \
  --extraAttributes gene_name \
  trimmed_aligned/*.bam
```
---
## Downstream Analysis

## Creating Seurat Object

```
counts <- read.table("gene_counts.txt", header = TRUE, row.names = 1, sep = "\t")
metadata <- readxl::read_excel("metadata.xlsx")

seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = "Tumor_vs_Normal", min.cells = 3, min.features = 200)
seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata)
```

## Quality control and filtering

```
seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
Seurat::VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)
```

## Normalization & Feature Selection

```
seurat_obj <- Seurat::SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
```

## Dimensionality Reduction & Clustering

```
seurat_obj <- Seurat::RunPCA(seurat_obj)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 1.2)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:10)
seurat_obj <- Seurat::RunTSNE(seurat_obj, dims = 1:10)
```

## Differential Expression Analysis

```
Seurat::Idents(seurat_obj) <- seurat_obj$Condition
deg_results <- Seurat::FindMarkers(seurat_obj, ident.1 = "Tumor", ident.2 = "Normal", min.pct = 0.1, logfc.threshold = 0.1)
deg_results$p_val_adj <- p.adjust(deg_results$p_val, method = "fdr")
write.csv(deg_results, "DEGs_Tumor_vs_Normal.csv")
```

## Visualization & Biomarker Discovery

```
EnhancedVolcano::EnhancedVolcano(deg_results, lab = rownames(deg_results), x = "avg_log2FC", y = "p_val_adj")

biomarkers <- dplyr::filter(deg_results, p_val_adj < 0.05, abs(avg_log2FC) > 1)
expr_data <- Seurat::GetAssayData(seurat_obj, assay = "SCT", slot = "scale.data")[rownames(biomarkers), ]
pheatmap::pheatmap(expr_data, annotation_col = seurat_obj@meta.data["Condition"])

```

## Enrichment Analysis

```
entrez_ids <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(biomarkers),
                                    column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

ego <- clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
clusterProfiler::dotplot(ego)

ekegg <- clusterProfiler::enrichKEGG(gene = entrez_ids, organism = 'hsa', pvalueCutoff = 0.05)
barplot(ekegg)
```

## ROC Analysis

```
library(pROC)
gene_expr <- Seurat::FetchData(seurat_obj, vars = "FOXA1")
roc_obj <- pROC::roc(seurat_obj$Condition, gene_expr$FOXA1)
plot(roc_obj, main = "ROC Curve for FOXA1")
pROC::auc(roc_obj)

```

## Results and Outputs

The pipeline generates the following outputs at different stages:

- **Quality reports** (`fastqc_reports/`)  
  - HTML reports from FastQC before and after trimming.
  
- **Trimmed reads** (`trimmed_data/`)  
  - FASTQ files after adapter and quality trimming using `fastp`.

- **Genome index** (`star_genome_index/`)  
  - STAR genome indices used for read alignment.

- **Aligned BAM files** (`alignments/`)  
  - Sorted and indexed BAM files generated by STAR .

- **Log files** (`alignments/*.log`)  
  - Log files containing alignment summaries and diagnostics.

- **Count matrix** (`featurecounts_output/`)  
  - Gene expression matrix generated using `featureCounts`.

- **Seurat object & metadata**  
  - Normalized and filtered single-cell data stored in R.

- **Clustering results**  
  - Cell clusters visualized with UMAP and t-SNE plots.

- **Marker gene identification**  
  - Top marker genes per cluster using `FindAllMarkers`.

- **DEG analysis results** (`DEGs_Tumor_vs_Normal.csv`)  
  - Differential expression results between tumor and normal cells.

- **Plots**  
  - Volcano plots, violin plots, and feature plots of DEGs and biomarkers.

- **Heatmaps**  
  - Expression heatmaps of top DEGs and validated biomarkers.

- **Enrichment results**  
  - GO and KEGG enrichment dotplots and barplots.

- **ROC curves & AUC values**  
  - ROC plots and AUC scores assessing biomarker diagnostic power.

Check alignment quality:

```bash
samtools flagstat alignments/SRRxxxxxxx_sorted.bam
```



## References

- [STAR Manual](https://github.com/alexdobin/STAR)  
  - Ultrafast universal RNA-seq aligner for spliced reads.

- [BWA Manual](http://bio-bwa.sourceforge.net/)  
  - Burrows-Wheeler Aligner for aligning short reads to a reference genome.

- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
  - Quality control tool for high-throughput sequence data.

- [fastp GitHub](https://github.com/OpenGene/fastp)  
  - A fast all-in-one preprocessing tool for FASTQ files.

- [Seurat Documentation](https://satijalab.org/seurat/)  
  - R package designed for QC, analysis, and exploration of single-cell RNA-seq data.

- [clusterProfiler Manual](https://yulab-smu.top/biomedical-knowledge-mining-book/)  
  - R package for comparing biological themes among gene clusters.

- [EnhancedVolcano GitHub](https://github.com/kevinblighe/EnhancedVolcano)  
  - Tool for creating publication-ready volcano plots.

- [DESeq2 Bioconductor](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)  
  - Differential gene expression analysis based on the negative binomial distribution.

- [org.Hs.eg.db Annotation Database](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)  
  - Genome-wide annotation for Human using Entrez Gene IDs.

- [featureCounts Manual](http://bioinf.wehi.edu.au/featureCounts/)  
  - Efficient program for assigning sequence reads to genomic features.

- [pROC Package](https://cran.r-project.org/web/packages/pROC/index.html)  
  - Tools for visualizing, smoothing and comparing ROC curves.


---

