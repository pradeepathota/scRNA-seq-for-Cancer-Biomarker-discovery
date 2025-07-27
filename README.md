# ScRNA-Seq Analysis Pipeline

This repository documents an **ScRNA-Seq analysis pipeline** including **FASTQ preprocessing, quality control, genome indexing, alignment**, and **downstream analysis** using **Seurat**, **DESeq2**, and **clusterProfiler**.

---
## 👨‍💻 Author

**Venkata Pradeep Kumar Athota**  
📧 pradeepathota3@gmail.com  
🔗 [LinkedIn](https://linkedin.com/in/pradeepathota)


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
  - [Creating Seurat Object](#creating--seurat-object)
  - [Quality Control & Filtering](#quality-control--filtering)
  - [Normalization & Feature Selection](#normalization--feature-selection)
  - [Dimensionality Reduction & Clustering](#dimensionality-reduction--clustering)
  - [Differential Expression Analysis](#differential-expression--analysis)
  - [Visualization & Biomarker Discovery](#visualization--biomarker-discovery)
  - [Enrichment Analysis](#enrichment-=analysis)
  - [ROC Analysis](#roc--analysis)
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
## Check alignment quality:

```bash
samtools flagstat alignments/SRRxxxxxxx_sorted.bam
```

---

## Job Submission example for alignment on SLURM
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
seurat_obj <- CreateSeuratObject(counts = read.table("counts.txt", header=TRUE, row.names=1))
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt")
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindClusters(seurat_obj)
```
## Differential Expression
```
deg <- FindMarkers(seurat_obj, ident.1 = "Tumor", ident.2 = "Normal")
```
## Enrichment Analysis
```
entrez_ids <- mapIds(org.Hs.eg.db, keys=rownames(deg), column="ENTREZID", keytype="SYMBOL", multiVals="first")
ego <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db)
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

