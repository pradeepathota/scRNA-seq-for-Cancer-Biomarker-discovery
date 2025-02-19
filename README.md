# ScRNA-Seq Analysis Pipeline

This repository documents an **ScRNA-Seq analysis pipeline** using **FASTQ preprocessing, quality control, genome indexing, and alignment** with **STAR and BWA**.

## Table of Contents
- [Overview](#overview)
- [Requirements](#requirements)
- [Data Preparation](#data-preparation)
- [Quality Control](#quality-control)
- [Trimming Reads](#trimming-reads)
- [Genome Indexing](#genome-indexing)
- [Read Alignment](#read-alignment)
- [Job Submission on SLURM](#job-submission-on-slurm)
- [Results and Outputs](#results-and-outputs)
- [Gene Quantification with featureCounts](#gene-quantification-with-featurecounts)
- [References](#references)

---

## Overview
This pipeline processes **ScRNA-Seq** data by performing:
1. **Data acquisition** (Downloading FASTQ files)
2. **Quality control** (FastQC)
3. **Read trimming** (fastp)
4. **Genome indexing** (STAR/BWA)
5. **Read alignment** (STAR/BWA)

---

## Requirements
Ensure the following tools are installed:

```bash
conda install -c bioconda fastqc fastp star bwa samtools
```

Or load them on **SLURM-based clusters**:

```bash
module load fastqc fastp star bwa samtools
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

### **BWA Alignment**

```bash
bwa mem -t 16 reference/Homo_sapiens.GRCh38.fa \
    trimmed_data/SRRxxxxxxx_1_trimmed.fastq trimmed_data/SRRxxxxxxx_2_trimmed.fastq > alignments/SRRxxxxxxx_aligned.sam
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
#SBATCH -A r00750

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

## Results and Outputs
The pipeline generates:
- **Quality reports** (`fastqc_reports/`)
- **Trimmed reads** (`trimmed_data/`)
- **Genome index** (`star_genome_index/`)
- **Aligned BAM files** (`alignments/`)
- **Log files** (`alignments/*.log`)

Check alignment quality:

```bash
samtools flagstat alignments/SRRxxxxxxx_sorted.bam
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

## References
- [STAR Manual](https://github.com/alexdobin/STAR)
- [BWA Manual](http://bio-bwa.sourceforge.net/)
- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [HISAT2 Documentation](https://daehwankimlab.github.io/hisat2/)

---

