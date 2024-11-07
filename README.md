# project-lab

terminal open command 
:- ssh username@quartz.uits.iu.edu

to create new environment in conda
:- "conda create -n new"

to install library
:- "conda install library name"

to dwonlaod the file
:- "prefetch"

'''
The command "zcat" is used to view the contents of a compressed file without uncompressing it permanently.


converting the files to fastq format
:- find . -name "*.sra" -exec fastq-dump {} \; this command to find all the sra files in your directory and convert them to fastQ format.

quality control
:- fastqc ".fastq".

trimming the data to remove the low quality reads
:- i did trimming with the "fastp" tool

:- here is the command "fastp -i *.fastq -o trimmed_*.fastq -h fastp_report.html -j fastp_report.json"

alignment
:- for the alignment w ehave to download the reference genome that is hg38 fro humann and to download that "wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz" and have to unzip that "gunzip Homo_sapiens.GRCh38.109.gtf.gz"
then the alignment using the bwa.

:- first we have to do indexing for the reference genome for that below is the command.

job submission in quartz
:-#!/bin/bash

#SBATCH -J bwa_index

#SBATCH -p general

#SBATCH -o bwa_index_output.log_%j.txt

#SBATCH -ebwa_index_error.log_%j.err

#SBATCH --mail-type=ALL

#SBATCH --mail-user=vathota@iu.edu

#SBATCH --nodes=1

#SBATCH --ntasks-per-node=1

#SBATCH --time=04:00:00

#SBATCH --mem=100G

#SBATCH -A r00750

#Load BWA module if available, or ensure BWA is accessible

module load bwa                       # Use this if BWA is available as a module

#or activate the environment where BWA is installed

#source activate bioenv               # Uncomment if using a Conda environment with BWA

# Run BWA indexing

bwa index GRCh38.primary_assembly.genome.fa

# alignment
the command for the alignment

:- bwa mem -t 4 GRCh38.primary_assembly.genome.fa trimmed_*.fastq > *_aligned_bwa.sam




