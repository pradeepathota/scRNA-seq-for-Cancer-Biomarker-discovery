# project-lab

#!/bin/bash

# Connect to Quartz Cluster
ssh username@quartz.uits.iu.edu

# Set up the environment it is not mandatory to create unless you require a specific environment to work in
conda create -n new_environment_name
conda activate new_environment_name
conda install library_name

# Download SRA files
prefetch <accession_id>

# Convert SRA files to FASTQ format
find . -name "*.sra" -exec fastq-dump {} ;

# View the contents of a compressed file
zcat file_name.gz

# Perform quality control using FastQC
fastqc *.fastq

# Trim low-quality reads using fastp
fastp -i input.fastq -o trimmed_output.fastq -h fastp_report.html -j fastp_report.json

# Download the human reference genome (GRCh38)
wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz

# Create a job script for genome indexing
cat <<EOT > bwa_index_job.sh
#!/bin/bash

#SBATCH -J bwa_index
#SBATCH -p general
#SBATCH -o bwa_index_output.log_%j.txt
#SBATCH -e bwa_index_error.log_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your_email>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:00:00
#SBATCH --mem=100G
#SBATCH -A r00750

module load bwa

# Indexing the reference genome
bwa index GRCh38.primary_assembly.genome.fa
EOT

# Submit the genome indexing job
sbatch bwa_index_job.sh

# Perform alignment using BWA
bwa mem -t 4 GRCh38.primary_assembly.genome.fa trimmed_*.fastq > output_aligned_bwa.sam





