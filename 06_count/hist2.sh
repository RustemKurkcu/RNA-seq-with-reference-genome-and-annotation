#!/bin/bash
#SBATCH --job-name=hisat2-index
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com
#SBATCH -o hisat2-index_%A.out
#SBATCH -e hisat2-index_%A.err
mkdir index
# download genome and annotation files
wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz

# unzip the files
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.104.gtf.gz

# load HISAT2 module
module load hisat2/2.2.1

# build HISAT2 index
hisat2-build -p 8 Homo_sapiens.GRCh38.dna.primary_assembly.fa index/Homo_sapiens_GRCh38

