#!/bin/bash
#SBATCH --job-name=get_human_genome
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Download genome and annotation from ENSEMBL
#################################################################

# load software
module load samtools/1.12

# output directory
GENOMEDIR=../human_genome
mkdir -p $GENOMEDIR

# we're using Homo sapiens from ensembl v105
    # we'll download the genome, GTF annotation, and transcript fasta
    # https://useast.ensembl.org/Homo_sapiens/Info/Index

# download the genome
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
# decompress it
gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz

# download the GTF annotation
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
# decompress it
gunzip Homo_sapiens.GRCh38.110.gtf.gz

# download the transcript fasta
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
# decompress it
gunzip Homo_sapiens.GRCh38.cds.all.fa.gz

# generate simple samtools fai indexes 
samtools faidx Homo_sapiens.GRCh38.dna.toplevel.fa
samtools faidx Homo_sapiens.GRCh38.cds.all.fa

# move everything to the genome directory
mv Homo_sapiens* $GENOMEDIR
