#!/bin/bash
#SBATCH --job-name=hisat2_index
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
# Index the Genome
#################################################################

# load software
module load hisat2/2.2.0

# input/output directories
OUTDIR=../genome/hisat2_index
mkdir -p $OUTDIR

GENOME=../genome/Homo_sapiens.GRCh38.dna.toplevel.fa

hisat2-build -p 16 $GENOME $OUTDIR/Fhet
