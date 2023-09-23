#!/bin/bash
#SBATCH --job-name=gtf_annotate
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#################################################################
# Annotate merged transcripts with ENSEMBL IDs
#################################################################

# load software
module load gffcompare/0.12.6

# input/output variables
INDIR=merged_transcripts
OUTDIR=annotated_transcripts
mkdir -p $OUTDIR

REFGTF=/home/FCAM/skurkcu/Marina/genome/hisat2_index2/Fhet.gtf
MERGEDGTF=$INDIR/merged.gtf

# run gffcompare
gffcompare -o gffcmp -r $REFGTF $MERGEDGTF

mv gffcmp* annotated_transcripts