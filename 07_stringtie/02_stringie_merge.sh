#!/bin/bash
#SBATCH --job-name=stringtie_merge
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

#################################################################
# Merge stringtie assemblies
#################################################################

# load necessary software
module load stringtie/2.1.5

# input/output directories and other resources

INDIR=transcripts
OUTDIR=merged_transcripts
mkdir -p $OUTDIR

# Generate a list of all GTF files from stringtie assembly
ls $INDIR/*gtf >gtflist.txt

# Use the GTF file path consistent with the one in the previous script
GTF=/home/FCAM/skurkcu/Marina/genome/hisat2_index2/Fhet.gtf

# run stringtie merge
stringtie --merge \
    -p 5 \
    -G $GTF \
    -o $OUTDIR/merged.gtf \
    gtflist.txt