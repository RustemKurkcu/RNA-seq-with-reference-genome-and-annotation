#!/bin/bash
#SBATCH --job-name=stringtie_assemble
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
# Assemble transcripts with stringtie
#################################################################

# load necessary software
module load stringtie/2.1.5
module load parallel/20180122

# input/output directories and other resources

INDIR=../04_align/alignment
OUTDIR=transcripts
mkdir -p $OUTDIR

# The correct GTF file path
GTF=/home/FCAM/skurkcu/Marina/genome/hisat2_index2/Fhet.gtf

# Use the accession list consistent with the one in the `qualimap` script
ACCLIST=./accessionlist11-12.txt

# run stringtie on all samples, up to 5 in parallel, 
# and capture the standard output and errors for troubleshooting.
cat $ACCLIST | \
parallel -j 5 \
    "stringtie \
        -o $OUTDIR/{}.gtf \
        -G $GTF \
        -p 2 \
        $INDIR/{}.bam"
