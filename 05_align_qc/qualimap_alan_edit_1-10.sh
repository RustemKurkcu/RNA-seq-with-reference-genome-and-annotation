#!/bin/bash 
#SBATCH --job-name=qualimap
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
# Considering memory requirements, I recommend using 25G if you run 5 processes with 5G each.
#SBATCH --mem=25G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

##################################
# calculate stats on alignments
##################################
# this time we'll use qualimap

# load software--------------------------------------------------------------------------
module load qualimap/2.2.1
module load parallel/20180122

# input, output directories--------------------------------------------------------------

INDIR=../04_align/alignment
OUTDIR=qualimap_reports
mkdir -p $OUTDIR

# gtf annotation is required here
# Ensure the GTF path is correct

GTF=/home/FCAM/skurkcu/Marina/genome/hisat2_index2/Fhet.gtf
# accession list
ACCLIST=./accessionlist1-10.txt

# run qualimap in parallel, capturing the standard output and errors for troubleshooting.
cat $ACCLIST | \
parallel -j 5 \
    "qualimap rnaseq \
        -bam $INDIR/{}.bam \
        -gtf $GTF \
        -outdir $OUTDIR/{} \
        --java-mem-size=5G > $OUTDIR/{}/qualimap.out 2> $OUTDIR/{}/qualimap.err"
