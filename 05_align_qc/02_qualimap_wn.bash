#!/bin/bash 
#SBATCH --job-name=qualimap
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
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

INDIR=/home/FCAM/skurkcu/CGiardina_mRNASeq_Aug_2023/AlanKuo_RNASeq_Aug2023/04_align/alignment/
OUTDIR=qualimap_reports
mkdir -p $OUTDIR

# gtf annotation is required here
GTF=/home/FCAM/skurkcu/CGiardina_mRNASeq_Feb2023/genome/hisat2_index/Homo_sapiens.GRCh38.106.chr.gtf

# accession list
ACCLIST=./
V3_1011A_S3_R.bam
# run qualimap in parallel
cat $ACCLIST | \
parallel -j 1 \
    qualimap \
        rnaseq \
        -bam $INDIR/{}.bam \
        -gtf $GTF \
        -outdir $OUTDIR/{} \
        --java-mem-size=6G  
