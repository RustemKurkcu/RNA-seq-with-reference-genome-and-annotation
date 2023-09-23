#!/bin/bash
#SBATCH --job-name=align_single
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH -o align_single.out
#SBATCH -e align_single.err

echo `hostname`

#################################################################
# Align reads to genome
#################################################################
module load hisat2/2.2.1
module load samtools/1.12

INDIR=/home/FCAM/skurkcu/CGiardina_mRNASeq_Aug_2023/AlanKuo_RNASeq_Aug2023/02_quality_control/trimmed_sequences
OUTDIR=alignments_single
mkdir -p $OUTDIR

INDEX=/home/FCAM/skurkcu/CGiardina_mRNASeq_Aug_2023/AlanKuo_RNASeq_Aug2023/genome/hisat2_index/Fhet

# Specify one sample for testing
SAMPLE=Kuo1_S131_L004

# run hisat2
hisat2 \
    -p 2 \
    -x $INDEX \
    -1 $INDIR/${SAMPLE}_R1_trim_1.fastq.gz \
    -2 $INDIR/${SAMPLE}_R2_trim_1.fastq.gz | \
samtools view -@ 1 -S -h -u - | \
samtools sort -@ 1 -T $SAMPLE - >$OUTDIR/$SAMPLE.bam

# index bam files
samtools index $OUTDIR/$SAMPLE.bam
