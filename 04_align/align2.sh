#!/bin/bash
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=[0-18]%5

echo `hostname`

#################################################################
# Align reads to genome
#################################################################
module load hisat2/2.2.1
module load samtools/1.12

INDIR=/home/FCAM/skurkcu/CGiardina_mRNASeq_Aug_2023/AlanKuo_RNASeq_Aug2023/02_quality_control/trimmed_sequences
OUTDIR=alignments
mkdir -p $OUTDIR

INDEX=/home/FCAM/skurkcu/CGiardina_mRNASeq_Aug_2023/AlanKuo_RNASeq_Aug2023/genome/hisat2_index/Fhet
ACCLIST=../01_raw_data/accessionlist.txt

NUM=$(expr ${SLURM_ARRAY_TASK_ID} + 1)
SAMPLE=$(sed -n ${NUM}p $ACCLIST)

# run hisat2
hisat2 \
	-p 2 \
	-x $INDEX \
	-1 $INDIR/${SAMPLE}R1_trim_1.fastq.gz \
	-2 $INDIR/${SAMPLE}R2_trim_1.fastq.gz | \
samtools view -@ 1 -S -h -u - | \
samtools sort -@ 1 -T $SAMPLE - >$OUTDIR/$SAMPLE.bam

# index bam files
samtools index $OUTDIR/$SAMPLE.bam