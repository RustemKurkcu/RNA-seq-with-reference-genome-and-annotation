#!/bin/bash
#SBATCH --job-name=align
#SBATCH -n 2
#SBATCH -N 2
#SBATCH -c 8
#SBATCH --mem=25G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=[1-18]%6

echo `hostname`

#################################################################
# Align reads to genome
#################################################################
module load hisat2/2.2.1
module load samtools/1.12

OUTDIR=alignment
mkdir -p $OUTDIR

# This assumes that accession_list.txt is in the same directory as the script
ACCESSION_LIST=$(dirname "${BASH_SOURCE[0]}")/accessionlist.txt

# Use sed to get the accession number for this task ID
ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $ACCESSION_LIST)

# Construct the input file names using the accession number
READ1=/home/FCAM/skurkcu/CGiardina_mRNASeq_Aug_2023/AlanKuo_RNASeq_Aug2023/02_quality_control/trimmed_sequences/${ACCESSION}R1_trim_1.fastq.gz
READ2=/home/FCAM/skurkcu/CGiardina_mRNASeq_Aug_2023/AlanKuo_RNASeq_Aug2023/02_quality_control/trimmed_sequences/${ACCESSION}R2_trim_1.fastq.gz

# Check if the input files exist
if [ ! -f $READ1 ] || [ ! -f $READ2 ]; then
    echo "Input files $READ1 and/or $READ2 do not exist. Exiting now ..."
    exit 1
fi

# run hisat2
hisat2 \
    -p 2 \
    -x /home/FCAM/skurkcu/CGiardina_mRNASeq_Feb2023/genome/hisat2_index/Fhet \
    -1 $READ1 \
    -2 $READ2 | \
samtools view -@ 1 -S -h -u - | \
samtools sort -@ 1 -T $ACCESSION - > $OUTDIR/$ACCESSION.bam

# index bam files
samtools index $OUTDIR/$ACCESSION.bam
