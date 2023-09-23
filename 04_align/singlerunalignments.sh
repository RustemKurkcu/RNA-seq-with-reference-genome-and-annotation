#!/bin/bash
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=13G
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

INDIR=../02_quality_control/trimmed_sequences
OUTDIR=alignments2
mkdir -p $OUTDIR

# this is an array job. 
	# one task will be spawned for each sample
	# for each task, we specify the sample as below
	# use the task ID to pull a single line, containing a single accession number from the accession list
	# then construct the file names in the call to hisat2 as below
	#SED is a text stream editor used on Unix systems to edit files quickly and efficiently. The tool searches through, 
	#replaces, adds, and deletes lines in a text file without opening the file in a text editor. Learn how to use the sed command and its options through easy-to-follow examples.

INDEX=../genome/hisat2_index/Fhet*

NUM=$(expr ${SLURM_ARRAY_TASK_ID} + 1)


# run hisat2
hisat2 \
	-p 12 \
	-x $INDEX \
	-1 $INDIR/1_R1_001.fastq.gz_trim.fastq.gz \
	-2 $INDIR/1_R2_001.fastq.gz_trim.fastq.gz| \
samtools view -@ 8 -S -h -u - | \
samtools sort -@ 8 -T 1_001.fastq.gz_trim.fastq.gz - >$OUTDIR/1_001.fastq.bam

# index bam files
samtools index $OUTDIR/$SAMPLE.bam
