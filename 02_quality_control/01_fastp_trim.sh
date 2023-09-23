#!/bin/bash
#SBATCH --job-name=fastp_trimming
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date
#Script Name: RNASeq Preprocessing using fastp
# 
# Description: This script preprocesses raw RNASeq data by trimming
# low-quality bases and Illumina adaptors using the fastp tool.
# The script processes multiple samples in parallel for efficiency,
# and generates both trimmed sequence files and quality reports.
# 
# Input:
# - Raw FASTQ files located in the 01_raw_data directory
# - An accession list (accessionlist.txt) indicating sample names
# 
# Output:
# - Trimmed sequence files placed in the trimmed_sequences directory
# - fastp quality reports (JSON and HTML) saved in fastp_reports directory
#
# Requirements:
# - fastp version 0.23.2
# - GNU parallel version 20180122
# 
# Usage: 
# Run this script in a directory containing the 01_raw_data directory 
# and the accessionlist.txt file.
# ====================
#################################################################
# Trimming/QC of reads using fastp
#################################################################
module load fastp/0.23.2
module load parallel/20180122

# set input/output directory variables
INDIR=../01_raw_data
REPORTDIR=fastp_reports
mkdir -p $REPORTDIR
TRIMDIR=trimmed_sequences
mkdir -p $TRIMDIR

ACCLIST=$INDIR/accessionlist.txt

# run fastp in parallel, 4 samples at a time
cat $ACCLIST | parallel -j 4 \
fastp \
    --in1 $INDIR/{}R1_1.fastq.gz \
    --in2 $INDIR/{}R2_1.fastq.gz \
    --out1 $TRIMDIR/{}R1_trim_1.fastq.gz \
    --out2 $TRIMDIR/{}R2_trim_1.fastq.gz \
    --json $REPORTDIR/{}R1_fastp.json \
    --html $REPORTDIR/{}R1_fastp.html
