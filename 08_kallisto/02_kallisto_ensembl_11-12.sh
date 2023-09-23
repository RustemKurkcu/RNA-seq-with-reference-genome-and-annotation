#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --mail-user=shan.kurkcu@uconn.edu
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
# Quantify transcript abundance with Kallisto
#################################################################

# load necessary software
module load kallisto/0.46.1
module load parallel/20180122
module load gffread/0.12.7
module load samtools/1.12

# input/output directories and other resources

INDIR=../04_align/alignment
ENSEMBL=quant_ensembl
mkdir -p $ENSEMBL

# Use the accession list consistent with the one in the `stringtie` script
ACCLIST=./accessionlist11-12.txt

# ensemble transcripts
CDS=/home/FCAM/skurkcu/Marina/genome/hisat2_index2/Homo_sapiens.GRCh38.cds.all.fa

# index transcript fasta file 
kallisto index -i ensembl_transcripts.idx $CDS

# quantify transcript abundance for each sample against ensembl sequences, up to 5 at a time
cat $ACCLIST | \
parallel -j 5 \
    "kallisto quant \
        -i ensembl_transcripts.idx \
        -o $ENSEMBL/{} \
        -b 100 \
        <(zcat $INDIR/{}_trim_1.fastq.gz) \
        <(zcat $INDIR/{}_trim_2.fastq.gz)"

# extract transcript to gene mapping for downstream analysis

# Use the GTF file path consistent with the one in the `stringtie` script
GTF=/home/FCAM/skurkcu/Marina/genome/hisat2_index2/Fhet.gtf

# use regular expressions to pull out IDs
paste \
    <(awk '$3 ~ /transcript/' $GTF | grep -oP "(?<=gene_id \")[A-Z0-9]+") \
    <(awk '$3 ~ /transcript/' $GTF | grep -oP "(?<=transcript_id \")[A-Z0-9]+") \
>ensembl_gene2tx.txt
