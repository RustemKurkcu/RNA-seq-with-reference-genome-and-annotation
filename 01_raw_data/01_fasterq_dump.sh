#!/bin/bash
#SBATCH --job-name=fastqer_dump_xanadu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=15G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Download fastq files from SRA 
#################################################################

# load parallel module
module load parallel/20180122
module load sratoolkit/3.0.1

# The data are 
    # PBMC_RNA_DATA
    


ACCLIST=accessionlist.txt

cat $ACCLIST | parallel -j 2 fasterq-dump

# compress the files 
ls *fastq | parallel -j 12 gzip

