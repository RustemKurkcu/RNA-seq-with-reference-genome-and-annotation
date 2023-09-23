#!/bin/bash
#SBATCH --job-name=htseq_count_sp
#SBATCH -n 2
#SBATCH -N 2
#SBATCH -c 6
#SBATCH --mem=35G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Generate Counts 
#################################################################
module load htseq/0.13.5
module load parallel/20180122

INDIR=/home/FCAM/skurkcu/CGiardina_mRNASeq_Aug_2023/AlanKuo_RNASeq_Aug2023/04_align/alignment
OUTDIR=counts_scratchspace
mkdir -p $OUTDIR

# accession list
ACCLIST=./accessionlist17-18.txt

# gtf formatted annotation file
GTF=/home/FCAM/skurkcu/Marina/genome/hisat2_index2/Fhet.gtf

# run htseq-count on each sample, up to 2 in parallel
cat $ACCLIST | \
parallel --compress -j 2 --tmpdir /scratch/ \
    "htseq-count \
        -s no \
        -r name \
        -f bam \
        $INDIR/{}.bam \
        $GTF \
        > $OUTDIR/{}.counts"

date
echo "Job finished!"
