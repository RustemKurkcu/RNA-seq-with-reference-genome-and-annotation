#!/bin/bash
#SBATCH --job-name=htseq_count
#SBATCH -n 2
#SBATCH -N 2
#SBATCH -c 6
#SBATCH --mem=30G
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

INDIR=/home/FCAM/skurkcu/CGiardina_mRNASeq_Feb2023/04.1_align/alignmentsv3
OUTDIR=counts_V3.2
mkdir -p $OUTDIR

# accession list
ACCLIST=/home/FCAM/skurkcu/CGiardina_mRNASeq_Feb2023/04.1_align/accession_listV3.2.txt

# gtf formatted annotation file
GTF=/home/FCAM/skurkcu/Marina/genome/hisat2_index2/Fhet.gtf

# run htseq-count on each sample, up to 5 in parallel
cat $ACCLIST | \
parallel --compress -j 1 --tmpdir /scratch/ \
    "htseq-count \
        -s no \
        -r name \
        -f bam \
        $INDIR/{}.bam \
        $GTF \
        > $OUTDIR/{}.counts"

date
echo "Job finished!"

