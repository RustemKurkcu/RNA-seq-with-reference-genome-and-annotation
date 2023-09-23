#!/bin/bash
#SBATCH --job-name=stringtie_merge
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
=
# Name: combine_counts.sh
# Description: Combine StringTie abundance.tsv counts from each sample directory into a single matrix

# Set the current directory to the directory where this script is located
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Navigate to the parent directory containing all the individual count folders
cd "$DIR/quant_stringtie"

# Step 1: Generate a list of all abundance.tsv files
find . -type d -mindepth 1 -maxdepth 1 | while read dir; do
    echo "$dir/abundance.tsv"
done > tsv_files.txt

# Step 2: Extract headers
# This takes the header from one of the abundance.tsv files 
# (assuming all files have the same order of genes)
awk 'NR==1{print $0}' $(head -1 tsv_files.txt) > counts_matrix.tsv

# Step 3: Combine all abundance.tsv files into a matrix
# This step loops through each file, extracts the counts column 
# and appends it to the counts_matrix.tsv
while read file; do
    # Extract the sample name from the directory name
    sample_name=$(echo $file | cut -d'/' -f2)
    
    # Extract counts from the current file and save to a temporary file
    awk -v name=$sample_name 'NR>1{print $5}' $file > tmp_counts.txt
    
    # Merge the counts into the matrix and replace the old matrix with the new one
    paste counts_matrix.tsv tmp_counts.txt > tmp && mv tmp counts_matrix.tsv
done < tsv_files.txt

# Clean up temporary files
rm tmp_counts.txt tsv_files.txt

echo "Counts combined successfully into counts_matrix.tsv"
