#!/bin/bash

# Description:
# This script aligns paired-end sequencing reads against of your H37Rv deletion mutant a given reference and produces a visualization using samplot.
#
# Usage:
# ./DeeDeeThree.sh <R1_fastq> <R2_fastq> <chromosome> <start_position> <end_position>
#
# Requirements:
# - FastQC
# - BWA
# - SAMtools
# - samplot

# Check if the required number of arguments are provided
if [ "$#" -ne 5 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <R1_fastq> <R2_fastq> <chromosome> <start_position> <end_position>"
    exit 1
fi

threads=6 

# Assign command-line arguments to variables
R1="$1"
R2="$2"
chromosome="$3"
start_position="$4"
end_position="$5"

# Name file
nm="${R1/_R1_001.fastq.gz/}"

# download fasta
curl --silent "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=448814763&&ncbi_phid=CE8DF97C4D5B3E710000000000300026" > ref.fasta
# static
ref="ref.fasta"

# Create a reference index for alignment using bwa
# since using only you own, this will already be done, but you can try make it look for the reference and then make it if it doesnt exist if you want. Otherwise it would remake it each time and overwrite files etc. (human genome take hours to index)
bwa index ${ref}


# Alignment pipeline
fastqc -t $threads "$R1" -o "${nm}_R1_fastqc"  # Quality control analysis on R1
# tee "${nm}_R1_fastqc.log" |  # Save FastQC log
fastqc -t $threads "$R2" -o "${nm}_R2_fastqc"  # Quality control analysis on R2

# Align R1 and R2 to the reference using bwa mem and save the output as a SAM file
bwa mem -t $threads "$ref" "$R1" "$R2" > "${nm}.sam"

# Convert and sort SAM to BAM
samtools view -bS "${nm}.sam" | 
    samtools sort -@ $threads -o "${nm}.bam"

# Index BAM file for visualization and analysis
samtools index "${nm}.bam"

# Cleanup intermediate files (uncomment if needed)
rm "${nm}.sam"

# Print completion message
echo "Alignment pipeline completed."

# Visualization command

# this should be pulled based on the gene name, so read in the gff file, and get the info

samplot plot -n "$nm" -b "${nm}.bam" -o "${nm}.png" -c <chromosome> -s <start_position> -e <end_position> -t "DEL"

# Exit the script
exit 0
