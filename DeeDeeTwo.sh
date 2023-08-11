#!/bin/bash

# Define input files
R1="$1"
R2="$2"
threads="${3-6}"

# Name file
nm="${R1/_R1_001.fastq.gz/}"

# Alignment pipeline
fastqc -t $threads "$R1" |  # Quality control analysis on R1
    tee "${nm}_R1_fastqc.log" |  # Save FastQC log
    fastqc -t $threads "$R2" -o "${nm}_R2_fastqc"  # Quality control analysis on R2

bwa mem -t $threads "$R1" "$R2" > "${nm}.sam"  # Align R1 and R2, save SAM output

# Convert and sort SAM to BAM
samtools view -bS "${nm}.sam" | 
    samtools sort -@ $threads -o "${nm}.bam"

# Index BAM file for visualization and analysis
samtools index "${nm}.bam"

# Cleanup intermediate files (uncomment if needed)
# rm "${nm}.sam"

# Print completion message
echo "Alignment pipeline completed."

# Exit the script
exit 0
