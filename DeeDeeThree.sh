#!/bin/bash

# Define input files
R1="$1"
R2="$2"
threads="${3-6}"

# Name file
nm="${R1/_R1_001.fastq.gz/}"

# static
ref="path to reference fasta"

# Create a reference index for alignment using bwa
# since using only you own, this will already be done, but you can try make it look for the reference and then make it if it doesnt exist if you want. Otherwise it would remake it each time and overwrite files etc. (human genome take hours to index)
#bwa index ${ref}


# Alignment pipeline
fastqc -t $threads "$R1" -o "${nm}_R1_fastqc"  # Quality control analysis on R1
# tee "${nm}_R1_fastqc.log" |  # Save FastQC log
fastqc -t $threads "$R2" -o "${nm}_R2_fastqc"  # Quality control analysis on R2

# Align R1 and R2 to the reference using bwa mem and save the output as a SAM file
bwa mem -t $threads 'name of reference' "$R1" "$R2" > "${nm}.sam"

# Convert and sort SAM to BAM
samtools view -bS "${nm}.sam" | 
    samtools sort -@ $threads -o "${nm}.bam"

# Index BAM file for visualization and analysis
samtools index "${nm}.bam"

# Cleanup intermediate files (uncomment if needed)
# rm "${nm}.sam"

# Print completion message
echo "Alignment pipeline completed."


# Visualization command
samplot plot -n "$nm" -b "${nm}.bam" -o "${nm}.png" -c <chromosome> -s <start_position> -e <end_position> -t <plot_type>

# Exit the script
exit 0
