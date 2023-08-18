#!/bin/bash
export PATH=$PATH:/home/semiquant/miniconda3/bin/
outdir="$1"
cd "$outdir"
ref="/var/www/html/scripts/deedee/refs/${3-Mycobacterium_tuberculosis_H37Rv_genome_v4.fasta}"
gene_in="${4-Rv0678}"
threads="4"
gff="${ref/_genome_v4.fasta/_gff_v4.gff}"

start="none"
line=$(grep -i -E "Locus=$gene_in|Name=$gene_in;" "$gff")
start=$(echo "$line" | awk -F '\t' '{print $4}')

if [[ "$start" =~ ^[0-9]+$ ]]; then
  end=$(echo "$line" | awk -F '\t' '{print $5}')
  chr=$(echo "$line" | awk -F '\t' '{print $1}')
else
  echo "Cant find input gene: $gene_in"
  echo "Available genes for reference are:"
  awk -F';Name=' '{print $1}' "$gff" | awk '{print $8}'
  exit 1
fi

# Create a reference index for alignment using bwa
# since using only you own, this will already be done, but you can try make it look for the reference and then make it if it doesnt exist if you want. Otherwise it would remake it each time and overwrite files etc. (human genome take hours to index)
#bwa index ${ref}


for R1 in $(ls -1 | grep -E '^.+_[Rr][1].*\.fastq\.gz$')
  do
    R2="${R1/_R1_/_R2_}"; R2="${R1/_r1_/_r2_}"
    nm="${R1/.*/}"
    # nm="$(basename $nm)"
   
    # Alignment pipeline
    fastqc -t $threads "$R1"
    fastqc -t $threads "$R2"
   
    # Align R1 and R2 to the reference using bwa mem and save the output as a SAM file
    bwa mem -t $threads "$ref" "$R1" "$R2" > "${nm}.sam"

    # Convert and sort SAM to BAM
    samtools view -bS "${nm}.sam" | samtools sort -@ $threads -o "${nm}.bam"

    # Index BAM file for visualization and analysis
    samtools index "${nm}.bam"

    samplot plot -n "$nm" -b "${nm}.bam" -o "${nm}.png" -c "$chr" -s $((start - 1000)) -e $((end + 1000)) -t "DEL"

    # Cleanup intermediate files
    rm "$R1" "$R2" "${nm}.sam"
done

email="$2"
/var/www/html/scripts/upload_to_user.sh "$outdir" "DeeDee" "$email"

exit 0




