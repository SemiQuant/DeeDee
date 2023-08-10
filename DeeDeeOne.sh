# Assign inputs 
R1="$1"
R2="$2"
threads="${3-6}"




 
#Name file, Align reference to mutant, convert sam to bam
fastqc -t $threads "$R1"
nm="${R1/_R1_001.fastq.gz/}"
fastqc -t $threads "$R2"
bwa mem -t $threads "$R1" "$R2" > "${nm}.sam"
samtools view -bS "${nm}.sam" |samtools sort -@ $threads -o "${nm}.bam"
# done

exit 0 