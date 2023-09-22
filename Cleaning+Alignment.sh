#!/bin/bash
FASTQ_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Data/PANCDB/hpapdata/HPAP-051/Islet-Studies/Islet-molecular-phenotyping-studies/Single-cell-RNAseq/Stanford_scRNAseq"
OUTPUT_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Data/PANCDB/Cleaned-N-Aligned-hpapdata"

#Loop over each FASTQ file

for fastq_file in "$FASTQ_DIR"/*-R1_fastq-data.fastq.gz; do
    base_name=$(basename "$fastq_file" -R1_fastq-data.fastq.gz) #Removing suffix to later find the R2 strand
    paired_file="$FASTQ_DIR/${base_name}-R2_fastq-data.fastq.gz"
    echo "base name: $base_name"s
    echo ""
    echo "$paired_file"
done

