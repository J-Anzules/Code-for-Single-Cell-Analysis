#!/bin/bash

FILE_LIST="/mnt/c/Users/jonan/Documents/1Work/scWork/Code-for-Single-Cell-Analysis/file_list_T2D.txt"
FASTQ_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Data/PANCDB/hpapdata/HPAP-051/Islet-Studies/Islet-molecular-phenotyping-studies/Single-cell-RNAseq/Stanford_scRNAseq"
OUTPUT_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Data/PANCDB/Cleaned-N-Aligned-hpapdata"
scMETADATA_NAME="metadata_T2D"

#Testing count


# Initialize the metadata file
METADATA_FILE="$OUTPUT_DIR/${scMETADATA_NAME}.csv"
echo "quant_file,donor,cell" > $METADATA_FILE


#Loop over each FASTQ file
for file_pwd in $(cat $FILE_LIST); do
    echo "Processing file: $file_pwd"

    for fastq_file in "$file_pwd"/*-R1_fastq-data.fastq.gz; do
        base_name=$(basename "$fastq_file" -R1_fastq-data.fastq.gz) #Removing suffix to later find the R2 strand
        paired_file="$FASTQ_DIR/${base_name}-R2_fastq-data.fastq.gz"
        
        echo "--------------------"
        echo "file name:"
        echo "$base_name"
        echo ""
        echo "Fastq file pwd:"
        echo "$fastq_file"
        echo ""
        echo "R2 strand pwd"
        echo "$paired_file"
        echo "--------------------"
        echo ""
        echo ""
        echo ""
        
        # Use fastp to clean the FASTQ files
        fastp   -i "$fastq_file" \
                -o "${base_name}_clean_R1.fastq.gz" \
                -I "$paired_file" \
                -O "${base_name}_clean_R2.fastq.gz"
        rm "${base_name}_clean_R1.fastq.gz" "${base_name}_clean_R2.fastq.gz"
        # I will be using the part of the file name to name the cell. If this doesn't work, the meta data has a similar cellID
        IFS='_' read -ra META <<< "$base_name"
                donor="${META[0]}"
                cell="${META[2]}"
                
        IFS='-' read -ra META2 <<< "$cell"
                cell="${META2[1]}"        
                # Append to the metadata file
                echo "$OUTPUT_DIR/${base_name}_quant/quant.sf,$donor,$cell" >> $METADATA_FILE
        # echo "--------------------"
    done
done