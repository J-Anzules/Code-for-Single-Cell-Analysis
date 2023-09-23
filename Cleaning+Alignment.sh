#!/bin/bash
############################################################################################
# Script Name: Single Cell RNA-seq Processing Pipeline
# Description: This script processes single-cell RNA-seq 
#              data by cleaning the FASTQ files with fastp 
#              and then aligning with Salmon. Additionally, 
#              metadata is extracted and saved.
# Author: Jonathan Anzules
# Contact: jonanzule@gmail.com
# Date Created: September 23, 2023
# Usage: ./Cleaning+Alignment.sh
# Dependencies: 
#   - fastp
#   - salmon
#   - Bash
############################################################################################


FILE_LIST="/mnt/c/Users/jonan/Documents/1Work/scWork/Code-for-Single-Cell-Analysis/file_list_T2D.txt"
FASTQ_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Data/PANCDB/hpapdata/HPAP-051/Islet-Studies/Islet-molecular-phenotyping-studies/Single-cell-RNAseq/Stanford_scRNAseq"
OUTPUT_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Data/PANCDB/Cleaned-N-Aligned-hpapdata"
SALM_INDEX="/mnt/c/Users/jonan/Documents/Genomes/Homo-Sapiens/GRCh38/GRCh38-Transcript-Salmon-index"

scMETADATA_NAME="metadata_T2D"

#Testing count
# echo "How many times should we test?"
# read count_end
# count=0

# Initialize the metadata file
METADATA_FILE="$OUTPUT_DIR/${scMETADATA_NAME}.csv"
echo "quant_file,donor,cell" > $METADATA_FILE

# Initializing the temporary folder
mkdir temp

#Loop over each FASTQ file
for file_pwd in $(cat $FILE_LIST); do
    echo "Processing file: $file_pwd"

    for fastq_file in "$file_pwd"/*-R1_fastq-data.fastq.gz; do
        # # Adding and testing the count
        # ((count++))
        # if [[ $count -ge $count_end ]]; then
        #     break
        # fi

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
        # File will be made in local directory, deleted later 
        fastp   -i "$fastq_file" \
                -o "./temp/${base_name}_clean_R1.fastq.gz" \
                -I "$paired_file" \
                -O "./temp/${base_name}_clean_R2.fastq.gz" \
                --json /dev/null \
                --html /dev/null

         # Check if fastp was successful
        if [ $? -eq 0 ]; then
            # Align with salmon
            salmon quant    -i "$SALM_INDEX" \
                            -l A \
                            -1 ./temp/"${base_name}_clean_R1.fastq.gz" \
                            -2 ./temp/"${base_name}_clean_R2.fastq.gz" \
                            -o "$OUTPUT_DIR/${base_name}_quant"

            if [ $? -eq 0 ]; then
                rm "./temp/${base_name}_clean_R1.fastq.gz" "./temp/${base_name}_clean_R2.fastq.gz"
            fi
        fi
      
        # I will be using the part of the file name to name the cell. OG metadata has a similar cellID
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
# removing temp folder
rm -r temp