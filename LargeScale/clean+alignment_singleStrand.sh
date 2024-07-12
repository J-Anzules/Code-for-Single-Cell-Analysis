#!/bin/bash
############################################################################################
# Script Name: Single Cell RNA-seq Processing Pipeline
# Description: This script processes single-cell RNA-seq 
#              data by cleaning the FASTQ files with fastp 
#              and then aligning with Salmon. Additionally, 
#              metadata is extracted and saved.
# Author: Jonathan Anzules
# Contact: jonanzule@gmail.com
# Date Created: July 10, 2023
# Usage: ./Cleaning+Alignment.sh
# Dependencies: 
#   - fastp
#   - salmon
#   - Bash
############################################################################################


OUTPUT_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta Cell Study/Data/LargeScale/quant_files_healthy/upenn1"
SALM_INDEX="/mnt/c/Users/jonan/Documents/Genomes/Homo-Sapiens/GRCh38/GRCh38-Transcript-Salmon-index"
DATA_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta Cell Study/Data/LargeScale/raw_download_healthy"
TEMP_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta Cell Study/Data/LargeScale"
scMETADATA_NAME="metadata_upenn1_healthy"

# Initialize the metadata file
METADATA_FILE="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta Cell Study/Data/LargeScale/${scMETADATA_NAME}.csv"
echo "quant_file,donor,cell" > "$METADATA_FILE"

# Get the first 5 FASTQ files and save them to an array
# IFS=$'\n' read -d '' -r -a FASTQ_FILES < <(ls -1 "${DATA_DIR}"/*.fastq* | head -n 10)
#For all files
IFS=$'\n' read -d '' -r -a FASTQ_FILES < <(ls -1 "${DATA_DIR}"/*.fastq*)


# Print the selected files (optional)
echo "Selected FASTQ files:"
for fastq_file in "${FASTQ_FILES[@]}"; do
    echo "${fastq_file}"
done

# Initializing the temporary folder
mkdir -p "$TEMP_DIR"/temp

# Loop over each FASTQ file
for fastq_file in "${FASTQ_FILES[@]}"; do
    echo "#############################################################"
    echo "#############################################################"
    echo "Processing file: $fastq_file"
    echo "#############################################################"
    echo "#############################################################"

    base_name=$(basename "$fastq_file" _fastq-data.fastq.gz) # Removing suffix

    echo "--------------------"
    echo "file name:"
    echo "$base_name"
    echo ""
    echo "Fastq file pwd:"
    echo "$fastq_file"
    echo "--------------------"
    echo ""
    echo ""
    echo ""
    
    # Use fastp to clean the FASTQ files
    # File will be made in local directory, deleted later 
    fastp   -i "$fastq_file" \
            -o "$TEMP_DIR/temp/${base_name}_clean.fastq.gz" \
            --json /dev/null \
            --html /dev/null

     # Check if fastp was successful
    if [ $? -eq 0 ]; then
        # Align with salmon
        salmon quant    -i "$SALM_INDEX" \
                        -l A \
                        -r "$TEMP_DIR/temp/${base_name}_clean.fastq.gz" \
                        -o "$OUTPUT_DIR/${base_name}_quant"
        # Check if salmon was successful
        if [ $? -eq 0 ]; then
            # rm "./temp/${base_name}_clean.fastq.gz"
            # I will be using the part of the file name to name the cell. OG metadata has a similar cellID
            IFS='_' read -ra META <<< "$base_name"
            donor="${META[0]}"
            cell="${META[2]}"
                    
            # IFS='-' read -ra META2 <<< "$cell"
            # cell="${META2[1]}"        
            # Append to the metadata file
            echo "$OUTPUT_DIR/${base_name}_quant/quant.sf,$donor,$cell" >> "$METADATA_FILE"
        fi
    fi
  
    
    echo "--------------------"
done

# Removing temp folder
rm -r "$TEMP_DIR"/temp
