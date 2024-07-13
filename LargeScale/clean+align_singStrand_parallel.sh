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
#   - GNU Parallel
############################################################################################

# OUTPUT_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/Testing"
OUTPUT_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/quant_files_healthy/upenn1"
SALM_INDEX="/mnt/c/Users/jonan/Documents/Genomes/Homo-Sapiens/GRCh38/GRCh38-Transcript-Salmon-index"
DATA_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/raw_download_healthy"
TEMP_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale"
scMETADATA_NAME="metadata_upenn1_healthy"
FASTQ_FILES_TXT="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Code-for-Single-Cell-Analysis/LargeScale/ftp_links/healthy/upenn_1_location.txt"

# Initialize the metadata file
METADATA_FILE="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta_Cell_Study/Data/LargeScale/${scMETADATA_NAME}.csv"
echo "quant_file,donor,cell" > "$METADATA_FILE"

# Read the file paths from the text file into an array
IFS=$'\n' read -d '' -r -a FASTQ_FILES < "$FASTQ_FILES_TXT"

# Print the selected files (optional)
echo "Selected FASTQ files:"
for fastq_file in "${FASTQ_FILES[@]}"; do
    echo "${fastq_file}"
done

if [ -d "$TEMP_DIR/temp" ]; then
  rm -r "$TEMP_DIR/temp"
fi

# Initializing the temporary folder
mkdir -p "$TEMP_DIR"/temp

# Function to process each FASTQ file
process_fastq() {
    fastq_file="$1"
    TEMP_DIR="$2"
    OUTPUT_DIR="$3"
    SALM_INDEX="$4"
    METADATA_FILE="$5"


    echo "#############################################################"
    echo "Processing file: $fastq_file"
    echo "#############################################################"

    base_name=$(basename "$fastq_file" _fastq-data.fastq.gz) # Removing suffix

    # Use fastp to clean the FASTQ files
    fastp -i "$fastq_file" -o "$TEMP_DIR/temp/${base_name}_clean.fastq.gz" --json /dev/null --html /dev/null

    # Check if fastp was successful
    if [ $? -eq 0 ]; then
        # Align with salmon
        salmon quant -i "$SALM_INDEX" -l A -r "$TEMP_DIR/temp/${base_name}_clean.fastq.gz" -o "$OUTPUT_DIR/${base_name}_quant"

        # Check if salmon was successful
        if [ $? -eq 0 ]; then
            IFS='_' read -ra META <<< "$base_name"
            donor="${META[0]}"
            cell="${META[2]}"
                    
            # Append to the metadata file
            echo "$OUTPUT_DIR/${base_name}_quant/quant.sf,$donor,$cell" >> "$METADATA_FILE"
        fi
    fi

    echo "--------------------"
}

export -f process_fastq

# Initialize the temporary folder
mkdir -p "$TEMP_DIR"/temp

# Run the processing in parallel, split into 3 sections
parallel -j 4 process_fastq {} "$TEMP_DIR" "$OUTPUT_DIR" "$SALM_INDEX" "$METADATA_FILE" ::: "${FASTQ_FILES[@]}"

# Removing temp folder
rm -r "$TEMP_DIR"/temp
