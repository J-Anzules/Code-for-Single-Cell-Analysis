#!/bin/bash

# -------------------------------------------------------------------------------------->
# File: process_fastq.sh
# Author: Jonathan Anzules
# Date: 11/1/23
# Description: This script processes FASTQ files for a large-scale analysis of beta cell
# heterogeneity in Type 2 Diabetes (T2D). It downloads data from specified FTP links, 
# cleans the FASTQ files using fastp, quantifies transcript abundance using Salmon,
# and compiles metadata information.
# -------------------------------------------------------------------------------------->

# Todo:  Createa a main that will run all of the other scripts
# Todo: Build a download and pass on script? Check with ben on best practices
#
#


# User defined directory containing your FASTQ files
FASTQ_DIR="/path/to/fastq/files"
# Todo this should be a static output directory. I could probably generalize the
# Output in a way that makes every step of the way reproducible.
OUTPUT_DIR="/path/to/output/directory"


# Initialize the metadata file
# TODO: I need to add a check for if the file exists
#           If there is a file, update name with current date.
#           check with Ben on best practice here
METADATA_FILE="$OUTPUT_DIR/metadata.csv"
echo "quant_file,donor,cell,batch" > $METADATA_FILE

# Loop over each FASTQ file
for fastq_file in "$FASTQ_DIR"/*_R1_fastq-data.fastq.gz; do
    base_name=$(basename "$fastq_file" _R1_fastq-data.fastq.gz)
    paired_file="$FASTQ_DIR/${base_name}_R2_fastq-data.fastq.gz"
    
    # Use fastp to clean the FASTQ files
    fastp -i "$fastq_file" -o "${base_name}_clean_R1.fastq.gz" -I "$paired_file" -O "${base_name}_clean_R2.fastq.gz"
    
    # Check if fastp was successful
    if [ $? -eq 0 ]; then
        # Align with salmon
        salmon quant -i /path/to/salmon/index -l A -1 "${base_name}_clean_R1.fastq.gz" -2 "${base_name}_clean_R2.fastq.gz" -o "$OUTPUT_DIR/${base_name}_quant"
        
        # Check if salmon was successful
        if [ $? -eq 0 ]; then
            # Delete the intermediate cleaned FASTQ files
            rm "${base_name}_clean_R1.fastq.gz" "${base_name}_clean_R2.fastq.gz"
            
            # 


            # Extract metadata information from the base_name
            # Here, I'm assuming the base_name is formatted as: DONOR_CELL_BATCH
            # Todo: I should quintuple check that this is the standard output of the data
            IFS='_' read -ra META <<< "$base_name"
            donor="${META[0]}"
            cell="${META[1]}"
            batch="${META[2]}"
            
            # Append to the metadata file
            echo "$OUTPUT_DIR/${base_name}_quant/quant.sf,$donor,$cell,$batch" >> $METADATA_FILE

            # Todo: Add a check for output and delete the raw data. Memory is limited.
        fi
    fi
done
