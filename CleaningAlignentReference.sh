#!/bin/bash

# Path to the directory containing your FASTQ files
FASTQ_DIR="/path/to/fastq/files"
OUTPUT_DIR="/path/to/output/directory"

# Initialize the metadata file
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
            
            # Extract metadata information from the base_name
            # Here, I'm assuming the base_name is formatted as: DONOR_CELL_BATCH
            IFS='_' read -ra META <<< "$base_name"
            donor="${META[0]}"
            cell="${META[1]}"
            batch="${META[2]}"
            
            # Append to the metadata file
            echo "$OUTPUT_DIR/${base_name}_quant/quant.sf,$donor,$cell,$batch" >> $METADATA_FILE
        fi
    fi
done
