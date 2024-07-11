#!/bin/bash

#Root ftp location
root_links="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta Cell Study/Code-for-Single-Cell-Analysis/LargeScale/ftp_links/Organized links/Healthy/"


# Input and output file paths
input_file="/mnt/data/ftpLinks.txt" # 6374
output_upenn_1="/mnt/data/upenn_1.txt"
output_upenn_2="/mnt/data/upenn_2.txt"
output_stanford="/mnt/data/stanford.txt"
temp_file="/mnt/data/temp.txt"

# Initialize output files
> "$output_upenn_1"
> "$output_upenn_2"
> "$output_stanford"
> "$temp_file"

# Read the input file line by line
while IFS= read -r line; do
    if [[ "$line" == *"Upenn_scRNAseq/fastq/HPAP-"*"_scRNA_"* ]]; then
        echo "$line" >> "$output_upenn_1"
    elif [[ "$line" == *"Upenn_scRNAseq/fastq/HPAP-"*"_10xscRNA_"* || "$line" == *"Upenn_scRNAseq/fastq/HPAP-"*"_scRNAseq-md5.txt"* ]]; then
        echo "$line" >> "$output_upenn_2"
    elif [[ "$line" == *"Stanford_scRNAseq/fastq/HPAP-"* ]]; then
        echo "$line" >> "$output_stanford"
    else
        echo "$line" >> "$temp_file"
    fi
done < "$input_file"

# Replace the original file with the temp file
mv "$temp_file" "$input_file"

echo "Separation complete. Check the following files for results:"
echo "Upenn data (type 1): $output_upenn_1"
echo "Upenn data (type 2): $output_upenn_2"
echo "Stanford data: $output_stanford"
