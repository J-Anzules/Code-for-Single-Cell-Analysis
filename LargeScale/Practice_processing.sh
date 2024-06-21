#!/bin/bash

# Input and output file paths
input_file="./ftp_links/ftpLinks.txt" # 6374
output_upenn_1="./ftp_links/upenn_1.txt"
output_upenn_2="./ftp_links/upenn_2.txt"
output_upenn_3="./ftp_links/upenn_3.txt"
output_upenn_4="./ftp_links/upenn_4.txt" # TODO: Do I have to change the HPAP079 to HPAP-079?
md5_file="./ftp_links/md5s.txt"
output_stanford="./ftp_links/stanford.txt"
temp_file="./ftp_links/temp.txt"

# Initialize output files
> "$output_upenn_1"
> "$output_upenn_2"
> "$output_upenn_3"
> "$md5_file"
> "$output_stanford"
> "$temp_file"

upenn_count=0
stanford_count=0
other_count=0
md5_count=0

# Read the input file line by line
# Upenn_scRNAseq/fastq/HPAP-120_FGC2468_111809_S3_L002_R1_001.fastq.gz
while IFS= read -r line; do
    if [[ $line == *"-md5"* ]]; then
        echo "$line" >> "md5_file"
        ((md5_count++))

    elif [[ "$line" == *"Upenn_scRNAseq/fastq/HPAP-"*"_scRNA_"* ]]; then
        echo "$line" >> "$output_upenn_1"
        ((upenn_count++))
        ((upenn_count_1++))

    elif [[ "$line" == *"Upenn_scRNAseq/fastq/HPAP"*"_L"* ]]; then
        echo "$line" >> "$output_upenn_2"
        ((upenn_count++))
        ((upenn_count_2++))

    elif [[ "$line" == *"Upenn_scRNAseq/fastq/HPAP-"* ]]; then
        echo "$line" >> "$output_upenn_3"
        ((upenn_count++))
        ((upenn_count_3++))
    
        ((upenn_count++))
    elif [[ "$line" == *"Stanford_scRNAseq/"* ]]; then
        echo "$line" >> "$output_stanford"
        ((stanford_count++))
    else
        echo "$line" >> "$temp_file"
        ((other_count++))
    fi
done < "$input_file"

# Adding numbers
total=$((upenn_count + stanford_count + other_count + md5_count))

echo "Upenn count = $upenn_count"
echo "Upenn count 1 = $upenn_count_1"
echo "Upenn count 2 = $upenn_count_2"
echo "Upenn count 3 = $upenn_count_3"
echo "-------------------------------"
echo " "
echo "Stanford count = $stanford_count"
echo "Other count = $other_count"
echo "md5 count = $md5_count"
echo "All counts = $total"