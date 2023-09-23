#!/bin/bash

count=0

file_pwd="/mnt/c/Users/jonan/Documents/1Work/scWork/Data/PANCDB/hpapdata/HPAP-051/Islet-Studies/Islet-molecular-phenotyping-studies/Single-cell-RNAseq/Stanford_scRNAseq/"
for fastq_file in "$file_pwd"/*-R1_fastq-data.fastq.gz; do
        # Adding and testing the count
        echo "$count"
        echo "$fastq_file"
        ((count++))
        if [[ $count -ge 2 ]]; then
            break
        fi
done