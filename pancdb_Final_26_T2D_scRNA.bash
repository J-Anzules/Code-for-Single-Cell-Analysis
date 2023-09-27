#!/usr/bin/env bash
#=========================================================================================
# I'm going to edit this script to download each file, process the data, then 
# , run fastp and salmon to clean and align the cells
#
# List of TO DO's
# - Make sure all data temp (Including FTP stuff) files are done outside of the code folder
#   - Somehow $p2 is what would need to be used to keep track of the output directory
#       directory that holds the data that needs to be processed
#   - Use $p2 to download all the riht data
#   - figure out what is wrong with the fastq folder and why isn't it being transfered over to the data_save_loc variable






# This script downloads HPAP experiments in your current working directory:
# * The data files (if requested) will be saved in "hpapdata" sub-folder;
# * Metadata files (if requested) will be saved in "metadata" sub-folder.
#
# IMPORTANT: Before starting to download experiment data, please:
# ==============================================================================
# (1) Register your SSH PUBLIC key at: https://hpap.pmacs.upenn.edu/explore/ftp
# (2) Modify "SSH_KEYFILE" on the next uncommented line to the path of your own
#     SSH PRIVATE key.
# Otherwise the script WON'T work.
# Read this article on SSH public key authentication:
# https://www.ssh.com/academy/ssh/public-key-authentication

SSH_KEYFILE="/home/jon/.ssh/id_rsa"

if [[ ! -r "${SSH_KEYFILE}" ]]; then
    echo "Error: SSH private key file not found"
    exit 1
fi

FILES="
hpapdata/HPAP-001/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq"
# hpapdata/HPAP-007/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-010/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-013/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-051/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-057/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-058/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-061/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-065/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-070/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-079/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-081/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-083/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-085/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-088/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-090/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-091/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-096/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-100/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-106/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-108/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-109/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-111/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-120/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-124/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# hpapdata/HPAP-126/Islet Studies/Islet molecular phenotyping studies/Single-cell RNAseq/Upenn_scRNAseq/fastq
# "

# Set IFS (Internal Field Separator)
IFS_BAK=$IFS
IFS=$'\n'

for f in $FILES; do
    p=$(dirname $f)
    echo "f variable:"
    echo "$f"
    echo "p variable:"
    echo "$p"
    p2=${p//' '/'-'} # escape space character in dir name
    data_save_loc="../Data/TestingDownload/PANCB/${p2}"
    echo  "File location: ${data_save_loc}"
    # mkdir -p "${data_save_loc}"
    f2=${f//' '/'\ '} # escape space character in file name
    # sftp -i ${SSH_KEYFILE} -r hpapsftp@hpap-test.pmacs.upenn.edu:"$f2" "$data_save_loc"
    
done

# Recover original IFS
IFS=$IFS_BAK
unset IFS_BAK


echo "------------------------"
echo "------------------------"
echo "------------------------"
ls "$data_save_loc"


# echo; echo "$(date): experiment data downloaded"; echo

# # Download metadata files into ./metadata
# META_URL="https://hpap.pmacs.upenn.edu/assets/metadata/latest"
# test -d ./metadata || mkdir ./metadata
# cd ./metadata

# curl --silent "$META_URL/PancDB_Donors.xlsx" -O
# curl --silent "$META_URL/PancDB_scRNA-seq_metadata_2022-06-15.xlsx" -O
# curl --silent "$META_URL/README.xlsx" -O

# cd ..
# echo "$(date): metadata downloaded"
