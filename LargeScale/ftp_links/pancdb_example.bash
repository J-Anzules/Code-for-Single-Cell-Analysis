#!/usr/bin/env bash
#
# This script downloads HPAP data/metadata into your working directory:
#   * Data (if requested) will be saved in "hpapdata" subdirectory;
#   * Metadata (if requested) will be saved in "metadata" subdirectory.


## ToDO
# 1. This file must take a number entry to test the pipeline for the first x number of files


DATA_SERVER="https://hpap.faryabilab.com"
DATA_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta Cell Study/Data/LargeScale/raw_download"
METADATA_DIR="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta Cell Study/Data/LargeScale"
FILES="/mnt/c/Users/jonan/Documents/1Work/scWork/Beta Cell Study/Code-for-Single-Cell-Analysis/LargeScale/ftp_links/upenn_1.txt"

# TODO: This must be chanced in a way that it takes in the value from another script. 

# Set IFS (Internal Field Separator)
IFS_BAK=$IFS
IFS=$'\n'

while IFS= read -r f; do
    echo "[$(date -Iseconds)] downloading $(basename $f) ..."
    # encoded_f="$(echo $f | sed 's/ /%20/g')"
    
    # This is to clean up the name of the file so it can work right
    f=$(echo "$f" | tr -d '\r\n')
    encoded_f=$(python3 -c "import urllib.parse; print(urllib.parse.quote('''$f'''))")
    # echo "Full file to download"
    # echo "${DATA_SERVER}/${encoded_f}"
    # echo " "
    # echo "File location"
    # echo "$DATA_DIR/$(basename $f)"
    curl --silent --create-dirs --output "$DATA_DIR/$(basename $f)" ${DATA_SERVER}/${encoded_f}

    if [ $? -ne 0 ]; then
        echo "Error downloading $f"
    fi

done < "$FILES"

# Recover original IFS
IFS=$IFS_BAK
unset IFS_BAK

# echo; echo "[$(date -Iseconds)] experiment data downloaded"; echo

# # Download metadata files into the metadata directory
# META_URL="https://hpap.pmacs.upenn.edu/assets/metadata/latest"
# test -d "$METADATA_DIR/metadata" || mkdir -p "$METADATA_DIR/metadata"
# cd "$METADATA_DIR/metadata"

# curl --verbose "$META_URL/PancDB_Donors.xlsx" -O
# curl --silent "$META_URL/PancDB_Stanford_scRNA-seq_metadata_2022-06-15.xlsx" -O
# curl --silent "$META_URL/PancDB_scRNA-seq_metadata_2023-12-22.xlsx" -O
# curl --silent "$META_URL/README.xlsx" -O

# cd ..
# echo "[$(date -Iseconds)] metadata downloaded"
