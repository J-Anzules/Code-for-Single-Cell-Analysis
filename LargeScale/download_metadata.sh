#!/bin/bash

# Set IFS (Internal Field Separator)
IFS_BAK=$IFS
IFS=$'\n'

for f in $FILES; do
    echo "[$(date -Iseconds)] downloading $(basename $f) ..."
    encoded_f="$(echo $f | sed 's/ /%20/g')"
    # curl --silent --create-dirs --output $f ${DATA_SERVER}/${encoded_f}
done

# Download metadata files into ./metadata
META_URL="https://hpap.pmacs.upenn.edu/assets/metadata/latest"
test -d ./metadata || mkdir ./metadata
cd ./metadata

curl --silent "$META_URL/PancDB_Donors.xlsx" -O
curl --silent "$META_URL/PancDB_Stanford_scRNA-seq_metadata_2022-06-15.xlsx" -O
curl --silent "$META_URL/PancDB_scRNA-seq_metadata_2023-12-22.xlsx" -O
curl --silent "$META_URL/README.xlsx" -O