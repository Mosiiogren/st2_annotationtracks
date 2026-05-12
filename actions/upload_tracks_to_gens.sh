#!/bin/bash

set -e

ANNOTATIONTRACKS=$1
OUTPUTFOLDER=$2

IFS="," read -a array <<< "${ANNOTATIONTRACKS}"

for file in "${array[@]}"
do
    echo "Copying $file to Gens"
    docker compose -f docker-compose.yaml cp "${OUTPUTFOLDER}$file" gens:/data/wgs/annotationtracks/
    echo "Uploading $file to Gens"
    docker compose -f docker-compose.yaml exec gens gens load annotations -b 38 -f /data/wgs/annotationtracks/$file
done
