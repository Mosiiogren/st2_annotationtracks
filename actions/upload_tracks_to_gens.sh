#!/bin/bash

set -e

ANNOTATIONTRACKS=$1
OUTPUTFOLDER=$2

IFS="," read -a array <<< "${ANNOTATIONTRACKS}"

for file in "${array[@]}"
do
    echo "Copying $file to Gens container"
    docker compose -f docker-compose.yaml cp "${OUTPUTFOLDER}$file" gens:/data/wgs/annotationtracks/
    
done
