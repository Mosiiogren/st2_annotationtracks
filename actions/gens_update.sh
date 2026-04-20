#!/bin/bash

set -e


ANNOTATIONTRACKS=$1



IFS="," read -a array <<< "${ANNOTATIONTRACKS}"

for file in "${array[@]}"
do
    echo "Uploading $file to Gens"
    docker compose -f docker-compose.yaml exec gens gens load annotations -b 38 -f /data/wgs/annotationtracks/$file
    
done





