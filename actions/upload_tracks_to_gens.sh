#!/bin/bash

set -e

ANNOTATIONTRACKSFILE=$1
GENSFOLDER=$2

docker compose -f docker-compose.yaml exec gens gens load annotations -b 38 -f ${GENSFOLDER}/${ANNOTATIONTRACKSFILE}

echo "Uploading to Gens is Working"
