#!/bin/bash

set -e

ANNOTATIONTRACKS=$1
GENSFOLDER=$2


docker compose -f docker-compose.yaml docker cp ${ANNOTATIONTRACKS} gens:${GENSFOLDER}

echo "Copy to Gens is Working"
