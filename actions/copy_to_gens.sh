#!/bin/bash

set -e

ANNOTATIONTRACKS=$1
GENSFOLDER=$2


docker -f docker-compose.yaml cp ${ANNOTATIONTRACKS} gens:${GENSFOLDER}

echo "Copy to Gens is Working"
