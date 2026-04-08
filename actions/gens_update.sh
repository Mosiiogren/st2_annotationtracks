#!/bin/bash

set -e

DOCKERCOMMAND=$1
MONGOCOMMAND=$2

docker compose -f docker-compose.yaml ${DOCKERCOMMAND}

echo "Working #1"

${MONGOCOMMAND}

echo "Working #2"
