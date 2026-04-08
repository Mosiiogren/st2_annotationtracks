#!/bin/bash

set -e

DOCKERCOMMAND=$1
MONGOCOMMAND=$2
USERNAME=$3
PASSWORD=$4

docker compose -f docker-compose.yaml ${DOCKERCOMMAND} -u ${USERNAME} -p ${PASSWORD}

echo "Working #1"

${MONGOCOMMAND}

echo "Working #2"
