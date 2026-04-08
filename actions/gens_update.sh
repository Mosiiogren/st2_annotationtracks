#!/bin/bash

set -e

DOCKERCOMMAND=$1
MONGOCOMMAND=$2
INPUTFIRST=$3
INPUTSECOND=$4

#docker compose -f docker-compose.yaml ${DOCKERCOMMAND}
docker compose -f docker-compose.yaml ${DOCKERCOMMAND} --username=${INPUTFIRST} --password=${INPUTSECOND} ${MONGOCOMMAND}

echo "Working #"

