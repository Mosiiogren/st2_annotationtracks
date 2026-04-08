#!/bin/bash

set -e

COMMAND = $1

docker compose -f docker-compose.yaml ${COMMAND}

echo "Working"
