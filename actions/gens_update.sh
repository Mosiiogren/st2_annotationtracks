#!/bin/bash

set -e

COMMAND = $0

docker compose -f docker-compose.yaml ${COMMAND}

echo "Working"
