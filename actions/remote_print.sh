#!/usr/bin/env bash

COUNT=${1:-"100"}
SLEEP_DELAY=${2:-"0.5"}

i="0"
while [ $i -lt ${COUNT} ];
do
    j=$[$i+1]
    if [ $(( $i % 2)) -eq 0 ]; then
        echo "stderr line ${j}" >&2
    else
        echo "stdout line ${j}"
    fi

    i=$[$i+1]

    sleep ${SLEEP_DELAY}
done
