#!/usr/bin/env bash

# while "$@"; do :; done
count=0
while "$@" 1> /dev/null; do
    (( count++ )) &&
    echo "-------------------------------------------"
    echo "[Info]: Run #${count} completed"
    echo "-------------------------------------------"
    sleep 2
done
echo "[Info]: Failed after ${count} runs"

