#!/usr/bin/env bash

# while "$@"; do :; done
count=0
while "$@" &> /dev/null; do
    echo "-------------------------------------------"
    echo "[Info]: Run #${count} currently in progress"
    echo "-------------------------------------------"
    (( count++ )) &&
    sleep 2
done
echo "[Info]: Failed after ${count} runs"

