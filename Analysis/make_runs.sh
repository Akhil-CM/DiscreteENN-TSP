#!/usr/bin/env bash

num_runs="10000"

log_file="output.log"

echo "Running Discrete ENN ${num_runs} times"
# for i in {0..10..2}
# for i in $(seq 1 2 20)
for i in $(seq 1 $num_runs)
do
    echo "Running ${i}" >> $log_file &&
   ./a.out --batch >> $log_file &&
    sleep 1
done
