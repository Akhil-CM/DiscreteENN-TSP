#!/usr/bin/env bash

# COMPILER="clang++"
COMPILER="g++"
"$COMPILER" -std=c++17 -march=native -fstrict-aliasing -Wall -Wextra -pedantic -Weffc++ "$1" &&
sleep 2 &&
time ./a.out
