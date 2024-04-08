#!/usr/bin/env bash

FLAGS_EXTRA="-lsfml-graphics -lsfml-window -lsfml-system"
# FLAGS_EXTRA=""

# COMPILER="clang++"
COMPILER="g++"
"$COMPILER" -std=c++17 -march=native -fstrict-aliasing -Wall -Wextra -pedantic -Weffc++ $FLAGS_EXTRA "$1" &&
sleep 2 &&
time ./a.out
