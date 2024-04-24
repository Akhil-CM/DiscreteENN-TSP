#!/usr/bin/env bash

# COMPILER="clang++"
COMPILER="g++"

CXX_STANDARD="-std=c++17"

FLAGS_CXX="-fstrict-aliasing -Wall -Wextra -pedantic"
FLAGS_CXX_EXTRA="-fsanitize=address,undefined"

FLAGS_RELEASE=""

FLAGS_EXTRA="-lsfml-graphics -lsfml-window -lsfml-system"
# FLAGS_EXTRA=""

if [[ "$1" == "run" ]]; then
    FLAGS_RELEASE="-march=native -O3 -DNDEBUG"
    FLAGS_CXX_EXTRA=""
    shift
fi

echo "Compiling the following:"
echo "$COMPILER $CXX_STANDARD $FLAGS_CXX $FLAGS_CXX_EXTRA $FLAGS_EXTRA $FLAGS_RELEASE $@"
"$COMPILER" $CXX_STANDARD $FLAGS_CXX $FLAGS_CXX_EXTRA $FLAGS_EXTRA $FLAGS_RELEASE "$@"
# sleep 2 &&
# time ./a.out
