#!/usr/bin/env bash

# COMPILER="clang++"
COMPILER="g++"

FLAGS_CXX="-march=native -fstrict-aliasing -Wall -Wextra -pedantic -fsanitize=address,undefined"
CXX_STANDARD="-std=c++17"

FLAGS_EXTRA="-lsfml-graphics -lsfml-window -lsfml-system"
# FLAGS_EXTRA=""

FLAGS_RELEASE=""

if [[ "$1" == "run" ]]; then
    FLAGS_RELEASE="-O3 -DNDEBUG"
    shift
fi

echo "Compiling the following:"
echo "$COMPILER $CXX_STANDARD $FLAGS_CXX $FLAGS_EXTRA $FLAGS_RELEASE $@"
"$COMPILER" $CXX_STANDARD $FLAGS_CXX $FLAGS_EXTRA $FLAGS_RELEASE "$@"
# sleep 2 &&
# time ./a.out
