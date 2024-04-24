#!/usr/bin/env bash

FLAGS_EXTRA="-lsfml-graphics -lsfml-window -lsfml-system"
# FLAGS_EXTRA=""

FLAGS_RELEASE="-O3 -DNDEBUG"

# COMPILER="clang++"
COMPILER="g++"

if [[ "$1" == "run" ]]; then
    shift
fi

echo "Compiling the following:"
echo "$COMPILER -std=c++17 -march=native -fstrict-aliasing -Wall -Wextra -pedantic $FLAGS_EXTRA $FLAGS_RELEASE $@"
"$COMPILER" -std=c++17 -march=native -fstrict-aliasing -Wall -Wextra -pedantic $FLAGS_EXTRA $FLAGS_RELEASE "$@"
# sleep 2 &&
# time ./a.out
