#!/usr/bin/env bash

EXE="tsp_enn"
CONFIG="Debug"

if [[ "$1" == "run" ]]; then
    CONFIG="Release"
fi

rm -rfv "${CONFIG}" &&
mkdir -p "${CONFIG}" &&
cmake -DCMAKE_CXX_COMPILER="/usr/bin/clang++" -DCMAKE_C_COMPILER="/usr/bin/clang" -DCMAKE_BUILD_TYPE="${CONFIG}" -S $PWD -B "${CONFIG}" &&
# cmake -DCMAKE_CXX_COMPILER="/usr/bin/g++" -DCMAKE_C_COMPILER="/usr/bin/gcc" -DCMAKE_BUILD_TYPE="${CONFIG}" -S $PWD -B "${CONFIG}" &&
cmake --build "${CONFIG}" &&
mkdir -p bin &&
cp "${CONFIG}/${EXE}" bin/"${EXE}"
# cp "${CONFIG}/compile_commands.json" "${PWD}/compile_commands.json"
