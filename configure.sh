#!/bin/bash

clear

CC=/usr/bin/gcc
CXX=/usr/bin/g++

DIR="build"
if [ ! -d "$DIR" ]; then
  mkdir ${DIR}
fi

cd ${DIR}

FILE="Makefile"
if [ -f "$FILE" ]; then
    make -k
else 
    cmake ..
    make -k
fi

#find -name '*.c' -o -name '*.h' | xargs clang-format -i --verbose