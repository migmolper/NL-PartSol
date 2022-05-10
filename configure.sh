#!/bin/bash

if [ "$(uname)" == "Darwin" ]; then
export PETSC_DIR=$HOME/petsc
export PETSC_ARCH=arch-darwin-c-debug    
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
export PETSC_DIR=$HOME/petsc
export PETSC_ARCH=arch-classic-docs
fi
#elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
#    # Do something under 32 bits Windows NT platform
#elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
#    # Do something under 64 bits Windows NT platform

clear

DIR="build"
if [ ! -d "$DIR" ]; then
  mkdir ${DIR}
fi

cd ${DIR}

FILE="Makefile" 
if [ -f "$FILE" ]; then
    make -k
else 
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DNO_KAHIP=True -DCMAKE_C_FLAGS="-O3"
#    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DNO_KAHIP=True -DCMAKE_C_FLAGS="-O3 -Wunused-variable"
#    cmake .. -DDEBUG_MODE=0 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DNO_KAHIP=True -DCMAKE_C_FLAGS="-O0 -g -Wall -Wpedantic -Wextra -Wunused-variable"
    make -j8

# decomment this to have it verbose
# make VERBOSE=1 -j4
# make -j8 
#    make -k
fi

#find -name '*.c' -o -name '*.h' | xargs clang-format -i --verbose