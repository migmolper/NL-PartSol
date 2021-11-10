#!/bin/bash

export CC=/usr/bin/mpicc
export CXX=/usr/bin/mpicxx

if [[ -d build ]]
then
    rm -r build
fi

mkdir build

cd build

cmake ..

make -k
