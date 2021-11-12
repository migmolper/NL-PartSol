#!/bin/bash

export CC=/usr/bin/mpicc
export CXX=/usr/bin/mpicxx

if [[ -d build ]]
then
    rm -r build
fi

mkdir build

cd build

cmake -D CMAKE_C_COMPILER=$CC CMAKE_CXX_CMPILER=$CXX ..

make -k
