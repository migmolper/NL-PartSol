#!/bin/zsh

export CC=$PETSC_DIR/$PETSC_ARCH/bin/mpicc
export CXX=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx

if [[ -d build ]]
then
    rm -r build
fi

mkdir build

cd build

cmake -D CMAKE_C_COMPILER=$CC CMAKE_CXX_COMPILER=$CXX ..

make -k
