#!/bin/bash

clear

## Configure environments
if [[ "$HOSTNAME" == "Hilbert" ]]
 then
 echo "We are in Hilbert (MM's machine)"
 export PETSC_DIR=$HOME/petsc
 export PETSC_ARCH=arch-linux-c-debug
 export PKG_CONFIG_PATH=$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig
# C_COMPILER=/usr/bin/gcc
 C_COMPILER=$PETSC_DIR/$PETSC_ARCH/bin/mpicc
# C_COMPILER=/usr/bin/clang-14

elif [[ "$HOSTNAME" == "istorage-00.hpc.cica.es" ]]
 then
 echo "We are in istorage (PA's cluster)"
 module load mpi/mpich-3.2-x86_64
 module load petsc-3.16.0
 export PKG_CONFIG_PATH=$PETSC_DIR/lib/pkgconfig
 C_COMPILER=mpicc
fi

## Check compiers
if ! command -v ${C_COMPILER} &> /dev/null
then
    echo -e ""${C_COMPILER}": "${RED}" False "${RESET}""
    echo "please, ask your administrator to install "${C_COMPILER}""
    exit   
else
    echo -e ""${C_COMPILER}": "${GREEN}" True "${RESET}""
fi

## Check if clang-format is installed
if ! command -v clang-format &> /dev/null
then
    echo -e "clang-format: "${RED}" False "${RESET}""
    echo "please, ask your administrator to run ''sudo apt-get install clang-tools''"
    exit   
else 
    echo -e "clang-format: "${GREEN}" True "${RESET}" "
fi

## Check if cmake is installed
if ! command -v cmake &> /dev/null
then
    echo "cmake could not be found"
    echo "please, ask your administrator to run sudo apt-get install make"
    exit   
fi

## Format the project to make it prettier :-)
#find -name '*.c' -o -name '*.h' | xargs clang-format -i --verbose

DIR="build"
if [ ! -d "$DIR" ]; then
  mkdir ${DIR}
fi

## Navigate inside of build
cd ${DIR}


## Configure and compile the project
FILE="Makefile" 
C_FLAGS="-O3"
# "-O3"
# "-O0 -g -Wall -Wpedantic -Wextra -Wunused-variable"
# "-O3 -Wunused-variable"

if [ -f "$FILE" ]; then
    make -k
else 
    cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=${C_COMPILER} \
    -DCMAKE_C_FLAGS="${C_FLAGS}"

    make -j8

# decomment this to have it verbose
# make VERBOSE=1 -j4
# make -j8 
#    make -k
fi
