#!/bin/bash

clear

# Uncomment your preferred IDE
PLATFORM="Unix Makefiles" # Generates standard UNIX makefiles.
#PLATFORM="Ninja" # Generates build.ninja files.
#PLATFORM="Ninja Multi-Config" # Generates build-<Config>.ninja files.
#PLATFORM="Watcom WMake" #  Generates Watcom WMake makefiles.
#PLATFORM="CodeBlocks - Ninja" # Generates CodeBlocks project files.
#PLATFORM="CodeBlocks - Unix Makefiles" # Generates CodeBlocks project files.
#PLATFORM="CodeLite - Ninja" # Generates CodeLite project files.
#PLATFORM="CodeLite - Unix Makefiles" # Generates CodeLite project files.
#PLATFORM="Eclipse CDT4 - Ninja" # Generates Eclipse CDT 4.0 project files.
#PLATFORM="Eclipse CDT4 - Unix Makefiles" # Generates Eclipse CDT 4.0 project files.
#PLATFORM="Kate - Ninja" # Generates Kate project files.
#PLATFORM="Kate - Unix Makefiles" # Generates Kate project files.
#PLATFORM="Sublime Text 2 - Ninja" # Generates Sublime Text 2 project files.
#PLATFORM="Sublime Text 2 - Unix Makefiles" # Generates Sublime Text 2 project files.

## Configure environments
if [[ "$HOSTNAME" == "Hilbert" ]]
 then
 echo "We are in Hilbert (MM's machine)"
 export PETSC_DIR=$HOME/petsc
 export PETSC_ARCH=arch-linux-c-debug
 export PKG_CONFIG_PATH=$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig
 C_COMPILER=/usr/local/openmpi-4.1.2/bin/mpicc

elif [[ "$HOSTNAME" == "MacBook-Pro-de-Miguel.local" ]]
 then
 echo "We are in MM's MacBool-Pro"
 export PETSC_DIR=$HOME/petsc
 export PETSC_ARCH=arch-darwin-c-debug
 export PKG_CONFIG_PATH=$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig
 C_COMPILER=$PETSC_DIR/$PETSC_ARCH/bin/mpicc

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

if [ -f "$FILE" ]; then
    make -k
else 
    cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=${C_COMPILER} \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
    -G "${PLATFORM}"
fi

if [[ "$PLATFORM" == "Unix Makefiles" ]]
then
cd build
make -k
elif [[ "$PLATFORM" == "Ninja" ]]
then
cd build
ninja
fi

