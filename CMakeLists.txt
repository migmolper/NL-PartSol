cmake_minimum_required(VERSION 3.10.2)
project("nl-partsol")

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# Add OpenMP
option(USE_OPENMP "Use OpenMP solver library" OFF)

# Add PETSC
option(USE_PETSC "Use PETSC solver library" ON)

# General compile settings
IF (NOT CMAKE_BUILD_TYPE)
 #SET(CMAKE_BUILD_TYPE "Debug")
 SET(CMAKE_BUILD_TYPE "Release")
ENDIF (NOT CMAKE_BUILD_TYPE)

add_subdirectory (nl-partsol) 
#add_subdirectory (test)
