cmake_minimum_required(VERSION 3.10.2)
project("nl-partsol")

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# Check OS
IF (APPLE)
    message("MACOSX system")
    add_definitions("-DSYS_MACOSX")
    SET(BLA_VENDOR Apple)
ELSEIF (UNIX)
    message("UNIX system")
    add_definitions("-DSYS_LINUX")
    SET(BLA_VENDOR OpenBLAS)
ELSE ()
    message(FATAL_ERROR "Your system has not been recongnised." )
ENDIF ()

# Plain strain conditions
option(USE_PLAINSTRAIN "Use Plain-Strain conditions" ON)

# Add LAPACK
option(USE_LAPACK "Use LAPACK solver library" ON)

# Add MPI
option(USE_MPI "Use MPI library" OFF)

# Add OpenMP
option(USE_OPENMP "Use OpenMP solver library" ON)

# Add PETSC
option(USE_PETSC "Use PETSC solver library" ON)

# Compile settings

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wpedantic -Wextra -Wunused-variable")
set(CMAKE_C_FLAGS_RELEASE "-Ofast")

add_subdirectory (nl-partsol) 
#add_subdirectory (test)
