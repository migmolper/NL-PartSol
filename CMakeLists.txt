cmake_minimum_required(VERSION 3.10.2)
project("nl-partsol")

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# Check OS
IF (UNIX)
    message("UNIX system")
    add_definitions("-DSYS_LINUX")
    SET(BLA_VENDOR OpenBLAS)
ELSEIF (APPLE)
    message("MACOSX system")
    add_definitions("-DSYS_MACOSX")
    SET(BLA_VENDOR Apple)
ELSE (UNIX)
    message(FATAL_ERROR "Your system has not been recongnised." )
ENDIF (UNIX)

# Add LAPACK
option(USE_LAPACK "Use LAPACK solver library" ON)

# Add OpenMP
option(USE_OPENMP "Use OpenMP solver library" ON)

# Add PETSC
option(USE_PETSC "Use PETSC solver library" ON)

# General compile settings
IF (NOT CMAKE_BUILD_TYPE)
 #SET(CMAKE_BUILD_TYPE "Debug")
 SET(CMAKE_BUILD_TYPE "Release")
ENDIF (NOT CMAKE_BUILD_TYPE)

add_subdirectory (nl-partsol) 
#add_subdirectory (test)
