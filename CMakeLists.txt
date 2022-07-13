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

# Add LAPACK
option(USE_LAPACK "Use LAPACK solver library" ON)

# Add OpenMP
option(USE_OPENMP "Use OpenMP solver library" ON)

# Add PETSC
option(USE_PETSC "Use PETSC solver library" ON)

# Add test
option(NL_PARTSOL_BUILD_TESTING "Enable testing" ON)


# Compile settings
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wpedantic -Wextra -Wunused-variable -fprofile-arcs -ftest-coverage")
set(CMAKE_C_FLAGS_RELEASE "-Ofast -fprofile-arcs -ftest-coverage")

###############################################################################
## dependencies ###############################################################
###############################################################################

# Add math library
link_libraries(m)

# Add PkgConfig
find_package(PkgConfig REQUIRED)

# Add LAPACK
if (USE_LAPACK)
  find_package(LAPACK REQUIRED)
  
  if(LAPACK_FOUND)
    link_libraries(-llapacke -llapack -lblas)
    add_definitions("-DUSE_LAPACK")
  endif()

endif()

# Add OpenMP
if (USE_OPENMP)
  find_package(OpenMP REQUIRED)

  if (OPENMP_FOUND)
    link_libraries(OpenMP::OpenMP_C)
    add_definitions("-DUSE_OPENMP")
  endif()

endif()

# Add PETSC
if (USE_PETSC)
  pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
  
  if (PETSC_FOUND)
    link_libraries(PkgConfig::PETSC)
    add_definitions("-DUSE_PETSC")
  endif()

endif()

# Add check library
link_libraries(check)
link_libraries(subunit)



###############################################################################
## file globbing ##############################################################
###############################################################################

# these instructions search the directory tree when cmake is
# invoked and put all files that match the pattern in the variables 
# `sources` and `data`

file(GLOB_RECURSE NL_PartSol_src
    nl-partsol/src/Constitutive/*.c 
    nl-partsol/src/Formulations/*.c 
    nl-partsol/src/InOutFun/*.c 
    nl-partsol/src/Linear-Solvers/*.c 
    nl-partsol/src/Matlib/*.c 
    nl-partsol/src/Nodes/*.c 
    nl-partsol/src/Particles/*.c)

file(GLOB_RECURSE DATA build/*)
# you can use set(NL_PartSol_src src/main.c) etc if you don't want to
# use globing to find files automatically


###############################################################################
## target definitions #########################################################
###############################################################################

# add the data to the target, so it becomes visible in some IDE
add_executable(${PROJECT_NAME} nl-partsol/src/driver-nl-partsol.c ${NL_PartSol_src})

# Add include folder
target_include_directories(${PROJECT_NAME} PUBLIC nl-partsol/include)

# this copies all resource files in the build directory
# we need this, because we want to work with paths relative to the executable
# file(COPY ${DATA} DESTINATION build)


###############################################################################
## Unit test ##################################################################
###############################################################################
if(NL_PARTSOL_BUILD_TESTING)
  
  enable_testing()

  file(GLOB_RECURSE DATA build/*)
  
  add_executable(test-nl-partsol nl-partsol/tests/test-nl-partsol.c ${NL_PartSol_src})
  target_include_directories(test-nl-partsol PUBLIC 
  nl-partsol/include
  nl-partsol/tests/Nodes)

  add_test(NAME test-nl-partsol COMMAND $<TARGET_FILE:test-nl-partsol>)


endif()


