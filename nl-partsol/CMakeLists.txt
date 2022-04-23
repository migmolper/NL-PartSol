
###############################################################################
## file globbing ##############################################################
###############################################################################

# these instructions search the directory tree when cmake is
# invoked and put all files that match the pattern in the variables 
# `sources` and `data`
file(GLOB_RECURSE NL_PartSol_src src/*.c include/*.h)
file(GLOB_RECURSE DATA ../build/*)
# you can use set(NL_PartSol_src src/main.c) etc if you don't want to
# use globing to find files automatically


###############################################################################
## target definitions #########################################################
###############################################################################

# add the data to the target, so it becomes visible in some IDE
add_executable(${PROJECT_NAME} ${NL_PartSol_src})

# Add include folder
target_include_directories(${PROJECT_NAME} PUBLIC include)

# this copies all resource files in the build directory
# we need this, because we want to work with paths relative to the executable
# file(COPY ${DATA} DESTINATION build)

###############################################################################
## dependencies ###############################################################
###############################################################################

# Add math library
target_link_libraries(${PROJECT_NAME} m)

# Add PkgConfig
find_package(PkgConfig REQUIRED)

# Add LAPACK
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
  message("Found Lapack")
  target_link_libraries(${PROJECT_NAME} ${LAPACK_LIBRARIES})
endif()

# Add OpenMP
option(USE_OpenMP "Use OpenMP solver library" OFF)
if (USE_OpenMP)
  find_package(OpenMP REQUIRED)
endif()
if (OPENMP_FOUND)
  message("Found OpenMP")
  target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_C)
  add_definitions("-DUSE_OPENMP")
endif()

# Add PETSC
option(USE_PETSC "Use PETSC solver library" OFF)
if (USE_PETSC)
  pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
endif()
if (PETSC_FOUND)
  message("Found PETSc")
  target_link_libraries(${PROJECT_NAME} PkgConfig::PETSC)
  add_definitions("-DUSE_PETSC")
endif()
 