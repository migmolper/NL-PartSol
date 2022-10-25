
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

# Plain strain
if(USE_PLAINSTRAIN)
  add_definitions("-DUSE_PLAINSTRAIN")
endif()

# MPI
if(USE_MPI)
find_package(MPI)
if (MPI_FOUND)
  include_directories(${MPI_CXX_INCLUDE_DIRS})
  link_libraries(${MPI_CXX_LIBRARIES})
endif()
endif()


# Add LAPACK
if (USE_LAPACK)
  find_package(LAPACK REQUIRED)
  
  if(LAPACK_FOUND)
    target_link_libraries(${PROJECT_NAME} -llapacke -llapack -lblas)
    add_definitions("-DUSE_LAPACK")
  endif()

endif()

# Add OpenMP
if (USE_OPENMP)
  find_package(OpenMP REQUIRED)

  if (OPENMP_FOUND)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_C)
    add_definitions("-DUSE_OPENMP")
  endif()

endif()

# Add PETSC
if (USE_PETSC)
  pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
  
  if (PETSC_FOUND)
    target_link_libraries(${PROJECT_NAME} PkgConfig::PETSC)
    add_definitions("-DUSE_PETSC")
  endif()

endif()
 
