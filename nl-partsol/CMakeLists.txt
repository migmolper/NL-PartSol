
###############################################################################
## file globbing ##############################################################
###############################################################################

# these instructions search the directory tree when cmake is
# invoked and put all files that match the pattern in the variables 
# `sources` and `data`
file(GLOB_RECURSE SOURCES src/*.c include/*.h)
file(GLOB_RECURSE DATA ../build/*)
# you can use set(SOURCES src/main.c) etc if you don't want to
# use globing to find files automatically


###############################################################################
## target definitions #########################################################
###############################################################################

# add the data to the target, so it becomes visible in some IDE
add_executable(${PROJECT_NAME} ${SOURCES})

# Add include folder
target_include_directories(${PROJECT_NAME} PUBLIC include)

# Add some compiler flags
set(default_flags "-Wall -fno-rtti")
set_target_properties("${PROJECT_NAME}" PROPERTIES
  debug "${default_flags} -O0 -g -Wall -Wpedantic -Wextra"
  release "${default_flags} -O2 "
  )

# this copies all resource files in the build directory
# we need this, because we want to work with paths relative to the executable
# file(COPY ${DATA} DESTINATION build)

###############################################################################
## dependencies ###############################################################
###############################################################################

# Add math library
target_link_libraries(${PROJECT_NAME} m)

# Use PkgConfig to find packages
find_package(PkgConfig REQUIRED)

# Add LAPACK
pkg_check_modules(LAPACK lapack IMPORTED_TARGET REQUIRED)
if(LAPACK_FOUND)
  target_link_libraries(${PROJECT_NAME} ${LAPACK_LIBRARIES})
endif()

# Add PETSc
pkg_check_modules(PETSC PETSc IMPORTED_TARGET REQUIRED)
if(PETSC_FOUND)
  
  list(APPEND COMPILE_OPTIONS ${PETSC_CFLAGS})

  include_directories(${PETSC_INCLUDE_DIRS})
  
  set(LINK_FLAGS "${LINK_FLAGS} ${PETSC_LDFLAGS}")   
  
  list(APPEND LIBRARIES ${PETSC_LINK_LIBRARIES})   
  
  set(CMAKE_REQUIRED_FLAGS ${PETSC_CFLAGS})   
  
  set(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDE_DIRS}")

  target_link_libraries(${PROJECT_NAME} ${PETSC_LIBRARIES})
endif()

# Add OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  message(STATUS "Found OpenMP ${OPENMP_VERSION}")
  target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_C)
endif()




