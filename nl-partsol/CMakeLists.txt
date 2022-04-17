
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

# Add PkgConfig
find_package(PkgConfig REQUIRED)

# Add LAPACK
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
  message("Found Lapack")
  target_link_libraries(${PROJECT_NAME} ${LAPACK_LIBRARIES})
endif()

# Add OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  message("Found OpenMP")
  target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_C)
endif()

# Add PETSC
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
if (PETSC_FOUND)
  message("Found PETSc")
  target_link_libraries(${PROJECT_NAME} PkgConfig::PETSC)
endif()


add_custom_target(debug
  COMMAND ${CMAKE_COMMAND} -DDEBUG_MODE=0 -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

add_custom_target(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )


  