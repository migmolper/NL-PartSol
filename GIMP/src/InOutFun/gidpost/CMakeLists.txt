cmake_minimum_required (VERSION 2.6)

project (gidpost)

# location of Zlib
message( "user = $ENV{USERNAME}")

# for Windows
if (WIN32 AND $ENV{USERNAME} STREQUAL "miguel")
  if ( EXISTS C:/GID64-git/libZ/zlib-1.2.11)
    set( ZLIB_ROOT "C:/GID64-git/libZ/zlib-1.2.11")
  elseif ( EXISTS C:/GID64-git/libZ/windows)
    set( ZLIB_ROOT "C:/GID64-git/libZ/windows")
  else( EXISTS C:/GID64-git/libZ/windows)
    message( "ZLIB: trying to guess Zlib location")
  endif( EXISTS C:/GID64-git/libZ/zlib-1.2.11)
endif (WIN32 AND $ENV{USERNAME} STREQUAL "miguel")

add_definitions ( -DHAVE_GIDPOST_CONFIG_H )

# The version number.
set (GP_VERSION_MAJOR 2)
set (GP_VERSION_MINOR 5)

set( GP_VERSION ${GP_VERSION_MAJOR}.${GP_VERSION_MINOR} )

option(BUILD_SHARED_LIBS "Build gidpost as a shared libraries." OFF)

if (BUILD_SHARED_LIBS)
  add_definitions (-DGIDPOST_SHARED)
endif (BUILD_SHARED_LIBS)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/gidpost_config.h.in"
  "${PROJECT_BINARY_DIR}/gidpost_config.h"
  )

install( FILES "${PROJECT_BINARY_DIR}/gidpost_config.h" DESTINATION "include" )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

include_directories ("${PROJECT_SOURCE_DIR}")
include_directories ("${PROJECT_SOURCE_DIR}/source")

add_subdirectory (source)

# Make sure the linker can find the gidpost library once it is
# built.
link_directories (${PROJECT_BINARY_DIR}/source) 

add_subdirectory (examples)
