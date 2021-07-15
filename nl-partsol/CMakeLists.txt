cmake_minimum_required(VERSION 3.10.2)
project("nl-partsol")

#Bring the headers
include_directories(include)

# Wildcard additions
# file(GLOB_RECURSE SOURCES "src/*.c")
file(GLOB_RECURSE SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "src/*.c")

add_executable(${PROJECT_NAME} ${SOURCES})

# Add compilation options
set(default_flags "-Wall -fno-rtti")
set_target_properties("${PROJECT_NAME}" PROPERTIES
  debug "${default_flags} -O0 -g -Wall -Wpedantic -Wextra"
  release "${default_flags} -O2 "
  )

# Add math library
target_link_libraries(${PROJECT_NAME} m)

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


add_custom_target(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

add_custom_target(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )
