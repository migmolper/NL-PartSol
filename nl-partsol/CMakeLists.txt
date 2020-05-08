cmake_minimum_required(VERSION 3.13.4)
project("nl-partsol")
 
#Bring the headers
include_directories(include)
 
# Wildcard additions
file(GLOB_RECURSE SOURCES "src/*.c" )

add_executable(nl-partsol ${SOURCES})

# Add libraries
target_link_libraries(nl-partsol m)
