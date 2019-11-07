option ( BUILD_FORTRAN_EXAMPLES "Build fortran examples, it is needed a fortran compiler" OFF )

if (BUILD_FORTRAN_EXAMPLES)
  enable_language ( Fortran )
endif (BUILD_FORTRAN_EXAMPLES)

add_executable ( testc_fd testpost_fd.c )
target_link_libraries ( testc_fd gidpost )

add_executable ( testc testpost.c )
target_link_libraries ( testc gidpost )

add_executable ( BugMeshGroup BugMeshGroup.c )
target_link_libraries ( BugMeshGroup gidpost )

add_executable ( BugMeshGroup02 BugMeshGroup02.c )
target_link_libraries ( BugMeshGroup02 gidpost )

add_executable ( BugMeshGroup03 BugMeshGroup03.c )
target_link_libraries ( BugMeshGroup03 gidpost )

if (CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -D_BSD_SOURCE")
endif (CMAKE_COMPILER_IS_GNUCXX)

if (CMAKE_Fortran_COMPILER)

  add_executable ( testf testpostfor.f )
  target_link_libraries ( testf gidpost )

  if (UNIX)
    target_link_libraries ( testf gfortran )
  endif (UNIX)

endif (CMAKE_Fortran_COMPILER)

add_executable( TestMultiplePartition TestMultiplePartition.cxx )
target_link_libraries( TestMultiplePartition gidpost )
