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

# Plain strain conditions
option(USE_PLAIN_STRAIN "Use Plain-Strain conditions" ON)

# Axial symmetry conditions
option(USE_AXIAL_SYMMETRY "Use Axial-Symmetry conditions" OFF)


# Add LAPACK
option(USE_LAPACK "Use LAPACK solver library" ON)

# Add MPI
option(USE_MPI "Use MPI library" OFF)

# Add OpenMP
option(USE_OPENMP "Use OpenMP solver library" ON)

# Add PETSC
option(USE_PETSC "Use PETSC solver library" ON)

# Compile settings
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wpedantic -Wextra -Wunused-variable")
set(CMAKE_C_FLAGS_RELEASE "-Ofast")

###############################################################################
## dependencies ###############################################################
###############################################################################

# Add math library
link_libraries(m)

# Add PkgConfig
find_package(PkgConfig REQUIRED)

# Plain strain
if(USE_PLAIN_STRAIN)
  add_definitions("-DUSE_PLAIN_STRAIN=1")
else(USE_PLAIN_STRAIN)
add_definitions("-DUSE_PLAIN_STRAIN=0")
endif(USE_PLAIN_STRAIN)

# Axial symmetry
if(USE_AXIAL_SYMMETRY)
  add_definitions("-DUSE_AXIAL_SYMMETRY=1")
else(USE_AXIAL_SYMMETRY)
add_definitions("-DUSE_AXIAL_SYMMETRY=0")
endif(USE_AXIAL_SYMMETRY)

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
    link_libraries(lapacke)
    link_libraries(lapack)
    link_libraries(blas)
    add_definitions("-DUSE_LAPACK")
  endif()
endif()

# Add OpenMP
if (USE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    add_definitions("-DUSE_OPENMP")
    if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    endif()
  endif()
endif()

# Add PETSc
if(USE_PETSC)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
if (PETSC_FOUND)
    link_libraries(PkgConfig::PETSC)
    add_definitions("-DUSE_PETSC")
    message ("PETSC linked")
else (PETSC_FOUND)
    message (FATAL_ERROR "${PROJECT_NAME}: Unable to find the PETSC folder.")
endif(PETSC_FOUND)
endif()


# Add MPI
find_package(MPI)
if (MPI_FOUND)
  add_definitions("-DUSE_MPI")
  include_directories(${MPI_C_INCLUDE_DIRS})
  link_libraries(${MPI_C_LIBRARIES})
endif()

###############################################################################
## Create executable ##########################################################
###############################################################################

file(GLOB_RECURSE NL_PartSol_src
    src/Constitutive/Fluid/Bingham-Fluid.c
    src/Constitutive/Fluid/Newtonian-Fluid.c
    src/Constitutive/Fracture/Beps.c
    src/Constitutive/Fracture/EigenErosion.c
    src/Constitutive/Fracture/EigenSoftening.c
    src/Constitutive/Hyperelastic/Hencky.c
    src/Constitutive/Hyperelastic/Mooney-Rivlin.c
    src/Constitutive/Hyperelastic/Neo-Hookean.c
    src/Constitutive/Hyperelastic/Saint-Venant-Kirchhoff.c
    src/Constitutive/Plasticity/Drucker-Prager.c
    src/Constitutive/Plasticity/Elastoplastic-Tangent-Matrix.c
    src/Constitutive/Plasticity/Lade-Duncan.c
    src/Constitutive/Plasticity/Matsuoka-Nakai.c
    src/Constitutive/Plasticity/Modified-Lade-Duncan.c
    src/Constitutive/Plasticity/Von-Mises.c
    src/Constitutive/Constitutive.c
    src/Formulations/Displacements/U-Analisys.c
    src/Formulations/Displacements/U-Discrete-Energy-Momentum.c
    src/Formulations/Displacements/U-Forward-Euler.c
    src/Formulations/Displacements/U-Generalized-Alpha.c
    src/Formulations/Displacements/U-Newmark-beta.c
    src/Formulations/Displacements/U-Static.c
    src/Formulations/Displacements/U-Verlet.c
    src/Formulations/Displacements-Pressure/U-p-Analisys.c
    src/Formulations/Displacements-Pressure/U-p-Newmark-beta.c
    src/Formulations/Displacements-WaterPressure/U-pw-Analisys.c
    src/Formulations/Displacements-WaterPressure/U-pw-Newmark-beta.c
    src/Formulations/Displacements-WaterPressure/U-pw-Verlet.c
    src/Formulations/Courant.c
    src/InOutFun/Analysis/Generate-One-Phase-Analysis.c
    src/InOutFun/Analysis/Generate-Soil-Water-Coupling.c
    src/InOutFun/Boundary-Conditions/NLPS-Read-u-Dirichlet-Boundary-Conditions.c
    src/InOutFun/Boundary-Conditions/NLPS-Read-u-Neumann-Boundary-Conditions.c
    src/InOutFun/Boundary-Conditions/NLPS-Read-upw-Dirichlet-Boundary-Conditions.c
    src/InOutFun/Boundary-Conditions/NLPS-Read-upw-Neumann-Boundary-Conditions.c
    src/InOutFun/Initial-Conditions/NLPS-Hidrostatic-Condition.c
    src/InOutFun/Initial-Conditions/NLPS-Initial-condition-nodes.c
    src/InOutFun/Initial-Conditions/Read_GramsInitials.c
    src/InOutFun/Material/Fluid/Compressible-Newtonian-Fluid.c
    src/InOutFun/Material/Fluid/Incompressible-Newtonian-Fluid.c
    src/InOutFun/Material/Linear/Linear-Elastic.c
    src/InOutFun/Material/Linear/Solid-Rigid.c
    src/InOutFun/Material/Hyperelastic/Hencky.c
    src/InOutFun/Material/Hyperelastic/Neo-Hookean.c
    src/InOutFun/Material/Hyperelastic/Saint-Venant-Kirchhoff.c
    src/InOutFun/Material/Plasticity/Drucker-Prager.c
    src/InOutFun/Material/Plasticity/Lade-Duncan.c
    src/InOutFun/Material/Plasticity/Matsuoka-Nakai.c
    src/InOutFun/Material/Plasticity/Von-Mises.c
    src/InOutFun/Material/NLPS-Read-Material-Mixtures.c
    src/InOutFun/Material/Read_GramsMaterials2.c
    src/InOutFun/Outputs/NLPS-Out-nodal-path-csv.c
    src/InOutFun/Outputs/NLPS-Out-particle-path-csv.c
    src/InOutFun/Outputs/Read_GramsOutputs.c
    src/InOutFun/Outputs/WriteCSV.c
    src/InOutFun/Outputs/WriteVtk.c
    src/InOutFun/File2Chain.c
    src/InOutFun/GnuPlotOutput.c
    src/InOutFun/NLPS-Read-Tools.c           
    src/InOutFun/OutputSimulation.c
    src/InOutFun/Parser.c
    src/InOutFun/print_ScreenMessage.c
    src/InOutFun/Read_Generate_Gravity_Field.c
    src/InOutFun/Read_GramsBox.c
    src/InOutFun/Read_GramsShapeFun.c
    src/InOutFun/Read_GramsTime.c
    src/InOutFun/Read-Delimited-File.c
    src/InOutFun/ReadCSV.c
    src/InOutFun/ReadCurve.c  
    src/Linear-Solvers/Conjugate-Gradient.c
    src/Linear-Solvers/dgetrs-LAPACK.c
    src/Linear-Solvers/Jacobi-Conjugate-Gradient.c
    src/Linear-Solvers/ksp-PETSC.c      
    src/Matlib/ChainOp.c                    
    src/Matlib/MathOp.c
    src/Matlib/MatrixOp.c
    src/Matlib/Solvers.c
    src/Matlib/TensorLib.c
    src/Nodes/aLME.c
    src/Nodes/GIMP.c
    src/Nodes/H8.c
    src/Nodes/LME.c
    src/Nodes/Nodes-Tools.c
    src/Nodes/Q4.c
    src/Nodes/Read-GID-Mesh.c
    src/Nodes/Shape-Functions.c
    src/Nodes/T3.c
    src/Nodes/T4.c
    src/Particles/compute-Strains.c
    src/Particles/Particles-Tools.c
    src/Globals.h
    src/InOutFun.h
    src/Macros.h
    src/Matlib.h 
    src/nl-partsol.h
    src/Particles.h
    src/Types.h)

file(GLOB_RECURSE DATA ../build/*)

# Create executables
add_executable(${PROJECT_NAME} driver-nl-partsol.c ${NL_PartSol_src})
target_include_directories(${PROJECT_NAME}  PUBLIC src)


