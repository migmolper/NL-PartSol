
#ifndef _U_NEWMARK_BETA_H_
#define _U_NEWMARK_BETA_H_

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Global libs
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"

// Shape functions and auxilar tools
#include "Nodes/Nodes-Tools.h"
#include "Nodes/Shape-Functions.h"

// Courant
#include "Formulations/Courant.h"

// Material lib
#include "Constitutive/Constitutive.h"

// HPC libs
#ifdef USE_OPENMP
#include <omp.h>
#endif

// Linear-Solver libs
#ifdef USE_PETSC
#include "Linear-Solvers/ksp-PETSC.h"
#else
#include "Linear-Solvers/dgetrs-LAPACK.h"
#endif

// 
#include "InOutFun.h"

#ifdef USE_PETSC
#include "petscviewerhdf5.h" 
#endif

/*
  Call global variables
*/
double Thickness_Plain_Stress;
Event *Out_nodal_path_csv;
Event *Out_particles_path_csv;
int Number_Out_nodal_path_csv;
int Number_Out_particles_path_csv;


// Global variuables
unsigned InitialStep;
unsigned TimeStep;
unsigned NumTimeStep;

/*!
  \brief Finite strains Newmark-beta
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
*/
int U_Newmark_Beta(
    Mesh FEM_Mesh, 
    Particle MPM_Mesh,
    Time_Int_Params Parameters_Solver);

#endif