
#ifndef _UPW_NEWMARK_BETA_H_
#define _UPW_NEWMARK_BETA_H_

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

// Linear-Solver lib
#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__ 
#include <Accelerate/Accelerate.h>
#endif

#include "Linear-Solvers/dgetrs-LAPACK.h"

#include <omp.h>

// 
#include "InOutFun.h"

/*
  Call global variables
*/
double Thickness_Plain_Stress;
Event *Out_nodal_path_csv;
Event *Out_particles_path_csv;
int Number_Out_nodal_path_csv;
int Number_Out_particles_path_csv;

/*!
  \fn void upw_Newmark_beta_Finite_Strains(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief u-pw formulation with one single set of material points. Implicit solver for finite strains based in the 
    Newmark-beta, see \cite Molinos_et_al_2021_CMAME.

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep  
*/
void upw_Newmark_beta_Finite_Strains(Mesh, Particle, Time_Int_Params);

#endif