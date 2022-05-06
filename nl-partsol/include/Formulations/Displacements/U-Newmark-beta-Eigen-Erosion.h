
#ifndef _U_NEWMARK_BETA_EIGEN_EROSION_H_
#define _U_NEWMARK_BETA_EIGEN_EROSION_H_

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

// Linear-Solver libs
#ifdef USE_PETSC
#include <petscksp.h>
//  #include "Linear-Solvers/"
#else
#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__ 
#include <Accelerate/Accelerate.h>
#endif
#include "Linear-Solvers/dgetrs-LAPACK.h"
#endif

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


// Global variuables
unsigned InitialStep;
unsigned TimeStep;
unsigned NumTimeStep;


#endif