#ifndef _U_VERLET_H_
#define _U_VERLET_H_

// clang-format off
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"
#include "Nodes/Nodes-Tools.h"
#include "Nodes/Shape-Functions.h"
#include "Formulations/Courant.h"
#include "Constitutive/Constitutive.h"
#include "InOutFun.h"
// clang-format on

// Linear-Solver libs
#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif
#include <omp.h>

// 
#include "InOutFun.h"

//  Call global variables
//double DeltaTimeStep;
//double Thickness_Plain_Stress;
//Event *Out_nodal_path_csv;
//Event *Out_particles_path_csv;
//int Number_Out_nodal_path_csv;
//int Number_Out_particles_path_csv;

/*!
  \brief Finite strains explicit predictor-corrector gamma = 0.5
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
int U_Verlet(
  Mesh FEM_Mesh, 
  Particle MPM_Mesh, 
  Time_Int_Params Parameters_Solver);

#endif