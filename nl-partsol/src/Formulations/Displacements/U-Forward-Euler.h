#ifndef _U_FORWARD_EULER_FINITE_STRAINS_H_
#define _U_FORWARD_EULER_FINITE_STRAINS_H_

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Global libs
#include "Macros.h"
#include "Types.h"
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
#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif
#include <omp.h>

// 
#include "InOutFun.h"

//  Call global variables
//double Thickness_Plain_Stress;

/*!
  \fn void U_Forward_Euler(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief Displacement formulation of the MPM with a Forward Euler as 
  time integrator scheme. The algorithm was taken from \cite Zhang_book_2016

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void U_Forward_Euler(Mesh, Particle, Time_Int_Params);

#endif