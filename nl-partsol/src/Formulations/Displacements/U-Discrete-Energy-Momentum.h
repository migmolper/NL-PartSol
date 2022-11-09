
#ifndef _U_DISCRETE_ENERGY_MOMENTUM_H_
#define _U_DISCRETE_ENERGY_MOMENTUM_H_

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
  \fn void U_Discrete_Energy_Momentum(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief Discrete energy-momentum method. Implicit iterative solver proposed 
  by \cite Simo_and_Tarnow_1992, preserves linear and angular momentum.

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep  
*/
void U_Discrete_Energy_Momentum(Mesh, Particle, Time_Int_Params);

#endif