#ifndef _U_GENERALIZED_ALPHA_H_
#define _U_GENERALIZED_ALPHA_H_

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


//  Global variables
//double Thickness_Plain_Stress;

/*!
  \brief The generalized-alpha algorithm here implemented is analogous
  to the one described in \cite Tran2019e

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void U_Generalized_alpha(
  Mesh FEM_Mesh, 
  Particle MPM_Mesh, 
  Time_Int_Params Parameters_Solver);

#endif