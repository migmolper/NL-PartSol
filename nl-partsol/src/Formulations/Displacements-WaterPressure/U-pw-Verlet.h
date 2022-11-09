
#ifndef _UPW_VERLET_BETA_H_
#define _UPW_VERLET_BETA_H_

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


/*!
  \fn void upw_Verlet(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief u-pw formulation with one single set of material points. Explicit solver for finite strains based in the 
    Newmark Predictor-Corrector, see \cite Molinos_et_al_2021_CMAME.

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep  
*/
void upw_Verlet(Mesh, Particle, Time_Int_Params);

#endif