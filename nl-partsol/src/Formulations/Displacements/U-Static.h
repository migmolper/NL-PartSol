
#ifndef _U_STATIC_H_
#define _U_STATIC_H_

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
#include "Particles/compute-Strains.h"

// Shape functions and auxilar tools
#include "Nodes/Nodes-Tools.h"
#include "Nodes/Shape-Functions.h"

// Courant
#include "Formulations/Courant.h"

// Material lib
#include "Constitutive/Constitutive.h"

// 
#include "InOutFun.h"

/**
 * @brief Run the static equations using a Finite strains Newmark-beta
 * 
 * @param FEM_Mesh Variable with the nodal information
 * @param MPM_Mesh Variable with the particle information
 * @param Parameters_Solver Time integration parameters
 * @return PetscErrorCode 
 */
PetscErrorCode U_Static(
  Mesh FEM_Mesh, 
  Particle MPM_Mesh, 
  Time_Int_Params Parameters_Solver);

#endif