/**
 * @file U-Newmark-beta.h
 * @author Miguel Molinos (@migmolper)
 * @brief Incremental Newmark-beta (finite strains)
 * @version 0.1
 * @date 2022-05-18
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef _U_NEWMARK_BETA_H_
#define _U_NEWMARK_BETA_H_

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

#ifdef USE_PETSC
#include "petscviewerhdf5.h"
#endif


/*!
  \brief Finite strains Newmark-beta

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
*/
int U_Newmark_Beta(Mesh FEM_Mesh, Particle MPM_Mesh,
                   Time_Int_Params Parameters_Solver);

#endif