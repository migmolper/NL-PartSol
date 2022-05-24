#ifndef _BEPS_CONSTITUTIVE_H_
#define _BEPS_CONSTITUTIVE_H_

// clang-format off
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
// clang-format on

// HPC libs
#ifdef USE_OPENMP
#include <omp.h>
#endif

/*******************************************************/

/**
 * @brief Compute the list of particles near to each particle 
 * to compute normalized variables in EigenErosion/EigenSoftening
 * 𝜖-neighbourhood
 * 
 * @param MPM_Mesh 
 * @param FEM_Mesh 
 * @param Initialize_Beps 
 */
void compute_Beps__Constitutive__(Particle MPM_Mesh, Mesh FEM_Mesh,
                                  bool Initialize_Beps);

/**************************************************************/

#endif