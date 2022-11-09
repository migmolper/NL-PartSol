#ifndef _EIGEN_EROSION_CONSTITUTIVE_H_
#define _EIGEN_EROSION_CONSTITUTIVE_H_

// clang-format off
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Matlib.h"
// clang-format on

// HPC libs
#ifdef USE_OPENMP
#include <omp.h>
#endif

/**
 * @brief Function to compute is a material point is or not eroded.
 * Here the notation is the same as in \cite Pandolfi_2012
 *
 * @param p Index of the particle
 * @param Damage_n Damage field (t = n)
 * @param Damage_n1 Damage field (t = n + 1)
 * @param W Strain-Energy field
 * @param J_n1 Jacobian field
 * @param Vol_0 Initial volume field
 * @param Stress_p Stress tensor in particle p
 * @param MatPro Define the material properties of the particle
 * @param Beps Define the particles close to each particle
 * @param DeltaX Mesh size
 * @return int
 */
int Eigenerosion__Constitutive__(unsigned p, const double *Damage_n,
                                 double *Damage_n1, const double *W,
                                 const double *J_n1, const double *Vol_0,
                                 const double *Stress_p, Material MatProp_p,
                                 const ChainPtr Beps_p, double DeltaX);

#endif