#ifndef _EIGEN_SOFTENING_CONSTITUTIVE_H_
#define _EIGEN_SOFTENING_CONSTITUTIVE_H_

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
 * Here the notation is the same as in \cite Navas_2017
 * 
 * @param p Index of the particle
 * @param Damage_n Damage field (t = n)
 * @param Damage_n1 Damage field (t = n + 1)
 * @param Strain_p Strain field
 * @param StrainF_n Crack opening field (t = n)
 * @param StrainF_n1 Crack opening field (t = n + 1)
 * @param Mass Mass field
 * @param Stress_p Stress field
 * @param MatPro Define the material properties of the particle
 * @param Beps Define the particles close to each particle
 * @return int
 */
int Eigensoftening__Constitutive__(unsigned p, const double *Damage_n,
                                    double *Damage_n1, const double *Strain_p,
                                    const double *StrainF_n, double *StrainF_n1,
                                    const double *Mass, const double *Stress,
                                    Material MatProp_p, const ChainPtr Beps_p);

#endif