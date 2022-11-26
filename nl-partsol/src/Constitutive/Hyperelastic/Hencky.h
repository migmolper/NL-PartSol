/**
 * @file Hencky.h
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-18
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef _HENCKY_CONSTITUTIVE_H_
#define _HENCKY_CONSTITUTIVE_H_

// clang-format off
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Matlib.h"
#include "Particles/compute-Strains.h"
// clang-format on

/**
 * @brief 
 * 
 * @param IO_State 
 * @param MatProp 
 * @return int 
 */
int compute_Kirchhoff_Stress_Hencky__Constitutive__(State_Parameters IO_State,
                                    Material MatProp);

/**
 * @brief 
 * 
 * @param Stiffness_density 
 * @param dN_alpha_n1 
 * @param dN_beta_n1 
 * @param IO_State 
 * @param MatProp 
 * @return int 
 */
int compute_stiffness_density_Hencky__Constitutive__(double *Stiffness_density,
                                             const double *dN_alpha_n1,
                                             const double *dN_beta_n1,
                                             const State_Parameters IO_State,
                                             Material MatProp);

#endif
