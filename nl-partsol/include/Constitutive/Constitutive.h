/**
 * @file Constitutive.h
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-18
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef _CONSTITUTIVE_H_
#define _CONSTITUTIVE_H_

// clang-format off
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Constitutive/Hyperelastic/Saint-Venant-Kirchhoff.h"
#include "Constitutive/Hyperelastic/Neo-Hookean.h"
#include "Constitutive/Plasticity/Matsuoka-Nakai.h"
#include "Constitutive/Plasticity/Lade-Duncan.h"
#include "Constitutive/Plasticity/Modified-Lade-Duncan.h"
#include "Constitutive/Plasticity/Drucker-Prager.h"
#include "Constitutive/Plasticity/Von-Mises.h"
#include "Constitutive/Fluid/Newtonian-Fluid.h"
#include "Constitutive/Plasticity/Elastoplastic-Tangent-Matrix.h"
#include "Constitutive/Fracture/Beps.h"
#include "Constitutive/Fracture/EigenErosion.h"
// clang-format on

/*******************************************************/



/*******************************************************/

/**
 * @brief 
 * 
 * @param p 
 * @param MPM_Mesh 
 * @param MatProp_p 
 * @return int 
 */
int Stress_integration__Constitutive__(
    int p, 
    Particle MPM_Mesh, 
    Material MatProp_p);
/*******************************************************/

/**
 * @brief 
 * 
 * @param p 
 * @param Stiffness_density 
 * @param dN_alpha_n1 
 * @param dN_beta_n1 
 * @param dN_alpha_n 
 * @param dN_beta_n 
 * @param alpha_4 
 * @param MPM_Mesh 
 * @param MatProp_p 
 * @return int 
 */
int stiffness_density__Constitutive__(
    int p, 
    double *Stiffness_density,
    const double *dN_alpha_n1,
    const double *dN_beta_n1,  
    const double *dN_alpha_n,
    const double *dN_beta_n,
    double alpha_4,
    Particle MPM_Mesh, 
    Material MatProp_p);

/*******************************************************/

#endif