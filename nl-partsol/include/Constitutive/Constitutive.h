#ifndef _CONSTITUTIVE_H_
#define _CONSTITUTIVE_H_

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

/*******************************************************/

/*!
    \param[in] p index of the particle  
    \param[in] MPM_Mesh  
    \param[in] MatProp_p Material properties
*/
int Stress_integration__Constitutive__(
    int p, 
    Particle MPM_Mesh, 
    Material MatProp_p);
/*******************************************************/

/*!
    \param[in] p 
    \param[in] Stiffness_density
    \param[in] dN_alpha_n1
    \param[in] dN_beta_n1  
    \param[in] dN_alpha_n
    \param[in] dN_beta_n
    \param[in] MPM_Mesh 
    \param[in] MatProp_p
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