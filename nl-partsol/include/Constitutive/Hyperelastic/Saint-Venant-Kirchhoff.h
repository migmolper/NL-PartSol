#ifndef _SAINT_VENANT_KIRCHHOFF_CONSTITUTIVE_H_
#define _SAINT_VENANT_KIRCHHOFF_CONSTITUTIVE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Matlib.h"
#include "Particles.h"

/*!

*/
double energy_Saint_Venant_Kirchhoff(Tensor, Material);
/*******************************************************/

/*
  \fn Tensor compute_1PK_Stress_Tensor_Saint_Venant_Kirchhoff(Tensor P, Tensor F, Material MatProp_p);
*/

State_Parameters compute_1PK_Stress_Tensor_Saint_Venant_Kirchhoff(State_Parameters,Material);
/*******************************************************/

/*!
  
*/
Tensor grad_energy_Saint_Venant_Kirchhoff(Tensor, Tensor, Material);
/*******************************************************/
/*!
  
*/
Tensor compute_stiffness_density_Saint_Venant_Kirchhoff(Tensor, Tensor, Material);
/*******************************************************/

#endif