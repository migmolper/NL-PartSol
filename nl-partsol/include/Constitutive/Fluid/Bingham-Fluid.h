
#ifndef _BINGHAM_FLUID_CONSTITUTIVE_H_
#define _BINGHAM_FLUID_CONSTITUTIVE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"

/*!
 * \fn State_Parameters compute_1PK_Stress_Tensor_Bingham_Fluid(State_Parameters Intput_SP,Material MatProp_p);
 * */
State_Parameters compute_1PK_Stress_Tensor_Bingham_Fluid(State_Parameters,Material);
/*******************************************************/  

/*!
 * \fn Tensor compute_stiffness_density_Bingham_Fluid(Tensor GRAD_I, Tensor GRAD_J, Tensor F, Tensor dFdt, double J, double alpha4, Material MatProp_p)
 * */ 
Tensor compute_stiffness_density_Bingham_Fluid(Tensor,Tensor,Tensor,Tensor,double, double,Material);
/*******************************************************/


#endif