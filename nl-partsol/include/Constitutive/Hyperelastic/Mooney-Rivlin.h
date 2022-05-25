#ifndef _MOONEY_RIVLIN_CONSTITUTIVE_H_
#define _MOONEY_RIVLIN_CONSTITUTIVE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Matlib.h"
#include "Particles.h"

/*******************************************************/

double energy_Mooney_Rivlin(
    Tensor C, 
    double J, 
    Material MatProp_p);
/*******************************************************/

State_Parameters compute_1PK_Stress_Tensor_Mooney_Rivlin(
    State_Parameters Intput_SP,
    Material MatProp_p);
/*******************************************************/

Tensor compute_stiffness_density_Mooney_Rivlin(Tensor GRAD_I, Tensor GRAD_J,
                                               Tensor F, double J,
                                               Material MatProp);
/*******************************************************/

#endif