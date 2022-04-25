#ifndef _VON_MISES_CONSTITUTIVE_H_
#define _VON_MISES_CONSTITUTIVE_H_


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/*******************************************************/

/*!
  \fn int compute_1PK_Von_Mises(State_Parameters IO_State, Material MatProp)
  
  \brief Compute Von Mises yield

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model

*/
int compute_1PK_Von_Mises(State_Parameters IO_State, Material MatProp);
/*******************************************************/



#endif