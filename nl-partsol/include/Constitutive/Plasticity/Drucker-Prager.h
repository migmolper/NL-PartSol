#ifndef _DRUCKER_PRAGER_CONSTITUTIVE_H_
#define _DRUCKER_PRAGER_CONSTITUTIVE_H_

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
  \fn int compute_1PK_Drucker_Prager(State_Parameters IO_State, Material MatProp);

  \brief Compute the behaviour of a Henky-Hyperelastic material \n
  using the Drucker-Prager yield criterium.

  \param[in] Input_SP : State parameters of the particle
  \param[in] MatProp : Material properties of the model
*/
int compute_1PK_Drucker_Prager(
    State_Parameters IO_State, 
    Material MatProp);

/**************************************************************/

#endif
