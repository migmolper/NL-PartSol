#ifndef _MODIFIED_LADE_DUNCAN_CONSTITUTIVE_H_
#define _MODIFIED_LADE_DUNCAN_CONSTITUTIVE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"


/*******************************************************/

/*!
  \fn int compute_1PK_Modified_Lade_Duncan(State_Parameters IO_State, Material MatProp)

  \brief Compute a family of smooth approximations of the Mohr-Coulomb model using a monolithic algorithm

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model
*/
int compute_1PK_Modified_Lade_Duncan(State_Parameters IO_State, Material MatProp);

/*******************************************************/

#endif