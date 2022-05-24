#ifndef _ELASTOPLASTIC_TANGENT_MATRIX_CONSTITUTIVE_H_
#define _ELASTOPLASTIC_TANGENT_MATRIX_CONSTITUTIVE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"

/*!
  \param[out] Stiffness_density Elastoplastic tangent matrix
  \param[in] dN_alpha_n1 Shape function gradient node A (t = n + 1)
  \param[in] dN_beta_n1 Shape function gradient node B (t = n + 1)
  \param[in] IO_State State variables
*/
int compute_1PK_elastoplastic_tangent_matrix(
  double *Stiffness_density,
  const double *dN_alpha_n1,
  const double *dN_beta_n1,  
  const State_Parameters IO_State);

/**************************************************************/

#endif
