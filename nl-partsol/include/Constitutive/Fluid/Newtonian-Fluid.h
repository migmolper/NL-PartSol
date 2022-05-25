#ifndef _NEWTONIAN_FLUID_CONSTITUTIVE_H_
#define _NEWTONIAN_FLUID_CONSTITUTIVE_H_

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

  \fn int compute_Kirchhoff_Stress_Tensor_Newtonian_Fluid(State_Parameters IO_State,Material MatProp_p);

  \brief Compute the First Piola-Kirchhoff stress tensor of a Newtonian fluid 

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model

*/
int compute_Kirchhoff_Stress_Tensor_Newtonian_Fluid(
  State_Parameters IO_State,
  Material MatProp_p);
/*******************************************************/

/*!
  \brief Assemble the factorised material tensor 

  \param[out] Stiffness_Density
  \param[in] dN_alpha_n1
  \param[in] dN_beta_n1
  \param[in] dN_alpha_n
  \param[in] dN_beta_n
  \param[in] IO_State
  \param[in] MatProp
*/
int compute_stiffness_density_Newtonian_Fluid(
  double * Stiffness_Density,
  const double * dN_alpha_n1,
  const double * dN_beta_n1,
  const double * dN_alpha_n,
  const double * dN_beta_n,   
  State_Parameters IO_State,
  Material MatProp);
/*******************************************************/

#endif