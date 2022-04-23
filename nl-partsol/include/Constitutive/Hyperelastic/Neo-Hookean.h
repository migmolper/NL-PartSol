#ifndef _NEOHOOK_CONSTITUTIVE_H_
#define _NEOHOOK_CONSTITUTIVE_H_

#include "Types.h"
#include "Globals.h"
#include "Matlib.h"

/*!

*/
double energy_Neo_Hookean_Wriggers(Tensor, double, Material);
/*******************************************************/

/*!
  \fn int compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(State_Parameters IO_State,Material MatProp);

  \brief Compute the First Piola-Kirchhoff stress tensor of a Neo-Hookean material 

  \param[in] IO_State : State parameters of the particle
  \param[in] MatProp : Material properties of the model
*/
int compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(
  State_Parameters IO_State,
  Material MatProp);
/*******************************************************/


Tensor compute_2PK_Stress_Tensor_Neo_Hookean_Wriggers(Tensor, Tensor, double, Material);
/*******************************************************/


/*!
  \fn Tensor compute_material_stiffness_density_Neo_Hookean_Wriggers(Tensor v, Tensor w, Tensor C, double J, Material MatProp)

  \brief Assemble the factorised material tensor 

  \param[in] v
  \param[in] w
  \param[in] C
  \param[in] J
  \param[in] MatProp
*/
Tensor compute_material_stiffness_density_Neo_Hookean_Wriggers(Tensor, Tensor, Tensor, double, Material);
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
int compute_stiffness_density_Neo_Hookean_Wriggers(
  double * Stiffness_Density,
  const double * dN_alpha_n1,
  const double * dN_beta_n1,
  const double * dN_alpha_n,
  const double * dN_beta_n,  
  State_Parameters IO_State,
  Material MatProp);
/*******************************************************/ 

#endif
