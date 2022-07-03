#ifndef _NEOHOOK_CONSTITUTIVE_H_
#define _NEOHOOK_CONSTITUTIVE_H_

// clang-format off
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Matlib.h"
#include "Particles.h"
// clang-format on

/**
 * @brief Compute the strain-energy of a Neo-Hookean material
 * 
 * @param b Left Cauchy Green tensor
 * @param J Jacobian
 * @param MatProp Material properties of the model
 * @return double 
 */
double compute_strain_energy_Neo_Hookean__Constitutive__(const double *b, double J,
                                                         Material MatProp);
/*******************************************************/

/**
 * @brief Compute the Kirchhoff stress tensor of a Neo-Hookean material
 * 
 * @param IO_State State parameters of the particle
 * @param MatProp Material properties of the model
 * @return int 
 */
int compute_Kirchhoff_Stress_Neo_Hookean__Constitutive__(
  State_Parameters IO_State,
  Material MatProp);
/*******************************************************/

/**
 * @brief 
 * 
 * @param Stiffness_Density 
 * @param dN_alpha_n1 
 * @param dN_beta_n1 
 * @param dN_alpha_n 
 * @param dN_beta_n 
 * @param IO_State 
 * @param MatProp 
 * @return int 
 */
int compute_stiffness_density_Neo_Hookean(
  double * Stiffness_Density,
  const double * dN_alpha_n1,
  const double * dN_beta_n1,
  const double * dN_alpha_n,
  const double * dN_beta_n,  
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

#endif
