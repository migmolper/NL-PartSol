/**
 * @file compute-Strains.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-11-26
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef _COMPUTE_STRAINS_H_
#define _COMPUTE_STRAINS_H_

// clang-format off
#include "Nodes/H8.h"
#include "Nodes/Q4.h"
#include "Nodes/T3.h"
#include "Nodes/T4.h"
// clang-format on

// clang-format off
#include "Nodes/H8.h"
#include "Nodes/Q4.h"
#include "Nodes/T3.h"
#include "Nodes/T4.h"
// clang-format on


typedef struct {

  unsigned Nnodes_p;

  const double * gradient_p;
  
  #if USE_AXIAL_SYMMETRY
  const double * shapefun_p;

  double R_p;
  #endif

} compute_strains_ctx;

/*******************************************************/

 /**
 * @brief Function to compute the increment of the deformation gradient
 * \f[
 * f_{n+1} = 1 + grad ( \Delta \phi )
 * \f]
 * 
 * @param f_n1 Increment of the deformation gradient
 * @param DeltaU Nodal displacement increment
 * @param ctx
 */
void update_increment_Deformation_Gradient__Particles__(
  double * f_n1, 
  const double * DeltaU, void* ctx);
/*******************************************************/

/**
 * @brief Function to update the deformation gradient rate  
 * as in \cite Love_and_Sulsky_2006_a. Takes the values from the previous step, 
 * and update it 
 * 
 * \f[ 
 * f_{n+1} = 1 + grad ( \Delta \phi )\\
 *  F_{n+1} = f_{n+1} F_{n}
 * \f]
 * 
 * @param dt_f_n1 Deformation gradient at t = n + 1 
 * @param DeltaV Nodal velocity increment
 * @param ctx
 */
void update_rate_increment_Deformation_Gradient__Particles__(
  double * dt_f_n1, 
  const double * DeltaV,
  void* ctx);
/*******************************************************/


/**
 * @brief Function to update the deformation gradient  
 * as in \cite Love_and_Sulsky_2006_a. Takes the values from the previous step, 
 * and update it
 * 
 * \f[
 *  F_{n+1} = f_{n+1} F_{n}
 * \f]
 * 
 * @param F_n1 Deformation gradient at t = n + 1
 * @param F_n Deformation gradient at t = n
 * @param f_n1 Increment of the deformation gadient
 */
void update_Deformation_Gradient_n1__Particles__(
  double * F_n1, 
  const double * F_n,
  const double * f_n1);
/*******************************************************/


/*!
  \brief Function to compute the F-bar using the approach proposed
  by Ortiz

*/
int get_locking_free_Deformation_Gradient_n1__Particles__(
  unsigned p,
  double DJ_patch,
  Particle MPM_Mesh);
/*******************************************************/

/**
 * @brief Function to update the deformation gradient rate  
 * as in \cite Love_and_Sulsky_2006_a. Takes the values from the previous step, 
 * and update it
 * 
 * \f[
 * f_{n+1} = 1 + grad ( \Delta \phi )\\
 * F_{n+1} = f_{n+1} F_{n}
 * \f]
 * 
 * @param dt_F_n1 Deformation gradient rate at t = n + 1
 * @param dt_f_n1 Rate of the increment of the deformation gradient 
 * @param F_n Deformation gradient at t = n
 * @param f_n1 Increment of the deformation gradient
 * @param dt_F_n Rate of the deformation gradient t = n
 */
void update_rate_Deformation_Gradient_n1__Particles__(
  double * dt_F_n1,
  const double * dt_f_n1,
  const double * F_n, 
  const double * f_n1,
  const double * dt_F_n);
/*******************************************************/



/**
 * @brief Function to compute the rate of the Jacobian, 
 * see Marsden and Hughes
 * \f[
 *  d_J_dt = J (FmT : dFdt)
 * \f]
 * 
 * @param d_J_dt Jacobian rate
 * @param J Jacobian
 * @param F Deformation gradient
 * @param d_F_dt Deformation gradient rate
 * @return int 
 */
int compute_Jacobian_Rate__Particles__(
  double * d_J_dt,
  double J, 
  const double * F,
  const double * d_F_dt);
/*******************************************************/

/**
 * @brief Compute the spatial velocity gradient
 * 
 * @param L Spatial velocity gradient
 * @param dFdt Rate of the deformation tensor
 * @param F Deformation gradient
 * @return int 
 */
int spatial_velocity_gradient__Particles__(double * L,const double * dFdt,const double * F);
/*******************************************************/



/*!
  \fn Tensor right_Cauchy_Green__Particles__(Tensor F)
  \brief Function to cumpute the right Cauchy-Green tensor defined as
  \f[
  C = (F)^T F
  \f]
  
  \param F Deformation gradient 

  \return The right Cauchy-Green tensor
*/
Tensor right_Cauchy_Green__Particles__(Tensor);
/*******************************************************/

/**
 * @brief Compute the left Cauchy-Green tensor
 * 
 * @param b Left Cauchy-Green tensor
 * @param F Deformation gradient
 */
void left_Cauchy_Green__Particles__(double * b, const double * F);
/*******************************************************/

/**
 * @brief Compute the Eulerian-Almansi strain tensor
 * 
 * @param e Eulerian-Almansi strain tensor
 * @param F Deformation gradient
 */
void eulerian_almansi__Particles__(double * e, const double * F);
/*******************************************************/

/*!
  \fn Tensor strain_Green_Lagrange__Particles__(Tensor C)
  \brief Function to cumpute the Lagrangian Strain tensor defined as
  \f[
  E = \frac{1}{2} [ C - I ]
  \f]
  
  \param C right Cauchy-Green tensor

  \return The Lagrangian Strain tensor 
*/
Tensor strain_Green_Lagrange__Particles__(Tensor);
/*******************************************************/

#endif // _COMPUTE_STRAINS_H_