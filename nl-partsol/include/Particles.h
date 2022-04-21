/*! \file Particles.h
    \brief File with the prototype of the material point functions/utilities
*/

#ifndef _PARTICLES_H_
#define _PARTICLES_H_

/*******************************************************/

/*!

*/
Tensor rate_inifinitesimal_Strain__Particles__(Matrix, Matrix);
/*******************************************************/

/*!

*/
Tensor infinitesimal_Strain__Particles__(Tensor, Tensor, double);
/*******************************************************/

/*!
  \brief Function to compute the increment of the deformation gradient

  \f[
  f_{n+1} = 1 + grad ( \Delta \phi )
  \f]

  \param[out] f_n1 : Increment of the deformation gradient
  \param[in] DeltaU : Nodal displacement increment
  \param[in] gradient_p : Shape function gradient
  \param[in] Nnodes_p : Number of nodes
*/
void update_increment_Deformation_Gradient__Particles__(
  double * f_n1, 
  const double * DeltaU, 
  const double * gradient_p, 
  unsigned Nnodes_p);
/*******************************************************/

/*!
  \brief Function to update the deformation gradient 
  as in \cite Love_and_Sulsky_2006_a. Takes the values from the previous step, 
  and update it

  \f[
  F_{n+1} = f_{n+1} F_{n}
  \f]

  \param[out] F_n1 Deformation gradient at t = n + 1
  \param[in] F_n Deformation gradient at t = n
  \param[in] f_n1 Increment of the deformation gadient
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

/*!
  \brief Function to update the deformation gradient rate 
  as in \cite Love_and_Sulsky_2006_a. Takes the values from the previous step, 
  and update it

  \f[
  f_{n+1} = 1 + grad ( \Delta \phi )\\
  F_{n+1} = f_{n+1} F_{n}
  \f]

  \param[out] dt_f_n1 Deformation gradient at t = n + 1
  \param[in] DeltaV Nodal velocity increment
  \param[in] gradient_p Shape function gradient 
  \param[in] Nnodes_p Number of nodes
*/
void update_rate_increment_Deformation_Gradient__Particles__(
  double * dt_f_n1, 
  const double * DeltaV, 
  const double * gradient_p,
  unsigned Nnodes_p);
/*******************************************************/

/*!
  \brief Function to update the deformation gradient rate 
  as in \cite Love_and_Sulsky_2006_a. Takes the values from the previous step, 
  and update it

  \f[
  f_{n+1} = 1 + grad ( \Delta \phi )\\
  F_{n+1} = f_{n+1} F_{n}
  \f]

  \param[out] dt_F_n1 Deformation gradient rate at t = n + 1
  \param[in] dt_f_n1 Rate of the increment of the deformation gradient 
  \param[in] F_n Deformation gradient at t = n
  \param[in] f_n1 Increment of the deformation gradient
  \param[in] dt_F_n Rate of the deformation gradient t = n
*/
void update_rate_Deformation_Gradient_n1__Particles__(
  double * dt_F_n1,
  const double * dt_f_n1,
  const double * F_n, 
  const double * f_n1,
  const double * dt_F_n);
/*******************************************************/

/*!
  \brief Function to compute the rate of the Jacobian,
  see Marsden and Hughes

  \f[
  d_J_dt = J (FmT : dFdt)
  \f]

  \param[out] d_J_dt Jacobian rate
  \param[in] J Jacobial 
  \param[in] F Deformation gradient
  \param[in] d_F_dt Deformation gradient rate
*/
int compute_Jacobian_Rate__Particles__(
  double * d_J_dt,
  double J, 
  const double * F,
  const double * d_F_dt);
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

/*!
  \fn void left_Cauchy_Green__Particles__(double * b, const double * F);

  \brief Compute the left Cauchy-Green tensor

  \param[out] b Left Cauchy-Green tensor
  \param[in] F Deformation gradient
*/
void left_Cauchy_Green__Particles__(double * b, const double * F);
/*******************************************************/

/*!
  \fn Tensor rate_of_deformation__Particles__(Tensor dFdt, Tensor Fm1)  

  \param dFdt Temporal derivative of the deformation gradient
  \param Fm1 Inverse of the deformation gradient  
*/
Tensor rate_of_deformation__Particles__(Tensor, Tensor);
/*******************************************************/

/*!
  \fn int spatial_velocity_gradient__Particles__(double * L,const double * dFdt,const double * F);

  \brief Compute the spatial velocity gradient

  \param[out] L Spatial velocity gradient
  \param[in] dFdt Rate of the deformation tensor
  \param[in] F Deformation gradient
*/
int spatial_velocity_gradient__Particles__(double * L,const double * dFdt,const double * F);
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

/*!

*/
double update_density__Particles__(double, double, Tensor);
/*******************************************************/

/*!

*/
Tensor explicit_integration_stress__Particles__(int, Particle, Material);

/*******************************************************/

int Stress_integration__Particles__(int p, Particle MPM_Mesh, Mesh FEM_Mesh,
                                     Material MatProp_p);
/*******************************************************/


Tensor tangent_matrix__Particles__(Tensor GRADIENT_pA,
 Tensor GRADIENT_pB, 
 Tensor F_n1_p,
 Tensor dFdt_n1_p,
 double J_p,
 double alpha_4,
 Material MatProp_p);

/*******************************************************/
/*!
  \fn Tensor configurational_midpoint_integration_Stress__Particles__(Tensor T_n1,Tensor T_n,double alpha)

  \brief Function to cumpute a midpoint tensor as in
  \cite Simo_and_Tarnow_1992.

  \f[
  T_{n+1} = \alpha T_{n+1} + (1 - \alpha) T_{n}
  \f]

  \param T_n1 : Tensor at t = n + 1
  \param T_n : Tensor at t = n
  \param alpha : Define the evaluation point t = n + \f$\alpha\f$
*/
Tensor configurational_midpoint_integration_Stress__Particles__(Tensor,Tensor, Tensor, Material);
/*******************************************************/

/*!
    \fn void average_strain_integration_Stress__Particles__(Tensor PK2,Tensor C_n1,Tensor C_n,Material Mat);

  \brief Compute the second Piola-Kirchhoff stress tensor
  with an average value of the right Cauchy-Green tensors in 
  t = n and t = n + 1
  \f[
  C^{n + 1/2} = \frac{C^{n+1} + C^{n}}{2}
  \f]
  And evaluating it the gradient of the internal energy function
  \f[
  S =  2 \cdot \Delta \hat{e}(C^{n + 1/2 )}) 
  \f]

  \param[in] PK2  Previous value of the stress tensor
  \param[in] C_n1 Right Cauchy-Green tensor at t = n + 1
  \param[in] C_n  Right Cauchy-Green tensor at t = n
  \param[in] Mat  Material properties.
*/
//Tensor average_strain_integration_Stress__Particles__(Tensor, Tensor, Tensor, Material);
/*******************************************************/

/*!

*/
double internal_energy__Particles__(Tensor, Tensor);
/*******************************************************/

/*!
  \brief update the internal energy of the particle using
  \f[
  W_n1 = W_n + dt * (P : dfdt)
  \f]

  \param[out] W Internal energy
  \param[in] P First Piola-Kichhoff stress tensor
  \param[in] dFdt Deformation gradient rate
  \param[in] dt Time step for the time integration
*/
void finite_strains_internal_energy__Particles__(
  double * W,
  const double * P,
  const double * dFdt,
  double dt);
/*******************************************************/
/*!

*/
void initial_position__Particles__(Matrix, Mesh, int);
/*******************************************************/

/*!

*/
Element nodal_set__Particles__(int, ChainPtr, int);
/*******************************************************/

/*!

*/
int search_particle_in_surrounding_elements__Particles__(int, Matrix, ChainPtr, Mesh);
/*******************************************************/

/*!

*/
void asign_to_nodes__Particles__(int, int, int, ChainPtr, Mesh);
/*******************************************************/


#endif
