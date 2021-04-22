/*! \file Particles.h
    \brief File with the prototype of the material point functions/utilities
*/

#ifndef _PARTICLES_H_
#define _PARTICLES_H_

/*******************************************************/

/*!

*/
Matrix body_loads__Particles__(Load *, int, int, int);
/*******************************************************/

/*!

*/
Matrix contact_loads__Particles__(Load *, int, int, int);
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
  \fn void update_DeformationGradient(Tensor F_n1,Tensor F_n,Matrix DeltaU,Matrix gradient_p)

  \brief Function to update the deformation gradient 
  as in \cite Love_and_Sulsky_2006_a. Takes the values from the previous step, 
  and update it

  \f[
  f_{n+1} = 1 + grad ( \Delta \phi )\\
  F_{n+1} = f_{n+1} F_{n}
  \f]

  \param F_n1 : Deformation gradient at t = n + 1
  \param F_n : Deformation gradient at t = n
  \param DeltaU : Nodal displacement increment (iterate)
  \param gradient_p : Shape function gradient for each node evaluated
  in the particle position
*/
void update_increment_Deformation_Gradient__Particles__(Tensor, Matrix, Matrix);
/*******************************************************/

/*!
  \fn void update_DeformationGradient(Tensor F_n1,Tensor F_n,Matrix DeltaU,Matrix gradient_p)

  \brief Function to update the deformation gradient 
  as in \cite Love_and_Sulsky_2006_a. Takes the values from the previous step, 
  and update it

  \f[
  f_{n+1} = 1 + grad ( \Delta \phi )\\
  F_{n+1} = f_{n+1} F_{n}
  \f]

  \param F_n1 : Deformation gradient at t = n + 1
  \param F_n : Deformation gradient at t = n
  \param DeltaU : Nodal displacement increment (iterate)
  \param gradient_p : Shape function gradient for each node evaluated
  in the particle position
*/
void   update_Deformation_Gradient_n1__Particles__(Tensor, Tensor, Tensor);
/*******************************************************/

/*!
  \fn double compute_Jacobian_Rate__Particles__(double J_p, Matrix Velocity, Matrix gradient_p)

*/
double compute_Jacobian_Rate__Particles__(double, Matrix, Matrix);
/*******************************************************/

/*!
  \fn Tensor right_Cauchy_Green__Particles__(Tensor F)
  \brief Function to cumpute the right Cauchy-Green tensor defined as
  \f[
  C = (F)^T F
  \f]
  
  \param F : Deformation gradient 

  \return The right Cauchy-Green tensor
*/
Tensor right_Cauchy_Green__Particles__(Tensor);
/*******************************************************/

/*!
  \fn Tensor logarithmic_strains__Particles__(Tensor C)
  \brief Function to cumpute the small strains countrepart of the 
  right Cauchy-Green tensor defined as
  \f[
  \varepsilon = 1/2*logC
  \f]
  
  \param C : right Cauchy-Green

  \return The logaritmic strain tensor
*/
Tensor logarithmic_strains__Particles__(Tensor);
/*******************************************************/
/*!
  \fn Tensor increment_Deformation_Gradient_exponential_strains__Particles__(Tensor D_E)

*/
Tensor increment_Deformation_Gradient_exponential_strains__Particles__(Tensor);
/*******************************************************/

/*!
  \fn Tensor strain_Green_Lagrange__Particles__(Tensor C)
  \brief Function to cumpute the Lagrangian Strain tensor defined as
  \f[
  E = \frac{1}{2} [ C - I ]
  \f]
  
  \param C : right Cauchy-Green tensor

  \return The Lagrangian Strain tensor 
*/
Tensor strain_Green_Lagrange__Particles__(Tensor);
/*******************************************************/

/*!

*/
Matrix compute_B_matrix__Particles__(Tensor, Tensor);
Matrix compute_BT_matrix__Particles__(Tensor, Tensor);
/*******************************************************/

/*!

*/
double update_density__Particles__(double, double, Tensor);
/*******************************************************/

/*!

*/
Tensor explicit_integration_stress__Particles__(Tensor, Tensor, Material);
/*******************************************************/
/*
\fn Tensor forward_integration_Stress__Particles__(Tensor S_p,Tensor F_n1_p,Material MatProp_p)
*/
Tensor forward_integration_Stress__Particles__(Tensor,Tensor,Material);
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

  \param PK2  : Previous value of the stress tensor
  \param C_n1 : Right Cauchy-Green tensor at t = n + 1
  \param C_n  : Right Cauchy-Green tensor at t = n
  \param Mat  : Material properties.
*/
Tensor average_strain_integration_Stress__Particles__(Tensor, Tensor, Tensor, Material);
/*******************************************************/

/*!
  \fn void average_itegration_Stress__Particles__(Tensor PK2,Tensor C_n1,Tensor C_n,Material Mat);

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

  \param PK2  : Previous value of the stress tensor
  \param C_n1 : Right Cauchy-Green tensor at t = n + 1
  \param C_n  : Right Cauchy-Green tensor at t = n
  \param Mat  : Material properties.
*/
Tensor average_itegration_Stress__Particles__(Tensor, Tensor, Tensor, Material);
/*******************************************************/


/*! \fn void compute_Piola_transformation__Particles__(Tensor S_p, Tensor sigma_k1, Tensor F_total, double J)

*/
void compute_Piola_transformation__Particles__(Tensor, Tensor, Tensor, double);


/*!

*/
double internal_energy__Particles__(Tensor, Tensor);
/*******************************************************/

/*!
  \fn double finite_strains_internal_energy__Particles__(Tensor F_p, Material MatProp_p)
*/
double finite_strains_internal_energy__Particles__(Tensor, Material, double);
/*******************************************************/
/*!

*/
void initial_position__Particles__(Matrix, Mesh, int);
/*******************************************************/

/*!

*/
void local_search__Particles__(Particle, Mesh);
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
void asign_to_nodes__Particles__(int, ChainPtr, Mesh);
/*******************************************************/

/*!

void update_plastic_deformation_gradient(Tensor D_E_plastic, Tensor F_plastic)

*/
void update_plastic_deformation_gradient__Particles__(Tensor, Tensor);
/*******************************************************/

#endif
