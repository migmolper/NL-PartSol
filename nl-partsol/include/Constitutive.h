/*! \file Constitutive.h
  \brief File with the prototype of the contitutive models
*/


#ifndef _CONSTITUTIVE_H_
#define _CONSTITUTIVE_H_

/*******************************************************/


/*!
  \fn Tensor SolidRigid(Tensor Strain)

  \brief This function is devoted to make a material point behaves as a
  solid rigid.

  \param Strain : Strain field of each particle 
*/
Tensor SolidRigid(Tensor Strain);
/*******************************************************/


/*!
  \fn State_Parameters LinearElastic(State_Parameters Intput_SP, Material MatProp_p)
  
  \brief This function is devoted to make a material point behaves as a
  linear elastic material.
  
  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model
*/
State_Parameters compute_kirchhoff_isotropic_linear_elasticity(State_Parameters, Material);
/*******************************************************/

/*!
  \fn void compute_particle_Damage(int p, Particle Particles, Mesh Nodes)

  \brief Compute if a particle is damaged and update it damage parameter chi

  \param p : Particle
  \param Particles : Information of the particle mesh
  \param Nodes : Informaction with the set of nodes
*/
void compute_particle_Damage(int, Particle, Mesh);
/*******************************************************/

/*!

*/
double energy_Saint_Venant_Kirchhoff(Tensor, Material);
/*******************************************************/

/*
  \fn Tensor compute_1PK_Stress_Tensor_Saint_Venant_Kirchhoff(Tensor P, Tensor F, Material MatProp_p);
*/

State_Parameters compute_1PK_Stress_Tensor_Saint_Venant_Kirchhoff(State_Parameters,Material);
/*******************************************************/

/*!
  
*/
Tensor grad_energy_Saint_Venant_Kirchhoff(Tensor, Tensor, Material);
/*******************************************************/
/*!
  
*/
Tensor compute_stiffness_density_Saint_Venant_Kirchhoff(Tensor, Tensor, Material);
/*******************************************************/


/*!

*/
double energy_Neo_Hookean_Wriggers(Tensor, double, Material);
/*******************************************************/


/*!

*/
State_Parameters compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(State_Parameters, Material);
Tensor compute_2PK_Stress_Tensor_Neo_Hookean_Wriggers(Tensor, Tensor, double, Material);
/*******************************************************/


/*!
  \fn Tensor compute_material_stiffness_density_Neo_Hookean_Wriggers(Tensor v, Tensor w, Tensor C, double J, Material MatProp)

  \brief Assemble the factorised material tensor 

  \param Tensor v
  \param Tensor w
  \param Tensor C
  \param double J
  \param Material MatProp
*/
Tensor compute_material_stiffness_density_Neo_Hookean_Wriggers(Tensor, Tensor, Tensor, double, Material);
/*******************************************************/

/*
Tensor compute_stiffness_density_Neo_Hookean_Wriggers(Tensor GRAD_I, Tensor GRAD_J, Tensor F, double J,Material MatProp)
*/
Tensor compute_stiffness_density_Neo_Hookean_Wriggers(Tensor, Tensor, Tensor, double, Material);
/*******************************************************/  

Matrix compute_D_matrix_Neo_Hookean_Wriggers(Tensor, double, Material);
/*******************************************************/

/*!
  \fn State_Parameters Von_Mises_backward_euler(State_Parameters Input_SP, Material MatProp)

  \brief Compute the plastic Von Mises model using an implicit backward euler radial returning 

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model
*/
State_Parameters Von_Mises_backward_euler(State_Parameters, Material);
/*******************************************************/

/*!
  \fn State_Parameters Von_Mises_forward_euler(State_Parameters Input_SP, Material MatProp)

  \brief Compute the plastic Von Mises model using an explicit forward euler radial returning 

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model
*/
State_Parameters Von_Mises_forward_euler(State_Parameters, Material);
/*******************************************************/


State_Parameters Drucker_Prager_backward_euler(State_Parameters, Material);
/*******************************************************/

/*!
  \fn int Frictional_Monolithic__Constitutive__(State_Parameters * Inputs_SP,Material MatProp)

  \brief Compute a family of smooth approximations of the Mohr-Coulomb model using a monolithic algorithm

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model
*/
int Frictional_Monolithic__Constitutive__(State_Parameters *,Material);
/*******************************************************/

/*!
  \fn int finite_strain_plasticity__Constitutive__(State_Parameters * Input_SP, Material MatProp)

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model
*/
int finite_strain_plasticity__Constitutive__(State_Parameters *,Material,int(* infinitesimal_plasticity)(State_Parameters *,Material));
/*******************************************************/

/*!

  \fn State_Parameters compute_1PK_Stress_Tensor_Newtonian_Fluid(Tensor P,State_Parameters Input_SP, Material MatProp_p)

*/
State_Parameters compute_1PK_Stress_Tensor_Newtonian_Fluid(State_Parameters,Material);
/*******************************************************/

/*!
\fn Tensor compute_stiffness_density_Newtonian_Fluid(Tensor GRAD_I,Tensor GRAD_J,Tensor F,Tensor dFdt,double J,double alpha4,Material MatProp_p)
*/
Tensor compute_stiffness_density_Newtonian_Fluid(Tensor,Tensor,Tensor,Tensor,double,double,Material);
/*******************************************************/

/*!
 * \fn State_Parameters compute_1PK_Stress_Tensor_Newtonian_Fluid_Incompressible(State_Parameters Intput_SP,Material MatProp_p)
 * */
 State_Parameters compute_1PK_Stress_Tensor_Newtonian_Fluid_Incompressible(State_Parameters,Material);
/*******************************************************/

/*!
 * \fn Tensor compute_stiffness_density_Newtonian_Fluid_Incompressible(Tensor GRAD_I,Tensor GRAD_J,Tensor F,Tensor dFdt,double J,double alpha4,Material MatProp_p);
 * */
Tensor compute_stiffness_density_Newtonian_Fluid_Incompressible(Tensor, Tensor,Tensor,Tensor,double,double,Material);
/*******************************************************/

/*!
 * \fn State_Parameters compute_1PK_Stress_Tensor_Bingham_Fluid(State_Parameters Intput_SP,Material MatProp_p);
 * */
State_Parameters compute_1PK_Stress_Tensor_Bingham_Fluid(State_Parameters,Material);
/*******************************************************/  

/*!
 * \fn Tensor compute_stiffness_density_Bingham_Fluid(Tensor GRAD_I, Tensor GRAD_J, Tensor F, Tensor dFdt, double J, double alpha4, Material MatProp_p)
 * */ 
Tensor compute_stiffness_density_Bingham_Fluid(Tensor,Tensor,Tensor,Tensor,double, double,Material);
/*******************************************************/

#endif

