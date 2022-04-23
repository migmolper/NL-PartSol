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
  \fn State_Parameters Von_Mises_forward_euler(State_Parameters Input_SP, Material MatProp)

  \brief Compute the plastic Von Mises model using an explicit forward euler radial returning 

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model
*/
State_Parameters Von_Mises_forward_euler(State_Parameters, Material);
/*******************************************************/

/*!
  \fn int compute_1PK_Von_Mises(State_Parameters IO_State, Material MatProp)
  
  \brief Compute Von Mises yield

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model

*/
int compute_1PK_Von_Mises(State_Parameters IO_State, Material MatProp);
/*******************************************************/




/*!
  \fn int compute_1PK_Matsuoka_Nakai(State_Parameters IO_State, Material MatProp)

  \brief Compute a family of smooth approximations of the Mohr-Coulomb model using a monolithic algorithm

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model
*/
int compute_1PK_Matsuoka_Nakai(State_Parameters IO_State, Material MatProp);
/*******************************************************/

/*!
  \fn int compute_1PK_Lade_Duncan(State_Parameters IO_State, Material MatProp)

  \brief Compute a family of smooth approximations of the Mohr-Coulomb model using a monolithic algorithm

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model
*/
int compute_1PK_Lade_Duncan(State_Parameters IO_State, Material MatProp);
/*******************************************************/

/*!
  \fn int compute_1PK_Modified_Lade_Duncan(State_Parameters IO_State, Material MatProp)

  \brief Compute a family of smooth approximations of the Mohr-Coulomb model using a monolithic algorithm

  \param Input_SP : State parameters of the particle
  \param MatProp : Material properties of the model
*/
int compute_1PK_Modified_Lade_Duncan(State_Parameters IO_State, Material MatProp);
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



int compute_1PK_elastoplastic_tangent_matrix(double *Stiffness_density, const double *dN_alpha, const double *dN_beta, const State_Parameters IO_State);

#endif

