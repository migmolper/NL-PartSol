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
  \fn Tensor LinearElastic(Tensor Strain, Tensor Stress, Material Properties)
  
  \brief This function is devoted to make a material point behaves as a
  linear elastic material.
  
  \param Strain : Strain field of each particle 
  \param Stress : Stress field of each particle
  \param Properties : Define the material properties of the particle
*/
Tensor LinearElastic(Tensor, Tensor, Material);
/*******************************************************/

/*!
  \fn void compute_particle_Damage(int p, GaussPoint Particles, Mesh Nodes)

  \brief Compute if a particle is damaged and update it damage parameter chi

  \param p : Particle
  \param Particles : Information of the particle mesh
  \param Nodes : Informaction with the set of nodes
*/
void compute_particle_Damage(int, GaussPoint, Mesh);
/*******************************************************/

/*!
  
*/
double energy_Saint_Venant_Kirchhoff(Tensor, Material);
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
Tensor Itegration_Stress_Configurational_Midpoint(Tensor,Tensor, Tensor, Material);
/*******************************************************/

/*!
  
*/
Tensor Itegration_Stress_Average_Strain(Tensor, Tensor, Tensor, Material);
/*******************************************************/

/*!
  
*/
Tensor Itegration_Stress_Average(Tensor, Tensor, Tensor, Material);
/*******************************************************/


#endif

