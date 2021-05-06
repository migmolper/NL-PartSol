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
Tensor volumetric_stress__LinearElastic__(double, Material);
Tensor deviatoric_stress__LinearElastic__(Tensor, Material);
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

Tensor compute_1PK_Stress_Tensor_Saint_Venant_Kirchhoff(Tensor,Tensor,Material);
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
Tensor compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(Tensor, Tensor, double, Material);
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


Plastic_status finite_strains_plasticity_Von_Mises(Tensor, Tensor, Tensor, Plastic_status, Material, double);
Plastic_status infinitesimal_strains_plasticity_Von_Mises(Tensor, Tensor, Plastic_status, Material);
/*******************************************************/  

Plastic_status finite_strains_plasticity_Drucker_Prager_Sanavia(Tensor, Tensor, Tensor, Plastic_status, Material, double);
Plastic_status infinitesimal_strains_plasticity_Drucker_Prager_Sanavia(Tensor, Tensor, Plastic_status, Material);
/*******************************************************/


// Compresible Newtonian fluid

/*!

  \fn Tensor compute_1PK_Stress_Tensor_Newtonian_Fluid(Tensor P,Tensor F,Tensor dFdt,double J, Material MatProp_p)

*/
Tensor compute_1PK_Stress_Tensor_Newtonian_Fluid(Tensor,Tensor,Tensor,double,Material);

/*!
\fn Tensor compute_stiffness_density_Newtonian_Fluid(Tensor GRAD_I,Tensor GRAD_J,Tensor F,Tensor dFdt,double J,double alpha4,Material MatProp_p)
*/
Tensor compute_stiffness_density_Newtonian_Fluid(Tensor,Tensor,Tensor,Tensor,double,double,Material);
/*******************************************************/  

#endif

