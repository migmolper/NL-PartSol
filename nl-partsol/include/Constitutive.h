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
  \fn void EigenerosionAlgorithm(int p,
  Matrix ji,
  Matrix W,
  Matrix Mass,
  Matrix Rho,
  Matrix Stress,
  Material MatPro,
  ChainPtr * Beps,
  double DeltaX)

  \brief Function to compute is a material point is or not eroded. 
  Here the notation is the same as in \cite Pandolfi_2012

  \param p : Index of the particle
  \param ji : Damage status 
  \param W : Internal work 
  \param Mass : Mass field
  \param Rho : Density field
  \param Stress : Stress field of each particle
  \param Properties : Define the material properties of the particle
  \param B_eps : Define the particles close to each particle
  \param DeltaX : Mesh size
*/

void EigenerosionAlgorithm(int, Matrix, Matrix,
			   Matrix, Matrix,
			   Matrix, Material,
			   ChainPtr *, double);

/*******************************************************/

/*!
  \fn void EigensofteningAlgorithm(int p,
  Matrix ji,
  Matrix Strain,
  Matrix StrainF,
  Matrix Mass,
  Matrix Stress,
  Material Properties,
  ChainPtr * Beps)

  \brief Function to compute is a material point is or not eroded. 
  Here the notation is the same as in \cite Navas_2017_ES

  \param p : Index of the particle
  \param ji : Matrix with the value of the damage parameter.
  \param Strain : Strain field of each particle. 
  \param StrainF : Value of the strain field at the failure init. 
  \param Mass : Matrix with the mass of the particle.
  \param Stress : Stress field of each particle
  \param Properties : Define the material properties of the particle
  \param Beps : Table with the list of neighbours per particle.
*/
void EigensofteningAlgorithm(int, Matrix, Matrix,
			     Matrix, Matrix, Matrix,
			     Material, ChainPtr *);

/*******************************************************/


#endif

