/*! \file Constitutive.h
    \brief File with the prototype of the contitutive models
*/


#ifndef _CONSTITUTIVE_H_
#define _CONSTITUTIVE_H_

/*******************************************************/


/*!
 * This function is devoted to make a material point behaves as a
 * solid rigid.
 * @param Strain 
 */
Tensor SolidRigid(Tensor Strain);

/*******************************************************/


/*!
 * This function is devoted to make a material point behaves as a
 * linear elastic material.
 * @param Strain 
 * @param Stress 
 * @param Mat 
 */
Tensor LinearElastic(Tensor Strain, Tensor Stress, Material Mat);

/*******************************************************/

/*! \fn void EigenerosionAlgorithm(int p,
  Matrix ji,
  Matrix W,
  Matrix Mass,
  Matrix Rho,
  Matrix Stress,
  Material MatPro,
  ChainPtr * Beps,
  double DeltaX)

  \brief Function to compute is a material point is or not eroded. 
  Here the notation is the same as in the paper : 
  A.Pandolfi & M.Ortiz.
  An eigenerosion approach to brittle fracture. 
  International Journal for Numerical Methods in Enginnering.
  92:694-714, 2012.

  Inputs :
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

void EigenerosionAlgorithm(int p, Matrix ji, Matrix W,
			   Matrix Mass, Matrix Rho,
			   Matrix Stress, Material MatPro,
			   ChainPtr * Beps, double DeltaX);

/*******************************************************/

/*!
 *  Pedro Navas, Rena C. Yu, Bo Li & Gonzalo Ruiz.
 *  Modeling the dynamic fracture in concrete: 
 *  an eigensoftening meshfree approach.
 *  International Journal of Impact Engineering.
 *  113 (2018) 9-20
 *  NOTE : Here notation is the same as in the paper.
 *  Inputs :
 *  @param Ji_k0 : Matrix with the value of the damage parameter.
 *  @param Mass : Matrix with the mass of the GP.
 *  @param StrainF : Value of the strain field at the failure init. 
 *  @param Beps : Table with the list of neighbours per GP.
 *  @param Neps : Number of neighbours per GP
 *  @param Num_GP : Number of GP of the mesh.
 */
void EigensofteningAlgorithm(int, Matrix, Matrix,
			     Matrix, Matrix, Matrix,
			     Material, ChainPtr *);

/*******************************************************/


#endif

