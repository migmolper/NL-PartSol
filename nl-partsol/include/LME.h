/*! \file LME.h
    \brief Local maximum entropy shape functions

    Shape functions based in :
    "" Local maximum-entropy approximation schemes : a seamless 
    bridge between finite elements and meshfree methods ""
    by \cite Arroyo2006.

    Here we employ the same nomenclature as in the paper. With the single
    different of the "l" variable wich represents the distances between the
    evaluation point and the neighborhood nodes.
*/

#ifndef _LME_H_
#define _LME_H_

/****************************************************************************/
/*!
  \fn void initialize__LME__(GaussPoint MPM_Mesh, Mesh FEM_Mesh);

  \brief  Initialize LME shape functions 

  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param Mesh FEM_Mesh : Nodes information
*/
void     initialize__LME__(GaussPoint, Mesh);
/****************************************************************************/
/*!
  \fn Matrix  beta_isotropic__LME__(Matrix Beta, Matrix l, double Gamma);

  \brief  Update the value of the thermalization parameter using a circular support.

  \param Matrix Beta : Previous value of the thermalization parameter 
  \param Matrix l : Matrix with the distances from nodes in the neiborghood to the particle
  \param double Gamma : Adimensional paramter to control the regularization parameter
*/
Matrix   beta__LME__(Matrix, Matrix, double);
/****************************************************************************/
/*!
  \fn Matrix beta_anisotropic__LME__(Matrix Beta, Matrix f);

  \brief Update the value of the termalization paramter with the increment
  of the deformation gradient

  \param Matrix Beta : Previous value of the thermalization parameter 
  \param Matrix f : Increment of the deformation gradient
*/
Matrix   lambda__LME__(Matrix, Matrix, Matrix);
Matrix   p__LME__(Matrix, Matrix, Matrix);
Matrix   dp__LME__(Matrix, Matrix);
ChainPtr tributary__LME__(Matrix, Matrix, int, Mesh);
#endif
