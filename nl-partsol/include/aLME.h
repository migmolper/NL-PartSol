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

#ifndef _aLME_H_
#define _aLME_H_

/****************************************************************************/

/*!
  \fn void initialize__aLME__(Particle MPM_Mesh, Mesh FEM_Mesh);

  \brief  Initialize LME shape functions 

  \param MPM_Mesh : Variable with the particle information
  \param FEM_Mesh : Variable wih information of the background set of nodes
*/
void initialize__aLME__(Particle, Mesh);
/****************************************************************************/

/*!
  \fn void initialize_beta__aLME__(Matrix Beta, double Gamma, double DeltaX);

  \brief  Compute the value of the thermalization parameter using a circular support.
  
  \param Beta : Termalization tensor
  \param Gamma : Adimensional paramter to control the regularization parameter.
  \param DeltaX : Minimum size in the all nodal set.
*/
void initialize_beta__aLME__(Matrix, double, double);
/****************************************************************************/

/*!
  \fn  Matrix metric__aLME__(Tensor F);

  \brief Return a general metric tensor to compute the locality parameter in the LME shape functions.

  \param F : Deformation gradient
*/
Matrix metric__aLME__(Tensor);
/****************************************************************************/

/*!
  \fn void update_lambda_Newton_Rapson__aLME__(Matrix l, Matrix lambda, Matrix Beta)

  \brief Function ot get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical methodis the Newton-Rapson.

  \param p : Current particle
  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Beta : Thermalization tensor.
*/
void update_lambda_Newton_Rapson__aLME__(int, Matrix, Matrix, Matrix);
/****************************************************************************/

/*!
  \fn void lambda_Nelder_Mead__aLME__(Matrix l, Matrix lambda, Matrix Beta)

  \brief Function ot get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical method is the Nelder-Mead.

  \param p : Current particle
  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Beta : Thermalization tensor.
*/
void update_lambda_Nelder_Mead__aLME__(int, Matrix, Matrix, Matrix);
/****************************************************************************/

/*!
  \fn Matrix p__aLME__(Matrix l, Matrix lambda, Matrix Beta)

  \brief Function to get the value of the shape function "pa" (1 x neighborhood) in the
  neighborhood nodes.

  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Beta : Thermalization tensor.

*/
Matrix p__aLME__(Matrix, Matrix, Matrix);
/****************************************************************************/

/*!
  \fn Matrix dp__aLME__(Matrix l, Matrix p) 

  \brief Compute the value of the shape function gradient "dp" (dim x neighborhood) in the neighborhood nodes

  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param p : Set with the evaluation of the shape function in the neighborhood nodes.
*/
Matrix dp__aLME__(Matrix, Matrix);
/****************************************************************************/

/*!

 \fn void local_search__aLME__(Particle MPM_Mesh, Mesh FEM_Mesh)

  \brief Compute the local search for the LME (update to reduce the number of computational nodes)

  \param MPM_Mesh : Variable with the particle information
  \param FEM_Mesh : Variable wih information of the background set of nodes
*/
void local_search__aLME__(Particle, Mesh);
/****************************************************************************/

#endif
