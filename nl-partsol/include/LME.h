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
  \fn void initialize__LME__(Particle MPM_Mesh, Mesh FEM_Mesh);

  \brief  Initialize LME shape functions 

  \param MPM_Mesh : Variable with the particle information
  \param FEM_Mesh : Nodes information
*/
void initialize__LME__(Particle, Mesh);
/****************************************************************************/

/*!
  \fn double beta__LME__(Matrix l, double Gamma, double DeltaX);

  \brief  Compute the value of the thermalization parameter using a circular support.
  
  \param Gamma : Adimensional paramter to control the regularization parameter.
  \param DeltaX : Minimum size in the all nodal set.
*/
double beta__LME__(double, double);
/****************************************************************************/

/*!
  \fn  Matrix metric__LME__(Tensor F);

  \brief Return a general metric tensor to compute the locality parameter in the LME shape functions.

  \param F : Deformation gradient
*/
Matrix metric__LME__(Tensor);
/****************************************************************************/

/*!
  \fn void update_lambda_Newton_Rapson__LME__(Matrix l, Matrix lambda, Matrix Metric, double Beta)

  \brief Function ot get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical methodis the Newton-Rapson.

  \param p : Current particle
  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Metric : Measure for the norm definition.
  \param Beta : Thermalization parameter.
*/
void update_lambda_Newton_Rapson__LME__(int, Matrix, Matrix, Matrix, double);
/****************************************************************************/

/*!
  \fn void lambda_Nelder_Mead__LME__(Matrix l, Matrix lambda, Matrix Metric, double Beta)

  \brief Function ot get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical method is the Nelder-Mead.

  \param p : Current particle
  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Metric : Measure for the norm definition.
  \param Beta : Thermalization parameter.
*/
void update_lambda_Nelder_Mead__LME__(int, Matrix, Matrix, Matrix, double);
/****************************************************************************/

/*!
  \fn Matrix p__LME__(Matrix l, Matrix lambda, Matrix Metric, double Beta)

  \brief Function to get the value of the shape function "pa" (1 x neighborhood) in the
  neighborhood nodes.

  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Metric : Measure for the norm definition.
  \param Beta : Thermalization parameter.

*/
Matrix p__LME__(Matrix, Matrix, Matrix, double);
/****************************************************************************/

/*!
  \fn Matrix dp__LME__(Matrix l, Matrix p) 

  \brief Compute the value of the shape function gradient "dp" (dim x neighborhood) in the neighborhood nodes

  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param p : Set with the evaluation of the shape function in the neighborhood nodes.
*/
Matrix dp__LME__(Matrix, Matrix);
/****************************************************************************/

/*!
  \fn Matrix tributary__LME__(Matrix X_p, Matrix Metric, double Beta_p, int I0, Mesh FEM_Mesh);

  \brief Compute a set with the sourrounding nodes of the particle

  \param X_p : Coordinates of the particle
  \param Metric : Measure for the norm definition
  \param Beta_p : Thermalization parameter of the particle
  \param I0 : Index of the closest node to the particle
  \param FEM_Mesh : Variable wih information of the background set of nodes
*/
ChainPtr tributary__LME__(int,Matrix, Matrix, double, int, Mesh);


#endif
