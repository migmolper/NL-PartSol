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

// Auxiliar functions to compute the shape functions
static double fa__LME__(Matrix, Matrix, double);
static double logZ__LME__(Matrix, Matrix, double);
static Matrix r__LME__(Matrix, Matrix);
static Matrix J__LME__(Matrix, Matrix, Matrix);

// Auxiliar functions for the Neldel Mead in the LME
static void initialise_lambda__LME__(int, Matrix, Matrix, Matrix, double);
static Matrix gravity_center_Nelder_Mead__LME__(Matrix);
static void order_logZ_simplex_Nelder_Mead__LME__(Matrix, Matrix);
static void expansion_Nelder_Mead__LME__(Matrix, Matrix, Matrix, Matrix, Matrix,
                                         double, double);
static void contraction_Nelder_Mead__LME__(Matrix, Matrix, Matrix, Matrix,
                                           Matrix, double, double);
static void shrinkage_Nelder_Mead__LME__(Matrix, Matrix, Matrix, double);
static ChainPtr tributary__LME__(int, Matrix, double, int, Mesh);

/****************************************************************************/

/*!
  \fn void initialize__LME__(Particle MPM_Mesh, Mesh FEM_Mesh);

  \brief  Initialize LME shape functions 

  \param MPM_Mesh : Variable with the particle information
  \param FEM_Mesh : Variable wih information of the background set of nodes
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
  \fn void update_lambda_Newton_Rapson__LME__(Matrix l, Matrix lambda, double Beta)

  \brief Function ot get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical methodis the Newton-Rapson.

  \param p : Current particle
  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Beta : Thermalization parameter.
*/
void update_lambda_Newton_Rapson__LME__(int, Matrix, Matrix, double);
/****************************************************************************/

/*!
  \fn void lambda_Nelder_Mead__LME__(Matrix l, Matrix lambda, double Beta)

  \brief Function ot get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical method is the Nelder-Mead.

  \param p : Current particle
  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Beta : Thermalization parameter.
*/
void update_lambda_Nelder_Mead__LME__(int, Matrix, Matrix, double);
/****************************************************************************/

/*!
  \fn Matrix p__LME__(Matrix l, Matrix lambda, double Beta)

  \brief Function to get the value of the shape function "pa" (1 x neighborhood) in the
  neighborhood nodes.

  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Beta : Thermalization parameter.

*/
Matrix p__LME__(Matrix, Matrix, double);
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

 \fn void local_search__LME__(Particle MPM_Mesh, Mesh FEM_Mesh)

  \brief Compute the local search for the LME (update to reduce the number of computational nodes)

  \param MPM_Mesh : Variable with the particle information
  \param FEM_Mesh : Variable wih information of the background set of nodes
*/
void local_search__LME__(Particle, Mesh);
/****************************************************************************/

#endif
