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


// Shape functions auxilar tools
#include "Nodes/Nodes-Tools.h"

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

 \fn int local_search__aLME__(Particle MPM_Mesh, Mesh FEM_Mesh)
  \brief Compute the local search for the LME (update to reduce the number of computational nodes)
  \param MPM_Mesh : Variable with the particle information
  \param FEM_Mesh : Variable wih information of the background set of nodes
*/
int local_search__aLME__(Particle, Mesh);
/****************************************************************************/

#endif
