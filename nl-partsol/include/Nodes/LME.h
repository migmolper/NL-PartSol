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

// Shape functions auxilar tools
#include "Nodes/Nodes-Tools.h"
#include <petsctao.h>

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

PetscErrorCode lambda__LME__(const PetscScalar *l_a, PetscScalar *lambda_tr,
                             PetscScalar Beta, PetscInt N_a);

/****************************************************************************/

/**
 * @brief Function to get the value of the shape function "pa" (1 x neighborhood) in the neighborhood nodes.
 * 
 * @param l_a Set than contanins vector form neighborhood nodes to particle.
 * @param lambda_tr Lagrange multiplier.
 * @param Beta Thermalization parameter.
 * @param N_a Number of nodes in the neighborhood
 * @return PetscScalar* 
 */
PetscScalar *p__LME__(const PetscScalar *l_a, const PetscScalar *lambda_tr,
                      PetscScalar Beta, PetscInt N_a);
/****************************************************************************/

/**
 * @brief Compute the value of the shape function gradient "dp" (dim x neighborhood) in the neighborhood nodes
 * 
 * @param l_a Set than contanins vector form neighborhood nodes to particle.
 * @param lambda_tr Lagrange multiplier.
 * @param Beta Thermalization parameter.
 * @param N_a Number of nodes in the neighborhood
 * @return PetscScalar* 
 */
PetscScalar *dp__LME__(const PetscScalar *l_a, PetscScalar *lambda_tr,
                              PetscScalar Beta, PetscInt N_a);

/****************************************************************************/

/**
 * @brief Compute the local search for the LME (update to reduce the number of computational nodes)
 * 
 * @param MPM_Mesh Variable with the particle information
 * @param FEM_Mesh Variable wih information of the background set of nodes
 * @return int 
 */
int local_search__LME__(Particle MPM_Mesh, Mesh FEM_Mesh);
/****************************************************************************/

#endif
