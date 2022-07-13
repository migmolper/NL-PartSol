#ifndef _LME_H_
#define _LME_H_

/**
 * @file LME.h
 * @author Miguel Molinos (@migmolper)
 * @brief Local maximum entropy shape functions
 *   Shape functions based in :
 *   "" Local maximum-entropy approximation schemes : a seamless
 *   bridge between finite elements and meshfree methods ""
 *   by \cite Arroyo2006. Here we employ the same nomenclature as in the paper.
 * With the single different of the "l" variable wich represents the distances
 * between the evaluation point and the neighborhood nodes.
 * @version 0.1
 * @date 2022-07-13
 *
 * @copyright Copyright (c) 2022
 *
 */

// Shape functions auxilar tools
#include "Nodes/Nodes-Tools.h"
#include <petsctao.h>

/****************************************************************************/

/**
 * @brief Initialize LME shape functions
 *
 * @param MPM_Mesh Structure with the particle set information
 * @param FEM_Mesh Structure with the nodal set information
 */
void initialize__LME__(Particle MPM_Mesh, Mesh FEM_Mesh);
/****************************************************************************/

/**
 * @brief Get the thermalization parameter beta using the global variable
 * gamma_LME.
 *
 * @param Gamma User-defined adimensional parameter used to control the value of
 * the thermalization parameter.
 * @param h_avg Average mesh size
 * @return double
 */
double beta__LME__(double, double);

/****************************************************************************/

/**
 * @brief Function to compute the lagrange multiplier required by the LME shape
 * functions. In this implementation we use TAO solver to solve the
 * minimization problem
 * @param l_a Distance from the nodes to the particle (xp - xa)
 * @param lambda_tr Value of the Lagrange multiplier, the input is a trial
 * @param Beta Thermalization parameter
 * @param N_a Number of nodes
 * @return PetscErrorCode
 */
PetscErrorCode lambda__LME__(const PetscScalar *l_a, PetscScalar *lambda_tr,
                             PetscScalar Beta, PetscInt N_a);

/****************************************************************************/

/**
 * @brief Function to get the value of the shape function "pa" (1 x
 * neighborhood) in the neighborhood nodes.
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
 * @brief Compute the value of the shape function gradient "dp" (dim x
 * neighborhood) in the neighborhood nodes
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
 * @brief Compute the local search for the LME (update to reduce the number of
 * computational nodes). Search the closest node to the particle based in its
 * previous position.
 *
 * @param MPM_Mesh Variable with the particle information
 * @param FEM_Mesh Variable wih information of the background set of nodes
 * @return int
 */
int local_search__LME__(Particle MPM_Mesh, Mesh FEM_Mesh);
/****************************************************************************/

#endif
