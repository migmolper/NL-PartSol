#ifndef _SHAPE_FUNCTIONS_H_
#define _SHAPE_FUNCTIONS_H_

#include <string.h>

// Global libs
#include "Nodes/aLME.h"
#include "Nodes/GIMP.h"
#include "Nodes/H8.h"
#include "Nodes/LME.h"
#include "Nodes/Q4.h"
#include "Nodes/T3.h"
#include "Nodes/T4.h"


void initialise_shapefun__MeshTools__(Particle,Mesh);
/*******************************************************/

/**
 * @brief Search particles in the neighbourhood of the closest node
 * 
 * @param MPM_Mesh Information of the particle set
 * @param FEM_Mesh Information of the nodal set
 * @return int 
 */
int local_search__MeshTools__(Particle MPM_Mesh, Mesh FEM_Mesh);
/*******************************************************/

Matrix   compute_N__MeshTools__(Element, Particle, Mesh);
/********************************************************************/

Matrix   compute_dN__MeshTools__(Element, Particle, Mesh);
/********************************************************************/

/*!
  \brief Extrapolate the shape function gradient of the n configuration to the 
  n+1 configuration using the incremental deformation gradient
  \param[in] Gradient_n_p Shape function gradient in the n configuration
  \param[in] d_phi Incremental deformation gradient
  \param[in] NumNodes Number of nodes of the shape function gradient
  \param[out] STATUS Status variables 
  \returns The shape function gradient in the n+1 configuration
*/
double * push_forward_dN__MeshTools__(
  const double * Gradient_n_p,
  const double * d_phi,
  unsigned NumNodes,
  int * STATUS);
/********************************************************************/

#endif
