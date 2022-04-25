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

void local_search__MeshTools__(Particle, Mesh);
/*******************************************************/

Matrix   compute_N__MeshTools__(Element, Particle, Mesh);
/********************************************************************/

Matrix   compute_dN__MeshTools__(Element, Particle, Mesh);
/********************************************************************/

int push_forward_dN__MeshTools__(
  double * Gradient_n1_p,
  const double * Gradient_n_p,
  const double * d_phi,
  unsigned NumNodes);
/********************************************************************/

#endif
