#ifndef _NODES_TOOLS_H_
#define _NODES_TOOLS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

// Global libs
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"

void generate_contour_nodes(Mesh);
/*******************************************************/

Mask     get_active_nodes__MeshTools__(Mesh);
/********************************************************************/

Mask get_active_dofs__MeshTools__(Mask, Mesh, int, int);
/********************************************************************/

Mask generate_Mask_for_static_condensation_upw__MeshTools__(Mask, Mesh, int, int);
/********************************************************************/

/**
 * @brief Get the set scalar field  MeshTools object
 * 
 * @param Field_Ap 
 * @param Field 
 * @param Nodes_p 
 * @param ActiveNodes 
 */
void get_set_scalar_field__MeshTools__(
  double * Field_Ap,
  const double * Field, 
  Element Nodes_p,
  Mask ActiveNodes);

/********************************************************************/

/**
 * @brief Get the set vectorial field  MeshTools object
 * 
 * @param Field_Ap 
 * @param Field 
 * @param Nodes_p 
 * @param ActiveNodes 
 */
void get_set_vectorial_field__MeshTools__(
  double * Field_Ap,
  const double * Field, 
  Element Nodes_p,
  Mask ActiveNodes);
/********************************************************************/

/*
	\df Matrix get_U_set_field_upw__MeshTools__(Matrix Field_upw, Element Nodes_p, Mask ActiveNodes)
*/
Matrix get_U_set_field_upw__MeshTools__(Matrix, Element, Mask);
/********************************************************************/

/*
	\df Matrix get_Pw_set_field_upw__MeshTools__(Matrix Field_upw, Element Nodes_p, Mask ActiveNodes)
*/
Matrix get_Pw_set_field_upw__MeshTools__(Matrix, Element, Mask);
/********************************************************************/

Matrix   get_nodes_coordinates__MeshTools__(ChainPtr, Matrix);
/********************************************************************/

double point_distance__MeshTools__(Matrix, Matrix);
/********************************************************************/

Matrix   compute_distance__MeshTools__(ChainPtr, Matrix, Matrix);
/********************************************************************/

Matrix   get_set_field_old__MeshTools__(Matrix, Element);
/********************************************************************/

int      get_closest_node__MeshTools__(Matrix, ChainPtr, Matrix);
/********************************************************************/

double interpolate_scalar_magnitude__MeshTools__(Matrix,Matrix);
/********************************************************************/

Tensor interpolate_vectorial_magnitude__MeshTools__(Matrix,Matrix);
/********************************************************************/

Tensor interpolate_vectorial_magnitude_gradient__MeshTools__(Matrix,Matrix);
/********************************************************************/

#endif
