/*! \file Nodes.h
    \brief File with the prototype of the interpolation techniques here adopted
*/

#ifndef _NODES_H_
#define _NODES_H_

Mesh     ReadGidMesh__MeshTools__(char *);
/********************************************************************/

void initialise_shapefun__MeshTools__(Particle,Mesh);
/*******************************************************/

void local_search__MeshTools__(Particle, Mesh);
/*******************************************************/

void generate_contour_nodes(Mesh);
/*******************************************************/

Mask     get_active_nodes__MeshTools__(Mesh);
/********************************************************************/

Mask get_active_dofs__MeshTools__(Mask, Mesh, int, int);
/********************************************************************/

Mask generate_Mask_for_static_condensation_upw__MeshTools__(Mask, Mesh, int, int);
/********************************************************************/

void get_set_field__MeshTools__(
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
