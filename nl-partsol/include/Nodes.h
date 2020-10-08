/*! \file Nodes.h
    \brief File with the prototype of the interpolation techniques here adopted
*/

#ifndef _NODES_H_
#define _NODES_H_

Mask     generate_NodalMask__MeshTools__(Mesh);
Mask     generate_Mask_for_static_condensation__MeshTools__(Mask, Mesh);
Matrix   get_set_field__MeshTools__(Matrix, Element, Mask);
Matrix   get_set_field_old__MeshTools__(Matrix, Element);
Matrix   compute_distance__MeshTools__(ChainPtr, Matrix, Matrix);
int      get_closest_node__MeshTools__(Matrix, ChainPtr, Matrix);
bool     inout_convex_set__MeshTools__(Matrix, ChainPtr, Matrix);
ChainPtr get_nodal_locality__MeshTools__(int, Mesh);
void     get_nodal_connectivity__MeshTools__(Mesh);
Matrix   get_nodes_coordinates__MeshTools__(ChainPtr, Matrix);
double   mesh_size__MeshTools__(Mesh);
Matrix   compute_N__MeshTools__(Element, GaussPoint, Mesh );
Matrix   compute_dN__MeshTools__(Element, GaussPoint, Mesh);

#endif
