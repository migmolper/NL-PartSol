/*! \file Nodes.h
    \brief File with the prototype of the interpolation techniques here adopted
*/

#ifndef _NODES_H_
#define _NODES_H_

Mesh     ReadGidMesh__MeshTools__(char *);
Mask     generate_NodalMask__MeshTools__(Mesh);
Mask     generate_Mask_for_static_condensation__MeshTools__(Mask, Mesh);
Matrix   get_set_field__MeshTools__(Matrix, Element, Mask);
Matrix   get_set_field_old__MeshTools__(Matrix, Element);
Matrix   compute_distance__MeshTools__(ChainPtr, Matrix, Matrix);
int      get_closest_node__MeshTools__(Matrix, ChainPtr, Matrix);
Matrix   get_nodes_coordinates__MeshTools__(ChainPtr, Matrix);
Matrix   compute_N__MeshTools__(Element, GaussPoint, Mesh );
Matrix   compute_dN__MeshTools__(Element, GaussPoint, Mesh);

#endif
