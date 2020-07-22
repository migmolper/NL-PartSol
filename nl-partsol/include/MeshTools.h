/*! \file ShapeFun.h
    \brief File with the prototype of the interpolation techniques here adopted
*/

#ifndef _MESHTOOLS_H_
#define _MESHTOOLS_H_

Matrix get_set_Coordinates(ChainPtr, Matrix, Matrix);
Matrix get_set_Field(Matrix, Element);
int get_closest_node_to(Matrix, ChainPtr, Matrix);
bool InOut_Element(Matrix, ChainPtr, Matrix);


ChainPtr DiscardElements(ChainPtr, Matrix, Matrix, Mesh);
ChainPtr get_locality_of_node(int, Mesh);
void GetNodalConnectivity(Mesh);
Matrix ElemCoordinates(ChainPtr, Matrix);
double mesh_size__MeshTools__(Mesh);
Matrix compute_N__MeshTools__(Element, GaussPoint, Mesh );
Matrix compute_dN__MeshTools__(Element, GaussPoint, Mesh);

#endif
