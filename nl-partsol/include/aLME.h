/*! \file aLME.h
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

void initialize__aLME__(GaussPoint, Mesh);

Tensor beta__aLME__(Tensor, Tensor);

Tensor cut_off__aLME__(Tensor, Tensor);

Tensor lambda__aLME__(Matrix, Tensor, Tensor);

Matrix p__aLME__(Matrix, Tensor, Tensor);

Matrix dp__aLME__(Matrix, Matrix);

ChainPtr tributary__aLME__(Matrix, Tensor, int, Mesh);

#endif