/*! \file LME.h
    \brief Local maximum entropy shape functions

    Shape functions based in :
    "" Local maximum-entropy approximation schemes : a seamless 
    bridge between finite elements and meshfree methods ""
    by M.Arroyo and M.Ortiz, 2006.

    Here we employ the same nomenclature as in the paper. With the single
    different of the "l" variable wich represents the distances between the
    evaluation point and the neighborhood nodes.
*/

#ifndef _LME_H_
#define _LME_H_

void     initialize__LME__(GaussPoint, Mesh);
Matrix   beta_isotropic__LME__(Matrix, Matrix, double);
Matrix   beta_anisotropic__LME__(Matrix, Matrix);
Matrix   lambda__LME__(Matrix, Matrix, Matrix);
Matrix   p__LME__(Matrix, Matrix, Matrix);
Matrix   dp__LME__(Matrix, Matrix);
ChainPtr tributary__LME__(Matrix, Matrix, int, Mesh);

#endif
