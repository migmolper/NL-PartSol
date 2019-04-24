#include "Utils.h"
#include "ShapeFunctions.h"

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

GaussPoint * AllocateGaussPoints(int);

GaussPoint Initialize_GP(int,
			 Matrix,
			 double,
			 double);

void Get_Lagrangian_CG_Tensor(GaussPoint *);

void Get_Eulerian_CG_Tensor(GaussPoint *); 

