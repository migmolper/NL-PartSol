#include "Utils.h"
#include "ShapeFunctions.h"

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

GaussPoint * AllocateGaussPoints(int);

void Initialize_GP(GaussPoint *,
		   int,
		   Matrix,
		   double,
		   double);

void Get_Lagrangian_CG_Tensor(GaussPoint *);

void Get_Eulerian_CG_Tensor(GaussPoint *); 

