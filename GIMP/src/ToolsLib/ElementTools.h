#include "Utils.h"
#include "ShapeFunctions.h"

#ifndef TypeDefinitions
#define TypeDefinitions
#endif


void Initialize_GP(GaussPoint *, Vector *);

void Initialize_Element(Element *,
			GaussPoint *,
			char *,
			double **,
			int *);

void Get_RefDeformation_Gradient(GaussPoint *,
				 Element *);

Matrix Get_dNdx_matrix(Element *,
		       GaussPoint *);

void Get_Lagrangian_CG_Tensor(GaussPoint *);

void Get_Eulerian_CG_Tensor(GaussPoint *); 
