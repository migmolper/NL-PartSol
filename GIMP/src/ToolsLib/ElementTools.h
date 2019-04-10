#include "Utils.h"
#include "ShapeFunctions.h"

#ifndef TypeDefinitions
#define TypeDefinitions
#endif


void Initialize_GP(GaussPoint , Matrix );

void Initialize_Element(Element ,
			GaussPoint *,
			char *,
			Matrix ,
			int *);

void Get_RefDeformation_Gradient(GaussPoint ,
				 Element );

Matrix Get_dNdx_matrix(Element *,
		       GaussPoint *);

void Get_Lagrangian_CG_Tensor(GaussPoint );

void Get_Eulerian_CG_Tensor(GaussPoint ); 

Matrix Get_Stiffness_Matrix(Element); 
