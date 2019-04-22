#include "Utils.h"
#include "ShapeFunctions.h"

#ifndef TypeDefinitions
#define TypeDefinitions
#endif


void Initialize_Element(Element *,
			GaussPoint *,
			char *,
			Matrix ,
			int *);

void Get_RefDeformation_Gradient(GaussPoint *,
				 Element *);

Matrix Get_dNdx_matrix(Element *,
		       GaussPoint *);

Matrix Get_Stiffness_Matrix(Element *); 
