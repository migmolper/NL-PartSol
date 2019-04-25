#include "../ToolsLib/Utils.h"
#include "ShapeFunctions.h"

#ifndef TypeDefinitions
#define TypeDefinitions
#endif


Element Initialize_Element(GaussPoint *,
			   char *,
			   Matrix ,
			   int *);

void Get_RefDeformation_Gradient(GaussPoint *,
				 Element *);

Matrix Get_dNdx_matrix(Element *,
		       GaussPoint *);

Matrix Get_Stiffness_Matrix(Element *); 
