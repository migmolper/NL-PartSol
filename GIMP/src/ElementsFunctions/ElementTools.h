#include "../ToolsLib/Utils.h"
#include "ShapeFunctions.h"

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* Functions definitions :*/
Matrix Get_RefDeformation_Gradient(Matrix,
				   Matrix);

Matrix Get_dNdx_matrix(Matrix,
		       Matrix,
		       int,int);

/* Matrix Get_Stiffness_Matrix(Element *);  */

Matrix Get_Geom_Mass_Matrix(GaussPoint,Element);

void ApplyBoundaryCondition(Matrix,int);
