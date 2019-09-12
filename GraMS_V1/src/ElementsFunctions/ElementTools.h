#include "../ToolsLib/Utils.h"
#include "ShapeFunctions.h"

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* Functions definitions :*/

/* Matrix Get_RefDeformation_Gradient(Matrix, */
/* 				   Matrix); */

int ** GetNodalConnectivity(Mesh FEM_Mesh);

Matrix Get_B_GP(Matrix, Matrix);

/* Matrix Get_Stiffness_Matrix(Element *);  */

/* Matrix Get_Geom_Mass_Matrix(GaussPoint,Mesh); */

Matrix GetNaturalCoordinates(Matrix,Matrix,Matrix);
