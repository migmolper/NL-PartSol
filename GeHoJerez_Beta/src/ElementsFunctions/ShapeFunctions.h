#include <stdio.h>
#include <stdlib.h>
#include "../ToolsLib/Utils.h"

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* Line of two nodes functions */
Matrix L2(Matrix);
Matrix dL2(Matrix);
Matrix Get_GlobalCoordinates_L2(Matrix,Matrix);
Matrix Get_RefDeformation_Gradient_L2(Matrix,Matrix);

/* Triangle of three nodes functions */
Matrix T3(Matrix);
Matrix dT3(Matrix);

/* Quadrilateral of four nodes functions */
Matrix Q4(Matrix);
Matrix dQ4(Matrix);
Matrix Get_F_Ref_Q4(Matrix,Matrix);
Matrix Get_dNdX_Q4(Matrix,Matrix);
Matrix Get_X_GC_Q4(Matrix,Matrix);
