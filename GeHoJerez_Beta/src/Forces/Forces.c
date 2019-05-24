
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"
#include "../ElementsFunctions/ShapeFunctions.h"

/*********************************************************************/


/* Matrix VolumeForces(){ */
/* } */

/* Matrix SurfaceForces(){ */
/* } */

/* Matrix LineForces(){ */
/* } */

/* Matrix PuntualForces(){ */
/* } */


/* Matrix GaussPointForces(GaussPoint * MaterialPoints) */
/* /\* */
/*   Forces defined in the Gauss Points : */
/*   Inputs  */
/* *\/ */
/* { */

/* } */

Matrix NodalForces(Element * ElementMesh,int NumberActiveNodes)
/*
  Forces defined in the nodes of the element:
  * Input :
  -> Mesh with the active elements
  -> Number of active nodes
  * Return : 
  -> Matrix (Nx1) with the nodal forces
*/
{
  Matrix F =  MatAlloc(NumberActiveNodes*NumberDimensions,1);
  
  return F;
}

