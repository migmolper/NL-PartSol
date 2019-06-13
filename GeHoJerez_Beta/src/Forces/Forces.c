
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"
#include "../ElementsFunctions/ShapeFunctions.h"

/*********************************************************************/


Matrix GaussPointForces(int Num_GP, int TimeStep)
/*
  Forces defined in the Gauss Points :
  Inputs
*/
{
  /* Definition of the local traction */
  Matrix t_GP;

  /* Allocate the local traction matrix for the gauss points */
  t_GP = MatAllocZ(Num_GP,NumberDimensions);

  /* Fill the matrix with the local forces */
  for(int i = 0 ; i<Num_GP ; i++){
    for(int j = 0 ; j<NumberDimensions ; j++){
      t_GP.nM[i][j] = ;
    }
  }
  
  return t_GP;
}

/*********************************************************************/

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

