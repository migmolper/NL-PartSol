#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"

Matrix GetBoundaryCondition(int TimeStep){
  Matrix BCC_val = MatAlloc(2,2);
  /* Fill the boundary conditions array */

  if(TimeStep == 0){
    BCC_val.nM[0][0] = 0;
    BCC_val.nM[1][0] = -1;

    BCC_val.nM[0][1] = 0;
    BCC_val.nM[1][1] = 0;
  }
  else if(TimeStep > 0){
    BCC_val.nM[0][0] = 0;
    BCC_val.nM[1][0] = 0;

    BCC_val.nM[0][1] = 0;
    BCC_val.nM[1][1] = 0;
  }
  
  return BCC_val;
}
