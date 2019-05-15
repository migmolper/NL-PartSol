#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"

double GetBoundaryCondition(int GP_i, int DOF_i,int TimeStep){
  double BCC_val;
  /* Fill the boundary conditions array */

  if( (TimeStep == 1) || (TimeStep == 2) ){
    
    if( (GP_i == 0 ) && (DOF_i == 1) ){
      BCC_val = -0.3;
    }
    else if( (GP_i == 20 ) && (DOF_i == 1) ){
      BCC_val = 0;
    }    
    else{
      BCC_val = NAN;
    }
   
  }
  else {

    if( (GP_i == 0 ) && (DOF_i == 0) ){
      BCC_val = 0;
    }
    else if( (GP_i == 9 ) && (DOF_i == 1) ){
      BCC_val = 0;
    }    
    else{
      BCC_val = NAN;
    }
    
  }
  
  return BCC_val;
}
