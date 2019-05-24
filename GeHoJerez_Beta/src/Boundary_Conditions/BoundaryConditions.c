#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"

/**********************************************************************/

void ApplyBoundaryCondition_Nod(Matrix Phi_n_Nod,int TimeStep){
  
  /* Fill the boundary conditions array */
  
  if( (TimeStep == 1) ){
    Phi_n_Nod.nM[1][0] = -1;
  }
  else {
    Phi_n_Nod.nM[0][0] = 0;
  }

  Phi_n_Nod.nM[1][8] = 0;
 
}

/*********************************************************************/

void ApplyBoundaryCondition_GP(Matrix Phi_n_GP,int TimeStep){
  
  /* Fill the boundary conditions array */
  
  if( (TimeStep == 1)){
    Phi_n_GP.nM[1][0] = -1;
  }
  else {
    Phi_n_GP.nM[0][0] = 0;
  }

  Phi_n_GP.nM[1][7] = 0;
 
}

/*********************************************************************/
