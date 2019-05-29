#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"

/**********************************************************************/

void ApplyBoundaryCondition_Nod(Mesh FEM_Mesh, Matrix Nod_Field, int Time)
/*
  Apply the boundary conditions
*/
{

  int Kind_BCC; /* Kind of the boundary condition 0/1 */
  int Id_BCC; /* Index of the node where we apply the BCC */
  
  for(int i=0 ; i<FEM_Mesh.NumNodesBound ; i++){
    Id_BCC = FEM_Mesh.NodesBound[i][0];
    for(int j=0 ; j<NumberDOF ; j++){
      Kind_BCC = FEM_Mesh.NodesBound[i][j+1];
      if(Kind_BCC == 1){
	if(NumberDOF>1)
	  Nod_Field.nM[j][Id_BCC] = FEM_Mesh.ValueBC.nM[i][j];
	if(NumberDOF==1)
	  Nod_Field.nV[Id_BCC] = FEM_Mesh.ValueBC.nV[i];	
      }
    }    
  } 
 
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
