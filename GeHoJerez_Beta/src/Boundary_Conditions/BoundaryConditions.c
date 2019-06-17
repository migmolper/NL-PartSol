#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"

/**********************************************************************/

void BCC_Nod_Momentum(Mesh FEM_Mesh,
		      BoundayConditions BCC_Momentum,
		      Matrix Nod_Field, int TimeStep)
/*
  Apply the boundary conditions
*/
{

  /* 0º  Check the time step */
  if( (TimeStep < 0) ||
      (TimeStep > BCC_Momentum.Value.Num)){
    puts("Error in BCC_Nod_Momentum() : The time step is out of the curve !!");
    exit(0);
  }
  
  int Kind_BCC; /* Kind of the boundary condition 0/1 */
  int Id_BCC; /* Index of the node where we apply the BCC */
  
  for(int i=0 ; i<BCC_Momentum.NumNodes ; i++){
    Id_BCC = FEM_Mesh.NodesBound[i][0];
    for(int j=0 ; j<Nod_Field.N_rows ; j++){
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

void BCC_GP_Forces(GaussPoint GP_Mesh, BoundayConditions Loads, int TimeStep)
/*
  Forces defined in the Gauss Points :
  Inputs
*/
{
  /* 0º  Check the time step */
  if( (TimeStep < 0) ||
      (TimeStep > Loads.Value.Num)){
    puts("Error in BCC_GP_Forces() : The time step is out of the curve !!");
    exit(0);
  }
  
  /* 1º Fill the matrix with the local forces */
  for(int i = 0 ; i<Loads.NumNodes ; i++){
    /* 2º Check if this GP has a force applied */
    if( (Loads.Nodes[i] > GP_Mesh.NumGP) ||
	(Loads.Nodes[i] < 0)){
      puts("Error in BCC_GP_Forces() : This GP does not exist !!");
      exit(0);
    }
    /* 3º Apply the force in the node */
    GP_Mesh.Phi.F.nM[Loads.Nodes[i]][Loads.Dim] =
      Loads.Value.Fx[TimeStep]; 
  }

}

/*********************************************************************/

/* void ApplyBoundaryCondition_GP(Matrix Phi_n_GP,int TimeStep){ */
  
/*   /\* Fill the boundary conditions array *\/ */
  
/*   if( (TimeStep == 1)){ */
/*     Phi_n_GP.nM[1][0] = -1; */
/*   } */
/*   else { */
/*     Phi_n_GP.nM[0][0] = 0; */
/*   } */

/*   Phi_n_GP.nM[1][7] = 0; */
 
/* } */

/*********************************************************************/
