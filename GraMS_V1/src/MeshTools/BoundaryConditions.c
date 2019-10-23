#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../MathTools/MathTools.h"
#include "../InOutFun/InOutFun.h"

/**********************************************************************/

void BCC_Nod_VALUE(Mesh FEM_Mesh, Matrix Nodal_VALUE, int TimeStep)
/*
  Apply the boundary conditions over the nodes 
*/
{

  /* 1º Define auxilar variables */
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Id_BCC; /* Index of the node where we apply the BCC */

  /* 2º Loop over the the boundaries */
  for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++){
    /* 3º Get the number of nodes of this boundarie */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;
    /* 4º Get the number of dimensions where the BCC it is applied */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;
    for(int j = 0 ; j<NumNodesBound ; j++){
      /* 5º Get the index of the node */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      /* 6º Loop over the dimensions of the boundary condition */
      for(int k = 0 ; k<NumDimBound ; k++){
	/* 7º Apply only if the direction is active (1) */
	if( (FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1) ||
	    (FEM_Mesh.Bounds.BCC_i[i].Dir[k] == -1)){
	  /* 8º Check if the curve it is on time */
	  if( (TimeStep < 0) ||
	      (TimeStep > FEM_Mesh.Bounds.BCC_i[i].Value[k].Num)){
	    puts("Error in BCC_Nodal_VALUE() : The time step is out of the curve !!");
	    exit(0);
	  }
	  /* 9º Assign the boundary condition */
	  Nodal_VALUE.nM[k][Id_BCC] =
	    FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*
	    (double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];
	}
      }
    }    
  }
  
  

}

/*********************************************************************/
