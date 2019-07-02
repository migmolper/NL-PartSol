#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "../InOutFun/InOutFun.h"

/**********************************************************************/

BoundaryConditions SetBCC(int NumNodes, int * Nodes,
			  char * Info, char * Curve_File){
  
  /* Define boundary conditions */
  BoundaryConditions BCC;

  /* Apply this boundary conditions */
  BCC.Nodes = Nodes;
  BCC.NumNodes = NumNodes;

  /* Asign curve of values */
  BCC.Value = ReadCurve(Curve_File);

  /* Copy information of the BCC */
  strcpy(BCC.Info,Info);

  /* Return boundary condition */
  return BCC;
}


/**********************************************************************/

void BCC_Nod_Momentum(Mesh FEM_Mesh, Matrix Nodal_MOMENTUM, int TimeStep)
/*
  Apply the momentum boundary conditions over the nodes 
*/
{

  /* 1º Define auxilar variables */
  int Id_BCC; /* Index of the node where we apply the BCC */
  

  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.TOP.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }
  /* 2º Loop over the nodes of the TOP boundary */
  for(int i = 0 ; i<FEM_Mesh.TOP.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = FEM_Mesh.TOP.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      Nodal_MOMENTUM.nM[j][Id_BCC] =
	FEM_Mesh.TOP.Value.Fx[TimeStep][j];
    }          
  }

  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.BOTTOM.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }
  /* 2º Loop over the nodes of the BOTTOM boundary */
  for(int i = 0 ; i<FEM_Mesh.BOTTOM.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = FEM_Mesh.BOTTOM.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      Nodal_MOMENTUM.nM[j][Id_BCC] =
	FEM_Mesh.BOTTOM.Value.Fx[TimeStep][j];
    }          
  }

  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.RIGHT.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }
  /* 2º Loop over the nodes of the RIGHT boundary */
  for(int i = 0 ; i<FEM_Mesh.RIGHT.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = FEM_Mesh.RIGHT.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      Nodal_MOMENTUM.nM[j][Id_BCC] =
	FEM_Mesh.RIGHT.Value.Fx[TimeStep][j];
    }          
  }

  
  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.LEFT.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }  
  /* 2º Loop over the nodes of the LEFT boundary */
  for(int i = 0 ; i<FEM_Mesh.LEFT.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = FEM_Mesh.LEFT.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      Nodal_MOMENTUM.nM[j][Id_BCC] =
	FEM_Mesh.LEFT.Value.Fx[TimeStep][j];
    }          
  } 
  

}

/*********************************************************************/



/*********************************************************************/

/* void BC_Nod_Springs(GaussPoint GP_Mesh, BoundaryConditions Loads, int TimeStep) */
/* /\* */
/*   Add to a surface springs to simulate the winkler */
/* *\/ */
/* { */
  
/* } */
