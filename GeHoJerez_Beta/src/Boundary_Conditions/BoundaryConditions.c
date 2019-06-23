#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"

/**********************************************************************/

void BCC_Nod_Momentum(Mesh FEM_Mesh,
		      BoundaryConditions BCC_Momentum,
		      Matrix Nodal_MOMENTUM, int TimeStep)
/*
  Apply the momentum boundary conditions over the nodes 
*/
{

  /* 0º  Check the time step */
  if( (TimeStep < 0) ||
      (TimeStep > BCC_Momentum.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }

  /* 1º Define auxilar variables */
  int Dim_BCC; /* Dimension where is applied the boundary condition */
  int Id_BCC; /* Index of the node where we apply the BCC */

  /* 2º Set the dimension */
  Dim_BCC = BCC_Momentum.Dim;
  
  /* 3º Loop over the nodes of the boundary */
  for(int i=0 ; i<BCC_Momentum.NumNodes ; i++){

    /* 4º Get the index of the element */
    Id_BCC = BCC_Momentum.Nodes[i];

    /* 5º Assign the BC to the nodes */
    Nodal_MOMENTUM.nM[Dim_BCC][Id_BCC] = BCC_Momentum.Value.Fx[TimeStep];
          
  } 

}

/*********************************************************************/

void BCC_GP_Forces(GaussPoint GP_Mesh, BoundaryConditions Loads, int TimeStep)
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

/* void BC_Nod_Springs(GaussPoint GP_Mesh, BoundaryConditions Loads, int TimeStep) */
/* /\* */
/*   Add to a surface springs to simulate the winkler */
/* *\/ */
/* { */
  
/* } */
