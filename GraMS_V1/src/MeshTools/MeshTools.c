#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../MathTools/MathTools.h"
#include "../Constitutive/Constitutive.h"
#include "MeshTools.h"

#define MAXNEIGHBOUR 10

/*********************************************************************/

int ** GetNodalConnectivity(Mesh FEM_Mesh){

  /* 0º Create an auxiliar table of pointer to store the information */
  int ** TableNeighbourNode;
  int ** NodeNeighbour;
  int * NumNeighbour;
  TableNeighbourNode =
    (int **)Allocate_MatrixZ(FEM_Mesh.NumNodesMesh,
			     MAXNEIGHBOUR,sizeof(int));
  NumNeighbour =
    (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh,sizeof(int));

  /* 1º Start the search of neighbour for each node */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    /* 2º Loop over all the elements in the mesh */
    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){
      /* 3º Loop over the all the node in an element */
      for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
	/* 4º If my node belong to the element */
	if(FEM_Mesh.Connectivity[j][k] == i){
	  if(NumNeighbour[i] == MAXNEIGHBOUR){
	    puts("Error in GetNodalConnectivity() : Max number of neighbour reached !!! ");
	    exit(0);
	  }
	  TableNeighbourNode[i][NumNeighbour[i]] = j;
	  NumNeighbour[i] += 1;
	}
      }
    }      
  }

  /* 5º Resize the table of pointer */
  NodeNeighbour = (int **)malloc((unsigned)FEM_Mesh.NumNodesMesh *
				 sizeof(int *));
  /* 6º Loop over the pointer table */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    /* 7º Save space for the each nodes neighbour */
    NodeNeighbour[i] = malloc((unsigned) (NumNeighbour[i] + 1) *
			      sizeof(int));
    /* 8º Check if it is not out of memory  */
    if (NodeNeighbour[i] == NULL){
      puts("Error in GetNodalConnectivity() : Out of memory !!! ");
      exit(0);
    }
    /* 9º Fill the new table */
    NodeNeighbour[i][0] = NumNeighbour[i];
    for(int j = 1 ; j<=NumNeighbour[i] ; j++){
      NodeNeighbour[i][j] = TableNeighbourNode[i][j-1];
    }    
  }
  
  /* 10º Free memory */
  free(TableNeighbourNode);
  free(NumNeighbour);

  /* 11º Return data */
  return NodeNeighbour; 
  
}

/*********************************************************************/

Matrix Get_B_GP(Matrix X_EC_GP,
		Matrix Element)
/*
   Get the B matrix (Usual in the classical formulation of 
   the finite element method )
   Inputs:
   - Matrix X_NC_GP : Element coordinates
   - Matrix Element : Coordinates of the element 
   (NumNodesElem x NumberDimensions)

   Outputs : Matrix B
*/
{

  /* 0º Define variables */
  /* Declaration of the output matrix (NdimVecStrain x Nnodes*Ndim) */
  Matrix B_GP;
  /* Derivatives of the SF evaluates in the GP (Ndim x Ndim) */
  Matrix dNdx_XG_GP; 

  /* 1º Select the case to solve */
  switch(NumberDimensions){
    
  case 1:  /* 1D stress tensor */
    
    puts("Error in Get_dNdi_matrix() : 1D cases not implemented yet");
    exit(0);
    
  case 2: /* 2D stress tensor */
    
    /* 2º Allocate the output */
    B_GP = MatAlloc(3,2*Element.N_rows);

    if(strcmp(Element.Info,"Quadrilateral") == 0){
      /* 3º Get the element gradient in natural coordinates*/
      dNdx_XG_GP = Get_dNdX_Q4(X_EC_GP,Element);
    }
        
    /* 4º Fill the array with the nodal partial derivation 
       of the reference element */    
    for(int i = 0 ; i<Element.N_rows ; i++){
      B_GP.nM[0][2*i] = dNdx_XG_GP.nM[0][i];
      B_GP.nM[1][2*i] = 0;
      B_GP.nM[2][2*i] = dNdx_XG_GP.nM[1][i];
      
      B_GP.nM[0][2*i + 1] = 0;
      B_GP.nM[1][2*i + 1] = dNdx_XG_GP.nM[1][i];
      B_GP.nM[2][2*i + 1] = dNdx_XG_GP.nM[0][i];      
    }

    /* 5º Free memory */
    FreeMat(dNdx_XG_GP);

    break;
    
  case 3: /* 3D stress tensor */
    puts("Error in Get_dNdi_matrix() : 3D cases not implemented yet");
    exit(0);
    
  default :
    puts("Error in Get_dNdi_matrix() : Wrong case select");
    exit(0);
  }
  
  return B_GP;
}

/*********************************************************************/

void GetNaturalCoordinates(Matrix X_EC_GP,
			   Matrix X_GC_GP,
			   Matrix Element_GC_Nod)
/* 
   The function return the natural coordinates of a point 
   inside of the element.
 
   Inputs :
   - Coordinates of the element nodes
   - Initial coordinate of the point to start the search 
   - Derivative function of the element

   Depending of the kind of element, we employ differents types
   of shape functions
*/
{
  
  X_EC_GP = Newton_Rapson(Get_X_GC_Q4,Element_GC_Nod,
			  Get_F_Ref_Q4,Element_GC_Nod,
			  X_GC_GP,X_EC_GP);
}


