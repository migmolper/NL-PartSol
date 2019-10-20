#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h> 
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

Matrix Get_B_GP(Matrix dNdX_GP)
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

  /* 1º Select the case to solve */
  switch(NumberDimensions){
    
  case 1:  /* 1D stress tensor */
    
    puts("Error in Get_dNdi_matrix() : 1D cases not implemented yet");
    exit(0);
    
  case 2: /* 2D stress tensor */
    
    /* 2º Allocate the output */
    B_GP = MatAlloc(3,2*dNdX_GP.N_cols);
    
    /* 4º Fill the array with the nodal partial derivation 
       of the reference element */    
    for(int i = 0 ; i<dNdX_GP.N_rows ; i++){
      B_GP.nM[0][2*i] = dNdX_GP.nM[0][i];
      B_GP.nM[1][2*i] = 0;
      B_GP.nM[2][2*i] = dNdX_GP.nM[1][i];
      
      B_GP.nM[0][2*i + 1] = 0;
      B_GP.nM[1][2*i + 1] = dNdX_GP.nM[1][i];
      B_GP.nM[2][2*i + 1] = dNdX_GP.nM[0][i];      
    }

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

/* A utility function to insert a node at the begining of a linked list*/
void push (Chain ** A, int new_I) 
{ 
    /* allocate node */
    Chain * new_Node = (Chain *)malloc(sizeof(Node)); 
  
    /* put in the data */
    new_Node->I = new_I; 
  
    /* link the old list off the new node */
    new_Node->next = (* A); 
  
    /* move the head to point to the new node */
    (*A) = new_Node; 
}

/*********************************************************************/

/* A utility function that returns true if data is  
   present in linked list else return false */
bool isPresent (Chain * A, int I) 
{ 
    Chain * t = A; 
    while (t != NULL){ 
        if (t->I == I) 
            return 1; 
        t = t->next; 
    } 
    return 0; 
} 

/*********************************************************************/

/* Function to get union of two linked lists
   A and B */
Chain * getUnion(Chain * A,
		Chain * B) 
{ 
    Chain * C = NULL; 
    Chain * t1 = A, * t2 = B; 
  
    // Insert all elements of A to the result list 
    while (t1 != NULL){ 
      push(&C, t1->I); 
      t1 = t1->next; 
    }
    
    // Insert those elements of B which are not 
    // present in result list 
    while (t2 != NULL){ 
        if (!isPresent(C, t2->I)) 
            push(&C, t2->I); 
        t2 = t2->next; 
    } 
  
    return C; 
}

/*********************************************************************/
  
/* Function to get intersection of two linked lists 
   A and B */
Chain * getIntersection(Chain * A,  
			Chain * B) 
{ 
    Chain * C = NULL; 
    Chain * t1 = A; 
  
    // Traverse A and search each element of it in 
    // B. If the element is present in B, then 
    // insert the element to result 
    while (t1 != NULL){ 
      if (isPresent(B, t1->I)) 
	push (&C, t1->I); 
      t1 = t1->next; 
    }
    
    return C; 
}

/*********************************************************************/
    
/* A utility function to print a linked list*/
void printList (Chain * A) 
{ 
    while (A != NULL) 
    { 
        printf ("%d ", A->I); 
        A = A->next; 
    } 
} 

/*********************************************************************/
