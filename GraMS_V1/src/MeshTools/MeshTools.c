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
	    puts("Error in GetNodalConnectivity(): Max number of neighbour reached !");
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
    for(int i = 0 ; i<dNdX_GP.N_cols ; i++){
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

ChainPtr ChainAlloc(int * NodesIndex, int NumNodes){

  /* Variable declaration */
  ChainPtr NewChain = NULL;

  /* Loop over the list to generate a List */
  for(int i = NumNodes-1; i > -1 ; i--){
    PushNode(&NewChain, NodesIndex[i]);    
  }

  /* Return the list */
  return NewChain;
}

/*********************************************************************/

void FreeChain(ChainPtr A){

  /* Loop index */
  ChainPtr INode = A;
  ChainPtr NextNode = NULL;
  
  while (INode != NULL){
    NextNode = INode->next;
    free(INode);
    INode = NextNode;    
   }

}

/*********************************************************************/

int LenghtChain(ChainPtr A){

  int NumElem = 0;
  ChainPtr INode = A;

  while(INode != NULL){
    NumElem ++;
    INode = INode->next;
  }

  return NumElem;
}

/*********************************************************************/

/* A utility function that returns true if data is  
   present in linked list else return false */
bool IsPresentNode (ChainPtr Node, int I) 
{
  /* Search index */
  ChainPtr INode = Node;
  
  /* Loop in the chain */
  while (INode != NULL){ 
    if (INode->I == I){
      return 1;
    }
    INode = INode->next; 
  }
  return 0; 
} 

/*********************************************************************/

/* A utility function to insert a node at the top of a linked list */
void PushNode (ChainPtr * TopNodePtr, int I_new) 
{
  /* Pointer to the new node of the chain */
  ChainPtr NewNodePtr;
  
  /* Allocate node */
  NewNodePtr = (ChainPtr)malloc(sizeof(Chain));

  /* Insert the node at the list top (stack) */
  if(NewNodePtr != NULL){
    /* put in the data */
    NewNodePtr->I = I_new;
    /* link the old list off the new node */
    NewNodePtr->next = (* TopNodePtr);
    /* move the head to point to the new node */
    (*TopNodePtr) = NewNodePtr; 
  }
  else{
    printf("Unable to insert node %i. No memory available.\n",
	   I_new);
  }
    
}

/*********************************************************************/

/* A utility function to extract a node of a linked list */
void PopNode (ChainPtr * TopNodePtr, int I_trash) 
{
  ChainPtr iPtr = (*TopNodePtr);
  ChainPtr PrevPtr = NULL;
  ChainPtr AuxPtr;

  while(iPtr != NULL){

    /* If the node is the one we want to extract */
    if(iPtr->I == I_trash){
      /* If the node is the first in the chain */
      if(PrevPtr == NULL){
	AuxPtr = iPtr->next;
	free(iPtr);
	(* TopNodePtr)->next = AuxPtr;
      }
      /* If the node is in the middle or at the end */
      else{
	PrevPtr->next = iPtr->next;
	free(iPtr);
      }
      /* Once the node is located, breack the loop */
      break;
    }

    /* The previous is the index */
    PrevPtr = iPtr;
    /* Update pointer index */
    iPtr = iPtr->next;
  }

}

/*********************************************************************/

/* Function to get union of two linked lists
   A and B */
ChainPtr GetUnion(ChainPtr A, ChainPtr B) 
{ 
    ChainPtr C = NULL; 
    ChainPtr iPtrA = A, iPtrB = B; 
  
    /* Insert all elements of A to the result list */
    while (iPtrA != NULL){
      /* Introduce a new element in the new chain */
      PushNode(&C, iPtrA->I);
      /* Updtate the interator index */
      iPtrA = iPtrA->next; 
    }
    
    /* Insert those elements of B which are not  */
    /* present in result list */ 
    while (iPtrB != NULL){
      /* Introduce a new element in the new chain */      
      if (!IsPresentNode(C, iPtrB->I)){
	PushNode(&C, iPtrB->I);
      }
      /* Updtate the interator index */
      iPtrB = iPtrB->next; 
    }

    return C; 
}

/*********************************************************************/
  
/* Function to get intersection of two linked lists 
   A and B */
ChainPtr GetIntersection(ChainPtr A,ChainPtr B) 
{ 
    ChainPtr C = NULL; 
    ChainPtr iPtrA = A; 
  
    /* Traverse A and search each element of it in  */
    /* B. If the element is present in B, then  */
    /* insert the element to result  */
    while (iPtrA != NULL){ 
      if (IsPresentNode(B, iPtrA->I)) 
	PushNode(&C, iPtrA->I); 
      iPtrA = iPtrA->next; 
    }
    
    return C; 
}

/*********************************************************************/

/* A utility function to print a linked list */
void printList (ChainPtr A) 
{ 
    while (A != NULL) 
    { 
        printf ("%d ", A->I); 
        A = A->next; 
    } 
} 

/*********************************************************************/
