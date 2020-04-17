#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "grams.h"

/*********************************************************************/

ChainPtr Pointer_to_Set(int * A_array, int NumNodes){

  /* Variable declaration */
  ChainPtr A_chain = NULL;

  /* Loop over the array to generate a chain */
  for(int i = NumNodes-1; i > -1 ; i--){
    push_to_Set(&A_chain, A_array[i]);    
  }

  /* Return the chain */
  return A_chain;
}

/*********************************************************************/

ChainPtr RangeChain(int Init, int End){

  if(Init > End){
    printf("%s : %s \n",
	   "Error in RangeChain",
	   "Init > End !!!");
    exit(0);
  }

  /* Variable declaration */
  ChainPtr A = NULL;
  int NumNodes = End-Init;

  /* Fill chain */
  for(int i = 0; i <= NumNodes ; i++){
    push_to_Set(&A,End-i);    
  }

  /* Return chain */
  return A;  
}

/*********************************************************************/

int * Set_to_Pointer(ChainPtr A_chain, int NumNodes){

  /* Variable declaration */
  ChainPtr iChain;
  int * A_Array = Allocate_Array(NumNodes, sizeof(int));

  if(A_chain == NULL){
    printf(" %s : %s \n",
	   "Error in Set_to_Pointer",
	   "The chain is empty");
    exit(0);
  }

  /* Loop over the chain to generate an array */
  iChain = A_chain;
  for(int i = 0;
      (i<NumNodes) || (iChain != NULL) ;
      i++, iChain = iChain->next){
    A_Array[i] = iChain->I;
  }

  /* Return the array */
  return A_Array;
}

/*********************************************************************/

void free_Set(ChainPtr A){

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

int get_Lenght_Set(ChainPtr A){
    
  int NumElem = 0;
  ChainPtr INode = A;

  if(A == NULL){
    return NumElem;
  }
  else{
    while(INode != NULL){
      NumElem ++;
      INode = INode->next;
    }
    return NumElem;
  }
}

/*********************************************************************/

/* A utility function that returns true if data is  
   present in linked list else return false */
bool is_in_Set(ChainPtr Node, int I) 
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
void push_to_Set (ChainPtr * TopNodePtr, int I_new) 
{
  /* Pointer to the new node of the chain */
  ChainPtr NewNodePtr;
  
  /* Allocate node */
  NewNodePtr = malloc(sizeof(ChainPtr));

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
    printf("%s %i : %s \n",
	   "Unable to insert node",
	   I_new,
	   "No memory available.");
  }
    
}

/*********************************************************************/

/* A utility function to extract a node of a linked list */
void pop_from_Set (ChainPtr * TopNodePtr, int I_trash) 
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
	(* TopNodePtr) = AuxPtr;
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

ChainPtr CopyChain(ChainPtr start1){
  
  ChainPtr start2=NULL;
  ChainPtr previous=NULL;

  while(start1!=NULL){
  
    ChainPtr temp = malloc (sizeof(ChainPtr));
    temp->I=start1->I;
    temp->next=NULL;

    if(start2==NULL)
      {
        start2=temp;
        previous=temp;
      }
    else
      {
        previous->next=temp;
        previous=temp;          
      }
    start1=start1->next;
  }
  return start2;
}

/*********************************************************************/

/* Function to get union of two linked lists
   A and B */
ChainPtr get_Union_Of(ChainPtr * Table, int NumTable)
{
    ChainPtr A = NULL;
    ChainPtr iTable;
    
    /* Loop in the table */
    for(int i = 0 ; i<NumTable ; i++ ){
      iTable = Table[i];
      while (iTable != NULL){
	/* Introduce a new element in the new chain */
	if (!is_in_Set(A, iTable->I)){
	  push_to_Set(&A, iTable->I);
	}
	/* Updtate the iterator index */
	iTable = iTable->next;
      }
      
    }

    return A;
}

/*********************************************************************/
  
/* Function to get intersection of two linked lists 
   A and B */
ChainPtr get_Intersection_Of(ChainPtr A,ChainPtr B) 
{ 
    ChainPtr C = NULL; 
    ChainPtr iPtrA = A; 
  
    /* Traverse A and search each element of it in  */
    /* B. If the element is present in B, then  */
    /* insert the element to result  */
    while (iPtrA != NULL){ 
      if (is_in_Set(B, iPtrA->I)) 
	push_to_Set(&C, iPtrA->I); 
      iPtrA = iPtrA->next; 
    }
    
    return C; 
}

/*********************************************************************/

/* A utility function to print a linked list */
void print_Set(ChainPtr A) 
{  
  while (A != NULL){ 
    printf ("%d \n", A->I); 
    A = A->next; 
  } 
} 

/*********************************************************************/

Matrix get_set_Coordinates(ChainPtr Set, Matrix X0, Matrix Coordinates)
{
  int Ndim = NumberDimensions;
  int SizeSet = get_Lenght_Set(Set);
  int I = 0;

  /* Allocate output */
  Matrix Set_Coordinates = MatAlloc(SizeSet,Ndim);

  /* Loop in the set */
  ChainPtr Aux_Set = Set;
  while (Aux_Set != NULL){ 
    /* Get coordinates local coodinates of each node in the set */
    for(int i = 0 ; i<Ndim ; i++){
      Set_Coordinates.nM[I][i] =  X0.nV[i] - Coordinates.nM[Aux_Set->I][i];
    }
    /* Update index */
    I++;	
    Aux_Set = Aux_Set->next; 
  }
  
  return Set_Coordinates;
}


/*********************************************************************/

void order_Set(ChainPtr * List1, ChainPtr * List0, Matrix Dist)
/*
  Ordenate recursively and array with distances and get a chain 
  with the positions in orden
*/
{
	
  /* Iterate while List0 is full of numbers */
  if((*List0) != NULL){
    
    ChainPtr INode = (* List0);
    double DistMax = 0.0;
    int I_DistMax;

    /* Loop over the chain */
    while(INode != NULL){
      /* Get the max distance of the matrix */
      if(Dist.nV[INode->I] > DistMax){
	DistMax = Dist.nV[INode->I];
	I_DistMax = INode->I;
      }
      /* Continue iterating */      
      INode = INode->next; 
    }
    
    /* Push and Pop node */
    push_to_Set(List1,I_DistMax);
    pop_from_Set(List0,I_DistMax);
    
    /* Recursive */
    order_Set(List1,List0,Dist);
  }
  
}

/*********************************************************************/
