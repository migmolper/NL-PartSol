#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "grams.h"

/*********************************************************************/

ChainPtr ArrayToChain(int * A_array, int NumNodes){

  /* Variable declaration */
  ChainPtr A_chain = NULL;

  /* Loop over the array to generate a chain */
  for(int i = NumNodes-1; i > -1 ; i--){
    PushNodeTop(&A_chain, A_array[i]);    
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
    PushNodeTop(&A,End-i);    
  }

  /* Return chain */
  return A;  
}

/*********************************************************************/

int * ChainToArray(ChainPtr A_chain, int NumNodes){

  /* Variable declaration */
  ChainPtr iChain;
  int * A_Array = Allocate_Array(NumNodes, sizeof(int));

  if(A_chain == NULL){
    printf(" %s : %s \n",
	   "Error in ChainToArray",
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
  
  if(A == NULL){
    printf(" %s : %s \n",
	   "Error in LenghtChain",
	   "The chain is empty");
    exit(0);
  }
  
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
void PushNodeTop (ChainPtr * TopNodePtr, int I_new) 
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
ChainPtr ChainUnion(ChainPtr * Table, int NumTable)
{
    ChainPtr A = NULL;
    ChainPtr iTable;
    
    /* Loop in the table */
    for(int i = 0 ; i<NumTable ; i++ ){
      iTable = Table[i];
      while (iTable != NULL){
	/* Introduce a new element in the new chain */
	if (!IsPresentNode(A, iTable->I)){
	  PushNodeTop(&A, iTable->I);
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
ChainPtr ChainIntersection(ChainPtr A,ChainPtr B) 
{ 
    ChainPtr C = NULL; 
    ChainPtr iPtrA = A; 
  
    /* Traverse A and search each element of it in  */
    /* B. If the element is present in B, then  */
    /* insert the element to result  */
    while (iPtrA != NULL){ 
      if (IsPresentNode(B, iPtrA->I)) 
	PushNodeTop(&C, iPtrA->I); 
      iPtrA = iPtrA->next; 
    }
    
    return C; 
}

/*********************************************************************/

/* A utility function to print a linked list */
void printList (ChainPtr A) 
{
  /* if(A == NULL){ */
  /*     printf(" %s : %s \n", */
  /* 	     "Error in printList", */
  /* 	     "The chain is empty"); */
  /*     exit(0); */
  /* } */
  
  while (A != NULL){ 
    printf ("%d \n", A->I); 
    A = A->next; 
  } 
} 

/*********************************************************************/
