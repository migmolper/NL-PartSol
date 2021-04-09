#include "nl-partsol.h"

/*********************************************************************/

ChainPtr memory_to_set__SetLib__(int * A_array, int NumNodes)
{

  /* Variable declaration */
  ChainPtr A_chain = NULL;

  /* Loop over the array to generate a chain */
  for(int i = NumNodes-1; i > -1 ; i--){
    push__SetLib__(&A_chain, A_array[i]);    
  }

  /* Return the chain */
  return A_chain;
}

/*********************************************************************/

ChainPtr range__SetLib__(int Init, int End){

  if(Init > End){
    printf("%s : %s \n",
	   "Error in range__SetLib__",
	   "Init > End !!!");
    exit(EXIT_FAILURE);
  }

  /* Variable declaration */
  ChainPtr A = NULL;
  int NumNodes = End-Init;

  /* Fill chain */
  for(int i = 0; i <= NumNodes ; i++){
    push__SetLib__(&A,End-i);    
  }

  /* Return chain */
  return A;  
}

/*********************************************************************/

int * set_to_memory__SetLib__(ChainPtr A_chain, int NumNodes){

  /* Variable declaration */
  ChainPtr iChain;
  int * A_Array = Allocate_Array(NumNodes, sizeof(int));

  if(A_chain == NULL){
    printf(" %s : %s \n",
	   "Error in set_to_memory__SetLib__",
	   "The chain is empty");
    exit(EXIT_FAILURE);
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

void free__SetLib__(ChainPtr * A)
{

  /* Loop index */
  ChainPtr INode = (*A);
  ChainPtr NextNode = NULL;
  
  while (INode != NULL)
  {
    NextNode = INode->next;
    free(INode);
    INode = NextNode;    
  }
  
  (*A) = NULL;
}

/*********************************************************************/

ChainPtr * alloc_table__SetLib__(int SizeTable)
{

  ChainPtr * SetTable;
  
  /* Tributary nodes for each particle */
  SetTable = (ChainPtr *)malloc(SizeTable*sizeof(ChainPtr));
  if(SetTable == NULL)
  {
    printf("%s : %s \n","alloc_table__SetLib__","Memory error");
    exit(EXIT_FAILURE);
  }
  for(int i = 0 ; i<SizeTable ; i++)
  {
    SetTable[i] = NULL;  
  }

  
  return SetTable;
}

/*********************************************************************/

void free_table__SetLib__(ChainPtr * A, int SizeTable)
{

  /* Loop in the table to free each set */
  for(int i = 0 ; i<SizeTable ; i++)
  {
    free__SetLib__(&A[i]);
  }

  /* free the table */
  free(A);  

}

/*********************************************************************/

int lenght__SetLib__(ChainPtr A)
{
  
  int NumElem = 0;
  ChainPtr INode = A;
  
  if(A == NULL)
  {
    return NumElem;
  }
  else
  {
    while(INode != NULL)
    {
      NumElem ++;
      INode = INode->next;
    }
    return NumElem;
  }
}

/*********************************************************************/

/* A utility function that returns true if data is  
   present in linked list else return false */
bool inout__SetLib__(ChainPtr Node, int I) 
{
  /* Search index */
  ChainPtr INode = Node;
  
  /* Loop in the chain */
  while (INode != NULL)
  { 
    if (INode->I == I)
    {
      return 1;
    }
    INode = INode->next; 
  }
  return 0; 
} 

/*********************************************************************/

/* A utility function to insert a node at the top of a linked list */
void push__SetLib__ (ChainPtr * TopNodePtr, int I_new) 
{
  /* Pointer to the new node of the chain */
  ChainPtr NewNodePtr;
  
  /* Allocate node */
  NewNodePtr = malloc(sizeof(ChainPtr));

  /* Insert the node at the list top (stack) */
  if(NewNodePtr != NULL)
  {
    /* put in the data */
    NewNodePtr->I = I_new;
    /* link the old list off the new node */
    NewNodePtr->next = (* TopNodePtr);
    /* move the head to point to the new node */
    (*TopNodePtr) = NewNodePtr; 
  }
  else
  {
    printf("%s %i : %s \n","Unable to insert node",I_new,"No memory available.");
  }
    
}

/*********************************************************************/

/* A utility function to extract a node of a linked list */
void pop__SetLib__(ChainPtr * TopNodePtr, int I_trash) 
{
  ChainPtr iPtr = (*TopNodePtr);
  ChainPtr PrevPtr = NULL;
  ChainPtr AuxPtr;

  while(iPtr != NULL){

    /* If the node is the one we want to extract */
    if(iPtr->I == I_trash)
    {
      /* If the node is the first in the chain */
      if(PrevPtr == NULL)
      {
        AuxPtr = iPtr->next;
        free(iPtr);
        (* TopNodePtr) = AuxPtr;
      }
      /* If the node is in the middle or at the end */
      else
      {
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

ChainPtr copy__SetLib__(ChainPtr start1)
{
  
  ChainPtr start2=NULL;
  ChainPtr previous=NULL;

  while(start1!=NULL)
  {
  
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

ChainPtr create_circular_set__SetLib__(ChainPtr A)
{
  ChainPtr Head = A; 
  ChainPtr A_circular = NULL;
  ChainPtr A_aux = NULL;
  ChainPtr New_Node = NULL;

  while(A!=NULL)
  {
  
    New_Node = malloc (sizeof(ChainPtr));
    New_Node->I=A->I;
    New_Node->next=Head;

    if(A_circular == NULL)
    {
      A_circular = New_Node;
      A_aux = New_Node;
    }
    else
    {
      A_aux->next=New_Node;
      A_aux=New_Node;          
    }

    A=A->next;

  }

  return A_circular;
}

/*********************************************************************/

/* Function to get union of two linked lists
   A and B */
ChainPtr union__SetLib__(ChainPtr * Table, int NumTable)
{
    ChainPtr A = NULL;
    ChainPtr iTable;
    
    /* Loop in the table */
    for(int i = 0 ; i<NumTable ; i++ ){
      iTable = Table[i];
      while (iTable != NULL){
	/* Introduce a new element in the new chain */
	if (!inout__SetLib__(A, iTable->I)){
	  push__SetLib__(&A, iTable->I);
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
ChainPtr intersection__SetLib__(ChainPtr A,ChainPtr B) 
{ 
    ChainPtr C = NULL; 
    ChainPtr iPtrA = A; 
  
    /* Traverse A and search each element of it in  */
    /* B. If the element is present in B, then  */
    /* insert the element to result  */
    while (iPtrA != NULL){ 
      if (inout__SetLib__(B, iPtrA->I)) 
	push__SetLib__(&C, iPtrA->I); 
      iPtrA = iPtrA->next; 
    }
    
    return C; 
}

/*********************************************************************/

/* A utility function to print a linked list */
void print__SetLib__(ChainPtr A) 
{  
  while (A != NULL){ 
    printf ("%d -> ", A->I); 
    A = A->next; 
  } 
  printf("NULL \n");
} 

/*********************************************************************/

void order__SetLib__(ChainPtr * List1, ChainPtr * List0, Matrix Dist)
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
    push__SetLib__(List1,I_DistMax);
    pop__SetLib__(List0,I_DistMax);
    
    /* Recursive */
    order__SetLib__(List1,List0,Dist);
  }
  
}

/*********************************************************************/
