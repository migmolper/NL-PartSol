#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 

/* Set of nodes */
typedef struct Set { 
  int I; /* Index of the node */
  struct Set * next;  /* Pointer to the next element */
} Set; 

/* Pointer to a set */
typedef Set * SetPtr;

/*********************************************************************/

void push_to_Set(SetPtr * TopNodePtr, int I_new)
/* 
   A utility function to insert a node at the top of a linked list 
*/
{
  /* Pointer to the new node of the set */
  SetPtr NewNodePtr;
  
  /* Allocate node */
  NewNodePtr = malloc(sizeof(SetPtr));

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


SetPtr ArrayToSet(int * A_array, int NumNodes)
{

  /* Variable declaration */
  SetPtr A_set = NULL;

  /* Loop over the array to generate a set */
  for(int i = NumNodes-1; i > -1 ; i--){
    push_to_Set(&A_set, A_array[i]);    
  }

  /* Return the set */
  return A_set;
}

/*********************************************************************/

SetPtr RangeSet(int Init, int End)
{

  if(Init > End){
    printf("%s : %s \n",
	   "Error in RangeSet",
	   "Init > End !!!");
    exit(0);
  }

  /* Variable declaration */
  SetPtr A = NULL;
  int NumNodes = End-Init;

  /* Fill set */
  for(int i = 0; i <= NumNodes ; i++){
    push_to_Set(&A,End-i);    
  }

  /* Return set */
  return A;  
}

/*********************************************************************/

/* int * SetToArray(SetPtr A_set, int NumNodes) */
/* { */

/*   /\* Variable declaration *\/ */
/*   SetPtr iSet; */
/*   int * A_Array = Allocate_Array(NumNodes, sizeof(int)); */

/*   if(A_set == NULL){ */
/*     printf(" %s : %s \n", */
/* 	   "Error in SetToArray", */
/* 	   "The set is empty"); */
/*     exit(0); */
/*   } */

/*   /\* Loop over the set to generate an array *\/ */
/*   iSet = A_set; */
/*   for(int i = 0; */
/*       (i<NumNodes) || (iSet != NULL) ; */
/*       i++, iSet = iSet->next){ */
/*     A_Array[i] = iSet->I; */
/*   } */

/*   /\* Return the array *\/ */
/*   return A_Array; */
/* } */

/*********************************************************************/

void free_Set(SetPtr A)
{

  /* Loop index */
  SetPtr INode = A;
  SetPtr NextNode = NULL;
  
  while (INode != NULL){
    NextNode = INode->next;
    free(INode);
    INode = NextNode;    
   }
}

/*********************************************************************/

int get_lenght_Set(SetPtr A)
{
  
  if(A == NULL){
    printf(" %s : %s \n",
	   "Error in get_lenght_Set",
	   "The set is empty");
    exit(0);
  }
  
  int NumElem = 0;
  SetPtr INode = A;

  while(INode != NULL){
    NumElem ++;
    INode = INode->next;
  }

  return NumElem;
}

/*********************************************************************/

bool is_in_Set(SetPtr Node, int I)
/*
  A utility function that returns true if data is  
  present in linked list else return false 
*/
{
  /* Search index */
  SetPtr INode = Node;
  
  /* Loop in the set */
  while (INode != NULL){ 
    if (INode->I == I){
      return 1;
    }
    INode = INode->next; 
  }
  return 0; 
} 

/*********************************************************************/


void pop_out_Set(SetPtr * TopNodePtr, int I_trash)
/* 
   A utility function to extract a node of a linked list
 */
{
  SetPtr iPtr = (*TopNodePtr);
  SetPtr PrevPtr = NULL;
  SetPtr AuxPtr;

  while(iPtr != NULL){

    /* If the node is the one we want to extract */
    if(iPtr->I == I_trash){
      /* If the node is the first in the set */
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

SetPtr copy_Set(SetPtr start1)
{
  
  SetPtr start2=NULL;
  SetPtr previous=NULL;

  while(start1!=NULL){
  
    SetPtr temp = malloc (sizeof(SetPtr));
    temp->I=start1->I;
    temp->next=NULL;

    if(start2==NULL){
      start2=temp;
      previous=temp;
    }
    else{
      previous->next=temp;
      previous=temp;          
    }
    start1=start1->next;
  }
  return start2;
}

/*********************************************************************/

SetPtr get_Union_Set(SetPtr * Table, int NumTable)
/* 
   Function to get union of n sets 
*/
{
    SetPtr A = NULL;
    SetPtr iTable;
    
    /* Loop in the table */
    for(int i = 0 ; i<NumTable ; i++ ){
      iTable = Table[i];
      while (iTable != NULL){
	/* Introduce a new element in the new set */
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
  
SetPtr get_Intersection_Set(SetPtr A, SetPtr B)
/* 
   Function to get intersection of two linked lists 
   A and B 
*/
{ 
    SetPtr C = NULL; 
    SetPtr iPtrA = A; 
  
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

void print_Set(SetPtr A)
/*
  A utility function to print a linked list 
*/
{  
  while (A != NULL){ 
    printf ("%d \n", A->I); 
    A = A->next; 
  } 
} 

/*********************************************************************/


void main(){
  SetPtr A = RangeSet(0, 10);
  SetPtr B = RangeSet(5, 15);
  SetPtr C = get_Intersection_Set(A,B);
  print_Set(C);
  free_Set(A);
  free_Set(B);
  free_Set(C);    
}
