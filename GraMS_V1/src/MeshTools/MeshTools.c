#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../MathTools/MathTools.h"
#include "../Constitutive/Constitutive.h"
#include "MeshTools.h"

/*********************************************************************/

void GetNodalConnectivity(Mesh FEM_Mesh){

  /* Variable declaration */
  int * Element_Connectivity;
  int NumNodesElem;
  
  /* 1º Start the search of neighbour for each node */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    /* 2º Loop over all the elements in the mesh */
    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){
      NumNodesElem = FEM_Mesh.NumNodesElem[j];
      Element_Connectivity = ChainToArray(FEM_Mesh.Connectivity[j],NumNodesElem);
      /* 3º Loop over the all the node in an element */
      for(int k = 0 ; k<NumNodesElem ; k++){
	/* 4º If my node belong to the element */
	if(Element_Connectivity[k] == i){
	  /* 5º Introduce the element in the chain */
	  PushNodeTop(&FEM_Mesh.NodeNeighbour[i], j);
	  /* 6º Update the counter */
	  FEM_Mesh.NumNeighbour[i] += 1;	  
	}
      }
      /* Free memory */
      free(Element_Connectivity);
    }      
  }
  
}

/*********************************************************************/

void GlobalSearchGaussPoints(GaussPoint MPM_Mesh, Mesh FEM_Mesh){

  /* Variables for the GP coordinates */
  Matrix X_GC_GP;
  X_GC_GP.N_rows = NumberDimensions;
  X_GC_GP.N_cols = 1;  
  X_GC_GP.n = NAN;
  Matrix X_EC_GP;
  X_EC_GP.N_rows = NumberDimensions;
  X_EC_GP.N_cols = 1;  
  X_EC_GP.n = NAN;

  /* Variables for the poligon */
  int NumVertex;
  int * Poligon_Connectivity;
  Matrix Poligon_Coordinates;
  ChainPtr ListNodes_I;
  Matrix lp; /* Just for GIMP -> Particle voxel */

  /* 1º Set to zero the active/non-active elements */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.ActiveNode[i] = 0;
  }

  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Assign the value to this auxiliar pointer */ 
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){

      /* 3º Connectivity of the Poligon */
      NumVertex = FEM_Mesh.NumNodesElem[j];
      Poligon_Connectivity =
	ChainToArray(FEM_Mesh.Connectivity[j],NumVertex);
     
      /* 4º Allocate the polligon Matrix and fill it */
      Poligon_Coordinates = MatAllocZ(NumVertex,NumberDimensions);
      for(int k = 0; k<NumVertex; k++){
	for(int l = 0 ; l<NumberDimensions ; l++){
	  Poligon_Coordinates.nM[k][l] =
	    FEM_Mesh.Coordinates.nM[Poligon_Connectivity[k]][l];
	}
      }

      /* 5º Free memory */
      free(Poligon_Connectivity);
      
      /* 6º Check out if the GP is in the Element */
      if(InOut_Poligon(X_GC_GP,Poligon_Coordinates) == 1){

	/* 7º Asign to the GP a element in the background mesh, just for 
	   searching porpuses */
	MPM_Mesh.Element_id[i] = j;
	
	/* 8º If the GP is in the element, get its natural coordinates */
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);

	/* 9º Assign the new connectivity of the GP */
	if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
	  MPM_Mesh.ListNodes[i] = CopyChain(FEM_Mesh.Connectivity[j]);
	}
	else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
	  lp.nV = MPM_Mesh.Phi.lp.nM[i];
	  MPM_Mesh.ListNodes[i] =
	    Tributary_Nodes_GIMP(X_EC_GP,MPM_Mesh.Element_id[i],
	  			 lp,FEM_Mesh);
	}
	  
	/* 10º Active those nodes that interact with the GP */
	ListNodes_I = MPM_Mesh.ListNodes[i];
	while(ListNodes_I != NULL){
	  FEM_Mesh.ActiveNode[ListNodes_I->I] += 1;
	  ListNodes_I = ListNodes_I->next; 
	}
	
      }
      
      /* 11º Free memory */
      FreeMat(Poligon_Coordinates);
      
    } 

  } 
  
}

/*********************************************************************/

void LocalSearchGaussPoints(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
/*
  Local search algorithm based on the velocity of the particle
*/
{

  /* Variables for the GP coordinates */
  Matrix X_GC_GP;
  X_GC_GP.N_rows = NumberDimensions;
  X_GC_GP.N_cols = 1;
  X_GC_GP.n = NAN;
  Matrix X_EC_GP;
  X_EC_GP.N_rows = NumberDimensions;
  X_EC_GP.N_cols = 1;
  X_EC_GP.n = NAN;

  /* Variables for the poligon description */
  Matrix Poligon_Coordinates;
  int * Poligon_Connectivity;
  int NumVertex;
  
  int Elem_i; /* Element of the GP i */
  Matrix V_GP; /* Velocity array */
  V_GP.N_rows = NumberDimensions;
  V_GP.N_cols = 1;
  V_GP.nM = NULL;
  V_GP.n = NAN;
  strcpy(V_GP.Info,"V_GP");
  double Search_Direction ; 
  Matrix V_GP_n;
  V_GP_n = MatAllocZ(1,NumberDimensions);
  int SearchVertex; /* Index to start the search */
  int * SearchList; /* Pointer to store the search list */
  ChainPtr ListNodes_I;
  Matrix lp;
  
  /* 1º Set to zero the active/non-active elements */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.ActiveNode[i] = 0;
  }

  /* 2º Loop over the GP */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 3º Get the global coordinate of the GP */
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    /* 4º Get the velocity vector of the GP and if the norm of the velocity 
     vector is zero, cicle */
    V_GP.nV = MPM_Mesh.Phi.vel.nM[i];
    if(Norm_Mat(V_GP,2) == 0){
      continue;
    }

    /* 4º Get the index of the initial element */
    Elem_i = MPM_Mesh.Element_id[i];

    /* 5º Connectivity of the Poligon */
    NumVertex = FEM_Mesh.NumNodesElem[Elem_i];
    
    Poligon_Connectivity =
      ChainToArray(FEM_Mesh.Connectivity[Elem_i],NumVertex);
      
    /* 6º Allocate the poligon Matrix and fill it */
    Poligon_Coordinates = MatAllocZ(NumVertex,NumberDimensions);
    for(int k = 0; k<NumVertex; k++){
      for(int l = 0 ; l<NumberDimensions ; l++){
	Poligon_Coordinates.nM[k][l] =
	  FEM_Mesh.Coordinates.nM[Poligon_Connectivity[k]][l];
      }
    }

    /* 7º Check if the GP is in the same element */
    if(InOut_Poligon(X_GC_GP,Poligon_Coordinates) == 1){
      
      /* If the GP is in the element, get its natural coordinates */
      X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];     
      Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);
      
      /* Active those nodes that interact with the GP */
      for(int k = 0 ; k<NumVertex ; k++){
	FEM_Mesh.ActiveNode[Poligon_Connectivity[k]] += 1;
      }

      /* Free memory */
      free(Poligon_Connectivity);
      FreeMat(Poligon_Coordinates);
      
    }
    /* 7º If the GP is not in the same element, search in the neighbour */
    else{
    
      /* 7aº As we are in a new element we set to zero the element coordinates
       and reset the chain with the nodal connectivity of the GP */
      for(int j = 0 ; j<NumberDimensions ; j++){
	MPM_Mesh.Phi.x_EC.nM[i][j] = 0;
      }
      FreeChain(MPM_Mesh.ListNodes[i]);
      
      /* 7bº Set to a NAN the SearchVertex in order to avoid bugs */
      SearchVertex = -999;

      /* 7dº Get the search direction for the vertex of the element
	 and check the search direction */
      for(int j = 0 ; j<NumVertex ; j++){
	/* First vertex */
	if(j == 0){ 
	  for(int k = 0 ; k<NumberDimensions ; k++){
	    Search_Direction =
	      2*Poligon_Coordinates.nM[0][k] -
	      Poligon_Coordinates.nM[1][k] -
	      Poligon_Coordinates.nM[NumVertex-1][k];
	    /* Get the velocity proyection in the search direction */
	    V_GP_n.nV[k] = V_GP.nV[k]*Search_Direction;
	  }
	}
	/* Last vertex */
	else if(j == (NumVertex - 1)){ 
	  for(int k = 0 ; k<NumberDimensions ; k++){
	    Search_Direction =
	      2*Poligon_Coordinates.nM[NumVertex-1][k] -
	      Poligon_Coordinates.nM[0][k] -
	      Poligon_Coordinates.nM[NumVertex-2][k];
	    /* Get the velocity proyection in the search direction */
	    V_GP_n.nV[k] = V_GP.nV[k]*Search_Direction;
	  }
	}
	/* The rest of the elements */
	else{ 
	  for(int k = 0 ; k<NumberDimensions ; k++){
	    Search_Direction =
	      2*Poligon_Coordinates.nM[j][k] -
	      Poligon_Coordinates.nM[j+1][k] -
	      Poligon_Coordinates.nM[j-1][k];
	    /* Get the velocity proyection in the search direction */
	    V_GP_n.nV[k] = V_GP.nV[k]*Search_Direction;
	  }
	}
	/* Check the projection of the velocity vector */
	if( (V_GP_n.nV[0] >= 0) && (V_GP_n.nV[1] >= 0)){
	  SearchVertex = Poligon_Connectivity[j];
	  break;
	}	  
      }


      /* Free memory */
      FreeMat(Poligon_Coordinates);
      free(Poligon_Connectivity);
     
      /* 7eº Check for errors */
      if (SearchVertex<0 || SearchVertex>=FEM_Mesh.NumNodesMesh){
	puts("Error in LocalSearchGaussPoints() : Search algorithm fails !!! ");
	exit(0);
      }

      /* 7fº Create the search list of this vertex */
      SearchList = ChainToArray(FEM_Mesh.NodeNeighbour[SearchVertex],
				FEM_Mesh.NumNeighbour[SearchVertex]);

      /* 7gº Search in the search list */
      for(int j = 0 ; j<FEM_Mesh.NumNeighbour[SearchVertex] ; j++){
	
	/* Discard the initial element for the search */
	if(SearchList[j] == Elem_i) continue;

	/* Connectivity of the Poligon */
	NumVertex = FEM_Mesh.NumNodesElem[SearchList[j]];
	Poligon_Connectivity = ChainToArray(FEM_Mesh.Connectivity[SearchList[j]],
					    NumVertex);
      
	/* Allocate the polligon Matrix and fill it */
	Poligon_Coordinates = MatAllocZ(NumVertex,NumberDimensions);
	for(int k = 0; k<NumVertex; k++){
	  for(int l = 0 ; l<NumberDimensions ; l++){
	    Poligon_Coordinates.nM[k][l] =
	      FEM_Mesh.Coordinates.nM[Poligon_Connectivity[k]][l];
	  }
	}


	/* Free poligon connectivity */
	free(Poligon_Connectivity);

	if(InOut_Poligon(X_GC_GP,Poligon_Coordinates) == 1){

	  /* Asign to the GP a element in the background mesh, just for 
	     searching porpuses */
	  MPM_Mesh.Element_id[i] = SearchList[j];

	  /* Get its natural coordinates */
	  X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	  Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);

	  /* Free memory */
	  FreeMat(Poligon_Coordinates);

	  /* Assign the new connectivity of the GP */
	  if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
	    MPM_Mesh.ListNodes[i] =
	      CopyChain(FEM_Mesh.Connectivity[SearchList[j]]);
	  }
	  else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
	    lp.nV = MPM_Mesh.Phi.lp.nM[i];
	    MPM_Mesh.ListNodes[i] =
	      Tributary_Nodes_GIMP(X_EC_GP,MPM_Mesh.Element_id[i],
	    			   lp,FEM_Mesh);
	  }
	  
	  /* Active those nodes that interact with the GP */
	  ListNodes_I = MPM_Mesh.ListNodes[i];
	  while(ListNodes_I != NULL){
	    FEM_Mesh.ActiveNode[ListNodes_I->I] += 1;
	    ListNodes_I = ListNodes_I->next; 
	  }

	  /* If this is true, stop the search */
	  break;
	}
	else{
	  /* Free memory */
	  FreeMat(Poligon_Coordinates);
	}
	
      }

      /* Free memory */
      free(SearchList);

      if(MPM_Mesh.Element_id[i] == Elem_i){
	printf(" %s %i %s %i !!! \n",
	       "Error in LocalSearchGaussPoints() : GP",i,
	       "is not in the neighbours of",Elem_i);
	exit(0);
      }
      
    }
  } /* Loop over the GP */
  
  /* 8º Free memory */
  FreeMat(V_GP_n);

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

/* /\* Function to copy the chain A in to the chain B *\/ */
/* ChainPtr CopyChain(ChainPtr A) */
/* { */

/*   ChainPtr B = NULL; */
/*   ChainPtr iPtrA = A; */

/*   /\* Insert all elements of A to the result list *\/ */
/*   while (iPtrA != NULL){ */
/*     /\* Introduce a new element in the new chain *\/ */
/*     PushNodeTop(&B, iPtrA->I); */
/*     /\* Updtate the interator index *\/ */
/*       iPtrA = iPtrA->next;  */
/*   } */
  
/*   return B; */
/* } */

ChainPtr CopyChain(ChainPtr start1){
  
  ChainPtr start2=NULL;
  ChainPtr previous=NULL;

  while(start1!=NULL){
  
    ChainPtr temp = (ChainPtr) malloc (sizeof(Chain));
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
ChainPtr ChainUnion(ChainPtr A, ChainPtr B) 
{ 
    ChainPtr C = NULL; 
    ChainPtr iPtrA = A, iPtrB = B; 
  
    /* Insert all elements of A to the result list */
    while (iPtrA != NULL){
      /* Introduce a new element in the new chain */
      PushNodeTop(&C, iPtrA->I);
      /* Updtate the interator index */
      iPtrA = iPtrA->next; 
    }
    
    /* Insert those elements of B which are not  */
    /* present in result list */ 
    while (iPtrB != NULL){
      /* Introduce a new element in the new chain */      
      if (!IsPresentNode(C, iPtrB->I)){
	PushNodeTop(&C, iPtrB->I);
      }
      /* Updtate the interator index */
      iPtrB = iPtrB->next; 
    }

    return C; 
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
    while (A != NULL) 
    { 
        printf ("%d \n", A->I); 
        A = A->next; 
    } 
} 

/*********************************************************************/
