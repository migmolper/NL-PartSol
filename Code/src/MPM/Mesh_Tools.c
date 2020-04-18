#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

/*********************************************************************/

Matrix GetInitialGaussPointPosition(Mesh FEM_Mesh, int GPxElement)
/*
 * 
 */
{

  int Ndim = NumberDimensions;
  int NumElemMesh = FEM_Mesh.NumElemMesh;
  Matrix N_GP;
  Matrix X_p = MatAllocZ(GPxElement*NumElemMesh,Ndim);
  strcpy(X_p.Info,"Global Coordinates");
  Matrix Xi_p = MatAllocZ(GPxElement,Ndim);
  Matrix Xi_p_j;
  Element Element;
  int Node;

  switch(GPxElement){
  case 1 :
    if(strcmp(FEM_Mesh.TypeElem,"Quadrilateral") == 0){
      /* Centred GP */
      Xi_p.nV[0] = 0.0;
      Xi_p.nV[1] = 0.0;
      /* Evaluate the shape function */
      N_GP = Q4_N(Xi_p);
      /* Get the coordinate of the center */
      for(int i = 0 ; i<NumElemMesh ; i++){
	Element = get_Element(i,FEM_Mesh.Connectivity[i],FEM_Mesh.NumNodesElem[i]);
	for(int k = 0 ; k<4 ; k++){
	  Node = Element.Connectivity[k];
	  for(int l = 0 ; l<2 ; l++){
	    X_p.nM[i][l] += N_GP.nV[k]*FEM_Mesh.Coordinates.nM[Node][l];
	  }
	}
	free(Element.Connectivity);
      }
      /* Free auxiliar matrix with the coordinates */
      FreeMat(Xi_p);
      /* Free value of the shape function in the GP */
      FreeMat(N_GP);
    }
    else if(strcmp(FEM_Mesh.TypeElem,"Triangle") == 0){
      Xi_p.nV[0] = (double)1/3;
      Xi_p.nV[1] = (double)1/3;
      /* Evaluate the shape function */
      N_GP = T3(Xi_p);
      /* Get the coordinate of the center */
      for(int i = 0 ; i<NumElemMesh ; i++){
	Element = get_Element(i, FEM_Mesh.Connectivity[i],
			      FEM_Mesh.NumNodesElem[i]);
	for(int k = 0 ; k<3 ; k++){
	  Node = Element.Connectivity[k];
	  for(int l = 0 ; l<2 ; l++){
	    X_p.nM[i][l] += N_GP.nV[k]*FEM_Mesh.Coordinates.nM[Node][l];
	  }
	}
	free(Element.Connectivity);
      }
      /* Free auxiliar matrix with the coordinates */
      FreeMat(Xi_p);
      /* Free value of the shape function in the GP */
      FreeMat(N_GP);
    }
    break;
  case 4:
    if(strcmp(FEM_Mesh.TypeElem,"Quadrilateral") == 0){
      /* Centred GP */
      Xi_p.nM[0][0] =  0.5;
      Xi_p.nM[0][1] =  0.5;
      Xi_p.nM[1][0] =  0.5;
      Xi_p.nM[1][1] = -0.5;
      Xi_p.nM[2][0] = -0.5;
      Xi_p.nM[2][1] =  0.5;
      Xi_p.nM[3][0] = -0.5;
      Xi_p.nM[3][1] = -0.5;
      /* Xi_p.nM[0][0] = (double)1/pow(3,0.5); */
      /* Xi_p.nM[0][1] = (double)1/pow(3,0.5); */
      /* Xi_p.nM[1][0] = (double)1/pow(3,0.5); */
      /* Xi_p.nM[1][1] = (double)-1/pow(3,0.5); */
      /* Xi_p.nM[2][0] = (double)-1/pow(3,0.5); */
      /* Xi_p.nM[2][1] = (double)1/pow(3,0.5); */
      /* Xi_p.nM[3][0] = (double)-1/pow(3,0.5); */
      /* Xi_p.nM[3][1] = (double)-1/pow(3,0.5); */
      /* Get the coordinate of the center */
      for(int i = 0 ; i<NumElemMesh ; i++){
	Element = get_Element(i, FEM_Mesh.Connectivity[i],
			      FEM_Mesh.NumNodesElem[i]);
	for(int j = 0 ; j<GPxElement ; j++){
	  /* Evaluate the shape function in the GP position */
	  Xi_p_j.nV = Xi_p.nM[j]; 
	  N_GP = Q4_N(Xi_p_j);
	  for(int k = 0 ; k<4 ; k++){
	    /* Connectivity of each element */
	    Node = Element.Connectivity[k];
	    for(int l = 0 ; l<NumberDimensions ; l++){
	      X_p.nM[i*GPxElement+j][l] +=
		N_GP.nV[k]*FEM_Mesh.Coordinates.nM[Node][l];
	    }
	  }
	  /* Free value of the shape function in the GP */
	  FreeMat(N_GP);
	}
	free(Element.Connectivity);
      }
      /* Free auxiliar matrix with the coordinates */
      FreeMat(Xi_p);
    }
    break;
  default :
    fprintf(stderr,"%s : %s \n",
	    "Error in GetInitialGaussPointPosition()",
	    "Wrong number of gauss point per element");
    exit(1);
  }
  
  return X_p;
}

/*********************************************************************/

double GetMinElementSize(Mesh FEM_Mesh)
/*
  Function to get the minimum mesh size.
*/
{

  /* Auxiliar variables of the function */
  int NumElemMesh = FEM_Mesh.NumElemMesh;
  int NumNodesElem; /* Number of nodes of each element */
  int * Connectivity; /* Connectivity of the element */
  Matrix Poligon; /* Element Poligon */
  Matrix X_eval = MatAllocZ(1,2); /* Where to evaluate the shape function */
  X_eval.nV[0] = 0;
  X_eval.nV[1] = 0;
  Matrix dNdx; /* Gradient of the shapefunction for each node */
  double MinElementSize_aux;
  double MinElementSize = 10e16;

  /* 1º Loop over the elements in the mesh */
  for(int i = 0 ; i<NumElemMesh ; i++){

    /* 2º Connectivity of the Poligon */
    NumNodesElem = FEM_Mesh.NumNodesElem[i];
    Connectivity = Set_to_Pointer(FEM_Mesh.Connectivity[i],NumNodesElem);
    
    /* 4º Get the gradient of the element for each node */
    if((NumNodesElem == 3) &&
       (NumberDimensions == 2)){ /* Triangular element */
      /* The poligon is a triangle */
      Poligon = MatAllocZ(3,2);
      /* Fill the triangle */
      for(int k = 0; k<3; k++){
	for(int l = 0 ; l<2 ; l++){
	  Poligon.nM[k][l] = FEM_Mesh.Coordinates.nM[Connectivity[k]][l];
	}
      }
      /* Get the gradient of the triangle */
      dNdx = Get_dNdX_T3(X_eval,Poligon);
      FreeMat(Poligon);
      
      /* Get the minimum minimum height of the triangle */
      for(int j = 0 ; j<3 ; j++){
	MinElementSize_aux =
	  1/pow(dNdx.nM[0][j]*dNdx.nM[0][j] +
		dNdx.nM[1][j]*dNdx.nM[1][j],0.5);
	MinElementSize = MINVAL(MinElementSize,MinElementSize_aux);
      }
      /* Free memory */
      FreeMat(dNdx);
      
    }
    else if((NumNodesElem == 4) &&
	    (NumberDimensions == 2)){ /* Quadrilateral element */
      /* The poligon is a quadrilateral */
      Poligon = MatAllocZ(4,2);

      /* Fill the poligon with vectors */
      for(int k = 0; k<3; k++){
	for(int l = 0 ; l<2 ; l++){
	  Poligon.nM[k][l] =
	    FEM_Mesh.Coordinates.nM[Connectivity[k+1]][l] -
	    FEM_Mesh.Coordinates.nM[Connectivity[k]][l];
	}
      }
      for(int l = 0 ; l<2 ; l++){
	Poligon.nM[3][l] = FEM_Mesh.Coordinates.nM[Connectivity[0]][l] -
	  FEM_Mesh.Coordinates.nM[Connectivity[3]][l];
      }
      
      /* Get the minimum minimum height of the triangle */
      for(int k = 0 ; k<4 ; k++){
	MinElementSize_aux = pow(Poligon.nM[k][0]*Poligon.nM[k][0] +
				 Poligon.nM[k][1]*Poligon.nM[k][1] , 0.5);
	MinElementSize = MINVAL(MinElementSize,MinElementSize_aux);
      }

      /* Free memory */
      FreeMat(Poligon);

    }
    else{
      printf("%s : %s %i %s \n",
	     "Error in GetMinElementSize",
	     "Element with ",
	     NumNodesElem,
	     "nodes is not implemented !!!" );
      exit(0); 
    }

    /* Free memory */
    free(Connectivity);
    
  }

  /* Free memory */
  FreeMat(X_eval);

  return MinElementSize;

}

/*********************************************************************/

Matrix ElemCoordinates(Mesh FEM_Mesh, ChainPtr Element_p)
/*
  Get the matrix with the coordinates of an element
*/
{

  int Ndim = NumberDimensions;
  int NumVertex = get_Lenght_Set(Element_p);
  Matrix Coordinates = MatAllocZ(NumVertex,Ndim);
  ChainPtr Idx = NULL;
  int I_Idx = 0;

  Idx = Element_p;

  while(Idx != NULL){

    /* Fill the elelemtn coordinates */
    for(int l = 0 ; l<Ndim ; l++){
      Coordinates.nM[I_Idx][l] = FEM_Mesh.Coordinates.nM[Idx->I][l];
    }
    
    /* Cycle */
    Ixd = Ixd->next;
    I_Idx++;
  }
  
  return Coordinates;
}

/*********************************************************************/

int search_particle_in(int p, Matrix X_p, ChainPtr ListElement, Mesh FEM_Mesh)
/*

*/
{
  ChainPtr Ixd = NULL;
  Element Element_i;
  int I_element = -999;
  int Nn; /* Numver of nodes of the element */
  ChainPtr Nodes;

  Ixd = ListElement;
  
  while(Ixd!=NULL){

    Nn = FEM_Mesh.NumNodesElem.nV[Ixd->I];
    Nodes = FEM_Mesh.Connectivity[Ixd->I];
    Element_i = get_Element(p, Nodes, Nn);

    /* Check if the particle is in the element */
    if(InOut_Element(X_p, Element_i, FEM_Mesh.Coordinates)){
      I_element = Ixd->I;
      break;
    }
    
    /* Cycle */
    Ixd = Ixd->next;

  }

  if(I_element == -999){
    fprintf(stderr,"%s : %s %i \n",
	    "Error in search_particle_in()",
	    "Not posible to find the particle",p);
    exit(EXIT_FAILURE);
  }

  
  return I_element;
}

/*********************************************************************/

void get_particle_tributary_nodes(GaussPoint MPM_Mesh, Mesh FEM_Mesh, int p){

  int Ndim = NumberDimensions;

  /* Indetify the node close to the particle */
  int I0 = MPM_Mesh.I0[p];

  /* Define auxiliar matrix for local/global coordinates */
  Matrix X_p = MatAssign(Ndim,1,NAN,MPM_Mesh.Phi.x_GC.nM[p],NULL);
  Matrix Xi_p = MatAssign(Ndim,1,NAN,MPM_Mesh.Phi.x_EC.nM[i],NULL);
  /* Number of nodes near the particle */
  int Num_Elements_Near_I0 = FEM_Mesh.NumNeighbour[I0];
  /* Lis of elements near the particle */
  ChainPtr Elements_Near_I0 = FEM_Mesh.NodeNeighbour[I0];
  /* Index of the element (Q4/uGIMP) */
  int IdxElement;
  /* Coordinates of the nodes of the Element */
  Matrix CoordElement;

  
  /* Free previous list of tributary nodes to the particle */
  free_Set(MPM_Mesh.ListNodes[p]);
  
  /* 6º Assign the new connectivity of the GP */
  if(strcmp(ShapeFunctionGP,"MPMQ4") == 0){

    /* Get the index of the element */
    IdxElement = search_particle_in_(p,X_p,Elements_Near_I0,FEM_Mesh);
    /* Asign connectivity */
    MPM_Mesh.ListNodes[p] = CopyChain(FEM_Mesh.Connectivity[IdxElement]);
    /* Get the coordinates of the element vertex */
    CoordElement = ElemCoordinates(FEM_Mesh,MPM_Mesh.ListNodes[p]);
    /* Compute local coordinates of the particle in this element */
    Q4_X_to_Xi(Xi_p,X_p,CoordElement);
    /* Free coordinates of the element */
    FreeMat(CoordElement);
    
  }
  else if(strcmp(ShapeFunctionGP,"uGIMP") == 0){

    /* Auxiliar variables for GIMP */
    Matrix lp = MatAssign(Ndim,1,NAN,MPM_Mesh.lp.nM[p],NULL);
    /* Get the index of the element */
    IdxElement = search_particle_in_(p,X_p,Elements_Near_I0,FEM_Mesh);
    /* Get the coordinates of the element vertex */
    CoordElement = ElemCoordinates(FEM_Mesh,MPM_Mesh.ListNodes[p]);
    /* Compute local coordinates of the particle in this element */
    Q4_X_to_Xi(Xi_p,X_p,CoordElement);
    /* Asign connectivity */
    MPM_Mesh.ListNodes[p] = uGIMP_Tributary_Nodes(Xi_p,IdxElement,lp,FEM_Mesh);  
    /* Calculate number of nodes */
    MPM_Mesh.NumberNodes[p] = get_Lenght_Set(MPM_Mesh.ListNodes[p]);
    
  }
  else if(strcmp(ShapeFunctionGP,"LME") == 0){
    
    /* Auxiliar variables for LME */
    Matrix lambda_p = MatAssign(Ndim,1,NAN,MPM_Mesh.lambda.nM[p],NULL);
    Matrix Delta_Xip; /* Distance from particles to the nodes */
    Matrix Beta_p = MatAssign(Ndim,1,NAN,MPM_Mesh.Beta.nM[p],NULL);
	
    /* Calculate connectivity with the previous value of beta */
    MPM_Mesh.ListNodes[p] = LME_Tributary_Nodes(X_p,Beta_p,I0,FEM_Mesh);
    
    /* Calculate number of nodes */
    MPM_Mesh.NumberNodes[p] = get_Lenght_Set(MPM_Mesh.ListNodes[p]);

    /* Generate nodal distance list */
    Delta_Xip = get_set_Coordinates(MPM_Mesh.ListNodes[p],
				    X_p, FEM_Mesh.Coordinates);
    	      
    /* Update Beta and Lambda for each particle */
    Beta_p = LME_Beta(Beta_p, Delta_Xip, gamma_LME);
    lambda_p = LME_lambda_NR(Delta_Xip, lambda_p, Beta_p);
    
    /* Free memory */
    FreeMat(Delta_Xip);
  }

}

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
      Element_Connectivity = Set_to_Pointer(FEM_Mesh.Connectivity[j],NumNodesElem);
      /* 3º Loop over the all the node in an element */
      for(int k = 0 ; k<NumNodesElem ; k++){
	/* 4º If my node belong to the element */
	if(Element_Connectivity[k] == i){
	  /* 5º Introduce the element in the chain */
	  push_to_Set(&FEM_Mesh.NodeNeighbour[i], j);
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

ChainPtr DiscardElements(ChainPtr SearchElem_0,
			 Matrix Corner_MAX,
			 Matrix Corner_MIN,
			 Mesh FEM_Mesh)
/*
  Auxiliary function to increase the computational accuracy in the
  GlobalSearchGaussPoints() function. We discard the elements where we can't have
  GPs because its vertex coordinates. 
*/
{   
  /* Found the tail */
  if (SearchElem_0 == NULL) {
    return NULL;
  }

  ChainPtr SearchElem_1;
  double Xmax_E, Ymax_E;
  double Xmin_E, Ymin_E;
  int Idx_Elem; /* Index of the element */
  int NumNodes; /* Num nodes of the element */
  int * Conect_Elem; /* Conectivity of the element */

  Idx_Elem = SearchElem_0->I;
  NumNodes = FEM_Mesh.NumNodesElem[Idx_Elem];
  Conect_Elem = Set_to_Pointer(FEM_Mesh.Connectivity[Idx_Elem],NumNodes);

  /* Init Corner_MAX and Corner_MIN */
  Xmax_E = FEM_Mesh.Coordinates.nM[Conect_Elem[0]][0];
  Ymax_E = FEM_Mesh.Coordinates.nM[Conect_Elem[0]][1];
  Xmin_E = FEM_Mesh.Coordinates.nM[Conect_Elem[0]][0];
  Ymin_E = FEM_Mesh.Coordinates.nM[Conect_Elem[0]][1];
  
  for(int i = 1 ; i<NumNodes ; i++){
    /* Corner_MAX */
    Xmax_E = MAXVAL(Xmax_E, FEM_Mesh.Coordinates.nM[Conect_Elem[i]][0]);
    Ymax_E = MAXVAL(Ymax_E, FEM_Mesh.Coordinates.nM[Conect_Elem[i]][1]);
    /* Corner_MIN */
    Xmin_E = MINVAL(Xmin_E, FEM_Mesh.Coordinates.nM[Conect_Elem[i]][0]);
    Ymin_E = MINVAL(Ymin_E, FEM_Mesh.Coordinates.nM[Conect_Elem[i]][1]);
  }

  /* Free */
  free(Conect_Elem);
  
  /* Test inferior limit */
  if ((Xmin_E<Corner_MIN.nV[0]) &&
      (Ymin_E<Corner_MIN.nV[1])){ 
    SearchElem_1 = SearchElem_0->next;
    free(SearchElem_0);
    SearchElem_1->next =
      DiscardElements(SearchElem_1->next,
		      Corner_MAX,Corner_MIN,
		      FEM_Mesh);
    return SearchElem_1;
  }
  /* Test superior limit */
  else if ((Xmax_E>Corner_MAX.nV[0]) &&
	   (Ymax_E>Corner_MAX.nV[1])){
    SearchElem_1 = SearchElem_0->next;
    free(SearchElem_0);
    SearchElem_1->next =
      DiscardElements(SearchElem_1->next,
		      Corner_MAX,Corner_MIN,
		      FEM_Mesh);
    return SearchElem_1;
  }
  /* Just keep going */
  else{
    SearchElem_0->next =
      DiscardElements(SearchElem_0->next,
		      Corner_MAX,Corner_MIN,
		      FEM_Mesh);
    return SearchElem_0;
  }
}

/*********************************************************************/

void LocalSearchGaussPoints(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
/*
  Local search algorithm based on the velocity of the particle.
*/
{

  /* Variables for the GP coordinates */
  int Ndim = NumberDimensions;
  Matrix X_p;
  Matrix V_p;

  /* Variables for the poligon description */
  Matrix Poligon_Coordinates;
  int * Poligon_Connectivity;
  int NumVertex;
  
  int Elem_i; /* Element of the GP i */

  int SearchVertex; /* Index to start the search */
  int * SearchList; /* Pointer to store the search list */
  ChainPtr ListNodes_I;


  /* */
  int I0_p;

  ChainPtr Neighbour_Nodes;

  /* Matrix with the distance from the particle to the nodes */
  Matrix X_Ip;
  
  /* Set to zero the active/non-active node, and the GPs in each element */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.NumParticles[i] = 0;
  }
  for(int i = 0 ; i<FEM_Mesh.NumElemMesh ; i++){
    free_Set(FEM_Mesh.I_particles[i]);
    FEM_Mesh.I_particles[i] = NULL;
  }

  /* Loop over the particles */
  for(int p = 0 ; p<MPM_Mesh.NumGP ; p++){

    /* Get the global coordinates and velocity of the particle */
    X_p = MatAssign(Ndim,1,NAN,MPM_Mesh.Phi.x_GC.nM[p],NULL);
    V_p = MatAssign(Ndim,1,NAN,MPM_Mesh.Phi.vel.nM[p],NULL);

    /* Check if the particle is static or is in movement */
    if(Norm_Mat(V_p,2) > 0){

      /* Get the index of the node close to the particle */
      I0_p = MPM_Mesh.I0[p];

      /* Get nodes close to the node I0_p */
      Neighbour_Nodes = get_NodesClose_toNode(I0_p, FEM_Mesh);

      /* Compute matrix with vectors from the particle to the nodes */
      X_Ip = get_set_Coordinates(Neighbour_Nodes,X_p,FEM_Mesh.Coordinates);

      /* Update the index of the node close to the particle */
      MPM_Mesh.I0[p] = get_node_near_to(X_Ip,);

      /* Update the tributary nodes of each particle */
      get_particle_tributary_nodes(MPM_Mesh,FEM_Mesh,p);

      /* Active those nodes that interact with the particle */
      ListNodes_p = MPM_Mesh.ListNodes[p];
      while(ListNodes_p != NULL){
	FEM_Mesh.NumParticles[ListNodes_p->I] += 1;
	push_to_Set(FEM_Mesh.I_particles[ListNodes_p->I],p);
	ListNodes_p = ListNodes_p->next; 
      }
      
    }
    else{
      /* Active those nodes that interact with the particle */
      ListNodes_p = MPM_Mesh.ListNodes[p];
      while(ListNodes_p != NULL){
	FEM_Mesh.NumParticles[ListNodes_p->I] += 1;
	push_to_Set(FEM_Mesh.I_particles[ListNodes_p->I],p);
	ListNodes_p = ListNodes_p->next; 
      }       
    }
}

/*********************************************************************/

void ComputeBeps(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
/*

*/
{

  int Ndim = NumberDimensions;
  
  /* Search radious */
  double epsilon;
  
  /* Number of GP and array with the element where it is each GP */
  int NumGP = MPM_Mesh.NumGP;
  int * I0 = MPM_Mesh.I0;
  int Mat_GP;
  
  /* List wit the nodes of the initial element */
  int Elem_GP; /* Element of the GP (search) */
  int NumNodesElem;
  int * NodesElem;  

  /* Chain with a list of elements near the initial element,
   and interator chain */
  ChainPtr * TableElements = NULL;
  ChainPtr ListElements;
  int NumElems;
  int * ArrayElements;

  /* Auxiliar chain */
  ChainPtr * TableGP = NULL;
  ChainPtr SearchGP;
  int * ArraySearchGP;
  int NumSearchGP;

  /* Distance */
  Matrix X0 = MatAssign(Ndim,1,NAN,NULL,NULL);
  Matrix X1 = MatAssign(Ndim,1,NAN,NULL,NULL);
  double Dist;

  /* Free the previous list and set to NULL */
  for(int i = 0 ; i<NumGP ; i++){
    free_Set(MPM_Mesh.Beps[i]);
    MPM_Mesh.Beps[i] = NULL;
  }

  /* Loop over the GP's and generate the list */
  for(int i = 0 ; i<NumGP ; i++){
    
    /* Get the search rad */
    Mat_GP = MPM_Mesh.MatIdx[i];
    epsilon = MPM_Mesh.Mat[Mat_GP].Ceps*FEM_Mesh.DeltaX;
        
    /* Number of nodes of the initial element and
       list of nodes */
    Elem_GP = I0[i];
    NumNodesElem = FEM_Mesh.NumNodesElem[Elem_GP];
    NodesElem = Set_to_Pointer(FEM_Mesh.Connectivity[Elem_GP],
			     NumNodesElem);
    
    /* Chain with the tributary elements, this is the list of element near the
       gauss point, including where it is */
    TableElements = malloc(NumNodesElem*sizeof(ChainPtr));
    ListElements = NULL;
    for(int j = 0 ; j<NumNodesElem ; j++){
      TableElements[j] = FEM_Mesh.NodeNeighbour[NodesElem[j]];
    }
    ListElements = get_Union_Of(TableElements,NumNodesElem);
    /* Free memory */
    free(NodesElem);
    free(TableElements);
    TableElements = NULL;
    
    /* First search : In elements near to each GP */
    NumElems = get_Lenght_Set(ListElements);
    ArrayElements = Set_to_Pointer(ListElements, NumElems);
    TableGP = malloc(NumElems*sizeof(ChainPtr));
    SearchGP = NULL;
    for(int j = 0 ; j<NumElems ; j++){
      TableGP[j] = FEM_Mesh.I_particles[ArrayElements[j]];
    }
    SearchGP = get_Union_Of(TableGP,NumElems);
    
    /* Free memory */
    free_Set(ListElements);
    ListElements = NULL;
    free(ArrayElements);
    free(TableGP);
    TableElements = NULL;
    
    /* Search in the cell and return the total search list modified  */
    /* without those GP inside of the cell and */
    /* include in Beps[i] the GP inside of the cell */
    NumSearchGP = get_Lenght_Set(SearchGP);
    ArraySearchGP = Set_to_Pointer(SearchGP,NumSearchGP);
    X0.nV = MPM_Mesh.Phi.x_GC.nM[i];
    for(int j = 0 ; j<NumSearchGP ; j++){
      X1.nV = MPM_Mesh.Phi.x_GC.nM[ArraySearchGP[j]];
      Dist = Distance(X1,X0);
       if (Dist < epsilon){
	 push_to_Set(&MPM_Mesh.Beps[i],ArraySearchGP[j]);
       }
    }
    /* Free chain */
    free_Set(SearchGP);
    SearchGP = NULL;
    free(ArraySearchGP);
  }
}

/*********************************************************************/

void GPinCell(ChainPtr * ListInCELL,
	      ChainPtr * GlobalList,
	      Matrix x_GC, int iGP,
	      double epsilon)
/*
  Search recursively the GP inside of a cell and if inside of a circle of
  radious epsilon  
*/
{

  int Ndim = NumberDimensions;
  ChainPtr iPtr = (* GlobalList);
  ChainPtr AuxPtr;
  Matrix X0, X1;
  double Dist;
  if(iPtr != NULL){
    X0 = MatAssign(Ndim,1,NAN,x_GC.nM[iGP],NULL);
    X1 = MatAssign(Ndim,1,NAN,x_GC.nM[iPtr->I],NULL);
    Dist = Distance(X1,X0);

    /* Found one to delete */
    if (Dist < epsilon){
      push_to_Set (ListInCELL,iPtr->I);
      AuxPtr = iPtr;
      iPtr = iPtr->next;
      free(AuxPtr);
      GPinCell(ListInCELL,&iPtr,
	       x_GC,iGP,epsilon);
    }
    /* Just keep going */
    else{
      iPtr = iPtr->next;
      GPinCell(ListInCELL,&iPtr,
	       x_GC,iGP,epsilon);
    }
  }
}

/*********************************************************************/

Element get_Element(int i_GP, ChainPtr ListNodes, int NumNodes){

  /* Define new element */
  Element GP_Element;

  /* Fill element */
  GP_Element.i_GP = i_GP;
  GP_Element.NumberNodes = NumNodes;
  GP_Element.Connectivity = Set_to_Pointer(ListNodes,NumNodes);

  return GP_Element;
}

/*********************************************************************/

bool InOut_Element(Matrix X_p, Element Elem_p, Matrix Coordinates){

  bool Is_In_Element = false;
  int * Element_Connectivity;
  Matrix Element_Coordinates;


  free(Element_Connectivity);
  FreeMat(Element_Coordinates);

  return Is_In_Element;
}

/*********************************************************************/

Matrix get_Element_Field(Matrix Nodal_Field, Element GP_Element){

  int * Element_Connectivity = GP_Element.Connectivity;
  int NumNodes = GP_Element.NumberNodes;
  Matrix Element_Field;
  int Ndim = NumberDimensions;
  int Ie;

  /* Allocate a matrix to store the nodal quatities in the element */
  Element_Field = MatAlloc(NumNodes,Ndim);

  /* Loop over the nodes of the element */
  for(int I = 0; I<NumNodes; I++){
    /* Get the node of the element */
    Ie = Element_Connectivity[I];
    /* Fill each dimension of the nodal quantitie */
    for(int i = 0 ; i<Ndim ; i++){
      Element_Field.nM[I][i] = Nodal_Field.nM[Ie][i];
    }
  }
  
  return Element_Field;
}

/*********************************************************************/

ChainPtr get_NodesClose_toNode(int I, Mesh FEM_Mesh){

  /* Define output */
  ChainPtr Nodes = NULL;

  /* Number of elements close to to the node */
  int NumNeighbour = FEM_Mesh.NumNeighbour[I];
  /* Index of the elements close to the node */
  int * NodeNeighbour = Set_to_Pointer(FEM_Mesh.NodeNeighbour[I],NumNeighbour);
  /* Table with the nodes of each element */
  ChainPtr * Table_ElemNodes = malloc(NumNeighbour*sizeof(ChainPtr));  
  for(int i = 0 ; i<NumNeighbour ; i++){
    Table_ElemNodes[i] = FEM_Mesh.Connectivity[NodeNeighbour[i]];
  }
  /* Free table with elements */
  free(NodeNeighbour);
  /* Get the union of this nodes */
  Nodes = get_Union_Of(Table_ElemNodes,NumNeighbour);
  /* Free table with the nodes of each elements */
  free(Table_ElemNodes);

  /* Return nodes close to the node I */
  return Nodes;
}

/*********************************************************************/

int get_node_near_to(Matrix Vector_to_Point){

  int Nnodes = Vector_to_Point.N_rows;
  ChainPtr List_Ord = NULL;
  ChainPtr List_Dis = RangeChain(0,Nnodes-1);
  Matrix Distance_to_Point;
  int Node;

  /* Get the modulus of the distance for each vector */
  Distance_to_Point = get_nurbs_distance(Vector_to_Point);
  
  /* Ordenate distances  */  
  order_Set(&List_Ord,&List_Dis,Distance_Xip);

  FreeMat(Distance_to_Point);

  /* Transform the list in to an array */
  List = Set_to_Pointer(List_Ord,NumNodes_GP);

  /* Free set in order */
  free_Set(List_Ord);

  /* Get the node near to the point */
  Node = List[0];

  /* Free the list */
  free(List);

  return Node;
}

/*********************************************************************/
