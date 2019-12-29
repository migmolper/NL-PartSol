#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "../GRAMS/grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

/*********************************************************************/

Boundaries GetBoundaryBox(Mesh BackgroundMesh){ 

  /* Auxiliary variable for BCCs */
  Boundaries BCC;
  /* Boundaries labels */
  char * BoundLabels [4] = {"top", "bottom", "left", "right"};

  /* Auxiliar to define the boundaries */
  int * NodesBound_aux;
  int * CounterNodesBound;
  int NumNodesBound;

  /* Count the number of elements that share this node */
  ChainPtr Elem_Conn; /* Loop over the connectivity chain */
  int Repeat_Nod; /* Counter */
  /* Variables that fills the boundaries nodes */
  int aux_RIGHT, aux_TOP, aux_LEFT, aux_BOTTOM;
  aux_RIGHT = 0;
  aux_TOP = 0;
  aux_LEFT = 0;
  aux_BOTTOM = 0;
  /* Define max and min value */
  double MAX_X, MAX_Y, MIN_X, MIN_Y;


  /* Kind of domain for the boundaries */
  strcpy(BCC.Info,"SQUARE");      
  
  /* Asign the number of boundaries for a square domain */
  BCC.NumBounds = 4;

  /* Generate boundaries for the domain */
  BCC.BCC_i = (Load *)Allocate_Array(BCC.NumBounds,sizeof(Load));

  /* Fill some parameters */
  for(int i = 0 ; i<BCC.NumBounds  ; i++){
    /* Set to zero the number of nodes in the boundary of the mesh */
    BCC.BCC_i[i].NumNodes = 0;
    /* Set the labels for a square domain */
    strcpy(BCC.BCC_i[i].Info,BoundLabels[i]);
  }

  /* Allocate an array of zeros to assign a 1
     to those nodes in the boundary */
  NodesBound_aux =
    (int *)Allocate_ArrayZ(BackgroundMesh.NumNodesMesh,sizeof(int));
  CounterNodesBound =
    (int *)Allocate_ArrayZ(4,sizeof(int));
  /* Set to zero the number of nodes in the boundary */
  NumNodesBound = 0;

  /* Set intial value */
  MAX_X = BackgroundMesh.Coordinates.nM[0][0];
  MAX_Y = BackgroundMesh.Coordinates.nM[0][1];
  MIN_X = BackgroundMesh.Coordinates.nM[0][0];
  MIN_Y = BackgroundMesh.Coordinates.nM[0][1];
	
  /* Iterate over the nodes to fin the nodes in the boundary */
  for(int i = 0 ; i<BackgroundMesh.NumNodesMesh ; i++){
	
    /* 4º Set the counter to zero */
    Repeat_Nod = 0;
	
    /* 5º Loop over the connectivity mesh to find the 
       repeted element */
    for(int j = 0 ; j<BackgroundMesh.NumElemMesh ; j++){
      Elem_Conn = BackgroundMesh.Connectivity[j];
      while(Elem_Conn != NULL){
	if((Elem_Conn->I) == i){
	  Repeat_Nod++;
	}
	Elem_Conn = Elem_Conn->next;
      }
    }
	
    /* Add this element to the boundary */
    if (Repeat_Nod < 4){
      NodesBound_aux[i] = 1;
      NumNodesBound++;
    }

    /* Get the max values of the boundary */
    MAX_X = MAXVAL(MAX_X,BackgroundMesh.Coordinates.nM[i][0]);
    MAX_Y = MAXVAL(MAX_Y,BackgroundMesh.Coordinates.nM[i][1]);
    MIN_X = MINVAL(MIN_X,BackgroundMesh.Coordinates.nM[i][0]);
    MIN_Y = MINVAL(MIN_Y,BackgroundMesh.Coordinates.nM[i][1]);
  }
      
  /* Count the number of nodes in each boundarie */
  for(int i = 0 ; i<BackgroundMesh.NumNodesMesh ; i++){
    if(NodesBound_aux[i] == 1){

      if(BackgroundMesh.Coordinates.nM[i][1] == MIN_Y){
	/* Count number of nodes in the bottom */
	CounterNodesBound[0]++;
      }
      if(BackgroundMesh.Coordinates.nM[i][0] == MAX_X){
	/* Count number of nodes in the right */
	CounterNodesBound[1]++;
      }
      if(BackgroundMesh.Coordinates.nM[i][1] == MAX_Y){
	/* Count the number of nodes in the top */
	CounterNodesBound[2]++;
      }
      if(BackgroundMesh.Coordinates.nM[i][0] == MIN_X){
	/* Count the number of nodes in the left */
	CounterNodesBound[3]++;
      }
    }
  }
      
  /* Allocate the arrays with the boundary nodes */
  for(int i = 0 ; i<BCC.NumBounds ; i++){
    BCC.BCC_i[i].NumNodes = CounterNodesBound[i];
    BCC.BCC_i[i].Nodes =
      (int *)Allocate_ArrayZ(CounterNodesBound[i],sizeof(int));
    BCC.BCC_i[i].Value =
      (Curve *)Allocate_Array(NumberDOF,sizeof(Curve));
  }

  /* Free data */
  free(CounterNodesBound);
      
  /* Fill the arrays with the boundary nodes */
  for(int i = 0 ; i<BackgroundMesh.NumNodesMesh ; i++){
    if(NodesBound_aux[i] == 1){

      if(BackgroundMesh.Coordinates.nM[i][1] == MIN_Y){
	BCC.BCC_i[0].Nodes[aux_BOTTOM] = i;
	aux_BOTTOM++;
      }
      if(BackgroundMesh.Coordinates.nM[i][0] == MAX_X){
	BCC.BCC_i[1].Nodes[aux_RIGHT] = i;
	aux_RIGHT++;
      }
      if(BackgroundMesh.Coordinates.nM[i][1] == MAX_Y){
	BCC.BCC_i[2].Nodes[aux_TOP] = i;
	aux_TOP++;
      }
      if(BackgroundMesh.Coordinates.nM[i][0] == MIN_X){
	BCC.BCC_i[3].Nodes[aux_LEFT] = i;
	aux_LEFT++;
      }
    }
  }
  
  return BCC;
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
    Connectivity =
      ChainToArray(FEM_Mesh.Connectivity[i],NumNodesElem);
    
    /* 4º Get the gradient of the element for each node */
    if((NumNodesElem == 3) &&
       (NumberDimensions == 2)){ /* Triangular element */
      /* The poligon is a triangle */
      Poligon = MatAllocZ(3,2);
      /* Fill the triangle */
      for(int k = 0; k<3; k++){
	for(int l = 0 ; l<2 ; l++){
	  Poligon.nM[k][l] =
	    FEM_Mesh.Coordinates.nM[Connectivity[k]][l];
	}
      }
      /* Get the gradient of the triangle */
      dNdx = Get_dNdX_T3(X_eval,Poligon);
      FreeMat(Poligon);
      
      /* Get the minimum minimum height of the triangle */
      for(int j = 0 ; j<3 ; j++){
	MinElementSize_aux =
	  1/pow(dNdx.nM[0][j]*dNdx.nM[0][j] + dNdx.nM[1][j]*dNdx.nM[1][j],0.5);
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

Matrix ElemCoordinates(Mesh FEM_Mesh, int * Connectivity, int NumVertex)
/*
  Get the matrix with the coordinates of an element
*/
{

  Matrix Coordinates;
   
  /* 2º Allocate the polligon Matrix and fill it */
  Coordinates = MatAllocZ(NumVertex,NumberDimensions);
  for(int k = 0; k<NumVertex; k++){
    for(int l = 0 ; l<NumberDimensions ; l++){
      Coordinates.nM[k][l] =
	FEM_Mesh.Coordinates.nM[Connectivity[k]][l];
    }
  }
  
  return Coordinates;
}

/*********************************************************************/

void GetListNodesGP(GaussPoint MPM_Mesh, Mesh FEM_Mesh, int iGP){

  int IdxElement = MPM_Mesh.Element_id[iGP];
  FreeChain(MPM_Mesh.ListNodes[iGP]);
  MPM_Mesh.ListNodes[iGP] = NULL;
  
  /* 6º Assign the new connectivity of the GP */
  if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
    /* Asign connectivity */
    MPM_Mesh.ListNodes[iGP] = CopyChain(FEM_Mesh.Connectivity[IdxElement]);
  }
  else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){

    /* Auxiliar variables for GIMP */
    Matrix lp; /* Particle voxel */
    Matrix X_EC_GP = /* Element coordinates */
      MatAssign(NumberDimensions,1,NAN,MPM_Mesh.Phi.x_EC.nM[iGP],NULL);
    
    /* Calculate connectivity */
    lp.nV = MPM_Mesh.lp.nM[iGP];
    MPM_Mesh.ListNodes[iGP] =
      Tributary_Nodes_GIMP(X_EC_GP,IdxElement,
			   lp,FEM_Mesh);
    
    /* Calculate number of nodes */
    MPM_Mesh.NumberNodes[iGP] = LenghtChain(MPM_Mesh.ListNodes[iGP]);
    
  }
  else if(strcmp(MPM_Mesh.ShapeFunctionGP,"LME") == 0){
    
    /* Auxiliar variables for LME */
    Matrix X_GC_GP = /* Global coordinates */
      MatAssign(NumberDimensions,1,NAN,MPM_Mesh.Phi.x_GC.nM[iGP],NULL);
    Matrix lambda_GP = /* Lagrange multipliers */
      MatAssign(NumberDimensions,1,NAN,NULL,NULL);
    Matrix Delta_Xip; /* Distance from GP to the nodes */
    double Beta; /* Tunning parameter */
    int NumNodes; /* Number of neibourghs */
    int * ListNodes; /* List of nodes */
    int I_iGP; /* Iterator for the neibourghs */
  
    /* Calculate connectivity */
    MPM_Mesh.ListNodes[iGP] =
      LME_Tributary_Nodes(X_GC_GP,IdxElement,
			  FEM_Mesh,MPM_Mesh.Gamma);
    
    /* Calculate number of nodes */
    MPM_Mesh.NumberNodes[iGP] = LenghtChain(MPM_Mesh.ListNodes[iGP]);

    /* Generate nodal distance list */
    NumNodes = MPM_Mesh.NumberNodes[iGP];
    ListNodes = ChainToArray(MPM_Mesh.ListNodes[iGP],NumNodes);
    Delta_Xip = MatAlloc(NumNodes,NumberDimensions);
    for(int k = 0 ; k<NumNodes ; k++){
      I_iGP = ListNodes[k];
      for(int l = 0 ; l<NumberDimensions ; l++){
	Delta_Xip.nM[k][l] =
	  X_GC_GP.nV[l]-
	  FEM_Mesh.Coordinates.nM[I_iGP][l];
      }
    }
    free(ListNodes);
	  
    /* Auxiliar lambda to update it */
    lambda_GP.nV = MPM_Mesh.lambda.nM[iGP];
    
    /* Calculate lagrange multipliers with Newton-Rapson */
    Beta = MPM_Mesh.Gamma/(FEM_Mesh.DeltaX*FEM_Mesh.DeltaX);
    lambda_GP = LME_lambda_NR(Delta_Xip, lambda_GP, Beta);
    
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
  Matrix X_GC_GP = MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  Matrix X_EC_GP = MatAssign(NumberDimensions,1,NAN,NULL,NULL);

  /* Variables for the poligon */
  int NumVertex;
  int * Poligon_Connectivity;
  Matrix Poligon_Coordinates;
  ChainPtr ListNodes_I;


  /* 1º Set to zero the active/non-active node, and the GPs in each 
   element */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.ActiveNode[i] = 0;
  }
  
  for(int i = 0 ; i<FEM_Mesh.NumElemMesh ; i++){
    FreeChain(FEM_Mesh.GPsElements[i]);
    FEM_Mesh.GPsElements[i] = NULL;
  }

  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Assign the value to this auxiliar pointer */ 
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){

      /* 3º Connectivity of the Poligon */
      NumVertex = FEM_Mesh.NumNodesElem[j];
      Poligon_Connectivity =
	ChainToArray(FEM_Mesh.Connectivity[j],NumVertex);

      /* 4º Get the coordinates of the element */
      Poligon_Coordinates =
	ElemCoordinates(FEM_Mesh,Poligon_Connectivity,NumVertex);
      
      /* 5º Check out if the GP is in the Element */
      if(InOut_Poligon(X_GC_GP,Poligon_Coordinates) == 1){

	/* 6º Asign to the GP a element in the background mesh, just for 
	   searching porpuses */
	MPM_Mesh.Element_id[i] = j;
	PushNodeTop(&FEM_Mesh.GPsElements[j],i);

	/* 7º If the GP is in the element, get its natural coordinates */
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);

	/* 8º Get list of nodes near to the GP */
	GetListNodesGP(MPM_Mesh,FEM_Mesh,i);
	
	/* 9º Active those nodes that interact with the GP */
	ListNodes_I = MPM_Mesh.ListNodes[i];
	while(ListNodes_I != NULL){
	  FEM_Mesh.ActiveNode[ListNodes_I->I] += 1;
	  ListNodes_I = ListNodes_I->next; 
	}
	
      }
      
      /* 10º Free memory */
      free(Poligon_Connectivity);
      FreeMat(Poligon_Coordinates);
      
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
  Conect_Elem = ChainToArray(FEM_Mesh.Connectivity[Idx_Elem],NumNodes);

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
  Matrix X_GC_GP = MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  Matrix X_EC_GP = MatAssign(NumberDimensions,1,NAN,NULL,NULL);

  /* Variables for the poligon description */
  Matrix Poligon_Coordinates;
  int * Poligon_Connectivity;
  int NumVertex;
  
  int Elem_i; /* Element of the GP i */
  Matrix V_GP = /* Velocity array */
    MatAssign(NumberDimensions,1,NAN,NULL,NULL); 
  double Search_Direction ; 
  Matrix V_GP_n;
  V_GP_n = MatAllocZ(1,NumberDimensions);
  int SearchVertex; /* Index to start the search */
  int * SearchList; /* Pointer to store the search list */
  ChainPtr ListNodes_I;
  
  /* 1º Set to zero the active/non-active node, and the GPs in each 
     element */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.ActiveNode[i] = 0;
  }
  for(int i = 0 ; i<FEM_Mesh.NumElemMesh ; i++){
    FreeChain(FEM_Mesh.GPsElements[i]);
    FEM_Mesh.GPsElements[i] = NULL;
  }

  /* 2º Loop over the GP */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* Get the velocity vector of the GP */
   
    /* 3º Get the global coordinates and velocity of the GP */
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];
    V_GP.nV = MPM_Mesh.Phi.vel.nM[i];

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

      /* Asign GP to the element */
      PushNodeTop(&FEM_Mesh.GPsElements[Elem_i],i);

      /* If the norm of the velocity
	 vector is zero get its natural coordinates */
      if(Norm_Mat(V_GP,2) > 0){
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);
      }
	      
      /* Get list of nodes near to the GP */
      GetListNodesGP(MPM_Mesh,FEM_Mesh,i);
            
      /* Active those nodes that interact with the GP */
      ListNodes_I = MPM_Mesh.ListNodes[i];
      while(ListNodes_I != NULL){
	FEM_Mesh.ActiveNode[ListNodes_I->I] += 1;
	ListNodes_I = ListNodes_I->next; 
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
	  PushNodeTop(&FEM_Mesh.GPsElements[SearchList[j]],i);

	  /* If the GP is in the element, get its natural coordinates */
	  X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	  Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);
	  
	  /* Get list of nodes near to the GP */
	  GetListNodesGP(MPM_Mesh,FEM_Mesh,i);

	  /* Free memory */
	  FreeMat(Poligon_Coordinates);
	  
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

void UpdateBeps(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
/*

*/
{
  
  /* Search radious */
  double epsilon;
  
  /* Number of GP and array with the element where it is each GP */
  int NumGP = MPM_Mesh.NumGP;
  int * Element_id = MPM_Mesh.Element_id;
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
  Matrix X0 = MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  Matrix X1 = MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  double Dist;

  /* Free the previous list and set to NULL */
  for(int i = 0 ; i<NumGP ; i++){
    FreeChain(MPM_Mesh.Beps[i]);
    MPM_Mesh.Beps[i] = NULL;
  }

  /* Loop over the GP's and generate the list */
  for(int i = 0 ; i<NumGP ; i++){
    
    /* Get the search radious */
    Mat_GP = MPM_Mesh.MatIdx[i];
    epsilon = MPM_Mesh.Mat[Mat_GP].Ceps*FEM_Mesh.DeltaX;
        
    /* Number of nodes of the initial element and
       list of nodes */
    Elem_GP = Element_id[i];
    NumNodesElem = FEM_Mesh.NumNodesElem[Elem_GP];
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],
			     NumNodesElem);
    
    /* Chain with the tributary elements, this is the list of element near the
       gauss point, including where it is */
    TableElements = malloc(NumNodesElem*sizeof(ChainPtr));
    ListElements = NULL;
    for(int j = 0 ; j<NumNodesElem ; j++){
      TableElements[j] = FEM_Mesh.NodeNeighbour[NodesElem[j]];
    }
    ListElements = ChainUnion(TableElements,NumNodesElem);
    /* Free memory */
    free(NodesElem);
    free(TableElements);
    TableElements = NULL;
    
    /* First search : In elements near to each GP */
    NumElems = LenghtChain(ListElements);
    ArrayElements = ChainToArray(ListElements, NumElems);
    TableGP = malloc(NumElems*sizeof(ChainPtr));
    SearchGP = NULL;
    for(int j = 0 ; j<NumElems ; j++){
      TableGP[j] = FEM_Mesh.GPsElements[ArrayElements[j]];
    }
    SearchGP = ChainUnion(TableGP,NumElems);
    
    /* Free memory */
    FreeChain(ListElements);
    ListElements = NULL;
    free(ArrayElements);
    free(TableGP);
    TableElements = NULL;
    
    /* Search in the cell and return the total search list modified  */
    /* without those GP inside of the cell and */
    /* include in Beps[i] the GP inside of the cell */
    NumSearchGP = LenghtChain(SearchGP);
    ArraySearchGP = ChainToArray(SearchGP,NumSearchGP);
    X0.nV = MPM_Mesh.Phi.x_GC.nM[i];
    for(int j = 0 ; j<NumSearchGP ; j++){
      X1.nV = MPM_Mesh.Phi.x_GC.nM[ArraySearchGP[j]];
      Dist = Distance(X1,X0);
       if (Dist < epsilon){
	 PushNodeTop(&MPM_Mesh.Beps[i],ArraySearchGP[j]);
       }
    }
    /* Free chain */
    FreeChain(SearchGP);
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

  ChainPtr iPtr = (* GlobalList);
  ChainPtr AuxPtr;
  Matrix X0, X1;
  double Dist;
  if(iPtr != NULL){
    X0 = MatAssign(NumberDimensions,1,NAN,x_GC.nM[iGP],NULL);
    X1 = MatAssign(NumberDimensions,1,NAN,x_GC.nM[iPtr->I],NULL);
    Dist = Distance(X1,X0);

    /* Found one to delete */
    if (Dist < epsilon){
      PushNodeTop (ListInCELL,iPtr->I);
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
