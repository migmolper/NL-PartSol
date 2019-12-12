#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "../GRAMS/grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))


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
  Matrix lp; /* Particle voxel (GIMP) */
  Matrix lambda_GP;   /* Lagrange multiplier for the GP */
  Matrix Delta_Xip; /* Distance from GP to the nodes (LME) */
  int NumNodes; /* Number of neibourghs (LME) */
  int * ListNodes; /* List of nodes (LME) */
  int GP_I; /* Iterator for the neibourghs (LME) */

  /* Usefull values */
  lambda_GP = MatAssign(NumberDimensions,1,NAN,NULL,NULL); 

  /* 1º Set to zero the active/non-active node */
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
	
	/* 9º Assign the new connectivity of the GP */
	if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
	  /* If the GP is in the element, get its natural coordinates */
	  X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	  Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);
	  /* Asign connectivity */
	  MPM_Mesh.ListNodes[i] = CopyChain(FEM_Mesh.Connectivity[j]);
	}
	else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
	  /* If the GP is in the element, get its natural coordinates */
	  X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	  Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);
	  /* Calculate connectivity */
	  lp.nV = MPM_Mesh.lp.nM[i];
	  MPM_Mesh.ListNodes[i] =
	    Tributary_Nodes_GIMP(X_EC_GP,MPM_Mesh.Element_id[i],
	  			 lp,FEM_Mesh);
	  /* Calculate number of nodes */
	  MPM_Mesh.NumberNodes[i] = LenghtChain(MPM_Mesh.ListNodes[i]);
	}
	else if(strcmp(MPM_Mesh.ShapeFunctionGP,"LME") == 0){
	  /* Calculate connectivity */
	  MPM_Mesh.ListNodes[i] =
	    LME_Tributary_Nodes(X_GC_GP,MPM_Mesh.Element_id[i],
				FEM_Mesh,MPM_Mesh.Gamma);
	  /* Calculate number of nodes */
	  MPM_Mesh.NumberNodes[i] = LenghtChain(MPM_Mesh.ListNodes[i]);

	  /* Generate nodal distance list */
	  NumNodes = MPM_Mesh.NumberNodes[i];
	  ListNodes = ChainToArray(MPM_Mesh.ListNodes[i],NumNodes);
	  Delta_Xip = MatAlloc(NumNodes,NumberDimensions);
	  for(int k = 0 ; k<NumNodes ; k++){
	    GP_I = ListNodes[k];
	    for(int l = 0 ; l<NumberDimensions ; l++){
	      Delta_Xip.nM[k][l] =
		X_GC_GP.nV[l]-
		FEM_Mesh.Coordinates.nM[GP_I][l];
	    }
	  }
	  free(ListNodes);
	  
	  /* Auxiliar lambda to update it */
	  lambda_GP.nV = MPM_Mesh.lambda.nM[i];
	  /* Calculate lagrange multipliers */
	  lambda_GP = LME_lambda(Delta_Xip, lambda_GP,
				 FEM_Mesh.DeltaX, MPM_Mesh.Gamma);
	  /* Free memory */
	  FreeMat(Delta_Xip);
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

  Matrix lp; /* Particle voxel (GIMP) */
  Matrix Delta_Xip; /* Distance from GP to the nodes (LME) */
  int NumNodes; /* Number of neibourghs (LME) */
  int * ListNodes; /* List of nodes (LME) */
  int GP_I; /* Iterator for the neibourghs (LME) */
  
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

      if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
	/* If the GP is in the element, get its natural coordinates */
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];     
	Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);
      }           
      /* Assign the new connectivity of the GP */
      if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
	/* If the GP is in the element, get its natural coordinates */
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];     
	Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);
	/* Calculate connectivity */
	lp.nV = MPM_Mesh.lp.nM[i];
	MPM_Mesh.ListNodes[i] =
	  Tributary_Nodes_GIMP(X_EC_GP,MPM_Mesh.Element_id[i],
			       lp,FEM_Mesh);
	/* Calculate number of nodes */
	MPM_Mesh.NumberNodes[i] = LenghtChain(MPM_Mesh.ListNodes[i]);
      }
      else if(strcmp(MPM_Mesh.ShapeFunctionGP,"LME") == 0){
	/* Calculate connectivity */
	MPM_Mesh.ListNodes[i] =
	  LME_Tributary_Nodes(X_GC_GP,MPM_Mesh.Element_id[i],
			      FEM_Mesh,MPM_Mesh.Gamma);
	/* Calculate number of nodes */
	MPM_Mesh.NumberNodes[i] = LenghtChain(MPM_Mesh.ListNodes[i]);
	/* Generate nodal distance list */
	NumNodes = MPM_Mesh.NumberNodes[i];
	ListNodes = ChainToArray(MPM_Mesh.ListNodes[i],NumNodes);
	Delta_Xip = MatAlloc(NumNodes,2);
	for(int k = 0 ; k<NumNodes ; k++){
	  GP_I = ListNodes[k];
	  for(int l = 0 ; l<NumberDimensions ; l++){
	    Delta_Xip.nM[k][l] =
	      MPM_Mesh.Phi.x_GC.nM[i][l]-
	      FEM_Mesh.Coordinates.nM[GP_I][l];
	  }
	}
	free(ListNodes);
	/* Calculate lagrange multipliers */	  
	MPM_Mesh.lambda =
	  LME_lambda(Delta_Xip, MPM_Mesh.lambda,
		     FEM_Mesh.DeltaX, MPM_Mesh.Gamma);
	/* Free memory */
	FreeMat(Delta_Xip);
      }
            
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

	  /* Assign the new connectivity of the GP */
	  if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
	    /* Get its natural coordinates */
	    X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	    Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);
	    /* Asign connectivity */
	    MPM_Mesh.ListNodes[i] =
	      CopyChain(FEM_Mesh.Connectivity[SearchList[j]]);
	  }
	  else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
	    /* Get its natural coordinates */
	    X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	    Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);
	    /* Calculate connectivity */
	    lp.nV = MPM_Mesh.lp.nM[i];
	    MPM_Mesh.ListNodes[i] =
	      Tributary_Nodes_GIMP(X_EC_GP,MPM_Mesh.Element_id[i],
	    			   lp,FEM_Mesh);
	    /* Calculate number of nodes */
	    MPM_Mesh.NumberNodes[i] = LenghtChain(MPM_Mesh.ListNodes[i]);
	  }
	  else if(strcmp(MPM_Mesh.ShapeFunctionGP,"LME") == 0){
	    /* Calculate connectivity */
	    MPM_Mesh.ListNodes[i] =
	      LME_Tributary_Nodes(X_GC_GP,MPM_Mesh.Element_id[i],
				  FEM_Mesh,MPM_Mesh.Gamma);
	    /* Calculate number of nodes */
	    MPM_Mesh.NumberNodes[i] = LenghtChain(MPM_Mesh.ListNodes[i]);
	    /* Generate nodal distance list */
	    NumNodes = MPM_Mesh.NumberNodes[i];
	    ListNodes = ChainToArray(MPM_Mesh.ListNodes[i],NumNodes);
	    Delta_Xip = MatAlloc(NumNodes,2);
	    for(int k = 0 ; k<NumNodes ; k++){
	      GP_I = ListNodes[k];
	      for(int l = 0 ; l<NumberDimensions ; l++){
		Delta_Xip.nM[k][l] =
		  MPM_Mesh.Phi.x_GC.nM[i][l]-
		  FEM_Mesh.Coordinates.nM[GP_I][l];
	      }
	    }
	    free(ListNodes);
	    /* Calculate lagrange multipliers */	  
	    MPM_Mesh.lambda =
	      LME_lambda(Delta_Xip, MPM_Mesh.lambda,
			 FEM_Mesh.DeltaX, MPM_Mesh.Gamma);
	    /* Free memory */
	    FreeMat(Delta_Xip);
	  }

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

ChainPtr * GP_Neighbours(GaussPoint MPM_Mesh, Mesh FEM_Mesh, double epsilon){

  int NumGP = MPM_Mesh.NumGP;
  int iElement;
  int * Element_id = MPM_Mesh.Element_id;
  int * ListElements;
  ChainPtr * ListNeighbours;
  ChainPtr SearchGP = RangeChain(0,NumGP-1); /* Not implemented */
  ChainPtr iSearchGP;
  Matrix X_GPi = MatAssign(1,NumberDimensions,NAN,NULL,NULL);
  Matrix X_GP  = MatAssign(1,NumberDimensions,NAN,NULL,NULL);
  
  /* List of GP near to each GP */
  ListNeighbours =
    (ChainPtr *)malloc(NumGP*sizeof(ChainPtr));
  if(ListNeighbours == NULL){
    printf("%s : %s \n",
	   "Error in GP_Neighbours",
	   "Memory error");
    exit(0);
  }

  /* Loop over the GP's and generate the list */
  for(int i = 0 ; i<NumGP ; i++){

    /* Get the coordenates of the search GP */
    X_GPi.nV = MPM_Mesh.Phi.x_GC.nM[i];

    /* Get the elements near to the GP */

    /* Set to NULL the list of each GP */
    ListNeighbours[i] = NULL;

    /* Search in the elements near to each GP */
    for(int k = 0 ; k<9 ; k++){

      /* Element index */
      iElement = ListElements[k];
      
      /* Search index */
      iSearchGP = SearchGP;
      
      while(iSearchGP != NULL){

	/* Search in the sourrounding elements */
	if(iElement == Element_id[iSearchGP->I]){
	  /* Check if is inside of the tolerance */
	  X_GP.nV = MPM_Mesh.Phi.x_GC.nM[iSearchGP->I];
	  if(PointDistance(X_GP,X_I) < epsilon){
	    PushNodeTop(&ListNeighbours[i],iSearchGP->I);
	  }
	}
      }
      
      /* Update search index */
      iSearchGP = iSearchGP->next;
    }    
    
  }


  return ListNeighbours;
  
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

