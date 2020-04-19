#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "grams.h"

/**************************************************/
/* Uniform Geralized Interpolation Material Point */
/**************************************************/

/*
  Shape functions based in : 
  "" The Generalized Interpolation Material Point Method ""
  by S.G.Bardenhagen and E.M.Kober, 2004
*/

/*            ^           */
/*          __|__         */
/*        _/  |  \_       */
/*      _/    |    \_     */
/*   __/      |      \__  */
/*  --o-------o-------o-- */
/*   (-1)    (0)     (1)  */

/*********************************************************************/

void uGIMP_Initialize(GaussPoint MPM_Mesh, Mesh FEM_Mesh){

  int Ndim = NumberDimensions;

  /* Variables for the GP coordinates */  
  Matrix X_GC_GP = MatAssign(Ndim,1,NAN,NULL,NULL);
  Matrix X_EC_GP = MatAssign(Ndim,1,NAN,NULL,NULL);

  Matrix lp; /* Particle voxel */

  
  double V_p; /* Volumen of the GP */
  double A_p;
  double rho_p;
  double th_p;
  double m_p;
  double l_p;
  int Mat_p;

  /* Variables for the poligon */
  int NumVertex;
  int * Poligon_Connectivity;
  Matrix Poligon_Coordinates;
  ChainPtr ListNodes_I;

  /* 1º Set to zero the active/non-active node, and the GPs in each 
     element */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.NumParticles[i] = 0;
  }
  
  for(int i = 0 ; i<FEM_Mesh.NumElemMesh ; i++){
    free_Set(FEM_Mesh.I_particles[i]);
    FEM_Mesh.I_particles[i] = NULL;
  }

  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 1º Get the voxel lenght */
    Mat_p = MPM_Mesh.MatIdx[i];
    rho_p = MPM_Mesh.Phi.rho.nV[i];
    m_p = MPM_Mesh.Phi.mass.nV[i];
    V_p = m_p/rho_p;
    th_p = MPM_Mesh.Mat[Mat_p].thickness;
    A_p = V_p/th_p;
    
    if(Ndim == 2){
      l_p = 0.5*pow(A_p,0.5);
    }
    for(int j =0 ; j<Ndim ; j++){
      MPM_Mesh.lp.nM[i][j] = l_p;
      lp.nV[j] = l_p;
    }

    /* 2º Assign the value to this auxiliar pointer */ 
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){

      /* 3º Connectivity of the Poligon */
      NumVertex = FEM_Mesh.NumNodesElem[j];
      Poligon_Connectivity =
	Set_to_Pointer(FEM_Mesh.Connectivity[j],NumVertex);

      /* 4º Get the coordinates of the element */
      Poligon_Coordinates = ElemCoordinates(Poligon_Connectivity,
	FEM_Mesh.Coordinates);
      
      /* 5º Check out if the GP is in the Element */
      if(InOut_Poligon(X_GC_GP,Poligon_Coordinates) == 1){

	/* 6º Asign to the GP a element in the background mesh, just for 
	   searching porpuses */
	MPM_Mesh.I0[i] = j;
	push_to_Set(&FEM_Mesh.I_particles[j],i);

	/* 7º If the GP is in the element, get its natural coordinates */
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	Q4_X_to_Xi(X_EC_GP,X_GC_GP,Poligon_Coordinates);

	/* 8º Get list of nodes near to the GP */
	free_Set(MPM_Mesh.ListNodes[i]);
	MPM_Mesh.ListNodes[i] = NULL;
	MPM_Mesh.ListNodes[i] =
	  uGIMP_Tributary_Nodes(X_GC_GP,MPM_Mesh.I0[i],lp,FEM_Mesh);
	MPM_Mesh.NumberNodes[i] = get_Lenght_Set(MPM_Mesh.ListNodes[i]);
	/* 9º Active those nodes that interact with the GP */
	ListNodes_I = MPM_Mesh.ListNodes[i];
	while(ListNodes_I != NULL){
	  FEM_Mesh.NumParticles[ListNodes_I->I] += 1;
	  ListNodes_I = ListNodes_I->next; 
	}

	/* 10º Free memory and go for the next GP */
	free(Poligon_Connectivity);
	FreeMat(Poligon_Coordinates);
	break;
      }
      
      /* 11º Free memory */
      free(Poligon_Connectivity);
      FreeMat(Poligon_Coordinates);
      
    } 

  }
  
}

/*********************************************************************/

/* Uniform GIMP shape function */
double uGIMP_Sip(double L, double lp, double Delta_xp){

  /* Evaluation of the shape function */

  if ((-lp < Delta_xp) && (Delta_xp <= lp)){
    return 1 - 0.5*(DSQR(Delta_xp) + lp*lp)*(double)1/(L*lp);
  }
  else if (((-L-lp) < Delta_xp) && (Delta_xp <= (-L+lp))){
    return (double)(0.25/(L*lp))*DSQR(L+lp+Delta_xp);
  }
  else if (((L-lp) < Delta_xp) && (Delta_xp <= (L+lp))){
    return (double)(0.25/(L*lp))*DSQR(L+lp-Delta_xp);
  }
  else if (((-L+lp) < Delta_xp) && (Delta_xp <= -lp)){
    return 1 + (double)Delta_xp/L;
  }
  else if ((lp < Delta_xp) && (Delta_xp <= (L-lp))){
    return 1 - (double)Delta_xp/L;
  }
  else {
    return (double)0.0;
  }
}

/*********************************************************************/

/* Uniform GIMP derivative shape function */
double uGIMP_dSip(double L, double lp, double Delta_xp){

  /* Evaluation of the shape function */
  if (((-L-lp) < Delta_xp) && (Delta_xp <= (-L+lp))){
    return (double)(0.5/(L*lp))*(L+lp+Delta_xp);
  }
  else if (((-L+lp) < Delta_xp) && (Delta_xp <= -lp)){
    return (double)1/L;
  }
  else if ((-lp < Delta_xp) && (Delta_xp <= lp)){
    return -(double)Delta_xp/(L*lp);
  }
  else if ((lp < Delta_xp) && (Delta_xp <= (L-lp))){
    return -(double)1/L;
  }
  else if (((L-lp) < Delta_xp) && (Delta_xp <= (L+lp))){
    return -(double)(0.5/(L*lp))*(L+lp-Delta_xp);
  }
  else {
    return (double)0.0;
  }
}

/*********************************************************************/

/* Uniform GIMP shape function 2D */
Matrix uGIMP_N(Matrix Delta_Xp, Matrix lp, double L){

  /* 1º Variable declaration */
  int Ndim = NumberDimensions;
  int Nnodes = Delta_Xp.N_rows;
  Matrix S_Ip = MatAllocZ(1,Nnodes);

  /* 2º Fill the shape function array */
  switch(Ndim){
  case 2 :
    for(int i = 0 ; i<Nnodes ; i++){
      /* 3º Shape function in this node */
      S_Ip.nV[i] =
	uGIMP_Sip(L, lp.nV[0], Delta_Xp.nM[i][0])*
	uGIMP_Sip(L, lp.nV[1], Delta_Xp.nM[i][1]);
    }
    break;
  case 3:
    for(int i = 0 ; i<Nnodes ; i++){
      /* 3º Shape function in this node */
      S_Ip.nV[i] =
	uGIMP_Sip(L, lp.nV[0], Delta_Xp.nM[i][0])*
	uGIMP_Sip(L, lp.nV[1], Delta_Xp.nM[i][1])*
	uGIMP_Sip(L, lp.nV[2], Delta_Xp.nM[i][2]);
    }
    break;
  }

  /* 4º Output */
  return S_Ip;
}

/*********************************************************************/

/* Uniform GIMP derivative shape function 2D */
Matrix uGIMP_dN(Matrix Delta_xp, Matrix lp, double L){

  /* 1º Variable declaration */
  int Ndim = NumberDimensions;
  int Nnodes = Delta_xp.N_rows;
  Matrix dS_Ip = MatAllocZ(Nnodes,Ndim); 
  
  /* 2º Fill the shape function array */
  switch(Ndim){    
  case 2 :
    for(int i = 0 ; i<Nnodes ; i++){
      /* 3º Gradient of the shape function for each node*/
      dS_Ip.nM[i][0] =
	uGIMP_dSip(L, lp.nV[0], Delta_xp.nM[i][0])*
	uGIMP_Sip(L, lp.nV[1], Delta_xp.nM[i][1]);
      dS_Ip.nM[i][1] =
	uGIMP_Sip(L, lp.nV[0], Delta_xp.nM[i][0])*
	uGIMP_dSip(L, lp.nV[1], Delta_xp.nM[i][1]);
    }
    break;
  case 3 :
    for(int i = 0 ; i<Nnodes ; i++){
      /* 3º Gradient of the shape function for each node*/
      dS_Ip.nM[i][0] =
	uGIMP_dSip(L, lp.nV[0], Delta_xp.nM[i][0])*
	uGIMP_Sip(L, lp.nV[1], Delta_xp.nM[i][1])*
	uGIMP_Sip(L, lp.nV[2], Delta_xp.nM[i][2]);
      dS_Ip.nM[i][1] =
	uGIMP_Sip(L, lp.nV[0], Delta_xp.nM[i][0])*
	uGIMP_dSip(L, lp.nV[1], Delta_xp.nM[i][1])*
	uGIMP_Sip(L, lp.nV[2], Delta_xp.nM[i][2]);
      dS_Ip.nM[i][2] =
	uGIMP_Sip(L, lp.nV[0], Delta_xp.nM[i][0])*
	uGIMP_Sip(L, lp.nV[1], Delta_xp.nM[i][1])*
	uGIMP_dSip(L, lp.nV[2], Delta_xp.nM[i][2]);      
    }
    break;
  }
  /* 4º Output */
  return dS_Ip;
}

/*********************************************************************/

ChainPtr uGIMP_Tributary_Nodes(Matrix X_EC_GP,
			      int Elem_GP,Matrix lp,
			      Mesh FEM_Mesh){

  ChainPtr * Table_Nodes = NULL;
  ChainPtr Triburary_Nodes = NULL;
  ChainPtr ChainElements = NULL;
  int * Tributary_Elements;
  int * NodesElem;
  int Elem_i;
  int Num_Elem;
  double Dist[2];

  /* Get the reference distance measured from the center of the element */
  for(int i = 0 ; i<2; i++){
    Dist[i] = 1 - 2*lp.nV[i]/FEM_Mesh.DeltaX;
  }

  /* Check if I am in the central area */
  if ((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
      (fabs(X_EC_GP.nV[1]) <= Dist[1])){
    Triburary_Nodes = CopyChain(FEM_Mesh.Connectivity[Elem_GP]);
  }    
  /* Check if I am in the 1º Quadrant */
  else if((X_EC_GP.nV[0]>=0) &&
	  (X_EC_GP.nV[1]>=0)){
    /* Create an array with the nodes of the element */
    NodesElem = Set_to_Pointer(FEM_Mesh.Connectivity[Elem_GP],4);
    if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */       
      ChainElements = get_Intersection_Of(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					  FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = get_Lenght_Set(ChainElements);
      Tributary_Elements = Set_to_Pointer(ChainElements,Num_Elem);
      free_Set(ChainElements);

       /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = get_Intersection_Of(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					  FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Num_Elem = get_Lenght_Set(ChainElements);
      Tributary_Elements = Set_to_Pointer(ChainElements,Num_Elem);
      free_Set(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[2]];
      Num_Elem = get_Lenght_Set(ChainElements);
      Tributary_Elements = Set_to_Pointer(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    /* Free memory */
    free(NodesElem);    
  }  
  /* Check if I am in the 2º Quadrant */
  else if((X_EC_GP.nV[0]<=0) &&
	  (X_EC_GP.nV[1]>=0)){
    /* Create an array with the nodes of the element */
    NodesElem = Set_to_Pointer(FEM_Mesh.Connectivity[Elem_GP],4);

    if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) >= Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = get_Intersection_Of(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Num_Elem = get_Lenght_Set(ChainElements);      
      Tributary_Elements =  Set_to_Pointer(ChainElements,Num_Elem);
      free_Set(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = get_Intersection_Of(FEM_Mesh.NodeNeighbour[NodesElem[3]],
					FEM_Mesh.NodeNeighbour[NodesElem[0]]);
      Num_Elem = get_Lenght_Set(ChainElements);
      Tributary_Elements = Set_to_Pointer(ChainElements,Num_Elem);
      free_Set(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[3]];
      Num_Elem = get_Lenght_Set(ChainElements);
      Tributary_Elements =  Set_to_Pointer(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }

    /* Free memory */
    free(NodesElem);
    
  }  
  /* Check if I am in the 3º Quadrant */
  else if((X_EC_GP.nV[0]<=0) &&
	  (X_EC_GP.nV[1]<=0)){
    /* Create an array with the nodes of the element */
    NodesElem = Set_to_Pointer(FEM_Mesh.Connectivity[Elem_GP],4);   
    if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = get_Intersection_Of(FEM_Mesh.NodeNeighbour[NodesElem[3]],
					FEM_Mesh.NodeNeighbour[NodesElem[0]]);
      Num_Elem = get_Lenght_Set(ChainElements);
      Tributary_Elements = Set_to_Pointer(ChainElements,Num_Elem);
      free_Set(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */
      ChainElements = get_Intersection_Of(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = get_Lenght_Set(ChainElements);
      Tributary_Elements = Set_to_Pointer(ChainElements,Num_Elem);
      free_Set(ChainElements);
      
      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(X_EC_GP.nV[1]) >= Dist[1]) &&
	    (fabs(X_EC_GP.nV[0]) >= Dist[0])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[0]];
      Num_Elem = get_Lenght_Set(ChainElements);
      Tributary_Elements = Set_to_Pointer(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    /* Free memory */
    free(NodesElem);    
  }  
  /* Check if it I am the 4º Quadrant */
  else if((X_EC_GP.nV[0]>=0) &&
	  (X_EC_GP.nV[1]<=0)){
    /* Create an array with the nodes of the element */
    NodesElem = Set_to_Pointer(FEM_Mesh.Connectivity[Elem_GP],4);

    if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = get_Intersection_Of(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = get_Lenght_Set(ChainElements);
      Tributary_Elements =  Set_to_Pointer(ChainElements,Num_Elem);
      free_Set(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = get_Intersection_Of(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					  FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = get_Lenght_Set(ChainElements);
      Tributary_Elements = Set_to_Pointer(ChainElements,Num_Elem);
      free_Set(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[1]];
      Num_Elem = get_Lenght_Set(ChainElements);
      Tributary_Elements = Set_to_Pointer(ChainElements,Num_Elem);
      
      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = get_Union_Of(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    /* Free memory */
    free(NodesElem);
  }
  else{
    printf("%s : %s \n",
	   "Error in Tributary_Nodes_GIMP",
	   "Unlocated GP in the element");
    printf("%s : (%f;%f) \n",
	   "Natural coordinates of the GP",
	   X_EC_GP.nV[0],X_EC_GP.nV[1]);
    exit(0);
  }
 
  return Triburary_Nodes;
}

  
/*   /\* Chain with the tributary elements, this is the list of element near the */
/*      gauss point, including where it is *\/ */

/*   /\* Iterate in the list and select the union of the sets of nodes *\/ */
/*   Table_Elem = malloc(NumNodesElem*sizeof(ChainPtr)); */
/*   for(int i = 0 ; i<NumNodesElem ; i++){ */
/*     Table_Elem[i] = FEM_Mesh.NodeNeighbour[NodesElem[i]]; */
/*   } */
/*   Triburary_Elements = get_Union_Of(Table_Elem,NumNodesElem); */
/*   /\* Free memory *\/ */
/*   free(NodesElem); */
/*   free(Table_Elem); */
/*   Table_Elem = NULL; */
  
/*   /\* List with the tributary nodes *\/ */
/*   Num_Elem = get_Lenght_Set(Triburary_Elements); */
/*   List_Elements = Set_to_Pointer(Triburary_Elements,Num_Elem); */

/*   /\* Free the chain with the tributary elements *\/ */
/*   free_Set(Triburary_Elements); */
  
/*   /\* Fill the chain with the preliminary tributary nodes *\/ */
/*   Table_ElemNodes = malloc(Num_Elem*sizeof(ChainPtr)); */
/*   for(int i = 0 ; i<Num_Elem ; i++){ */
/*     Table_ElemNodes[i] = FEM_Mesh.Connectivity[List_Elements[i]]; */
/*   } */
  
/*   List_Nodes = get_Union_Of(Table_ElemNodes,Num_Elem); */
  
/*   /\* Free the array wit the list of tributary elements *\/ */
/*   free(List_Elements); */
/*   free(Table_ElemNodes); */
/*   Table_ElemNodes = NULL; */
  
/*   /\* Initialize the iterator to iterate over the list of tributary nodes *\/ */
/*   iPtr = List_Nodes; */

/*   /\* Loop over the chain with the tributary nodes *\/ */
/*   while(iPtr != NULL){ */

/*     /\* Assign to a pointer the coordinates of the nodes *\/ */
/*     X_I.nV = FEM_Mesh.Coordinates.nM[iPtr->I]; */

/*     /\* Get a vector from the GP to the node *\/ */
/*     Distance = Sub_Mat(X_GP,X_I); */

/*     /\* If the node is near the GP push in the chain *\/ */
/*     if(Norm_Mat(Distance,2) <= Ra){ */
/*       PushNodeTop(&Triburary_Nodes,iPtr->I); */
/*     } */

/*     /\* Free memory of the distrance vector *\/ */
/*     FreeMat(Distance); */

/*     /\* Update pointer index *\/ */
/*     iPtr = iPtr->next; */
/*   } */
/*   /\* Free memory *\/ */
/*   free_Set(List_Nodes); */
  
/*   return Triburary_Nodes; */
/* } */


