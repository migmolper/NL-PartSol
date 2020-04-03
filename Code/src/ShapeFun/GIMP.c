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

void GIMP_Initialize(GaussPoint MPM_Mesh, Mesh FEM_Mesh){

  /* Variables for the GP coordinates */  
  Matrix X_GC_GP = MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  Matrix X_EC_GP = MatAssign(NumberDimensions,1,NAN,NULL,NULL);

  Matrix lp; /* Particle voxel */
  double Vol_GP; /* Volumen of the GP */

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
    free_Set(FEM_Mesh.GPsElements[i]);
    FEM_Mesh.GPsElements[i] = NULL;
  }

  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 1º Get the voxel lenght */
    Vol_GP = MPM_Mesh.Phi.mass.nV[i]/MPM_Mesh.Phi.rho.nV[i];
    for(int j =0 ; j<NumberDimensions ; j++){
      MPM_Mesh.lp.nM[i][j] = 0.5*pow(Vol_GP,(double)1/NumberDimensions);;
    }
    lp.nV = MPM_Mesh.lp.nM[i];

    /* 2º Assign the value to this auxiliar pointer */ 
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){

      /* 3º Connectivity of the Poligon */
      NumVertex = FEM_Mesh.NumNodesElem[j];
      Poligon_Connectivity =
	Set_to_Pointer(FEM_Mesh.Connectivity[j],NumVertex);

      /* 4º Get the coordinates of the element */
      Poligon_Coordinates =
	ElemCoordinates(FEM_Mesh,Poligon_Connectivity,NumVertex);
      
      /* 5º Check out if the GP is in the Element */
      if(InOut_Poligon(X_GC_GP,Poligon_Coordinates) == 1){

	/* 6º Asign to the GP a element in the background mesh, just for 
	   searching porpuses */
	MPM_Mesh.Element_id[i] = j;
	push_to_Set(&FEM_Mesh.GPsElements[j],i);

	/* 7º If the GP is in the element, get its natural coordinates */
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	Q4_X_to_Xi(X_EC_GP,X_GC_GP,Poligon_Coordinates);

	/* 8º Get list of nodes near to the GP */
	free_Set(MPM_Mesh.ListNodes[i]);
	MPM_Mesh.ListNodes[i] = NULL;
	MPM_Mesh.ListNodes[i] =
	  Tributary_Nodes_GIMP(X_GC_GP,MPM_Mesh.Element_id[i],lp,FEM_Mesh);
	/* MPM_Mesh.ListNodes[i] = */
	/*   Tributary_Nodes_GIMP(X_EC_GP,MPM_Mesh.Element_id[i],lp,FEM_Mesh); */
	MPM_Mesh.NumberNodes[i] = get_Lenght_Set(MPM_Mesh.ListNodes[i]);
	/* 9º Active those nodes that interact with the GP */
	ListNodes_I = MPM_Mesh.ListNodes[i];
	while(ListNodes_I != NULL){
	  FEM_Mesh.ActiveNode[ListNodes_I->I] += 1;
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
double uGIMP(double L, double lp, double Delta_xp){

  /* Evaluation of the shape function */

  if ((-lp < Delta_xp) && (Delta_xp <= lp)){
    return 1 - 0.5*(Delta_xp*Delta_xp + lp*lp)*(double)1/(L*lp);
  }
  else if (((-L-lp) < Delta_xp) && (Delta_xp <= (-L+lp))){
    return (double)(0.25/(L*lp))*(L+lp+Delta_xp)*(L+lp+Delta_xp);
  }
  else if (((L-lp) < Delta_xp) && (Delta_xp <= (L+lp))){
    return (double)(0.25/(L*lp))*(L+lp-Delta_xp)*(L+lp-Delta_xp);
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
double d_uGIMP(double L, double lp, double Delta_xp){

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
Matrix GIMP_2D(Matrix Delta_Xp, Matrix lp, double L){

  /* 1º Variable declaration */
  Matrix S_Ip = MatAlloc(1,Delta_Xp.N_rows);

  /* 2º Fill the shape function array */
  for(int i = 0 ; i<Delta_Xp.N_rows ; i++){
    /* 3º Shape function in this node */
    S_Ip.nV[i] =
      uGIMP(L, lp.nV[0], Delta_Xp.nM[i][0])*
      uGIMP(L, lp.nV[1], Delta_Xp.nM[i][1]);
  }

  /* 4º Output */
  return S_Ip;
}

/*********************************************************************/

/* Uniform GIMP derivative shape function 2D */
Matrix dGIMP_2D(Matrix Delta_xp, Matrix lp, double L){

  /* 1º Variable declaration */
  Matrix dS_Ip = MatAlloc(2,Delta_xp.N_rows); 
  
  /* 2º Fill the shape function array */
  for(int i = 0 ; i<Delta_xp.N_rows ; i++){
    /* 3º Gradient of the shape function for each node*/
    for(int j = 0, k = 1 ; j<=1 && k>=0 ; j++, k--){
      dS_Ip.nM[j][i] =
    	d_uGIMP(L, lp.nV[j], Delta_xp.nM[i][j]) *
    	uGIMP(L, lp.nV[k], Delta_xp.nM[i][k]);
    }
  }

  /* 4º Output */
  return dS_Ip;
}

/*********************************************************************/

ChainPtr Tributary_Nodes_GIMP(Matrix X_GP, int Elem_GP,
			      Matrix lp, Mesh FEM_Mesh){

  Matrix Distance; /* Distance between node and GP */
  Matrix X_I = MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  ChainPtr * Table_Elem = NULL;
  ChainPtr Triburary_Nodes = NULL;
  ChainPtr List_Nodes = NULL;
  ChainPtr * Table_ElemNodes = NULL;
  ChainPtr Triburary_Elements = NULL;
  ChainPtr iPtr = NULL;
  int NumNodesElem = /* Number of nodes of the element */
    FEM_Mesh.NumNodesElem[Elem_GP];
  int Num_Elem;
  int * List_Elements;
  int * NodesElem = /* List of nodes of the element */
    Set_to_Pointer(FEM_Mesh.Connectivity[Elem_GP],NumNodesElem);
  double Ra = /* Get the search radius */
    FEM_Mesh.DeltaX + lp.nV[0];
  
  /* Chain with the tributary elements, this is the list of element near the
     gauss point, including where it is */

  /* Iterate in the list and select the union of the sets of nodes */
  Table_Elem = malloc(NumNodesElem*sizeof(ChainPtr));
  for(int i = 0 ; i<NumNodesElem ; i++){
    Table_Elem[i] = FEM_Mesh.NodeNeighbour[NodesElem[i]];
  }
  Triburary_Elements = get_Union_Of(Table_Elem,NumNodesElem);
  /* Free memory */
  free(NodesElem);
  free(Table_Elem);
  Table_Elem = NULL;
  
  /* List with the tributary nodes */
  Num_Elem = get_Lenght_Set(Triburary_Elements);
  List_Elements = Set_to_Pointer(Triburary_Elements,Num_Elem);

  /* Free the chain with the tributary elements */
  free_Set(Triburary_Elements);
  
  /* Fill the chain with the preliminary tributary nodes */
  Table_ElemNodes = malloc(Num_Elem*sizeof(ChainPtr));
  for(int i = 0 ; i<Num_Elem ; i++){
    Table_ElemNodes[i] = FEM_Mesh.Connectivity[List_Elements[i]];
  }
  
  List_Nodes = get_Union_Of(Table_ElemNodes,Num_Elem);
  
  /* Free the array wit the list of tributary elements */
  free(List_Elements);
  free(Table_ElemNodes);
  Table_ElemNodes = NULL;
  
  /* Initialize the iterator to iterate over the list of tributary nodes */
  iPtr = List_Nodes;

  /* Loop over the chain with the tributary nodes */
  while(iPtr != NULL){

    /* Assign to a pointer the coordinates of the nodes */
    X_I.nV = FEM_Mesh.Coordinates.nM[iPtr->I];

    /* Get a vector from the GP to the node */
    Distance = Sub_Mat(X_GP,X_I);

    /* If the node is near the GP push in the chain */
    if(Norm_Mat(Distance,2) <= Ra){
      push_to_Set(&Triburary_Nodes,iPtr->I);
    }

    /* Free memory of the distrance vector */
    FreeMat(Distance);

    /* Update pointer index */
    iPtr = iPtr->next;
  }
  /* Free memory */
  free_Set(List_Nodes);
  
  return Triburary_Nodes;
}

