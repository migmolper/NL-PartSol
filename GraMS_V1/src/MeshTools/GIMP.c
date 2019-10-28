#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../MathTools/MathTools.h"
#include "MeshTools.h"

/***********************************************/
/******* 2D cuadrilateral linear element *******/
/***********************************************/

/*            ^           */
/*          __|__         */
/*        _/  |  \_       */
/*      _/    |    \_     */
/*   __/      |      \__  */
/*  --o-------o-------o-- */
/*   (-1)    (0)     (1)  */

/* Uniform GIMP shape function */
double uGIMP(double L, double lp, double Delta_xp){

  /* Evaluation of the shape function */
  if ((0.5*lp <= fabs(Delta_xp)) && (fabs(Delta_xp) <= L-0.5*lp)){
    return 1 - fabs(Delta_xp)/L;
  }
  else if((L-0.5*lp <= fabs(Delta_xp)) && (fabs(Delta_xp) < L + 0.5*lp)){
    return 0.5*(L+0.5*lp-fabs(Delta_xp))*(L+0.5*lp-fabs(Delta_xp))/(L*lp);
  }
  else if(fabs(Delta_xp) < 0.5*lp){
    return 1 - 0.25*(4*Delta_xp*Delta_xp+lp*lp)/(L*lp);
  }
  else{
    return 0;
  }
  
}

/*********************************************************************/

/* Uniform GIMP derivative shape function */
double d_uGIMP(double L, double lp, double Delta_xp){

  /* Evaluation of the shape function */
  if((0.5*lp <= fabs(Delta_xp)) && (fabs(Delta_xp) <= L-0.5*lp)){
    return -(1/L)*SignumFunct(Delta_xp);
  }
  else if((L-0.5*lp <= fabs(Delta_xp)) && (fabs(Delta_xp) < L+0.5*lp)){
    return -SignumFunct(Delta_xp)*(L+0.5*lp-fabs(Delta_xp))/(L*lp);
  }
  else if(fabs(Delta_xp) < 0.5*lp){
    return -2*Delta_xp/(L*lp);
  }
  else{
    return 0;
  }
  
}

/*********************************************************************/

/* Uniform GIMP shape function 2D */
Matrix GIMP_2D(Matrix Delta_xp, Matrix lp, double L){

  /* 1º Variable declaration */
  Matrix S_Ip = MatAlloc(1,Delta_xp.N_rows);

  /* 2º Fill the shape function array */
  for(int i = 0 ; i<Delta_xp.N_rows ; i++){
    /* 3º Shape function in this node */
    S_Ip.nV[i] =
      uGIMP(L, lp.nV[0], Delta_xp.nM[i][0])*
      uGIMP(L, lp.nV[1], Delta_xp.nM[i][1]);
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

ChainPtr Tributary_Nodes_GIMP(Matrix X_EC_GP,
			      int Elem_GP,Matrix lp,
			      Mesh FEM_Mesh){

  ChainPtr Triburary_Nodes;
  ChainPtr ChainElements;
  int * Tributary_Elements;
  int * NodesElem;
  int Elem_i;
  double Dist[2];

  /* Get the reference distance measured from the center of the element */
  for(int i = 0 ; i<2; i++){
    Dist[i] = 1 - lp.nV[i]/FEM_Mesh.DeltaX;
  }

  /* Check if I am in the central area */
  if ((fabs(X_EC_GP.nV[0]) < Dist[0]) &&
      (fabs(X_EC_GP.nV[1]) < Dist[1])){    
    Triburary_Nodes = FEM_Mesh.Connectivity[Elem_GP];
  }
  
  /* Check if I am in the 1º Quadrant */
  else if((X_EC_GP.nV[0]>0) &&
	  (X_EC_GP.nV[1]>0)){

    /* Create an array with the nodes of the element */
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);

    if((X_EC_GP.nV[0] > Dist[0]) &&
       (X_EC_GP.nV[1] < Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */       
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Tributary_Elements =  ChainToArray(ChainElements,2);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<2 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }
    else if((X_EC_GP.nV[0] < Dist[0]) &&
	    (X_EC_GP.nV[1] > Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Tributary_Elements =  ChainToArray(ChainElements,2);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<2 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }
    else if((X_EC_GP.nV[0] > Dist[0]) &&
	    (X_EC_GP.nV[1] > Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[0]];
      Tributary_Elements =  ChainToArray(ChainElements,4);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<4 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }

    /* Free memory */
    free(NodesElem);
    
  }
  
  /* Check if I am in the 2º Quadrant */
  else if((X_EC_GP.nV[0]<0) &&
	  (X_EC_GP.nV[1]>0)){

    /* Create an array with the nodes of the element */
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);

    if((X_EC_GP.nV[0] > -Dist[0]) &&
       (X_EC_GP.nV[1] > Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Tributary_Elements =  ChainToArray(ChainElements,2);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<2 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }
    else if((X_EC_GP.nV[0] < -Dist[0]) &&
	    (X_EC_GP.nV[1] < Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[1]],
					FEM_Mesh.NodeNeighbour[NodesElem[2]]);
      Tributary_Elements =  ChainToArray(ChainElements,2);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<2 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }
    else if((X_EC_GP.nV[0] < -Dist[0]) &&
	    (X_EC_GP.nV[1] > Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[1]];
      Tributary_Elements =  ChainToArray(ChainElements,4);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<4 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }

    /* Free memory */
    free(NodesElem);
    
  }
  
  /* Check if I am in the 3º Quadrant */
  else if((X_EC_GP.nV[0]<0) &&
	  (X_EC_GP.nV[1]<0)){
    
    /* Create an array with the nodes of the element */
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);
    
    if((X_EC_GP.nV[0] < -Dist[0]) &&
       (X_EC_GP.nV[1] > -Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[1]],
					FEM_Mesh.NodeNeighbour[NodesElem[2]]);
      Tributary_Elements = ChainToArray(ChainElements,2);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<2 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }
    else if((X_EC_GP.nV[0] > -Dist[0]) &&
	    (X_EC_GP.nV[1] < -Dist[1])){
      
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Tributary_Elements = ChainToArray(ChainElements,2);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<2 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }
    else if((X_EC_GP.nV[1] < -Dist[1]) &&
	    (X_EC_GP.nV[0] < -Dist[0])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[2]];
      Tributary_Elements =  ChainToArray(ChainElements,4);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<4 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }

    /* Free memory */
    free(NodesElem);
    
  }
  
  /* Check if it I am the 4º Quadrant */
  else if((X_EC_GP.nV[0]>0) &&
	  (X_EC_GP.nV[1]<0)){

    /* Create an array with the nodes of the element */
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);

    if((X_EC_GP.nV[0] < Dist[0]) &&
       (X_EC_GP.nV[1] < -Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Tributary_Elements =  ChainToArray(ChainElements,2);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<2 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }
    else if((X_EC_GP.nV[0] > Dist[0]) &&
	    (X_EC_GP.nV[1] > -Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Tributary_Elements =  ChainToArray(ChainElements,2);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<2 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }
    else if((X_EC_GP.nV[0] > Dist[0]) &&
	    (X_EC_GP.nV[1] < -Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[3]];
      Tributary_Elements =  ChainToArray(ChainElements,4);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<4 ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes = ChainUnion(Triburary_Nodes,
				     FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }

    /* Free memory */
    free(NodesElem);
    
  }
 
  return Triburary_Nodes;
}

/*********************************************************************/
