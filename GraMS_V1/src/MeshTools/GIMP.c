#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../MathTools/MathTools.h"
#include "MeshTools.h"

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

/* Uniform GIMP shape function */
double uGIMP(double L, double lp, double Delta_xp){

  /* Evaluation of the shape function */

  /* if(fabs(Delta_xp) < lp){ */
  /*   return 1 - 0.5*(Delta_xp*Delta_xp+lp*lp)*(1/(L*lp)); */
  /* } */
  /* else if ((lp <= fabs(Delta_xp)) && (fabs(Delta_xp) <= L-lp)){ */
  /*   return 1 - fabs(Delta_xp)*(1/L); */
  /* } */
  /* else if((L-lp <= fabs(Delta_xp)) && (fabs(Delta_xp) < L+lp)){ */
  /*   return 0.25*(L+lp-fabs(Delta_xp))*(L+lp-fabs(Delta_xp))*(1/(L*lp)); */
  /* } */
  /* else{ */
  /*   return 0; */
  /* } */

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

  /* if((lp <= fabs(Delta_xp)) && (fabs(Delta_xp) <= L-lp)){ */
  /*   return -SignumFunct(Delta_xp)*(1/L); */
  /* } */
  /* else if((L-lp <= fabs(Delta_xp)) && (fabs(Delta_xp) <= L+lp)){ */
  /*   return -SignumFunct(Delta_xp)*(L+lp-fabs(Delta_xp))/(2*L*lp); */
  /* } */
  /* else if(fabs(Delta_xp) <= lp){ */
  /*   return -Delta_xp/(L*lp); */
  /* } */
  /* else{ */
  /*   return 0; */
  /* } */

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

ChainPtr Tributary_Nodes_GIMP(Matrix X_EC_GP,
			      int Elem_GP,Matrix lp,
			      Mesh FEM_Mesh){

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
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);
    if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */       
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }
    else if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[2]];
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);
      
    }
    /* Free memory */
    free(NodesElem);    
  }  
  /* Check if I am in the 2º Quadrant */
  else if((X_EC_GP.nV[0]<=0) &&
	  (X_EC_GP.nV[1]>=0)){
    /* Create an array with the nodes of the element */
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);

    if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) >= Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Num_Elem = LenghtChain(ChainElements);      
      Tributary_Elements =  ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[3]],
					FEM_Mesh.NodeNeighbour[NodesElem[0]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[3]];
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements =  ChainToArray(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);      
    }

    /* Free memory */
    free(NodesElem);
    
  }  
  /* Check if I am in the 3º Quadrant */
  else if((X_EC_GP.nV[0]<=0) &&
	  (X_EC_GP.nV[1]<=0)){
    /* Create an array with the nodes of the element */
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);   
    if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[3]],
					FEM_Mesh.NodeNeighbour[NodesElem[0]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);
      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */      
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[1]) >= Dist[1]) &&
	    (fabs(X_EC_GP.nV[0]) >= Dist[0])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[0]];
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);      
    }
    /* Free memory */
    free(NodesElem);    
  }  
  /* Check if it I am the 4º Quadrant */
  else if((X_EC_GP.nV[0]>=0) &&
	  (X_EC_GP.nV[1]<=0)){
    /* Create an array with the nodes of the element */
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);

    if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements =  ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */      
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */      
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[1]];
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);      
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

/*********************************************************************/
