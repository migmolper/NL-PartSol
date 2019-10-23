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
double uGIMP(double L, double lp, double Xp, double Xi){

  /* Variable definition */
  double S_Ip;
  double Delta_xp;

  /* Calcule the distance from the GP to the node */
  Delta_xp = Xp-Xi;

  /* Evaluation of the shape function */
  if ((-L+lp < Delta_xp) && (Delta_xp <= -lp)){
    S_Ip = 1 + (Xp-Xi)/L;
  }
  else if ((lp < Delta_xp) && (Delta_xp <= L-lp)){
    S_Ip = 1 - (Xp-Xi)/L;
  }
  else if ((-lp < Delta_xp) && (Delta_xp <= lp)){
    S_Ip = 1 - ((Xp-Xi)*(Xp-Xi) + lp*lp)/(2*L*lp);
  }
  else if ((-L-lp < Delta_xp) && (Delta_xp <= -L+lp)){
    S_Ip = ((L+lp+Xp-Xi)*(L+lp+Xp-Xi))/(1/(4*L*lp));
  }
  else if ((L-lp < Delta_xp) && (Delta_xp <= L+lp)){
    S_Ip = ((L+lp-Xp+Xi)*(L+lp-Xp+Xi))/(4*L*lp);
  }
  else if (fabs(Delta_xp) >= L+lp){
    S_Ip = 0;
  }
  
  return S_Ip;
}

/*********************************************************************/

/* Uniform GIMP derivative shape function */
double d_uGIMP(double L, double lp, double Xp, double Xi){

  /* Variable definition */
  double dS_Ip;
  double Delta_xp;

  /* Calcule the distance from the GP to the node */
  Delta_xp = Xp-Xi;

  /* Evaluation of the shape function */
  if ((-L+lp < Delta_xp) && (Delta_xp <= -lp)){
    dS_Ip = 1/L;
  }
  else if ((lp < Delta_xp) && (Delta_xp <= L-lp)){
    dS_Ip = 1/L;
  }
  else if ((-lp < Delta_xp) && (Delta_xp <= lp)){
    dS_Ip = - (Xp-Xi)/(L*lp);
  }
  else if ((-L-lp < Delta_xp) && (Delta_xp <= -L+lp)){
    dS_Ip = (L+lp+Xp-Xi)/(2*L*lp);
  }
  else if ((L-lp < Delta_xp) && (Delta_xp <= L+lp)){
    dS_Ip = (L+lp-Xp+Xi)/(4*L*lp);
  }
  else if (fabs(Delta_xp) >= L+lp){
    dS_Ip = 0;
  }
  
  return dS_Ip;
}

/*********************************************************************/

/* Uniform GIMP shape function 2D */
Matrix GIMP_2D(Matrix X_GC_GP, Matrix lp, Matrix Element, double L){

  /* 1º Variable declaration */
  Matrix S_Ip = MatAlloc(1,Element.N_rows);
  Matrix X_GC_I;

  /* 2º Fill the shape function array */
  for(int i = 0 ; i<Element.N_rows ; i++){

    /* 3º Coordinate of the node */
    X_GC_I.nV = Element.nM[i];

    /* 4º Shape function in this node */
    S_Ip.nV[i] =
      uGIMP(L, lp.nV[0], X_GC_GP.nV[0], X_GC_I.nV[0])*
      uGIMP(L, lp.nV[1], X_GC_GP.nV[1], X_GC_I.nV[1]);
  }

  /* 5º Output */
  return S_Ip;
}

/*********************************************************************/

/* Uniform GIMP derivative shape function 2D */
Matrix dGIMP_2D(Matrix X_GC_GP, Matrix lp, Matrix Nodes, double L){

  /* 1º Variable declaration */
  Matrix dS_Ip = MatAlloc(2,Nodes.N_rows);
  Matrix X_GC_I;
  
  /* 2º Fill the shape function array */
  for(int i = 0 ; i<Nodes.N_rows ; i++){
    
    /* 3º Coordinate of the node */
    X_GC_I.nV = Nodes.nM[i];

    /* 4º Gradient of the shape function for each node*/
    for(int j = 0, k = 2 ; j<2 && k>0 ; j++, k--){
      dS_Ip.nM[i][j] =
	d_uGIMP(L, lp.nV[j], X_GC_GP.nV[j], X_GC_I.nV[j]) *
	uGIMP(L, lp.nV[k], X_GC_GP.nV[k], X_GC_I.nV[k]);    
    }
  }

  /* 5º Output */
  return dS_Ip;
}

/*********************************************************************/

ChainPtr Tributary_Nodes_GIMP(Matrix X_EC_GP,
			      ChainPtr NodesElem,
			      Matrix lp,double L,
			      int ** NodeNeighbour){

  ChainPtr Triburary_Nodes_n1;
  double Dist[2];

  /* Get the reference distance measured from the center of the element */
  for(int i = 0 ; i<2; i++){
    Dist[i] = 1 - lp.nV[i]/L;
  }

  /* Check if I am in the central area */
  if ((fabs(X_EC_GP.nV[0]) < Dist[0]) &&
      (fabs(X_EC_GP.nV[1]) < Dist[1])){    
    Triburary_Nodes_n1 = NodesElem;
  }
  
  /* Check if I am in the 1º Quadrant */
  else if((X_EC_GP.nV[0]>0) &&
	  (X_EC_GP.nV[1]>0)){

    if((X_EC_GP.nV[0] > Dist[0]) &&
       (X_EC_GP.nV[1] < Dist[1])){
      
    }
    else if((X_EC_GP.nV[0] < Dist[0]) &&
	    (X_EC_GP.nV[1] > Dist[1])){

    }
    else if((X_EC_GP.nV[0] > Dist[0]) &&
	    (X_EC_GP.nV[1] > Dist[1])){

    }
    
  }
  
  /* Check if I am in the 2º Quadrant */
  else if((X_EC_GP.nV[0]<0) &&
	  (X_EC_GP.nV[1]>0)){

    if((X_EC_GP.nV[0] > -Dist[0]) &&
       (X_EC_GP.nV[1] > Dist[1])){

    }
    else if((X_EC_GP.nV[0] < -Dist[0]) &&
	    (X_EC_GP.nV[1] < Dist[1])){

    }
    else if((X_EC_GP.nV[0] < -Dist[0]) &&
	    (X_EC_GP.nV[1] > Dist[1])){

    }
    
  }
  
  /* Check if I am in the 3º Quadrant */
  else if((X_EC_GP.nV[0]<0) &&
	  (X_EC_GP.nV[1]<0)){
    
    if((X_EC_GP.nV[0] < -Dist[0]) &&
       (X_EC_GP.nV[1] > -Dist[1])){

    }
    else if((X_EC_GP.nV[0] > -Dist[0]) &&
	    (X_EC_GP.nV[1] < -Dist[1])){

    }
    else if((X_EC_GP.nV[1] < -Dist[1]) &&
	    (X_EC_GP.nV[0] < -Dist[0])){

    }
    
  }
  
  /* Check if it I am the 4º Quadrant */
  else if((X_EC_GP.nV[0]>0) &&
	  (X_EC_GP.nV[1]<0)){

    if((X_EC_GP.nV[0] < Dist[0]) &&
       (X_EC_GP.nV[1] < -Dist[1])){

    }
    else if((X_EC_GP.nV[0] > Dist[0]) &&
	    (X_EC_GP.nV[1] > -Dist[1])){

    }
    else if((X_EC_GP.nV[0] > Dist[0]) &&
	    (X_EC_GP.nV[1] < -Dist[1])){

    }
    
  }
 
  return Triburary_Nodes_n1;
}

/*********************************************************************/
