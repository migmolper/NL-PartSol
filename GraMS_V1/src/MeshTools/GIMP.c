#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../MathTools/MathTools.h"


double uGIMP(double L, double lp, double Xp, double Xi){
  
  double S_Ip;
  double Delta_xp;

  Delta_xp = Xp-Xi;
  
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


double d_uGIMP(double L, double lp, double Xp, double Xi){
  
  double dS_Ip;
  double Delta_xp;

  Delta_xp = Xp-Xi;
  
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

Matrix GIMP_2D(Matrix X_GC_GP, Matrix lp, Matrix Element, double L){

  /* 1º Variable declaration */
  Matrix S_Ip = MatAlloc(1,Element.n);
  Matrix X_GC_I;

  /* 2º Fill the shape function array */
  for(int i = 0 ; i<Element.n ; i++){

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

Matrix dGIMP_2D(Matrix X_GC_GP, Matrix lp, Matrix Nodes, double L){

  /* 1º Variable declaration */
  Matrix dS_Ip = MatAlloc(2,Nodes.Num);
  Matrix X_GC_I;
  
  /* 2º Fill the shape function array */
  for(int i = 0 ; i<Nodes.Num ; i++){
    
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
