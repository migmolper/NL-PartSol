#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../MathTools/MathTools.h"

/***********************************************/
/************** 1D linear element **************/
/***********************************************/

/*              */
/*  o-------o   */
/* (0)     (1)  */

Matrix L2(Matrix X_e){

  if( (X_e.n > 1 ) || (X_e.n < 0 ) ){
    printf("Error in L2() : Out of the element bounds !!! \n");
    exit(0);
  }    

  /* Definition and allocation */
  Matrix N_ref =  MatAlloc(1,2);

  N_ref.nV[0] = 1-X_e.n;
  N_ref.nV[1] = X_e.n;
  
  return N_ref;
}

/* Derivatives of the shape functions */
Matrix dL2(Matrix X_e){

  if( (X_e.n > 1 ) || (X_e.n < 0 ) ){
    printf("Error in dL2() : Out of the element bounds !!! \n");
    exit(0);
  }    

  /* Definition and allocation */
  Matrix dNdX_ref = MatAlloc(1,2);

  /* Fill the matrix */
  dNdX_ref.nV[0] = -1;
  dNdX_ref.nV[1] = +1;

  return dNdX_ref;
}

/* Global coordinates of the four nodes quadrilateral */
Matrix Get_X_GC_L2(Matrix X_NC_GP,Matrix X_GC_Nodes)
/*
This function evaluate the position of the GP in the element, and get it global coordiantes    
 */
{
  /* 0º Variable declaration */
  Matrix N_ref;
  Matrix X_GC_GP;
  
  N_ref = L2(X_NC_GP);

  X_GC_GP.n = N_ref.nV[0]*X_GC_Nodes.nM[0][0] +
    N_ref.nV[1]*X_GC_Nodes.nM[1][0];

  FreeMat(N_ref);

  return X_GC_GP;
 
}

/* Deformation gradient of the two-nodes linear element */
Matrix Get_F_Ref_L2(Matrix X_NC_GP,Matrix X_GC_Nodes)
/*
  Get the deformation gradient of the reference element:

  F_ref = grad(N_{\alpha}) \otimes x_{\alpha}
  
  Inputs :
  - X_g -> This are the coordinates of the nodes
  - dNdX_Ref_GP -> Derivative gradient evaluated in the GP

  Output :
  - F_Ref -> Deformation gradient of the reference element
  evaluated in the GP
*/
{
  /* Variable declaration */
  Matrix F_Ref;
  Matrix dNdX_Ref_GP;


  /* 1º Evaluate the derivarive of the shape function in the GP */
  dNdX_Ref_GP = dL2(X_NC_GP);

  /* 2º Get the F_ref */
  F_Ref.n = dNdX_Ref_GP.nV[0]*X_GC_Nodes.nM[0][0] +
    dNdX_Ref_GP.nV[1]*X_GC_Nodes.nM[1][0];
    
  /* 3º Free memory */
  FreeMat(dNdX_Ref_GP);
    

  return F_Ref;
}


/*********************************************************************/






