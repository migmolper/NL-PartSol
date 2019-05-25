#include <stdio.h>
#include <stdlib.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"

/***********************************************/
/************** 1D linear element **************/
/***********************************************/

/*              */
/*  o-------o   */
/* (0)     (1)  */

Matrix L2(Matrix X_e){

  /* Definition and allocation */
  Matrix N_ref =  MatAlloc(1,2);

  N_ref.nV[0] = 1-X_e.n;
  N_ref.nV[1] = X_e.n;
  
  return N_ref;
}

/* Derivatives of the shape functions */
Matrix dL2(Matrix X_e){

  /* Definition and allocation */
  Matrix dNdX_ref = MatAlloc(1,2);

  /* Fill the matrix */
  dNdX_ref.nV[0] = -1;
  dNdX_ref.nV[1] = +1;

  return dNdX_ref;
}

/***********************************************/
/********* 2D triangle linear element **********/
/***********************************************/

/* (2)     */
/*  o      */
/*  |\     */
/*  | \    */
/*  o--o   */
/* (0) (1) */

/* Shape functions */
Matrix T3(Matrix X_e){
  
  /* Definition and allocation */
  Matrix N_ref =  MatAlloc(1,3);

  /* Fill the array */
  N_ref.nV[0] = 1 - X_e.nV[0] - X_e.nV[1]; 
  N_ref.nV[1] = X_e.nV[0];
  N_ref.nV[2] = X_e.nV[1];
  
  return N_ref;
}

/* Derivatives of the shape functions */
Matrix dT3(Matrix X_e){
  
  /* Definition and allocation */
  Matrix dNdX_ref = MatAlloc(2,3);
  
  /* Fill the matrix */
  /* Node 0 */
  dNdX_ref.nM[0][0] = - 1; /* \frac{\partial N0}{\partial \xi} */
  dNdX_ref.nM[1][0] = - 1; /* \frac{\partial N0}{\partial \eta} */
  /* Node 1 */
  dNdX_ref.nM[0][1] = + 1; /* \frac{\partial N1}{\partial \xi} */
  dNdX_ref.nM[1][1] = + 0; /* \frac{\partial N1}{\partial \eta} */
  /* Node 2 */
  dNdX_ref.nM[0][2] = + 0; /* \frac{\partial N2}{\partial \xi} */
  dNdX_ref.nM[1][2] = + 1; /* \frac{\partial N2}{\partial \eta} */
  
  return dNdX_ref;
}





/***********************************************/
/******* 2D cuadrilateral linear element *******/
/***********************************************/

/* (3)     (2)  */
/*  o-------o   */
/*  |       |   */
/*  |       |   */
/*  o-------o   */
/* (0)     (1)  */

/* Shape functions */
Matrix Q4(Matrix X_e){
  
  /* Definition and allocation */
  Matrix N_ref =  MatAlloc(1,4);

  /* Fill the array */
  N_ref.nV[0] = 0.25 - 0.25*X_e.nV[0] - 0.25*X_e.nV[1] + 0.25*(X_e.nV[0])*(X_e.nV[1]);
  N_ref.nV[1] = 0.25 + 0.25*X_e.nV[0] - 0.25*X_e.nV[1] - 0.25*(X_e.nV[0])*(X_e.nV[1]);
  N_ref.nV[2] = 0.25 + 0.25*X_e.nV[0] + 0.25*X_e.nV[1] + 0.25*(X_e.nV[0])*(X_e.nV[1]);
  N_ref.nV[3] = 0.25 - 0.25*X_e.nV[0] + 0.25*X_e.nV[1] - 0.25*(X_e.nV[0])*(X_e.nV[1]);
  
  return N_ref;
}

/* Derivatives of the shape functions */
Matrix dQ4(Matrix X_e){
  
  /* Definition and allocation */
  Matrix dNdX_ref = MatAlloc(2,4);
  
  /* Fill the matrix */
  /* Node 1 */
  dNdX_ref.nM[0][0] = - 0.25 + 0.25*(X_e.nV[1]); /* \frac{\partial N0}{\partial \xi} */
  dNdX_ref.nM[1][0] = - 0.25 + 0.25*(X_e.nV[0]); /* \frac{\partial N0}{\partial \eta} */
  /* Node 2 */
  dNdX_ref.nM[0][1] = + 0.25 - 0.25*(X_e.nV[1]); /* \frac{\partial N1}{\partial \xi} */
  dNdX_ref.nM[1][1] = - 0.25 - 0.25*(X_e.nV[0]); /* \frac{\partial N1}{\partial \eta} */
  /* Node 3 */
  dNdX_ref.nM[0][2] = + 0.25 + 0.25*(X_e.nV[1]); /* \frac{\partial N2}{\partial \xi} */
  dNdX_ref.nM[1][2] = + 0.25 + 0.25*(X_e.nV[0]); /* \frac{\partial N2}{\partial \eta} */
  /* Node 4 */
  dNdX_ref.nM[0][3] = - 0.25 - 0.25*(X_e.nV[1]); /* \frac{\partial N3}{\partial \xi} */
  dNdX_ref.nM[1][3] = + 0.25 - 0.25*(X_e.nV[0]); /* \frac{\partial N3}{\partial \eta} */
  
  return dNdX_ref;
}


/*********************************************************************/


Matrix Get_RefDeformation_Gradient_Q4(Matrix X_NC_GP,Matrix X_GC_Nodes)
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
  dNdX_Ref_GP = dQ4(X_NC_GP);

  /* 2º Get the F_ref */
  F_Ref = Scalar_prod(dNdX_Ref_GP, /* grad(N_{\alpha}) */
		      X_GC_Nodes);  /* x_{\alpha} */

  /* 3º Free memory */
  free(dNdX_Ref_GP.nM);
    

  return F_Ref;
}
