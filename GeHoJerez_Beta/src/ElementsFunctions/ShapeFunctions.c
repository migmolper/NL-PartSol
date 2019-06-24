#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

/* Global coordinates of the four nodes quadrilateral */
Matrix Get_GlobalCoordinates_L2(Matrix X_NC_GP,Matrix X_GC_Nodes)
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

  free(N_ref.nV);

  return X_GC_GP;
 
}

/* Deformation gradient of the two-nodes linear element */
Matrix Get_RefDeformation_Gradient_L2(Matrix X_NC_GP,Matrix X_GC_Nodes)
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
  free(dNdX_Ref_GP.nM);
    

  return F_Ref;
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

  if( (fabs(X_e.nV[0]) > 1 ) || (fabs(X_e.nV[1]) > 1 ) ){
    printf("Error in Q4() : Out of the element bounds !!! \n");
    exit(0);
  }    
  
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

  if( (fabs(X_e.nV[0]) > 1 ) || (fabs(X_e.nV[1]) > 1 ) ){
    printf("Error in dQ4() : Out of the element bounds !!! \n");
    exit(0);
  }    
  
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

/* Global coordinates of the four nodes quadrilateral */
Matrix Get_GlobalCoordinates_Q4(Matrix X_NC_GP,Matrix X_GC_Nodes)
/*
This function evaluate the position of the GP in the element, and get it global coordiantes    
 */
{
  /* 0º Variable declaration */
  Matrix N_ref;
  Matrix X_GC_GP;

  /* 1º Evaluate the Q4 element in the element coordinates */
  N_ref = Q4(X_NC_GP);

  /* 2º Allocate the output coordinates */
  X_GC_GP = MatAlloc(2,1);

  /* 3º Get the global coordinates for this element coordiantes in this element */
  X_GC_GP.nV[0] =
    N_ref.nV[0]*X_GC_Nodes.nM[0][0] +
    N_ref.nV[1]*X_GC_Nodes.nM[1][0] +
    N_ref.nV[2]*X_GC_Nodes.nM[2][0] +
    N_ref.nV[3]*X_GC_Nodes.nM[3][0];
  
  X_GC_GP.nV[1] =
    N_ref.nV[0]*X_GC_Nodes.nM[0][1] +
    N_ref.nV[1]*X_GC_Nodes.nM[1][1] +
    N_ref.nV[2]*X_GC_Nodes.nM[2][1] +
    N_ref.nV[3]*X_GC_Nodes.nM[3][1];

  /* 4º Free memory */
  free(N_ref.nV);

  /* 5º Output */
  return X_GC_GP;
 
}

/* Deformation gradient of the four-nodes quadrilateral */
Matrix Get_Jacobian_Q4(Matrix X_NC_GP,Matrix X_GC_Nodes)
/*
  Get the jacobian of the transformation of the reference element :

  J = grad(N_{\alpha}) \cdot x_{\alpha}
  
  Inputs :
  - X_g -> This are the coordinates of the nodes
  - dNdX_Ref_GP -> Derivative gradient evaluated in the GP

  Output :
  - Jacobian -> Jacobian of the reference element
  evaluated in the GP
*/
{
  /* Variable declaration */
  Matrix Jacobian = MatAlloc(2,2);
  Matrix dNdX_Ref_GP;

  /* 1º Evaluate the derivarive of the shape function in the GP */
  dNdX_Ref_GP = dQ4(X_NC_GP); 

  /* 2º Get the F_ref */
  Jacobian = Scalar_prod(dNdX_Ref_GP,X_GC_Nodes);
  
  /* 3º Free memory */
  free(dNdX_Ref_GP.nM);
  
  return Jacobian;
}


/* Element gradient in the real element */
Matrix Get_dNdX_Q4(Matrix X_EC_GP,Matrix Element)
/*
  - Matrix X_EC_GP : Element coordinates of the gauss point
  - Matrix Element : Coordinates of the element (NumNodesElem x NumberDimensions)
*/
{
  
  /* 0º Definition and allocation */
  Matrix dNdX_Ref_GP; /* Derivative of the shape function evaluated in the GP (Ndim x Nnodes) */
  Matrix J_GP; /* Jacobian of the transformation evaluated in the GP (Ndim x Ndim) */
  Matrix J_GP_m1; /* Inverse of the Jacobian */
  Matrix J_GP_Tm1; /* Transpose of the inverse Jacobian */
  Matrix dNdx_GP; /* Derivatives of the shape function evaluates in the GP (Ndim x Ndim) */

  /* 1º Evaluate the gradient of the shape function in the GP */
  dNdX_Ref_GP = dQ4(X_EC_GP);
  
  /* 2º Get the Jacobian of the transformation evaluated in the GP */
  J_GP = Get_Jacobian_Q4(X_EC_GP,Element);
    
  /* 3º Get the inverse of the deformation gradient */
  J_GP_m1 = Get_Inverse(J_GP), 
    free(J_GP.nM);
  /* 4º Get the transpose of the inverse of the Jacobian */
  J_GP_Tm1 = Transpose_Mat(J_GP_m1),
    free(J_GP_m1.nM);
  
  /* 5º Get the gradient of the shape functions in global coordinates */
  dNdx_GP = Scalar_prod(J_GP_Tm1,dNdX_Ref_GP),
    free(J_GP_Tm1.nM),
    free(dNdX_Ref_GP.nM);

  /* 6º Return result */
  return dNdx_GP;
}
