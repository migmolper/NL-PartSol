#include "nl-partsol.h"

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

  /* Error check */
  if( (fabs(X_e.nV[0]) > 1 ) ||
      (fabs(X_e.nV[1]) > 1 ) ||
      (1 - X_e.nV[0] - X_e.nV[1] < 0)){
    printf("Error in T3() : Out of the element bounds !!! \n");
    exit(EXIT_FAILURE);
  }    
  
  /* Definition and allocation */
  Matrix N_ref =  MatAlloc(1,3);

  /* Fill the array */
  N_ref.nV[0] = 1 - X_e.nV[0] - X_e.nV[1]; 
  N_ref.nV[1] = X_e.nV[0];
  N_ref.nV[2] = X_e.nV[1];
  
  return N_ref;
}

/*********************************************************************/

/* Derivatives of the shape functions */
Matrix dT3(Matrix X_e){

  /* Error check */
  if( (fabs(X_e.nV[0]) > 1 ) ||
      (fabs(X_e.nV[1]) > 1 ) ||
      (1 - X_e.nV[0] - X_e.nV[1] < 0)){
    printf("Error in T3() : Out of the element bounds !!! \n");
    exit(EXIT_FAILURE);
  }
  
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


/*********************************************************************/

/* Deformation gradient of the reference element for the three-nodes triangle */
Matrix Get_F_Ref_T3(Matrix X_NC_GP,Matrix X_GC_Nodes)
/*
  Get the deformation gradient of the transformation of the reference element :

  F = grad(N_{\alpha}) \cdot x_{\alpha}
  
  Inputs :
  - X_NC_GP -> This are the natural coordinates of the GP
  - X_GC_Nodes -> Global coordinates of the element nodes

  Output :
  - F -> Deformation gradient of the reference element
  evaluated in the GP
*/
{
  /* Variable declaration */
  Matrix dNdX_Ref_GP;
  Matrix X_alpha = MatAlloc(2,1);
  Matrix dNdx_alpha = MatAlloc(1,2);
  Matrix F_Ref_alpha;
  Matrix F_Ref = MatAllocZ(2,2);

  /* 1º Evaluate the derivarive of the shape function in the GP */
  dNdX_Ref_GP = dT3(X_NC_GP); 

  /* 2º Get the F_Ref doing a loop over the nodes of the element */
  for(int i = 0 ; i<3 ; i++){

    /* 3º Fill arrays for the tensorial product */
    for(int j = 0 ; j<2 ; j++){
      X_alpha.nV[j] = X_GC_Nodes.nM[i][j];
      dNdx_alpha.nV[j] = dNdX_Ref_GP.nM[j][i];
    }

    /* 4º Get the nodal contribution */
    F_Ref_alpha = Tensorial_prod(X_alpha,dNdx_alpha);

    /* 5º Increment the reference deformation gradient */
    F_Ref = Incr_Mat(F_Ref, F_Ref_alpha);

    /* 6º Free data of the nodal contribution */
    FreeMat(F_Ref_alpha);
    
  }
  
  /* 7º Free memory */
  FreeMat(dNdX_Ref_GP);
  FreeMat(X_alpha);
  FreeMat(dNdx_alpha);

  /* 8º Output */
  return F_Ref;
}

/*********************************************************************/

/* Element gradient in the real element */
Matrix Get_dNdX_T3(Matrix X_EC_GP,Matrix Element)
/*
  - Matrix X_EC_GP : Element coordinates of the gauss point
  - Matrix Element : Coordinates of the element (NumNodesElem x NumberDimensions)
*/
{
  
  /* 0º Definition and allocation */
  Matrix dNdX_Ref_GP; /* Derivative of the shape function evaluated in the GP (Ndim x Nnodes) */
  Matrix F_GP; /* Deformation gradient of the transformation evaluated in the GP (Ndim x Ndim) */
  Matrix F_GP_m1; /* Inverse of the deformation gradient */
  Matrix F_GP_Tm1; /* Transpose of the deformation gradient */
  Matrix dNdx_GP; /* Derivatives of the shape function evaluates in the GP (Ndim x Ndim) */

  /* 1º Evaluate the gradient of the shape function in the GP */
  dNdX_Ref_GP = dT3(X_EC_GP);
  
  /* 2º Get the deformation gradient of the transformation evaluated in the GP */
  F_GP = Get_F_Ref_T3(X_EC_GP,Element);
    
  /* 3º Get the inverse of the deformation gradient */
  F_GP_m1 = Get_Inverse(F_GP);
  FreeMat(F_GP);
  
  /* 4º Get the transpose of the inverse of the deformation gradient */
  F_GP_Tm1 = Transpose_Mat(F_GP_m1);
  FreeMat(F_GP_m1);
  
  /* 5º Get the gradient of the shape functions in global coordinates */
  dNdx_GP = Scalar_prod(F_GP_Tm1,dNdX_Ref_GP);
  FreeMat(F_GP_Tm1);
  FreeMat(dNdX_Ref_GP);

  /* 6º Return result */
  return dNdx_GP;
}

/*********************************************************************/

/* Global coordinates of the four nodes quadrilateral */
Matrix Get_X_GC_T3(Matrix X_NC_GP,Matrix X_GC_Nodes)
/*
This function evaluate the position of the GP in the element, and get it global coordiantes    
 */
{
  /* 0º Variable declaration */
  Matrix N_ref;
  Matrix X_GC_GP;

  /* 1º Evaluate the Q4 element in the element coordinates */
  N_ref = T3(X_NC_GP);

  /* 2º Allocate the output coordinates */
  X_GC_GP = MatAlloc(2,1);

  /* 3º Get the global coordinates for this element coordiantes in this element */
  X_GC_GP.nV[0] =
    N_ref.nV[0]*X_GC_Nodes.nM[0][0] +
    N_ref.nV[1]*X_GC_Nodes.nM[1][0] +
    N_ref.nV[2]*X_GC_Nodes.nM[2][0];
  
  X_GC_GP.nV[1] =
    N_ref.nV[0]*X_GC_Nodes.nM[0][1] +
    N_ref.nV[1]*X_GC_Nodes.nM[1][1] +
    N_ref.nV[2]*X_GC_Nodes.nM[2][1];

  /* 4º Free memory */
  FreeMat(N_ref);

  /* 5º Output */
  return X_GC_GP;
 
}

/*********************************************************************/

void Get_X_EC_T3(Matrix X_EC_GP,
		 Matrix X_GC_GP,
		 Matrix Element_GC_Nod)
/* 
   The function return the natural coordinates of a point 
   inside of the element.
 
   Inputs :
   - Coordinates of the element nodes
   - Initial coordinate of the point to start the search 
   - Derivative function of the element

   Depending of the kind of element, we employ differents types
   of shape functions
*/
{
  
  X_EC_GP = Newton_Rapson(Get_X_GC_T3,Element_GC_Nod,
			  Get_F_Ref_T3,Element_GC_Nod,
			  X_GC_GP,X_EC_GP);
}

/*********************************************************************/
