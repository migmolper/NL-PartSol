#include "nl-partsol.h"

/*
  Auxiliar functions
 */
static Matrix F_Ref__T3__(Matrix,Matrix);
static Matrix dN_Ref__T3__(Matrix);
static Matrix Xi_to_X__T3__(Matrix,Matrix);
static void   X_to_Xi__T3__(Matrix,Matrix,Matrix);

/*********************************************************************/

/*
  Shape functions 
*/
Matrix N__T3__(Matrix X_e){

  /* Error check */
  if( (fabs(X_e.nV[0]) > 1 ) ||
      (fabs(X_e.nV[1]) > 1 ) ||
      (1 - X_e.nV[0] - X_e.nV[1] < 0)){
    printf("Error in N__T3__() : Out of the element bounds !!! \n");
    exit(EXIT_FAILURE);
  }    
  
  /* Definition and allocation */
  Matrix N_ref =  alloc__MatrixLib__(1,3);

  /* Fill the array */
  N_ref.nV[0] = 1 - X_e.nV[0] - X_e.nV[1]; 
  N_ref.nV[1] = X_e.nV[0];
  N_ref.nV[2] = X_e.nV[1];
  
  return N_ref;
}

/*********************************************************************/

static Matrix dN_Ref__T3__(Matrix X_e)
{
  
  /* Error check */
  if( (fabs(X_e.nV[0]) > 1 ) ||
      (fabs(X_e.nV[1]) > 1 ) ||
      (1 - X_e.nV[0] - X_e.nV[1] < 0)){
    printf("Error in dN_Ref__T3__() : Out of the element bounds !!! \n");
    exit(EXIT_FAILURE);
  }
  
  /* Definition and allocation */
  Matrix dNdX_ref = alloc__MatrixLib__(2,3);
  
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
Matrix F_Ref__T3__(Matrix X_NC_GP,Matrix X_GC_Nodes)
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
  Matrix X_alpha = alloc__MatrixLib__(2,1);
  Matrix dNdx_alpha = alloc__MatrixLib__(1,2);
  Matrix F_Ref_alpha;
  Matrix F_Ref = allocZ__MatrixLib__(2,2);

  /* 1º Evaluate the derivarive of the shape function in the GP */
  dNdX_Ref_GP = dN_Ref__T3__(X_NC_GP); 

  /* 2º Get the F_Ref doing a loop over the nodes of the element */
  for(int i = 0 ; i<3 ; i++){

    /* 3º Fill arrays for the tensorial product */
    for(int j = 0 ; j<2 ; j++){
      X_alpha.nV[j] = X_GC_Nodes.nM[i][j];
      dNdx_alpha.nV[j] = dNdX_Ref_GP.nM[j][i];
    }

    /* 4º Get the nodal contribution */
    F_Ref_alpha = dyadic_product__MatrixLib__(X_alpha,dNdx_alpha);

    /* 5º Increment the reference deformation gradient */
    F_Ref = increment__MatrixLib__(F_Ref, F_Ref_alpha);

    /* 6º Free data of the nodal contribution */
    free__MatrixLib__(F_Ref_alpha);
    
  }
  
  /* 7º Free memory */
  free__MatrixLib__(dNdX_Ref_GP);
  free__MatrixLib__(X_alpha);
  free__MatrixLib__(dNdx_alpha);

  /* 8º Output */
  return F_Ref;
}

/*********************************************************************/

/* Element gradient in the real element */
Matrix dN__T3__(Matrix X_EC_GP,Matrix Element)
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
  dNdX_Ref_GP = dN_Ref__T3__(X_EC_GP);
  
  /* 2º Get the deformation gradient of the transformation evaluated in the GP */
  F_GP = F_Ref__T3__(X_EC_GP,Element);
    
  /* 3º Get the inverse of the deformation gradient */
  F_GP_m1 = inverse__MatrixLib__(F_GP);
  free__MatrixLib__(F_GP);
  
  /* 4º Get the transpose of the inverse of the deformation gradient */
  F_GP_Tm1 = transpose__MatrixLib__(F_GP_m1);
  free__MatrixLib__(F_GP_m1);
  
  /* 5º Get the gradient of the shape functions in global coordinates */
  dNdx_GP = matrix_product__MatrixLib__(F_GP_Tm1,dNdX_Ref_GP);
  free__MatrixLib__(F_GP_Tm1);
  free__MatrixLib__(dNdX_Ref_GP);

  /* 6º Return result */
  return dNdx_GP;
}

/*********************************************************************/

/* Global coordinates of the four nodes quadrilateral */
static Matrix Xi_to_X__T3__(Matrix X_NC_GP,Matrix X_GC_Nodes)
/*
This function evaluate the position of the GP in the element, and get it global coordiantes    
 */
{
  /* 0º Variable declaration */
  Matrix N_ref;
  Matrix X_GC_GP;

  /* 1º Evaluate the Q4 element in the element coordinates */
  N_ref = N__T3__(X_NC_GP);

  /* 2º Allocate the output coordinates */
  X_GC_GP = alloc__MatrixLib__(2,1);

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
  free__MatrixLib__(N_ref);

  /* 5º Output */
  return X_GC_GP;
 
}

/*********************************************************************/

void X_to_Xi__T3__(Matrix X_EC_GP,
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
  
  X_EC_GP = Newton_Rapson__MatrixSolvers__(Xi_to_X__T3__,Element_GC_Nod,
			  F_Ref__T3__,Element_GC_Nod,
			  X_GC_GP,X_EC_GP);
}

/*********************************************************************/
