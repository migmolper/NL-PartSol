#include <stdio.h>
#include <stdlib.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"


/***********************************************/
/*** Four nodes cuadrilateral lineal element ***/
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
  N_ref.nV[0] = 1 - X_e.nV[0] - X_e.nV[1] + (X_e.nV[0])*(X_e.nV[1]);
  N_ref.nV[1] = X_e.nV[0] - (X_e.nV[0])*(X_e.nV[1]);
  N_ref.nV[2] = (X_e.nV[0])*(X_e.nV[1]);
  N_ref.nV[3] = X_e.nV[1] - X_e.nV[1]*(X_e.nV[0]);
  
  return N_ref;
}

/* Derivatives of the shape functions */
Matrix dQ4(Matrix X_e){
  
  /* Definition and allocation */
  Matrix dNdX_ref = MatAlloc(2,4);
  
  /* Fill the matrix */
  /* $\xi$ */
  dNdX_ref.nM[0][0] = X_e.nV[1] - 1; /* \frac{\partial N0}{\partial \xi} */
  dNdX_ref.nM[0][1] = 1 - X_e.nV[1]; /* \frac{\partial N1}{\partial \xi} */
  dNdX_ref.nM[0][2] = X_e.nV[1]; /* \frac{\partial N2}{\partial \xi} */
  dNdX_ref.nM[0][3] = - X_e.nV[1]; /* \frac{\partial N3}{\partial \xi} */
  /* $\eta$ */
  dNdX_ref.nM[1][0] = X_e.nV[0] - 1; /* \frac{\partial N0}{\partial \eta} */
  dNdX_ref.nM[1][1] = - X_e.nV[0]; /* \frac{\partial N1}{\partial \eta} */
  dNdX_ref.nM[1][2] = X_e.nV[0]; /* \frac{\partial N2}{\partial \eta} */
  dNdX_ref.nM[1][3] = 1 - X_e.nV[0]; /* \frac{\partial N3}{\partial \eta} */
  
  return dNdX_ref;
}

