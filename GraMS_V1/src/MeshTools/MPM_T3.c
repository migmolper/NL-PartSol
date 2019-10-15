#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../MathTools/MathTools.h"

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


/*********************************************************************/
