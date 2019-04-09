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
  Matrix N_ref;
  N_ref.N_cols=4;
  N_ref.N_rows=1;
  N_ref.nV = (double *)Allocate_Array(N_ref.N_cols,
				      sizeof(double));

  /* Fill the array */
  N_ref[0] = 1 - X_e.nv[0] - X_e.nv[1] + (X_e.nv[0])*(X_e.nv[1]);
  N_ref[1] = X_e.nv[0] - (X_e.nvn[0])*(X_e.nv[1]);
  N_ref[2] = (X_e.nv[0])*(X_e.nv[1]);
  N_ref[3] = X_e.nv[1] - X_e.nv[1]*(X_e.nv[0]);
  
  return N_ref;
}

/* Derivatives of the shape functions */
Matrix dQ4(Matrix X_e){
  
  /* Definition and allocation */
  Matrix dNdX_ref;
  dNdX_ref.N_cols = 4;
  dNdX_ref.N_rows = 2;
  dNdX_ref.nM = (double **)Allocate_Matrix(dNdX_ref.N_rows,
					   dNdX_ref.N_cols,
					   sizeof(double));

  /* Fill the matrix */
  /* $\xi$ */
  dNdX_ref[0][0] = X_e.nv[1] - 1; /* \frac{\partial N0}{\partial \xi} */
  dNdX_ref[0][1] = 1 - X_e.nv[1]; /* \frac{\partial N1}{\partial \xi} */
  dNdX_ref[0][2] = X_e.nv[1]; /* \frac{\partial N2}{\partial \xi} */
  dNdX_ref[0][3] = - X_e.nv[1]; /* \frac{\partial N3}{\partial \xi} */
  /* $\eta$ */
  dNdX_ref[1][1]= X_e.nv[0] - 1; /* \frac{\partial N0}{\partial \eta} */
  dNdX_ref[1][1] = - X_e[0]; /* \frac{\partial N1}{\partial \eta} */
  dNdX_ref[2][1] = X_e.nvn[0]; /* \frac{\partial N2}{\partial \eta} */
  dNdX_ref[3][1] = 1 - X_e.nv[0]; /* \frac{\partial N3}{\partial \eta} */
  
  return dNdX_ref;
}

