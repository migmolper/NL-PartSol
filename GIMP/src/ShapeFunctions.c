#include <stdio.h>
#include <stdlib.h>
#include "include/Utils.h"
#include "include/TypeDefinitions.h"

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
double * Q4(Vector * X_e){
  
  double * n = (double *)Allocate_Array(4,sizeof(double));
  
  n[0] = 1 - X_e->n[0] - X_e->n[1] + (X_e->n[0])*(X_e->n[1]);
  n[1] = X_e->n[0] - (X_e->n[0])*(X_e->n[1]);
  n[2] = (X_e->n[0])*(X_e->n[1]);
  n[3] = X_e->n[1] - X_e->n[1]*(X_e->n[0]);
  
  return n;
}

/* Derivatives of the shape functions */
double ** dQ4(Vector * X_e){
  
  double ** dn = (double **)Allocate_Matrix(4,2,sizeof(double));

  /* $\xi$ */
  dn[0][0] = X_e->n[1] - 1; /* \frac{\partial N0}{\partial \xi} */
  dn[1][0] = 1 - X_e->n[1]; /* \frac{\partial N1}{\partial \xi} */
  dn[2][0] = X_e->n[1]; /* \frac{\partial N2}{\partial \xi} */
  dn[3][0] = - X_e->n[1]; /* \frac{\partial N3}{\partial \xi} */
  /* $\eta$ */
  dn[0][1] = X_e->n[0] - 1; /* \frac{\partial N0}{\partial \eta} */
  dn[1][1] = - X_e->n[0]; /* \frac{\partial N1}{\partial \eta} */
  dn[2][1] = X_e->n[0]; /* \frac{\partial N2}{\partial \eta} */
  dn[3][1] = 1 - X_e->n[0]; /* \frac{\partial N3}{\partial \eta} */
  
  return dn;
}



