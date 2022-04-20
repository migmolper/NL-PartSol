#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Macros.h"

/*********************************************************************/

int Conjugate_Gradient_Method(
    double * Tangent_Stiffness, 
    double * Residual, 
    double * U,
    unsigned Nactivedofs)
{

  int STATUS = EXIT_SUCCESS;
  unsigned Order = Nactivedofs;
  double aux;
  double Tol_r;
  double alpha_k, beta_k;
  double dividend, divisor;
  double Norm_r;
  double * r_k;
  double * r_k1;
  double * p;
  unsigned Num_Iter, Num_Iter_Max;

  Num_Iter = 0;
  Num_Iter_Max = 25;
  Tol_r = 0.0001;

  /* Allocate the residual and the p basis arrays */
  r_k = (double *)calloc(Order, __SIZEOF_DOUBLE__);
  r_k1 = (double *)calloc(Order, __SIZEOF_DOUBLE__);
  p = (double *)calloc(Order, __SIZEOF_DOUBLE__);
  if((r_k == NULL)
  || (r_k1 == NULL)
  || (p == NULL)){
        fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
        return EXIT_FAILURE;
  } 
  
  // First trial
  for (unsigned i = 0; i < Order; i++) {

    aux = 0;

    for (unsigned j = 0; j < Order; j++) {
      aux += Tangent_Stiffness[i*Order + j] * U[j];
    }

    r_k[i] = Residual[i] - aux;

    p[i] = r_k[i];
  }

  // Compute the norm of the residual
  Norm_r = 0.0;
  for (unsigned i = 0; i < Order; i++) {
        Norm_r += DSQR(r_k[i]);
  }
  Norm_r = pow(Norm_r, 0.5);

  if (Norm_r < Tol_r) {
    free(r_k);
    free(r_k1);
    free(p);     
    return STATUS;
  }

  while (Norm_r > Tol_r) {

    // Calcule the stopping criteria
    if (Num_Iter > Num_Iter_Max) { 
      printf("%s : %s \n \t %s : %f \n",
               "Warning in Conjugate_Gradient_Method",
               "Maximum number of iterations", "Norm of r", Norm_r);
      break;
    }

    /* 1th step : Get alpha */
    alpha_k = 0;
    dividend = 0;
    divisor = 0;
    for (unsigned i = 0; i < Order; i++) {
      dividend += r_k[i] * r_k[i];
      aux = 0;
      for (unsigned j = 0; j < Order; j++) {
        aux += Tangent_Stiffness[i*Order + j] * p[j];
      }
      divisor += p[i] * aux;
    }
    alpha_k = dividend / divisor;

    for (unsigned i = 0; i < Order; i++) {
      U[i] += alpha_k * p[i];
      aux = 0;
      for (unsigned j = 0; j < Order; j++) {
        aux += Tangent_Stiffness[i*Order + j] * p[j];
      }
      r_k1[i] = r_k[i] - alpha_k * aux;
    }


    // Compute the norm of the residual (k iteration)
    Norm_r = 0.0;
    for (unsigned i = 0; i < Order; i++) {
        Norm_r += DSQR(r_k1[i]);
    }
    Norm_r = pow(Norm_r, 0.5);

    // 
    beta_k = 0;
    dividend = 0;
    divisor = 0;
    for (unsigned i = 0; i < Order; i++) {
      dividend += r_k1[i] * r_k1[i];
      divisor += r_k[i] * r_k[i];
    }
    beta_k = dividend / divisor;

    for (unsigned i = 0; i < Order; i++) {
      p[i] = r_k1[i] + beta_k * p[i];
      r_k[i] = r_k1[i];
    }

  }

  /* Free memory */
  free(r_k);
  free(r_k1);
  free(p);

  return STATUS;
}

/*********************************************************************/