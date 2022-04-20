#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Macros.h"

/*********************************************************************/

int Jacobi_Conjugate_Gradient_Method(
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
  double * lumped_Tangent_Stiffness;
  double Norm_r;
  double * r_k;
  double * r_k1;
  double * z_k;
  double * z_k1;
  double * p;
  int STOP_CRITERIA = 0;
  int Num_Iter, Num_Iter_Max;

  Num_Iter = 0;
  Num_Iter_Max = 25;
  Tol_r = 0.0000001;

  /* Allocate the residual arrays and the basis p array */
  r_k = (double *)calloc(Order, __SIZEOF_DOUBLE__);
  r_k1 = (double *)calloc(Order, __SIZEOF_DOUBLE__);
  z_k = (double *)calloc(Order, __SIZEOF_DOUBLE__);
  z_k1 = (double *)calloc(Order, __SIZEOF_DOUBLE__);
  p = (double *)calloc(Order, __SIZEOF_DOUBLE__);
  if((r_k == NULL)
  || (r_k1 == NULL)
  || (z_k == NULL)
  || (z_k1 == NULL)
  || (p == NULL)){
        fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
        return EXIT_FAILURE;
  } 

  // First trial
  for (unsigned i = 0; i < Order; i++) {

    /* Calculate the residual $r^1 = Residual - A x^1$ */
    aux = 0;
    for (unsigned j = 0; j < Order; j++) {
      aux += Tangent_Stiffness[i*Order + j] * U[j];
    }
    r_k[i] = Residual[i] - aux;

    z_k[i] = (1.0 / Tangent_Stiffness[i*Order + i]) * r_k[i];

    p[i] = z_k[i];
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
    free(z_k);
    free(z_k1);
    free(p);   
    return STATUS;
  }

  /* Get the Lumped-Mass matrix */
  lumped_Tangent_Stiffness = (double *)calloc(Order*Order, __SIZEOF_DOUBLE__);
  if(lumped_Tangent_Stiffness == NULL){
        fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
        return EXIT_FAILURE;
  } 

//  lumped_Tangent_Stiffness = lumped__MatrixLib__(Tangent_Stiffness);

  while (Norm_r > Tol_r) {

    // Calcule the stopping criteria
    if (Num_Iter > Num_Iter_Max) {
      printf("%s : %s \n \t %s : %f \n",
               "Warning in Jacobi_Conjugate_Gradient_Method",
               "Maximum number of iterations", "Norm of r", Norm_r);
      break;
    }      

    /* 1th step : Get alpha */
    alpha_k = 0;
    dividend = 0;
    divisor = 0;
    for (int i = 0; i < Order; i++) {
      /* Scalar product -> dividend = r_k^T \cdot r_k */
      dividend += r_k[i] * z_k[i];
      /* Scalar product -> divisor = r_k^T \cdot Tangent_Stiffness \cdot r_k */
      aux = 0;
      for (int j = 0; j < Order; j++) {
        aux += Tangent_Stiffness[i*Order + j] * p[j];
      }
      divisor += p[i] * aux;
    }
    alpha_k = dividend / divisor;

    for (int i = 0; i < Order; i++) {
      /* 2th Step : Get the solution array (Update) */
      U[i] += alpha_k * p[i];
      /* 3th Step : Calcule the residual */
      aux = 0;
      for (int j = 0; j < Order; j++) {
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


    beta_k = 0;
    dividend = 0;
    divisor = 0;
    for (int i = 0; i < Order; i++) {
      /* 4th step : */
      z_k1[i] = ((double)1 / lumped_Tangent_Stiffness[i]) * r_k1[i];
      /* 5th step :*/
      /* Scalar product -> dividend = r_{k+1}^T \cdot r_{k+1} */
      dividend += z_k1[i] * r_k1[i];
      /* Scalar product -> divisor = r_k^T \cdot r_k */
      divisor += z_k[i] * r_k[i];
    }
    beta_k = dividend / divisor;

    for (int i = 0; i < Order; i++) {
      /* 6th step : Update the basis vector $p$ */
      p[i] = z_k1[i] + beta_k * p[i];

      /* 7th step : Update the residual $r_{k}$ and the $z_k$ vector  */
      r_k[i] = r_k1[i];
      z_k[i] = z_k1[i];
    }
  }

  /* Free memory */
  free(r_k);
  free(r_k1);
  free(z_k);
  free(z_k1);
  free(p);
  free(lumped_Tangent_Stiffness);

  return STATUS;
}

/*********************************************************************/