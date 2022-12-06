
/**
 * @file Elastoplastic-Tangent-Matrix.c
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Constitutive/Plasticity/Elastoplastic-Tangent-Matrix.h"
#include "Globals.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/**************************************************************/

/*!
    \param[out] eigval_b_e Eigenvalues of b elastic trial.
    \param[out] eigvec_b_e Eigenvector of b elastic trial.
    \param[in] b_e (n) Elastic left Cauchy-Green.
*/
static int __spectral_decomposition_b_e(double *eigval_b_e, double *eigvec_b_e,
                                        const double *b_e);
/**************************************************************/

/*!
   \param[out] eigval_T Eigenvalues of the Kirchhoff stress tensor.
   \param[in] T Kirchhoff stress tensor
*/
static int __eigenvalues_kirchhoff(double *eigval_T, const double *T);
/**************************************************************/

int compute_stiffness_elastoplastic__Constitutive__(double *Stiffness_density,
                                             const double *dN_alpha_n1,
                                             const double *dN_beta_n1,
                                             const State_Parameters IO_State) {

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;

#if NumberDimensions == 2
  Stiffness_density[0] = 0.0;
  Stiffness_density[1] = 0.0;
  Stiffness_density[2] = 0.0;
  Stiffness_density[3] = 0.0;
#else
  Stiffness_density[0] = 0.0;
  Stiffness_density[1] = 0.0;
  Stiffness_density[2] = 0.0;
  Stiffness_density[3] = 0.0;
  Stiffness_density[4] = 0.0;
  Stiffness_density[5] = 0.0;
  Stiffness_density[6] = 0.0;
  Stiffness_density[7] = 0.0;
  Stiffness_density[8] = 0.0;
#endif

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  printf("u: [%e, %e] \n", u[0], u[1]);
  printf("v: [%e, %e] \n", v[0], v[1]);
#endif
#endif

#if NumberDimensions == 2
  double eigval_b_e[2] = {0.0, 0.0};
  double eigvec_b_e[4] = {0.0, 0.0, 0.0, 0.0};
  double u__o__v[2][2];
#else
  double eigval_b_e[3] = {0.0, 0.0, 0.0};
  double eigvec_b_e[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double u__o__v[3][3];
#endif

  STATUS = __spectral_decomposition_b_e(eigval_b_e, eigvec_b_e, IO_State.b_e);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __spectral_decomposition_b_e" RESET "\n");
    return EXIT_FAILURE;
  }

  double eigval_T[3] = {0.0, 0.0, 0.0};

  STATUS = __eigenvalues_kirchhoff(eigval_T, IO_State.Stress);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __eigenvalues_kirchhoff" RESET "\n");
    return EXIT_FAILURE;
  }


  // Do the diadic product of gradient directions
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      u__o__v[i][j] = dN_beta_n1[i] * dN_alpha_n1[j];
    }
  }


  for (unsigned A = 0; A < Ndim; A++) {

    double u_A = 0.0;
    double v_A = 0.0;
    for (unsigned i = 0; i < Ndim; i++) 
    {
      u_A += dN_alpha_n1[i]*eigvec_b_e[A + i * Ndim];
      v_A += dN_beta_n1[i]*eigvec_b_e[A + i * Ndim];
    }

    for (unsigned B = 0; B < Ndim; B++) {

      double u_B = 0.0;
      double v_B = 0.0;
      for (unsigned i = 0; i < Ndim; i++) 
      {
        u_B += dN_alpha_n1[i]*eigvec_b_e[B + i * Ndim];        
        v_B += dN_beta_n1[i]*eigvec_b_e[B + i * Ndim];
      }

      double C_ep_AB = IO_State.C_ep[A * Ndim + B];
      double v_A__dot__u_B = u_B*v_A;
      double u_A__dot__v_B = u_A*v_B;    
      double u_B__dot__v_B = u_B*v_B;
      
      for (unsigned i = 0; i < Ndim; i++) {
        for (unsigned j = 0; j < Ndim; j++) {
          Stiffness_density[i * Ndim + j] +=
              C_ep_AB * u_A__dot__v_B * eigvec_b_e[A + i * Ndim] * eigvec_b_e[B + j * Ndim];

          if (A != B) {
            if (fabs(eigval_b_e[B] - eigval_b_e[A]) > 1E-14) {
              Stiffness_density[i * Ndim + j] +=
                  0.5 *
                  ((eigval_T[B] - eigval_T[A]) /
                   (eigval_b_e[B] - eigval_b_e[A])) *
                  (eigval_b_e[B] * u_B__dot__v_B*(eigvec_b_e[A + i * Ndim] * eigvec_b_e[A + j * Ndim]) +
                   eigval_b_e[A] * v_A__dot__u_B*(eigvec_b_e[A + i * Ndim] * eigvec_b_e[B + j * Ndim]));
            } 
          }          
        }
      }
    }
  }

  // Assemble the geometrical contribution to the tanget matrix
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      for (unsigned k = 0; k < Ndim; k++) {
        Stiffness_density[i * Ndim + j] +=
            -IO_State.Stress[i * Ndim + k] * u__o__v[k][j];
      }
    }
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __spectral_decomposition_b_e(double *eigval_b_e, double *eigvec_b_e,
                                        const double *b_e) {

  lapack_int n = NumberDimensions;
  lapack_int lda = NumberDimensions;

#if NumberDimensions == 2
  eigvec_b_e[0] = b_e[0];
  eigvec_b_e[1] = b_e[1];
  eigvec_b_e[2] = b_e[2];
  eigvec_b_e[3] = b_e[3];
#else
  eigvec_b_e[0] = b_e[0];
  eigvec_b_e[1] = b_e[1];
  eigvec_b_e[2] = b_e[2];
  eigvec_b_e[3] = b_e[3];
  eigvec_b_e[4] = b_e[4];
  eigvec_b_e[5] = b_e[5];
  eigvec_b_e[6] = b_e[6];
  eigvec_b_e[7] = b_e[7];
  eigvec_b_e[8] = b_e[8];
#endif
 
  lapack_int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec_b_e, lda, eigval_b_e);

  if (info > 0) {
    fprintf(stderr,
            "" RED "Error in dsyev_(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    fprintf(stderr,
            "" RED "Error in dsyev_(): the %i-th argument had an "
            "illegal value." RESET "\n",
            abs(info));
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __eigenvalues_kirchhoff(double *eigval_T, const double *T) {

  unsigned Ndim = NumberDimensions;
  lapack_int n = NumberDimensions;
  lapack_int lda = NumberDimensions;

#if NumberDimensions == 2
  double T_aux[4] = {T[0], T[1], T[2], T[3]};
#else
  double T_aux[9] = {T[0], T[1], T[2], T[3], T[4], T[5], T[6], T[7], T[8]};
#endif

  lapack_int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, T_aux, lda, eigval_T);

  if (info > 0) {
    fprintf(stderr,
            "" RED "Error in dsyev_(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    fprintf(stderr,
            "" RED "Error in dsyev_(): the %i-th argument had an "
            "illegal value." RESET "\n",
            abs(info));
    return EXIT_FAILURE;
  }


#if NumberDimensions == 2
  eigval_T[2] = T[4];
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/
