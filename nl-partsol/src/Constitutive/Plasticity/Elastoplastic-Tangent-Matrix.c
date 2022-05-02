

#include "Constitutive/Plasticity/Elastoplastic-Tangent-Matrix.h"

/**************************************************************/ 

/*!
    \param[out] eigval_b_e Eigenvalues of b elastic trial.
    \param[out] eigvec_b_e Eigenvector of b elastic trial.
    \param[in] b_e (n) Elastic left Cauchy-Green.
*/
static int __spectral_decomposition_b_e(
    double *eigval_b_e /**< [out] Eigenvalues of b elastic trial. */,
    double *eigvec_b_e /**< [out] Eigenvector of b elastic trial. */,
    const double *b_e /**< [in] (n) Elastic left Cauchy-Green.*/);
/**************************************************************/

/*!
   \param[out] eigval_T Eigenvalues of the Kirchhoff stress tensor.
   \param[in] T Kirchhoff stress tensor
*/
static int __eigenvalues_kirchhoff(
    double *eigval_T,
    const double *T);
/**************************************************************/

int compute_1PK_elastoplastic_tangent_matrix(
  double *Stiffness_density,
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
#else
  double eigval_b_e[3] = {0.0, 0.0, 0.0};
  double eigvec_b_e[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif

  STATUS = __spectral_decomposition_b_e(eigval_b_e, eigvec_b_e, IO_State.b_e);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __spectral_decomposition_b_e" RESET "\n");
    return EXIT_FAILURE;
  }

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  printf("eigval_b_e: [%e, %e] \n", eigval_b_e[0], eigval_b_e[1]);
  puts("eigvec_b_e: ");
  printf("%e, %e \n", eigvec_b_e[0], eigvec_b_e[1]);
  printf("%e, %e \n", eigvec_b_e[2], eigvec_b_e[3]);
#endif
#endif


  double eigval_T[3] = {0.0, 0.0, 0.0};

  STATUS = __eigenvalues_kirchhoff(eigval_T, IO_State.Stress);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __eigenvalues_kirchhoff" RESET "\n");
    return EXIT_FAILURE;
  }


#if NumberDimensions == 2

  double n1[2] = {0.0, 0.0};
  double n2[2] = {0.0, 0.0};

  double m[4][4] = {
      {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},
  };

  double mu[4][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};

  double mv[4][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};

  double u__o__v[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

#else
  No esta implementado
#endif

  // Generate the matrix with the eigenbases
  // keeping in mind that LAPACK returns column-wise matrix
  for (unsigned A = 0; A < Ndim; A++) {
    for (unsigned B = 0; B < Ndim; B++) {
      for (unsigned i = 0; i < Ndim; i++) {
        for (unsigned j = 0; j < Ndim; j++) {
          m[A * Ndim + B][i * Ndim + j] =
              eigvec_b_e[A + i * Ndim] * eigvec_b_e[B + j * Ndim];
        }
      }
    }
  }

  // Do the projection mu and mv
  for (unsigned A = 0; A < Ndim; A++) {
    for (unsigned B = 0; B < Ndim; B++) {
      for (unsigned i = 0; i < Ndim; i++) {
        for (unsigned j = 0; j < Ndim; j++) {
          mv[A * Ndim + B][i] += m[A * Ndim + B][i * Ndim + j] * dN_alpha_n1[j];
          mu[A * Ndim + B][i] += m[A * Ndim + B][i * Ndim + j] * dN_beta_n1[j];
        }
      }
    }
  }

  // Do the diadic product of gradient directions
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      u__o__v[i][j] = dN_beta_n1[i] * dN_alpha_n1[j];
    }
  }

  // Assemble the material contribution to the tanget matrix
  for (unsigned A = 0; A < Ndim; A++) {
    for (unsigned B = 0; B < Ndim; B++) {

      for (unsigned i = 0; i < Ndim; i++) {
        for (unsigned j = 0; j < Ndim; j++) {
          Stiffness_density[i * Ndim + j] += IO_State.C_ep[A * Ndim + B] *
                                (mv[A * Ndim + A][i] * mu[B * Ndim + B][j]);

          if (A != B) {
            if(fabs(eigval_b_e[B] - eigval_b_e[A]) > 1E-10)
            {
              Stiffness_density[i * Ndim + j] +=
                  0.5 *
                  ((eigval_T[B] - eigval_T[A]) /
                  (eigval_b_e[B] - eigval_b_e[A])) *
                  (eigval_b_e[B] * (mv[A * Ndim + B][i] * mu[A * Ndim + B][j]) +
                  eigval_b_e[A] * (mv[A * Ndim + B][i] * mu[B * Ndim + A][j]));
            }
            else
            {
              Stiffness_density[i * Ndim + j] += 0.0;
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
        Stiffness_density[i * Ndim + j] += -IO_State.Stress[i * Ndim + k] * u__o__v[k][j];
      }
    }
  }

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  puts("Stiffness_density_p: ");
  printf("%e, %e\n", Stiffness_density[0], Stiffness_density[1]);
  printf("%e, %e\n", Stiffness_density[2], Stiffness_density[3]);
#endif
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __spectral_decomposition_b_e(double *eigval_b_e, double *eigvec_b_e,
                                        const double *b_e) {

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

  /* Locals */
  int n = NumberDimensions;
  int lda = NumberDimensions;
  int ldvl = NumberDimensions;
  int ldvr = NumberDimensions;
  int info;
  int lwork;
  double wkopt;
  double *work;

  /*
    Query and allocate the optimal workspace
  */
  lwork = -1;
  dsyev_("V", "L", &n, eigvec_b_e, &lda, eigval_b_e, &wkopt, &lwork, &info);

  /* Check for convergence */
  if (info > 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): the %i-th argument had an "
            "illegal value." RESET "\n",
            abs(info));
    return EXIT_FAILURE;
  }

  lwork = (int)wkopt;
  work = (double *)malloc(lwork * sizeof(double));

  dsyev_("V", "L", &n, eigvec_b_e, &lda, eigval_b_e, work, &lwork, &info);
  /* Check for convergence */
  if (info > 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): the %i-th argument had an "
            "illegal value." RESET "\n",
            abs(info));
    return EXIT_FAILURE;
  }

  free(work);

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __eigenvalues_kirchhoff(double *eigval_T, const double *T) {

unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  double T_aux[4] = {
    T[0], T[1],
    T[2], T[3]};
#else
  double T_aux[9] = {
    T[0], T[1], T[2],
    T[3], T[4], T[5],
    T[6], T[7], T[8]};
#endif


  /* Locals */
  int n = NumberDimensions;
  int lda = NumberDimensions;
  int ldvl = NumberDimensions;
  int ldvr = NumberDimensions;
  int info;
  int lwork;
  double wkopt;
  double *work;

  /* Local arrays */
#if NumberDimensions == 2
  int IPIV[2] = {0, 0};
  double wi[2];
  double vl[4];
#else
  int IPIV[3] = {0, 0, 0};
  double wi[3];
  double vl[9];
#endif

  /*
    Query and allocate the optimal workspace
  */
  lwork = -1;
  dsyev_("N", "L", &n, T_aux, &lda, eigval_T, &wkopt, &lwork, &info);
  lwork = (int)wkopt;
  work = (double *)malloc(lwork * sizeof(double));

  /* Check for convergence */
  if (info > 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): the %i-th argument had an "
            "illegal value." RESET "\n",
            abs(info));
    return EXIT_FAILURE;
  }

  dsyev_("N", "L", &n, T_aux, &lda, eigval_T, work, &lwork, &info);
  /* Check for convergence */
  if (info > 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): the %i-th argument had an "
            "illegal value." RESET "\n",
            abs(info));
    return EXIT_FAILURE;
  }

  free(work);

#if NumberDimensions == 2
  eigval_T[2] = T[4];
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/

