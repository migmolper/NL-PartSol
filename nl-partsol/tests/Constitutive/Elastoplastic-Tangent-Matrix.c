#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

#define RESET "\033[0m"
#define RED "\033[31m"

/**************************************************************/
/******************* Material Parameters **********************/
/**************************************************************/
#define NumberDimensions 2

/**************************************************************/

typedef struct {
  /*!
   * Particle identifier
   * */
  int Particle_Idx;

  /*!
   * Stress/strain parameters
   * */
  double *Stress;
  double *Strain;
  double Pressure;

  /*!
   * Finite strain kinematic parameters
   * */
  double *d_phi;
  double *D_phi;
  double *Fbar;
  double *rate_D_phi;
  double J;

  /*!
   * Plasticity parameters
   * */
  double *Back_stress;
  double *b_e;

  double Cohesion;
  double Yield_stress;

  // Internal hardening variables
  double *Equiv_Plast_Str; // Equivalent plastic strain
  double *Kappa;           // Hardening Parameter
  double *a_ep;            // Elastoplastic tangent matrix
  //
  bool *Failure;

} State_Parameters;

static int
__compute_u_v(double *u /**< [out] component-free auxiliar variable */,
              double *v /**< [out] component-free auxiliar variable */,
              const double *dN_alpha /**< [in] grad shape function (alpha) */,
              const double *dN_beta /**< [in] grad shape function (beta) */,
              const double *D_phi /**< [in] Total deformation gradient. */);

static int __spectral_decomposition_b_e(
    double *eigval_b_e /**< [out] Eigenvalues of b elastic trial. */,
    double *eigvec_b_e /**< [out] Eigenvector of b elastic trial. */,
    const double *b_e /**< [in] (n) Elastic left Cauchy-Green.*/);

static int __eigenvalues_kirchhoff(
    double *eigval_T /**< [out] Eigenvalues of the Kirchhoff stress tensor. */,
    const double *P /**< [in] Nominal stress tensor */,
    const double *D_phi /**< [in] Total deformation gradient. */);

int compute_1PK_elastoplastic_tangent_matrix(double *A_ep,
                                             const double *dN_alpha,
                                             const double *dN_beta,
                                             const State_Parameters IO_State);

/**************************************************************/

int main() {

  int STATUS = EXIT_SUCCESS;

  double a_ep[4] = {10.11, 35.3, 35.3, 47.05};
  double A_ep[4] = {0.0, 0.0, 0.0, 0.0};
  double dN_alpha[2] = {0.35, 0.42};
  double dN_beta[2] = {0.1, 0.8};
  double tau[4] = {200, 345, 345, 652};
  double FPK[4];
  double D_phi[4] = {1.2, 0.5, 0.7, 0.8};
  double b_e[4] = {0.8, 0.6, 0.6, 0.1};

  double D_phi_mT[4] = {1.2, 0.7, 0.5, 0.8};

  int INFO;
  int N = 2;
  int LDA = 2;
  int LWORK = 2;
  int IPIV[2] = {0, 0};
  double WORK[2] = {0, 0};

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N, &N, D_phi_mT, &LDA, IPIV, &INFO);
  // Check output of dgetrf
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(
          stderr,
          "" RED
          "Error in dgetrf_(): the %i-th argument had an illegal value" RESET
          "",
          abs(INFO));
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetrf_(): D_phi_mT(%i,%i) %s \n %s \n %s \n %s" RESET
              " \n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor D_phi_mT is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  dgetri_(&N, D_phi_mT, &LDA, IPIV, WORK, &LWORK, &INFO);
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(stderr,
              "" RED "Error in dgetri_(): the %i-th argument of dgetrf_ had an "
              "illegal value" RESET "\n",
              abs(INFO));
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetri_(): D_phi_mT(%i,%i) %s \n %s \n %s \n %s " RESET
              "\n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor D_phi_mT is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  FPK[0] = tau[0] * D_phi_mT[0] + tau[1] * D_phi_mT[2];
  FPK[1] = tau[0] * D_phi_mT[1] + tau[1] * D_phi_mT[3];
  FPK[2] = tau[2] * D_phi_mT[0] + tau[3] * D_phi_mT[2];
  FPK[3] = tau[2] * D_phi_mT[1] + tau[3] * D_phi_mT[3];

  State_Parameters IO_State;
  IO_State.Stress = FPK;
  IO_State.D_phi = D_phi;
  IO_State.b_e = b_e;
  IO_State.a_ep = a_ep;

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  puts("FPK: ");
  printf("%e, %e \n", FPK[0], FPK[1]);
  printf("%e, %e \n", FPK[2], FPK[3]);

  puts("D_phi_mT: ");
  printf("%e, %e \n", D_phi_mT[0], D_phi_mT[1]);
  printf("%e, %e \n", D_phi_mT[2], D_phi_mT[3]);
#endif
#endif

  STATUS = compute_1PK_elastoplastic_tangent_matrix(A_ep, dN_alpha, dN_beta,
                                                    IO_State);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,
            "" RED "Error in compute_1PK_elastoplastic_tangent_matrix" RESET
            "\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**************************************************************/
int compute_1PK_elastoplastic_tangent_matrix(double *A_ep,
                                             const double *dN_alpha,
                                             const double *dN_beta,
                                             const State_Parameters IO_State) {

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;

#if NumberDimensions == 2
  double u[2] = {0.0, 0.0};
  double v[2] = {0.0, 0.0};
#else
  No esta implementado
#endif

  STATUS = __compute_u_v(u, v, dN_alpha, dN_beta, IO_State.D_phi);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __compute_u_v" RESET "\n");
    return EXIT_FAILURE;
  }

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

#if NumberDimensions == 2
  double eigval_T[2] = {0.0, 0.0};
#else
  double eigval_T[3] = {0.0, 0.0, 0.0};
#endif

  STATUS = __eigenvalues_kirchhoff(eigval_T, IO_State.Stress, IO_State.D_phi);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __eigenvalues_kirchhoff" RESET "\n");
    return EXIT_FAILURE;
  }

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  printf("eigval_T: [%e, %e] \n", eigval_T[0], eigval_T[1]);
#endif
#endif

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
          mv[A * Ndim + B][i] += m[A * Ndim + B][i * Ndim + j] * v[j];
          mu[A * Ndim + B][i] += m[A * Ndim + B][i * Ndim + j] * u[j];
        }
      }
    }
  }

  // Do the diadic product of gradient directions
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      u__o__v[i][j] = u[i] * v[j];
    }
  }

  // Assemble the material contribution to the tanget matrix
  for (unsigned A = 0; A < Ndim; A++) {
    for (unsigned B = 0; B < Ndim; B++) {

      for (unsigned i = 0; i < Ndim; i++) {
        for (unsigned j = 0; j < Ndim; j++) {
          A_ep[i * Ndim + j] += IO_State.a_ep[A * Ndim + B] *
                                (mv[A * Ndim + A][i] * mu[B * Ndim + B][j]);

          if (A != B) {
            A_ep[i * Ndim + j] +=
                0.5 *
                ((eigval_T[B] - eigval_T[A]) /
                 (eigval_b_e[B] - eigval_b_e[A])) *
                (eigval_b_e[B] * (mv[A * Ndim + B][i] * mu[A * Ndim + B][j]) +
                 eigval_b_e[A] * (mv[A * Ndim + B][i] * mu[B * Ndim + A][j]));
          }
        }
      }
    }
  }

  // Assemble the geometrical contribution to the tanget matrix
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      for (unsigned k = 0; k < Ndim; k++) {
        A_ep[i * Ndim + j] += -IO_State.Stress[i * Ndim + k] * u__o__v[k][j];
      }
    }
  }

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  puts("A_ep: ");
  printf("%e, %e\n", A_ep[0], A_ep[1]);
  printf("%e, %e\n", A_ep[2], A_ep[3]);
#endif
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __compute_u_v(double *u, double *v, const double *dN_alpha,
                         const double *dN_beta, const double *D_phi) {

#if NumberDimensions == 2

  double D_phi_mT[4] = {0.0, 0.0, 0.0, 0.0};

  D_phi_mT[0] = D_phi[0];
  D_phi_mT[1] = D_phi[2];
  D_phi_mT[2] = D_phi[1];
  D_phi_mT[3] = D_phi[3];

  // compute the inverse of D_phi
  int INFO;
  int N = 2;
  int LDA = 2;
  int LWORK = 2;
  int IPIV[2] = {0, 0};
  double WORK[2] = {0, 0};

#else
  No esta implementado
#endif

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N, &N, D_phi_mT, &LDA, IPIV, &INFO);
  // Check output of dgetrf
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(
          stderr,
          "" RED
          "Error in dgetrf_(): the %i-th argument had an illegal value" RESET
          "",
          abs(INFO));
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetrf_(): D_phi_mT(%i,%i) %s \n %s \n %s \n %s" RESET
              " \n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor D_phi_mT is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  dgetri_(&N, D_phi_mT, &LDA, IPIV, WORK, &LWORK, &INFO);
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(stderr,
              "" RED "Error in dgetri_(): the %i-th argument of dgetrf_ had an "
              "illegal value" RESET "\n",
              abs(INFO));
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetri_(): D_phi_mT(%i,%i) %s \n %s \n %s \n %s " RESET
              "\n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor D_phi_mT is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

#if NumberDimensions == 2

  u[0] = dN_beta[0];
  u[1] = dN_beta[1];

  v[0] = D_phi_mT[0] * dN_alpha[0] + D_phi_mT[1] * dN_alpha[1];
  v[1] = D_phi_mT[2] * dN_alpha[0] + D_phi_mT[3] * dN_alpha[1];

#else
  No esta implementado
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
  No esta implementado
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

static int __eigenvalues_kirchhoff(double *eigval_T, const double *P,
                                   const double *D_phi) {

#if NumberDimensions == 2
  double T[4] = {0.0, 0.0, 0.0, 0.0};
  T[0] = P[0] * D_phi[0] + P[1] * D_phi[1];
  T[1] = P[0] * D_phi[2] + P[1] * D_phi[3];
  T[2] = P[2] * D_phi[0] + P[3] * D_phi[1];
  T[3] = P[2] * D_phi[2] + P[3] * D_phi[3];
#else
  double T[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  No esta implementado
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
  dsyev_("N", "L", &n, T, &lda, eigval_T, &wkopt, &lwork, &info);
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

  dsyev_("N", "L", &n, T, &lda, eigval_T, work, &lwork, &info);
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