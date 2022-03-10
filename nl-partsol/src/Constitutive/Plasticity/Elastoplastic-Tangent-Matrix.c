#include <math.h>
#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

static int __compute_u_v(
    double * u /**< [out] component-free auxiliar variable */,
    double * v /**< [out] component-free auxiliar variable */,
    const double * dN_alpha /**< [in] grad shape function (alpha) */,
    const double * dN_beta /**< [in] grad shape function (beta) */, 
    const double * D_phi /**< [in] Total deformation gradient. */);

static int __spectral_decomposition_b_e(
    double *eigval_b_e /**< [out] Eigenvalues of b elastic trial. */,
    double *eigvec_b_e /**< [out] Eigenvector of b elastic trial. */,
    const double *b_e /**< [in] (n) Elastic left Cauchy-Green.*/);

static int __eigenvalues_kirchhoff(
    double *eigval_T/**< [out] Eigenvalues of the Kirchhoff stress tensor. */,
    const double *P /**< [in] Nominal stress tensor */,
    const double *D_phi /**< [in] Total deformation gradient. */);

/**************************************************************/ 
int compute_1PK_elastoplastic_tangent_matrix(
    const double * dN_alpha,
    const double * dN_beta,
    State_Parameters IO_State)
{

  int STATUS = EXIT_SUCCESS;

#if NumberDimensions == 2

double u[2] = {0.0,0.0};
double v[2] = {0.0,0.0};

#else
  No esta implementado
#endif

  STATUS = __compute_u_v(u, v, dN_alpha,dN_beta, IO_State.D_phi);
  if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __compute_u_v" RESET "\n");
      return EXIT_FAILURE;
  }
  
  double eigval_T[3] = {0.0, 0.0, 0.0};
  double eigval_b_e[3] = {0.0, 0.0, 0.0};
  double eigvec_b_e[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  

  STATUS = __spectral_decomposition_b_e(eigval_b_e, eigvec_b_e, IO_State.b_e);
  if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __spectral_decomposition_b_e" RESET "\n");
      return EXIT_FAILURE;
  }

  STATUS = __eigenvalues_kirchhoff(eigval_T,IO_State.Stress,IO_State.D_phi);
  if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __eigenvalues_kirchhoff" RESET "\n");
      return EXIT_FAILURE;
  }

#if NumberDimensions == 2

 double nI[2] = {0.0,0.0};  
 double nII[2] = {0.0,0.0};
 double nI_nI[4] = {0.0,0.0,0.0};
 double nI_nII[4] = {0.0,0.0,0.0};
 double nII_nI[4] = {0.0,0.0,0.0};
 double nII_nII[4] = {0.0,0.0,0.0};
 
 nI[0] = eigvec_b_e[0];
 nI[1] = eigvec_b_e[1];
 nII[0] = eigvec_b_e[3];
 nII[1] = eigvec_b_e[4];

for(unsigned i = 0 ; i<2; i++)
{
    for(unsigned j = 0 ; j<2; j++)
    {
        nI_nI[i*2 + j] = nI[i]*nI[j];
        nI_nII[i*2 + j] = nI[i]*nII[j];
        nII_nI[i*2 + j] = nII[i]*nI[j];
        nII_nII[i*2 + j] = nII[i]*nII[j];
    }
}

#else
  No esta implementado
#endif


return EXIT_SUCCESS;
}

/**************************************************************/

static int __compute_u_v(double * u,double * v,const double * dN_alpha,
                         const double * dN_beta,const double * D_phi){


double D_phi_mT[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

#if NumberDimensions == 2

  D_phi_mT[0] = D_phi[0];
  D_phi_mT[1] = D_phi[2];
  D_phi_mT[3] = D_phi[1];
  D_phi_mT[4] = D_phi[3];
  D_phi_mT[8] = D_phi[4];


#else
  No esta implementado
#endif


  // compute the inverse of D_phi
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};

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

v[0] = D_phi_mT[0]*dN_alpha[0] + D_phi_mT[1]*dN_alpha[1];
v[1] = D_phi_mT[3]*dN_alpha[0] + D_phi_mT[4]*dN_alpha[1];

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
  eigvec_b_e[3] = b_e[2];
  eigvec_b_e[4] = b_e[3];
  eigvec_b_e[8] = b_e[4];

#else
  No esta implementado
#endif

  /* Locals */
  int n = 3;
  int lda = 3;
  int ldvl = 3;
  int ldvr = 3;
  int info;
  int lwork;
  double wkopt;
  double *work;
    
  /* Local arrays */
  int IPIV[3] = {0, 0, 0};
  double wi[3];
  double vl[9];

  /*
    Query and allocate the optimal workspace
  */
  lwork = -1;
  dsyev_("V","L",&n, eigvec_b_e, &lda, eigval_b_e,&wkopt, &lwork,&info);
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

  dsyev_("V","L",&n, eigvec_b_e, &lda, eigval_b_e,work, &lwork,&info);
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

double T[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

#if NumberDimensions == 2

  T[0] = P[0]*D_phi[0] + P[1]*D_phi[1];
  T[1] = P[0]*D_phi[2] + P[1]*D_phi[3];
  T[3] = P[2]*D_phi[0] + P[3]*D_phi[1];
  T[4] = P[2]*D_phi[2] + P[3]*D_phi[3];
  T[8] = P[4]*D_phi[4];

#else
  No esta implementado
#endif

  /* Locals */
  int n = 3;
  int lda = 3;
  int ldvl = 3;
  int ldvr = 3;
  int info;
  int lwork;
  double wkopt;
  double *work;
    
  /* Local arrays */
  int IPIV[3] = {0, 0, 0};
  double wi[3];
  double vl[9];

  /*
    Query and allocate the optimal workspace
  */
  lwork = -1;
  dsyev_("N","L",&n, T, &lda, eigval_T,&wkopt, &lwork,&info);
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

  dsyev_("N","L",&n, T, &lda, eigval_T,work, &lwork,&info);
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