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
  double * A_ep,
  const double * dN_alpha,
  const double * dN_beta,
  const State_Parameters IO_State)
{

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;

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

  double n1[2] = {0.0,0.0};  
  double n2[2] = {0.0,0.0};
 
  int m[4][4] =
  {
      {0.0,0.0,0.0},
      {0.0,0.0,0.0},
      {0.0,0.0,0.0},
      {0.0,0.0,0.0},
  };

  double mu[4][2] = 
  {
    {0.0,0.0},
    {0.0,0.0},
    {0.0,0.0},
    {0.0,0.0}
  };

  double mv[4][2] =
  {
    {0.0,0.0},
    {0.0,0.0},
    {0.0,0.0},
    {0.0,0.0}
  }; 

  n1[0] = eigvec_b_e[0];
  n1[1] = eigvec_b_e[1];
  n2[0] = eigvec_b_e[3];
  n2[1] = eigvec_b_e[4];

#else
  No esta implementado
#endif

  // Generate the 
  for(unsigned A = 0 ; A<Ndim; A++)
  {
    for(unsigned B = 0 ; B<Ndim; B++)
    {
      for (unsigned i = 0; i < Ndim; i++)
      {
        for (unsigned j = 0; j < Ndim; j++)
        {      
          m[A*Ndim + B][i*Ndim + j] = n1[i]*n1[j];
        }
      }
    }
  }

  for (unsigned A = 0; A < Ndim; A++)
  {
    for (unsigned B = 0; B < Ndim; B++)
    {
      for (unsigned i = 0; i < Ndim; i++)
      {
        for (unsigned j = 0; j < Ndim; j++)
        {
          mv[A*Ndim + B][i] += m[A*Ndim + B][i*Ndim+j]*v[j];
          mu[A*Ndim + B][i] += m[A*Ndim + B][i*Ndim+j]*u[j];
        }
      }  
    }  
  }
  

  for(unsigned A = 0 ; A<Ndim; A++)
  {
    for (unsigned B = 0; B < Ndim; B++)
    {
      for (unsigned i = 0; i < Ndim; i++)
      {
        for (unsigned j = 0; j < Ndim; j++)
        {      
          A_ep[A*Ndim + B] += IO_State.e_ep[A*Ndim + B]*(mv[A*Ndim + A][i] * mv[B*Ndim + B][j]) 
          + (A != B)*0.5*((eigval_T[B] - eigval_T[A])/(eigval_b_e[B] - eigval_b_e[A]))*
          (eigval_T[B]*(mv[A*Ndim + B][i] * mv[A*Ndim + B][j]) + eigval_T[A]*(mv[A*Ndim + B][i] * mv[B*Ndim + A][j]));
        }
      }
    }
  }


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