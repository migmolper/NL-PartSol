#include <math.h>
#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/*******************************************************/

void update_increment_Deformation_Gradient__Particles__(
  double * DF_p,
  const double * DeltaU,
  const double * gradient_p,
  unsigned Nnodes_p) {

  unsigned Ndim = NumberDimensions;
  
  // Initialise with the identity tensor
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      DF_p[i*Ndim + j] = 1.0 * (i == j);
    }
  }

  // Ad the nodal contribution to the train tensor
  for (unsigned A = 0; A < Nnodes_p; A++) {
    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        DF_p[i*Ndim + j] += DeltaU[A*Ndim + i] * gradient_p[A*Ndim + j];
      }
    }

  }
}

/*******************************************************/

void update_rate_increment_Deformation_Gradient__Particles__(
  double * dt_DF_p, 
  const double * DeltaV, 
  const double * gradient_p,
  unsigned Nnodes_p) {

  /* Variable definition */
  unsigned Ndim = NumberDimensions;

  // Initialise with the identity tensor
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      dt_DF_p[i*Ndim + j] = 0.0;
    }
  }

  // Ad the nodal contribution to the train tensor
  for (unsigned A = 0; A < Nnodes_p; A++) {
    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        dt_DF_p[i*Ndim + j] += DeltaV[A*Ndim + i]* gradient_p[A*Ndim + j];
      }
    }
  }
}

/*******************************************************/

void update_Deformation_Gradient_n1__Particles__(
  double * F_n1, 
  const double * F_n,
  const double * f_n1) {

  int Ndim = NumberDimensions;
  double aux;

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      /*
              Set to zero the deformation gradient at t = n + 1
      */
      F_n1[i*Ndim + j] = 0;

      /*
        Compute row-column multiplication
      */
      aux = 0;
      for (int k = 0; k < Ndim; k++) {
        aux += f_n1[i*Ndim + k] * F_n[k*Ndim + j];
      }

      /*
        New value
      */
      F_n1[i*Ndim + j] = aux;
    }
  }
}

/*******************************************************/

int get_locking_free_Deformation_Gradient_n1__Particles__(
  unsigned p,
  double DJ_patch,
  Particle MPM_Mesh) {

  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned MatIndx_p = MPM_Mesh.MatIdx[p];

  double DJ_p;
  double DJ_averaged;
  double averaged_DF_vol;
  double alpha = MPM_Mesh.Mat[MatIndx_p].alpha_Fbar;

  double * F_total = MPM_Mesh.Phi.F_n1.nM[p];
  double * DF_p = MPM_Mesh.Phi.DF.nM[p];
  double * Fbar_n = MPM_Mesh.Phi.Fbar.nM[p];

#if NumberDimensions == 2
  double Fbar_n1[4] = {
    0.0,0.0,
    0.0,0.0};
  double DF_bar_p[4] = {
    0.0,0.0,
    0.0,0.0};
#else
  double Fbar_n1[9] = {
    0.0,0.0,0.0,
    0.0,0.0,0.0,
    0.0,0.0,0.0};
  double DF_bar_p[9] = {
    0.0,0.0,0.0,
    0.0,0.0,0.0,
    0.0,0.0,0.0};
#endif


  // Compute the averaged jacobian of the deformation gradient
  DJ_p = I3__TensorLib__(DF_p);
  if (DJ_p <= 0.0) {
    fprintf(stderr, ""RED"Negative jacobian in particle %i"RESET" \n",p);
    return EXIT_FAILURE;
  }

  DJ_averaged = DJ_patch / DJ_p;

  // Compute the averaged volume of the deformation gradient
  averaged_DF_vol = pow(DJ_averaged, (double)1.0 / Ndim);

  // Compute the incremental F-bar
  for (unsigned i = 0; i < Ndim*Ndim; i++) {
      DF_bar_p[i] = averaged_DF_vol * DF_p[i];
  }

  // Compute the new F-bar
  matrix_product__TensorLib__(Fbar_n1, DF_bar_p, Fbar_n);  

  // Update the deformation gradient to avoid locking (F-bar)
  for (unsigned i = 0; i < Ndim*Ndim; i++) {
    Fbar_n[i] = alpha * F_total[i] + (1 - alpha) * Fbar_n1[i];
  }

  return STATUS;
}

/*******************************************************/

void update_rate_Deformation_Gradient_n1__Particles__(
  double * dt_F_n1,
  const double * dt_f_n1,
  const double * F_n, 
  const double * f_n1,
  const double * dt_F_n) {
  
  unsigned Ndim = NumberDimensions;
  double aux;

  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      /*
        Set to zero the deformation gradient at t = n + 1
      */
      dt_F_n1[i*Ndim + j] = 0.0;

      /*
        Compute row-column multiplication
      */
      aux = 0.0;
      for (unsigned k = 0; k < Ndim; k++) {
        aux += dt_f_n1[i*Ndim + k] * F_n[k*Ndim + j] + f_n1[i*Ndim + k] * dt_F_n[k*Ndim + j];
      }

      /*
        New value
      */
      dt_F_n1[i*Ndim + j] = aux;
    }
  }
}

/*******************************************************/

int compute_Jacobian_Rate__Particles__(
  double * d_J_dt,
  double J, 
  const double * F,
  const double * d_F_dt) {
  
  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  double F_mT[4] = {
    0.0,0.0,
    0.0,0.0};
#else
  double F_mT[9] = {
    0.0,0.0,0.0,
    0.0,0.0,0.0,
    0.0,0.0,0.0};
#endif

  STATUS = compute_adjunt__TensorLib__(F_mT, F);
  if(STATUS == EXIT_FAILURE){
    fprintf(stderr, ""RED"Error in compute_adjunt__TensorLib__()"RESET" \n");
    return EXIT_FAILURE;
  }

  double F_mT__x__F = 0.0;
  for(unsigned i = 0; i<Ndim*Ndim; i++){
    F_mT__x__F += F_mT[i] * d_F_dt[i];
  }

  *d_J_dt = J * F_mT__x__F;

  return STATUS;
}

/*******************************************************/

int spatial_velocity_gradient__Particles__(
  double * L, 
  const double * dFdt,
  const double * F) {
    
  int STATUS = EXIT_SUCCESS;

  int Ndim = NumberDimensions;
     
#if NumberDimensions == 2

double F_m1[4] = {
  F[0], F[1],
  F[2], F[3]};
  
  int INFO;
  int N = 2;
  int LDA = 2;
  int LWORK = 2;
  int IPIV[2] = {0, 0};
  double WORK[2] = {0, 0};
#else

double F_m1[9] = {
  F[0], F[1], F[2],
  F[3], F[4], F[5],
  F[6], F[7], F[8]};
 
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};
#endif

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N, &N, F_m1, &LDA, IPIV, &INFO);
  // Check output of dgetrf
  if (INFO != 0) {
    if (INFO < 0) {
      printf(
          "" RED
          "Error in dgetrf_(): the %i-th argument had an illegal value " RESET
          "\n",
          abs(INFO));
    } else if (INFO > 0) {

      printf("" RED
             "Error in dgetrf_(): F_m1(%i,%i) %s \n %s \n %s \n %s " RESET
             "\n",
             INFO, INFO, "is exactly zero. The factorization",
             "has been completed, but the factor F_m1 is exactly",
             "singular, and division by zero will occur if it is used",
             "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  dgetri_(&N, F_m1, &LDA, IPIV, WORK, &LWORK, &INFO);
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(stderr, "" RED "%s: the %i-th argument %s" RESET "\n",
              "Error in dgetri_()", abs(INFO), "had an illegal value");
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetri_(): F_m1(%i,%i) %s \n %s \n %s \n %s " RESET
              "\n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor F_m1 is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }


  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
      L[i*Ndim + j] = 0.0;
      for (unsigned k = 0; k < Ndim; k++)
      {
        L[i*Ndim + j] += dFdt[i*Ndim + k]*F_m1[k*Ndim + j];
      }
    }
  }
  
    
  return STATUS;
}

/*******************************************************/

Tensor right_Cauchy_Green__Particles__(Tensor F) {
  /* Define output */
  Tensor C = alloc__TensorLib__(2);
  /* Define the number of dimensions */
  int Ndim = NumberDimensions;

  /* Compute C = F^T F */
  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      for (int k = 0; k < Ndim; k++) {
        C.N[i][j] += F.N[k][i] * F.N[k][j];
      }
    }
  }

  return C;
}

/*******************************************************/

void left_Cauchy_Green__Particles__(double * b, const double * F) {

#if NumberDimensions == 2
  b[0] = F[0]*F[0] + F[1]*F[1];
  b[1] = F[0]*F[2] + F[1]*F[3];
  b[2] = b[1];
  b[3] = F[2]*F[2] + F[3]*F[3];
#else  
  b[0] = F[0]*F[0] + F[1]*F[1] + F[2]*F[2];
  b[1] = F[0]*F[3] + F[1]*F[4] + F[2]*F[5];
  b[2] = F[0]*F[6] + F[1]*F[7] + F[2]*F[8];
  b[3] = b[1];
  b[4] = F[3]*F[3] + F[4]*F[4] + F[5]*F[5];
  b[5] = F[3]*F[6] + F[4]*F[7] + F[5]*F[8];
  b[6] = b[2];
  b[7] = b[5];
  b[8] = F[6]*F[6] + F[7]*F[7] + F[8]*F[8];
#endif

}

/*******************************************************/

Tensor strain_Green_Lagrange__Particles__(Tensor C) {
  /* Define output */
  Tensor E = alloc__TensorLib__(2);
  /* Define eye tensor */
  Tensor Identity = Identity__TensorLib__();
  /* Define the number of dimensions */
  int Ndim = NumberDimensions;

  /* Compute E = 1/2 * [ C - Identity]  */
  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      E.N[i][j] = 0.5 * (C.N[i][j] - Identity.N[i][j]);
    }
  }

  free__TensorLib__(Identity);

  return E;
}

/**************************************************************/
