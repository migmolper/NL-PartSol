
// clang-format off
#include <math.h>
#include <stdbool.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"
#include "Particles/compute-Strains.h"
// clang-format on

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/*******************************************************/

void update_increment_Deformation_Gradient__Particles__(
  double * DF_p,
  const double * DeltaU,
  void* ctx) {

  unsigned Ndim = NumberDimensions;
  
  const double * gradient_p = ((compute_strains_ctx*)ctx)->gradient_p;
  unsigned Nnodes_p = ((compute_strains_ctx*)ctx)->Nnodes_p;

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

#if USE_AXIAL_SYMMETRY
  const double * shapefun_p = ((compute_strains_ctx*)ctx)->shapefun_p;
  double R_p = ((compute_strains_ctx*)ctx)->R_p; 
  for (unsigned A = 0; A < Nnodes_p; A++) {
    DF_p[4] += DeltaU[A*Ndim] * shapefun_p[A]/R_p;
  }
#endif
}

/*******************************************************/

void update_rate_increment_Deformation_Gradient__Particles__(
  double * dt_DF_p, 
  const double * DeltaV, 
  void* ctx) {

  /* Variable definition */
  unsigned Ndim = NumberDimensions;
  const double * gradient_p = ((compute_strains_ctx*)ctx)->gradient_p;
  unsigned Nnodes_p = ((compute_strains_ctx*)ctx)->Nnodes_p;

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

#if USE_AXIAL_SYMMETRY
  const double * shapefun_p = ((compute_strains_ctx*)ctx)->shapefun_p;
  double R_p = ((compute_strains_ctx*)ctx)->R_p; 
  for (unsigned A = 0; A < Nnodes_p; A++) {
    dt_DF_p[4] += DeltaV[A*Ndim] * shapefun_p[A]/R_p;
  }
#endif

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

#if USE_AXIAL_SYMMETRY
  F_n1[4] = f_n1[4] * F_n[4];
#endif

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

#if USE_AXIAL_SYMMETRY
  Fbar_n[4] = alpha * F_total[4] + (1 - alpha) * Fbar_n1[4];
#endif

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

#if USE_AXIAL_SYMMETRY
  dt_F_n1[4] = dt_f_n1[4] * F_n[4] + f_n1[4] * dt_F_n[4];
#endif

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

#if USE_AXIAL_SYMMETRY
  F_mT__x__F += F_mT[4] * d_F_dt[4];
#endif


  *d_J_dt = J * F_mT__x__F;

  return STATUS;
}

/*******************************************************/

int spatial_velocity_gradient__Particles__(
  double * L, 
  const double * dFdt,
  const double * F) {
    
  int STATUS = EXIT_SUCCESS;

  unsigned int Ndim = NumberDimensions;
     
#if NumberDimensions == 2
  double F_m1[5] ;  
  int INFO;
  int N = 2;
  int LDA = 2;
  int LWORK = 2;
  int IPIV[2] = {0, 0};
  double WORK[2] = {0, 0};
#else
  double F_m1[9];
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};
#endif


for (unsigned int i = 0; i < Ndim*Ndim; i++) {
  F_m1[i] = F[i];
}

#if USE_AXIAL_SYMMETRY
  F_m1[4] = F[4];
#endif

  compute_inverse__TensorLib__(F_m1,F);


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
  
#if USE_AXIAL_SYMMETRY
  L[4] = dFdt[4]*F_m1[4];
#endif
    
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
#if USE_AXIAL_SYMMETRY
  b[4] = F[4]*F[4];  
#endif
  
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

void eulerian_almansi__Particles__(double * e, const double * F) {

  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  double b[4] = {
    0.0,0.0,
    0.0,0.0};
  double b_m1[4] = {
    0.0,0.0,
    0.0,0.0};
  double Identity[4] = {
    1.0,0.0,
    0.0,1.0};
#else  
  double b[9] = {
    0.0,0.0,0.0,
    0.0,0.0,0.0,
    0.0,0.0,0.0};
  double b_m1[9] = {
    0.0,0.0,0.0,
    0.0,0.0,0.0,
    0.0,0.0,0.0};
  double Identity[9] = {
    1.0,0.0,0.0,
    0.0,1.0,0.0,
    0.0,0.0,1.0};
#endif

  left_Cauchy_Green__Particles__(b, F);

  compute_inverse__TensorLib__(b_m1,b);

  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
      e[i*Ndim + j] = 0.5*(Identity[i*Ndim + j] - b_m1[i*Ndim + j]);
    }
  }
  
#if USE_AXIAL_SYMMETRY
  e[4] = 0.5*(Identity[4] - b_m1[4]);
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
