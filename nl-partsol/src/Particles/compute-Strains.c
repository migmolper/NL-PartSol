#include <math.h>
#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/*************************************************************/

Tensor rate_inifinitesimal_Strain__Particles__(Matrix Velocity,
                                               Matrix Gradient) {
  unsigned Ndim = NumberDimensions;
  Tensor Rate_Strain = alloc__TensorLib__(2);
  Tensor Velocity_A;
  Tensor Gradient_A;
  Tensor VoG_A;

  unsigned NodesElem = Gradient.N_rows;

  /* Compute rate of strain */
  for (unsigned A = 0; A < NodesElem; A++) {
    /* Assign from matrix to tensor */
    Velocity_A = memory_to_tensor__TensorLib__(Velocity.nM[A], 1);
    Gradient_A = memory_to_tensor__TensorLib__(Gradient.nM[A], 1);

    /* Compute the dyadic product of the nodal velocity and the
       gradient of the shape functions */
    VoG_A = dyadic_Product__TensorLib__(Velocity_A, Gradient_A);

    /* Ad the nodal contribution to the train tensor */
    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        Rate_Strain.N[i][j] += 0.5 * (VoG_A.N[i][j] + VoG_A.N[j][i]);
      }
    }
    /* Free memory */
    free__TensorLib__(VoG_A);
  }

  return Rate_Strain;
}

/*******************************************************/

Tensor infinitesimal_Strain__Particles__(Tensor Strain, Tensor Rate_Strain,
                                         double TimeStep) {
  int Ndim = NumberDimensions;
  /* Check in the input its is ok */
  if ((Strain.Order == 2) && (Rate_Strain.Order == 2)) {
    /* Update strain tensor with the rate of strain tensor */
    for (int i = 0; i < Ndim; i++) {
      for (int j = 0; j < Ndim; j++) {
        Strain.N[i][j] += TimeStep * Rate_Strain.N[i][j];
      }
    }
  } else {
    fprintf(stderr, "%s : %s %s !!! \n",
            "Error in infinitesimal_Strain__Particles__()",
            "The input should be", "two tensors of 2nd order and a scalar");
    exit(EXIT_FAILURE);
  }

  return Strain;
}

/*******************************************************/

void update_increment_Deformation_Gradient__Particles__(Tensor DF_p,
                                                        Matrix DeltaU,
                                                        Matrix gradient_p) {

  /* Variable definition */
  unsigned Ndim = NumberDimensions;
  unsigned Nnodes_p = DeltaU.N_rows;
  Tensor f_n1;
  Tensor DeltaU_A;
  Tensor gradient_A;
  Tensor gradient_DeltaU_A;

  /*
    Compute increment of the deformation gradient
    f_n1 = I + Delta_u 0 gradient_N
  */

  /* Initialise with the identity tensor */
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      DF_p.N[i][j] = 1 * (i == j);
    }
  }

  for (unsigned A = 0; A < Nnodes_p; A++) {

    /* Assign from matrix to tensor */
    DeltaU_A = memory_to_tensor__TensorLib__(DeltaU.nM[A], 1);
    gradient_A = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);

    /* Compute the dyadic product of the nodal velocity and the
    gradient of the shape functions */
    gradient_DeltaU_A = dyadic_Product__TensorLib__(DeltaU_A, gradient_A);

    /* Ad the nodal contribution to the train tensor */
    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        DF_p.N[i][j] += gradient_DeltaU_A.N[i][j];
      }
    }

    /* Free memory */
    free__TensorLib__(gradient_DeltaU_A);
  }
}

/*******************************************************/

void update_rate_increment_Deformation_Gradient__Particles__(
    Tensor dt_DF_p, Matrix DeltaV, Matrix gradient_p) {

  /* Variable definition */
  unsigned Ndim = NumberDimensions;
  unsigned Nnodes_p = DeltaV.N_rows;
  Tensor f_n1;
  Tensor DeltaV_A;
  Tensor gradient_A;
  Tensor gradient_DeltaV_A;

  /*
    Compute increment of the deformation gradient
    dt_f_n1 = I + (Delta_V o gradient_N)
  */

  /* Initialise with the identity tensor */
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      dt_DF_p.N[i][j] = 0.0;
    }
  }

  for (unsigned A = 0; A < Nnodes_p; A++) {

    /* Assign from matrix to tensor */
    DeltaV_A = memory_to_tensor__TensorLib__(DeltaV.nM[A], 1);
    gradient_A = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);

    /* Compute the dyadic product of the nodal velocity and the
        gradient of the shape functions */
    gradient_DeltaV_A = dyadic_Product__TensorLib__(DeltaV_A, gradient_A);

    /* Ad the nodal contribution to the train tensor */
    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        dt_DF_p.N[i][j] += gradient_DeltaV_A.N[i][j];
      }
    }

    /* Free memory */
    free__TensorLib__(gradient_DeltaV_A);
  }
}

/*******************************************************/

void update_Deformation_Gradient_n1__Particles__(Tensor F_n1, Tensor F_n,
                                                 Tensor f_n1) {
  int Ndim = NumberDimensions;
  double aux;

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      /*
              Set to zero the deformation gradient at t = n + 1
      */
      F_n1.N[i][j] = 0;

      /*
        Compute row-column multiplication
      */
      aux = 0;
      for (int k = 0; k < Ndim; k++) {
        aux += f_n1.N[i][k] * F_n.N[k][j];
      }

      /*
        New value
      */
      F_n1.N[i][j] = aux;
    }
  }
}

/*******************************************************/

void get_locking_free_Deformation_Gradient_n1__Particles__(int p,
                                                           double DJ_patch,
                                                           Particle MPM_Mesh) {

  int Ndim = NumberDimensions;
  int MatIndx_p = MPM_Mesh.MatIdx[p];

  double DJ_p;
  double DJ_averaged;
  double averaged_DF_vol;

  double alpha = MPM_Mesh.Mat[MatIndx_p].alpha_Fbar;

  Tensor F_total = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p], 2);
  Tensor DF_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.DF.nM[p], 2);
  Tensor Fbar_n = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Fbar.nM[p], 2);
  Tensor DF_bar_p = alloc__TensorLib__(2);

  // Compute the averaged jacobian of the deformation gradient
  DJ_p = I3__TensorLib__(DF_p);
  DJ_averaged = DJ_patch / DJ_p;

  // Compute the averaged volume of the deformation gradient
  averaged_DF_vol = pow(DJ_averaged, (double)1 / Ndim);

  // Compute the incremental F-bar
  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      DF_bar_p.N[i][j] = averaged_DF_vol * DF_p.N[i][j];
    }
  }

  // Compute the new F-bar
  Tensor Fbar_n1 = matrix_product__TensorLib__(DF_bar_p, Fbar_n);

  // Update the deformation gradient to avoid locking (F-bar)
  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      Fbar_n.N[i][j] = alpha * F_total.N[i][j] + (1 - alpha) * Fbar_n1.N[i][j];
    }
  }

  free__TensorLib__(DF_bar_p);
  free__TensorLib__(Fbar_n1);
}

/*******************************************************/

void update_rate_Deformation_Gradient_n1__Particles__(Tensor dt_F_n1,
                                                      Tensor dt_f_n1,
                                                      Tensor F_n, Tensor f_n1,
                                                      Tensor dt_F_n) {
  int Ndim = NumberDimensions;
  double aux;

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      /*
        Set to zero the deformation gradient at t = n + 1
      */
      dt_F_n1.N[i][j] = 0;

      /*
        Compute row-column multiplication
      */
      aux = 0;
      for (int k = 0; k < Ndim; k++) {
        aux += dt_f_n1.N[i][k] * F_n.N[k][j] + f_n1.N[i][k] * dt_F_n.N[k][j];
      }

      /*
        New value
      */
      dt_F_n1.N[i][j] = aux;
    }
  }
}

/*******************************************************/

double compute_Jacobian_Rate__Particles__(double J_p, Tensor F_p,
                                          Tensor dt_F_p) {
  /* Variable definition */
  double dt_J_p = 0;
  Tensor inverse_F_p = Inverse__TensorLib__(F_p);
  Tensor transpose_inverse_F_p = transpose__TensorLib__(inverse_F_p);

  dt_J_p = J_p * inner_product__TensorLib__(transpose_inverse_F_p, dt_F_p);

  free__TensorLib__(inverse_F_p);
  free__TensorLib__(transpose_inverse_F_p);

  return dt_J_p;
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
