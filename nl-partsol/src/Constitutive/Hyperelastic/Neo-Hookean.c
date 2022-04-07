#include <math.h>
#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/**************************************************************/

double energy_Neo_Hookean_Wriggers(Tensor C, double J, Material MatProp_p) {
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /* Material parameters */
  double ElasticModulus = MatProp_p.E;
  double nu = MatProp_p.nu;
  double G = ElasticModulus / (2 * (1 + nu));
  double lambda = nu * ElasticModulus / ((1 - nu * 2) * (1 + nu));
  double I1_C = I1__TensorLib__(C);
  double f_J = 0.25 * lambda * (J * J - 1) - 0.5 * lambda * log(J) - G * log(J);

  double W = f_J + 0.5 * G * (I1_C - Ndim);

  return W;
}

/**************************************************************/

State_Parameters
compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(State_Parameters Intput_SP,
                                               Material MatProp_p) {
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /*
    Output state parameter
  */
  State_Parameters Output_SP;

  /* Get information from the state parameter */
  Tensor F = memory_to_tensor__TensorLib__(Intput_SP.D_phi, 2);
  Tensor P = memory_to_tensor__TensorLib__(Intput_SP.Stress, 2);
  double J = Intput_SP.J;

  /* Material parameters */
  double ElasticModulus = MatProp_p.E;
  double nu = MatProp_p.nu;
  double G = ElasticModulus / (2 * (1 + nu));
  double lambda = nu * ElasticModulus / ((1 - nu * 2) * (1 + nu));
  double J2 = J * J;

  /*
    Auxiliar tensors
  */
  Tensor Fm1 = Inverse__TensorLib__(F);

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      P.N[i][j] =
          lambda * 0.5 * (J2 - 1) * Fm1.N[j][i] + G * (F.N[i][j] - Fm1.N[j][i]);
    }
  }

  /*
    Plane strain conditions
  */
  if (Ndim == 2) {
    Intput_SP.Stress[4] = lambda * 0.5 * (J2 - 1);
  }

  /*
    Free tensors
  */
  free__TensorLib__(Fm1);

  return Output_SP;
}

/**************************************************************/

int compute_stiffness_density_Neo_Hookean_Wriggers(
  double * Stiffness_Density,
  const double * dN_alpha_n,
  const double * dN_beta_n, 
  State_Parameters IO_State_p,
  Material MatProp) {

  /*
    Number of dimensions
  */
  int Ndim = NumberDimensions;

  /*
    Material parameters
  */
  double ElasticModulus = MatProp.E;
  double nu = MatProp.nu;
  double G = ElasticModulus / (2 * (1 + nu));
  double lambda = nu * ElasticModulus / ((1 - nu * 2) * (1 + nu));
  double J = IO_State_p.J;
  double sqr_J = J*J;
  double alpha = lambda * sqr_J;
  double beta = G - 0.5 * lambda * (sqr_J - 1);

  double * D_phi_n = IO_State_p.D_phi;
  double * d_phi = IO_State_p.d_phi;

#if NumberDimensions == 2
  double d_phi_mT[4] = {
    d_phi[0], d_phi[2], 
    d_phi[1], d_phi[3]};

  // Parameters for dgetrf_ and dgetri_
  int INFO;
  int N = 2;
  int LDA = 2;
  int LWORK = 2;
  int IPIV[2] = {0, 0};
  double WORK[2] = {0, 0};
#else
  double d_phi_mT[9] = {
    d_phi[0], d_phi[3], d_phi[6],
    d_phi[1], d_phi[4], d_phi[7], 
    d_phi[2], d_phi[5], d_phi[8]};
  
  // Parameters for dgetrf_ and dgetri_
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};
#endif

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N, &N, d_phi_mT, &LDA, IPIV, &INFO);
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
             "Error in dgetrf_(): d_phi_mT(%i,%i) %s \n %s \n %s \n %s " RESET
             "\n",
             INFO, INFO, "is exactly zero. The factorization",
             "has been completed, but the factor d_phi_mT is exactly",
             "singular, and division by zero will occur if it is used",
             "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  dgetri_(&N, d_phi_mT, &LDA, IPIV, WORK, &LWORK, &INFO);
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(stderr, "" RED "%s: the %i-th argument %s" RESET "\n",
              "Error in dgetri_()", abs(INFO), "had an illegal value");
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetri_(): d_phi_mT(%i,%i) %s \n %s \n %s \n %s " RESET
              "\n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor d_phi_mT is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  // Do the projection of the shape function gradient to the n + 1 configuration
#if NumberDimensions == 2
  double dN_alpha_n1[2] = {0.0,0.0};
  double dN_beta_n1[2] = {0.0,0.0};
#else  
  double dN_alpha_n1[3] = {0.0,0.0,0.0};
  double dN_beta_n1[3] = {0.0,0.0,0.0};
#endif

  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
      dN_alpha_n1[i] += d_phi_mT[i*Ndim + j]*dN_alpha_n[j];
      dN_beta_n1[i] += d_phi_mT[i*Ndim + j]*dN_beta_n[j];
    }
  }
  
  // Compute the left Cauchy-Green (t = n)
#if NumberDimensions == 2
  double b_n[4] = {
    0.0,0.0,
    0.0,0.0};
#else  
  double b_n[9] = {
    0.0,0.0,0.0,
    0.0,0.0,0.0,
    0.0,0.0,0.0};
#endif
  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
      for (unsigned k = 0; k < Ndim; k++)
      {
        b_n[i*Ndim + j] += D_phi_n[i*Ndim + k]*D_phi_n[j*Ndim + k];
      }      
    }
  }

  double lenght_0 = 0.0;
  double lenght_0_aux = 0.0;
  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
      lenght_0_aux += b_n[i*Ndim + j]*dN_alpha_n[j];
    }
    lenght_0 += dN_beta_n[i]*lenght_0;  
    lenght_0_aux = 0.0;
  }


  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      Stiffness_Density[i*Ndim + j] += 
      alpha * dN_alpha_n1[i]*dN_beta_n1[j] 
      + G * lenght_0 * (i == j) 
      + beta * dN_alpha_n1[j]*dN_beta_n1[i];
    }
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

Tensor compute_2PK_Stress_Tensor_Neo_Hookean_Wriggers(Tensor grad_e, Tensor C,
                                                      double J,
                                                      Material MatProp_p) {
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /* Material parameters */
  double ElasticModulus = MatProp_p.E;
  double nu = MatProp_p.nu;
  double G = ElasticModulus / (2 * (1 + nu));
  double lambda = nu * ElasticModulus / ((1 - nu * 2) * (1 + nu));
  double J2 = J * J;

  /*
    Auxiliar tensors
  */
  Tensor Identity = Identity__TensorLib__();
  Tensor C_m1 = Inverse__TensorLib__(C);

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      grad_e.N[i][j] = lambda * 0.5 * (J2 - 1) * C_m1.N[i][j] +
                       G * (Identity.N[i][j] - C_m1.N[i][j]);
    }
  }

  /*
    Free tensors
  */
  free__TensorLib__(Identity);
  free__TensorLib__(C_m1);

  return grad_e;
}

/**************************************************************/

Tensor compute_material_stiffness_density_Neo_Hookean_Wriggers(
    Tensor v, Tensor w, Tensor C, double J, Material MatProp) {

  /*
    Number of dimensions
  */
  int Ndim = NumberDimensions;

  /*
    Material parameters
  */
  double ElasticModulus = MatProp.E;
  double nu = MatProp.nu;
  double G = ElasticModulus / (2 * (1 + nu));
  double lambda = nu * ElasticModulus / ((1 - nu * 2) * (1 + nu));
  double J2 = J * J;
  double alpha = lambda * J2;
  double beta = 0.5 * lambda * (J2 - 1) - G;

  /*
    Stifness density tensor
  */
  Tensor C_mat = alloc__TensorLib__(2);

  /*
    Auxiliar variables
  */
  Tensor Cm1;
  Tensor Cm1_dot_v;
  Tensor Cm1_dot_w;
  Tensor Cm1_dot_v_o_Cm1_dot_w;
  double Cm1_dot_w_dot_v;

  Cm1 = Inverse__TensorLib__(C);
  Cm1_dot_v = vector_linear_mapping__TensorLib__(Cm1, v);
  Cm1_dot_w = vector_linear_mapping__TensorLib__(Cm1, w);
  Cm1_dot_w_dot_v = inner_product__TensorLib__(Cm1_dot_w, v);
  Cm1_dot_v_o_Cm1_dot_w = dyadic_Product__TensorLib__(Cm1_dot_v, Cm1_dot_w);

  for (int A = 0; A < Ndim; A++) {
    for (int B = 0; B < Ndim; B++) {
      C_mat.N[A][B] += alpha * Cm1_dot_v_o_Cm1_dot_w.N[A][B] -
                       beta * Cm1_dot_w_dot_v * Cm1.N[A][B] -
                       beta * Cm1_dot_v_o_Cm1_dot_w.N[B][A];
    }
  }

  /*
    Free memory
   */
  free__TensorLib__(Cm1);
  free__TensorLib__(Cm1_dot_w);
  free__TensorLib__(Cm1_dot_v);
  free__TensorLib__(Cm1_dot_v_o_Cm1_dot_w);

  return C_mat;
}

/**************************************************************/

Matrix compute_D_matrix_Neo_Hookean_Wriggers(Tensor C, double J,
                                             Material MatProp)
/*

*/
{

  /*
    Material parameters
  */
  double ElasticModulus = MatProp.E;
  double nu = MatProp.nu;
  double G = ElasticModulus / (2 * (1 + nu));
  double lambda = nu * ElasticModulus / ((1 - nu * 2) * (1 + nu));
  double J2 = J * J;

  Tensor Cm1 = Inverse__TensorLib__(C);

  Matrix D = allocZ__MatrixLib__(3, 3);

  D.nM[0][0] = lambda * J2 * Cm1.N[0][0] * Cm1.N[0][0] -
               (lambda * (J2 - 1) - 2 * G) * Cm1.N[0][0] * Cm1.N[0][0];
  D.nM[0][1] = lambda * J2 * Cm1.N[0][0] * Cm1.N[1][1] -
               (lambda * (J2 - 1) - 2 * G) * Cm1.N[0][1] * Cm1.N[1][0];
  D.nM[1][0] = D.nM[0][1];
  D.nM[0][2] = lambda * J2 * Cm1.N[0][0] * Cm1.N[0][1] -
               (lambda * (J2 - 1) - 2 * G) * Cm1.N[0][0] * Cm1.N[1][0];
  D.nM[2][0] = D.nM[0][2];
  D.nM[1][1] = lambda * J2 * Cm1.N[1][1] * Cm1.N[1][1] -
               (lambda * (J2 - 1) - 2 * G) * Cm1.N[1][1] * Cm1.N[1][1];
  D.nM[1][2] = lambda * J2 * Cm1.N[1][1] * Cm1.N[0][1] -
               (lambda * (J2 - 1) - 2 * G) * Cm1.N[1][0] * Cm1.N[1][1];
  D.nM[2][1] = D.nM[1][2];
  D.nM[2][2] = lambda * J2 * Cm1.N[0][1] * Cm1.N[0][1] -
               0.5 * (lambda * (J2 - 1) - 2 * G) *
                   (Cm1.N[0][0] * Cm1.N[1][1] + Cm1.N[0][1] * Cm1.N[0][1]);

  free__TensorLib__(Cm1);

  return D;
}

/**************************************************************/
