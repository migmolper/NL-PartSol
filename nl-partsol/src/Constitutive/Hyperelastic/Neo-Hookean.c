/**
 * @file Neo-Hookean.c
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Constitutive/Hyperelastic/Neo-Hookean.h"
#include "Globals.h"
#include <stdlib.h>

/**************************************************************/

double compute_strain_energy_Neo_Hookean__Constitutive__(const double *b, double J,
                                                         Material MatProp) {
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /* Material parameters */
  double ElasticModulus = MatProp.E;
  double nu = MatProp.nu;
  double G = ElasticModulus / (2 * (1 + nu));
  double lambda = nu * ElasticModulus / ((1 - nu * 2) * (1 + nu));
  double I1_b = I1__TensorLib__(b);
  double f_J = 0.25 * lambda * (J * J - 1) - 0.5 * lambda * log(J) - G * log(J);

  double W = f_J + 0.5 * G * (I1_b - Ndim);

  return W;
}

/**************************************************************/

int compute_Kirchhoff_Stress_Neo_Hookean__Constitutive__(State_Parameters IO_State,
                                         Material MatProp) {

  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;

  // Get information from the state parameter
  double *T = IO_State.Stress;
  double *D_phi_n1 = IO_State.D_phi_n1;
  double J = IO_State.J;

  // Material parameters
  double ElasticModulus = MatProp.E;
  double nu = MatProp.nu;
  double G = ElasticModulus / (2 * (1 + nu));
  double lambda = nu * ElasticModulus / ((1 - nu * 2) * (1 + nu));
  double J2 = J * J;
  double c0 = lambda * 0.5 * (J2 - 1.0);

  // Define and compute auxiliar tensor
#if NumberDimensions == 2
  double b_n1[4];
  double Identity[4] = {1.0, 0.0, 0.0, 1.0};
#else
  double b_n1[9];
  double Identity[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
#endif

  // Compute the left Cauchy-Green (t = n)
  left_Cauchy_Green__Particles__(b_n1, D_phi_n1);

  // Compute stress
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      T[i * Ndim + j] = c0 * Identity[i * Ndim + j] +
                        G * (b_n1[i * Ndim + j] - Identity[i * Ndim + j]);
    }
  }

#if NumberDimensions == 2
  T[4] = c0;
#endif

//! Compute deformation energy
 *(IO_State.W) = compute_strain_energy_Neo_Hookean__Constitutive__(b_n1, J, MatProp);

  return STATUS;
}

/**************************************************************/

int compute_stiffness_density_Neo_Hookean(
    double *Stiffness_Density, const double *dN_alpha_n1,
    const double *dN_beta_n1, const double *dN_alpha_n, const double *dN_beta_n,
    State_Parameters IO_State, Material MatProp) {

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;

  // State parameters
  double *D_phi_n = IO_State.D_phi_n;
  double J = IO_State.J;

  //  Material parameters
  double ElasticModulus = MatProp.E;
  double nu = MatProp.nu;
  double G = ElasticModulus / (2 * (1 + nu));
  double lambda = nu * ElasticModulus / ((1 - nu * 2) * (1 + nu));
  double sqr_J = J * J;
  double c0 = lambda * sqr_J;
  double c1 = G - 0.5 * lambda * (sqr_J - 1);

  // Define auxiliar variables
#if NumberDimensions == 2
  double b_n[4];
#else
  double b_n[9];
#endif

  // Compute the left Cauchy-Green (t = n)
  left_Cauchy_Green__Particles__(b_n, D_phi_n);

  // Compute
  double lenght_0 = 0.0;
  double lenght_0_aux = 0.0;
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      lenght_0_aux += b_n[i * Ndim + j] * dN_alpha_n[j];
    }
    lenght_0 += dN_beta_n[i] * lenght_0_aux;
    lenght_0_aux = 0.0;
  }

  // Assemble Stiffness Density matrix
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      Stiffness_Density[i * Ndim + j] = c0 * dN_alpha_n1[i] * dN_beta_n1[j] +
                                        G * lenght_0 * (i == j) +
                                        c1 * dN_alpha_n1[j] * dN_beta_n1[i];
    }
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

Tensor compute_2PK_Stress_Tensor_Neo_Hookean_Wriggers(Tensor grad_e, Tensor C,
                                                      double J,
                                                      Material MatProp) {
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /* Material parameters */
  double ElasticModulus = MatProp.E;
  double nu = MatProp.nu;
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