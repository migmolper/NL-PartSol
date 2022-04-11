#include <math.h>
#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/**************************************************************/

double energy_Neo_Hookean_Wriggers(Tensor C, double J, Material MatProp) {
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /* Material parameters */
  double ElasticModulus = MatProp.E;
  double nu = MatProp.nu;
  double G = ElasticModulus / (2 * (1 + nu));
  double lambda = nu * ElasticModulus / ((1 - nu * 2) * (1 + nu));
  double I1_C = I1__TensorLib__(C);
  double f_J = 0.25 * lambda * (J * J - 1) - 0.5 * lambda * log(J) - G * log(J);

  double W = f_J + 0.5 * G * (I1_C - Ndim);

  return W;
}

/**************************************************************/

int compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(
  State_Parameters IO_State,
  Material MatProp) {

  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;

  // Get information from the state parameter
  double * P = IO_State.Stress;
  double * D_phi_n1 = IO_State.D_phi_n1;
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
  double D_phi_mT[4];
#else
  double D_phi_mT[9];
#endif

  compute_adjunt__TensorLib__(D_phi_mT,D_phi_n1);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,"" RED "Error in compute_adjunt__TensorLib__" RESET "\n");
    return EXIT_FAILURE;
  }

  // Compute stress
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      P[i*Ndim + j] =
          c0 * D_phi_mT[i*Ndim + j] 
          + G * (D_phi_n1[i*Ndim + j] - D_phi_mT[i*Ndim + j]);
    }
  }

#if NumberDimensions == 2
    P[4] = c0;
#endif

  return STATUS;
}

/**************************************************************/

int compute_stiffness_density_Neo_Hookean_Wriggers(
  double * Stiffness_Density,
  const double * dN_alpha_n,
  const double * dN_beta_n, 
  State_Parameters IO_State_p,
  Material MatProp) {

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;

  // State parameters
  double * D_phi_n = IO_State_p.D_phi_n;
  double * d_phi = IO_State_p.d_phi;  
  double J = IO_State_p.J;

  //  Material parameters
  double ElasticModulus = MatProp.E;
  double nu = MatProp.nu;
  double G = ElasticModulus / (2 * (1 + nu));
  double lambda = nu * ElasticModulus / ((1 - nu * 2) * (1 + nu));
  double sqr_J = J*J;
  double c0 = lambda * sqr_J;
  double c1 = G - 0.5 * lambda * (sqr_J - 1);

  // Define auxiliar variables
#if NumberDimensions == 2
  double d_phi_mT[4];
  double dN_alpha_n1[2];
  double dN_beta_n1[2];
  double b_n[4];
#else
  double d_phi_mT[9];
  double dN_alpha_n1[3];
  double dN_beta_n1[3];  
  double b_n[9];
#endif

  // Compute the adjunt of the incremental deformation gradient
  STATUS = compute_adjunt__TensorLib__(d_phi_mT, d_phi);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,"" RED "Error in compute_adjunt__TensorLib__" RESET "\n");
    return EXIT_FAILURE;
  }

  // Do the projection of the shape function gradient to the n + 1 configuration
  for (unsigned i = 0; i < Ndim; i++)
  {
    dN_alpha_n1[i] = 0.0;
    dN_beta_n1[i] = 0.0;
    for (unsigned j = 0; j < Ndim; j++)
    {
      dN_alpha_n1[i] += d_phi_mT[i*Ndim + j]*dN_alpha_n[j];
      dN_beta_n1[i] += d_phi_mT[i*Ndim + j]*dN_beta_n[j];
    }
  }
  
  // Compute the left Cauchy-Green (t = n)
  left_Cauchy_Green__Particles__(b_n, D_phi_n);

  // Compute
  double lenght_0 = 0.0;
  double lenght_0_aux = 0.0;
  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
      lenght_0_aux += b_n[i*Ndim + j]*dN_alpha_n[j];
    }
    lenght_0 += dN_beta_n[i]*lenght_0_aux;  
    lenght_0_aux = 0.0;
  }


  // Assemble Stiffness Density matrix
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      Stiffness_Density[i*Ndim + j] = 
      c0 * dN_alpha_n1[i]*dN_beta_n1[j] 
      + G * lenght_0 * (i == j) 
      + c1 * dN_alpha_n1[j]*dN_beta_n1[i];
    }
  }

  return STATUS;
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