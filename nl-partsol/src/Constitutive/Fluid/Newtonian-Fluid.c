#include <math.h>
#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/**************************************************************/

int compute_1PK_Stress_Tensor_Newtonian_Fluid(
  State_Parameters IO_State,
  Material MatProp_p) {

  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;

  //  Take information from input state parameters
  double * P = IO_State.Stress;
  double * D_phi_n1 = IO_State.D_phi_n1;
  double * dFdt = IO_State.dFdt;
  double J = IO_State.J;

  //  Material parameters
  double p0 = MatProp_p.ReferencePressure;
  double mu = MatProp_p.Viscosity;
  double n = MatProp_p.n_Macdonald_model;
  double K = MatProp_p.Compressibility;
  double c0 = -J * (p0 + (K / n) * (pow(J, -n) - 1.0));
  double c1 = J * mu;

  // Define and compute auxiliar tensor
#if NumberDimensions == 2  
  double D_phi_mT[4];
  double L[4];
  double E[4];  
  double E__x__D_phi_mT[4];
#else
  double D_phi_mT[9];
  double L[9];
  double E[9];  
  double E__x__D_phi_mT[9];
#endif

  compute_adjunt__TensorLib__(D_phi_mT,D_phi_n1);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,"" RED "Error in compute_adjunt__TensorLib__" RESET "\n");
    return EXIT_FAILURE;
  }

  STATUS = spatial_velocity_gradient__Particles__(L, dFdt, D_phi_n1);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,"" RED "Error in spatial_velocity_gradient__Particles__" RESET "\n");
    return EXIT_FAILURE;
  }
  symmetrise__TensorLib__(E, L);  

#if NumberDimensions == 2 
  double trace_E = E[0] + E[2];
#else 
  double trace_E = E[0] + E[4] + E[8];
#endif

  matrix_product__TensorLib__(E__x__D_phi_mT,E,D_phi_mT); 


  // Compute stress
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      P[i*Ndim + j] = c0 * D_phi_mT[i*Ndim + j] +
                  2.0 * c1 * E__x__D_phi_mT[i*Ndim + j] -
                  (2.0 / 3.0) * c1 * trace_E * D_phi_mT[i*Ndim + j];
    }
  }

#if NumberDimensions == 2
  IO_State.Stress[4] = c0;
#endif


  return STATUS;
}

/**************************************************************/

int compute_stiffness_density_Newtonian_Fluid(
  double * Stiffness_Density,
  const double * dN_alpha_n,
  const double * dN_beta_n, 
  State_Parameters IO_State,
  Material MatProp) {

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;

  // State parameters
  double * D_phi_n1 = IO_State.D_phi_n1;
  double * D_phi_n = IO_State.D_phi;
  double * d_phi = IO_State.d_phi;
  double * dFdt = IO_State.dFdt;
  double J = IO_State.J;
  double alpha4 = IO_State.alpha_4;

  // Material parameters
  double p0 = MatProp.ReferencePressure;
  double mu = MatProp.Viscosity;
  double n = MatProp.n_Macdonald_model;
  double K = MatProp.Compressibility;
  double c0 = J * mu;
  double c1 = J * p0 + J * (K / n) * (pow(J, -n) - 1.0) - K * pow(J, 1.0 - n) + (2.0 / 3.0) * alpha4 * c0;
  double c2 = J * p0 + J * (K / n) * (pow(J, -n) - 1) + alpha4 * c0;

  // Define auxiliar variables
#if NumberDimensions == 2
  double d_phi_mT[4];
  double L[4];
  double E[4];
  double b_n[4];
  double dN_alpha_n1[2] = {0.0,0.0};
  double dN_beta_n1[2] = {0.0,0.0};
  double E_dN_alpha_n1[2] = {0.0,0.0};
  double E_dN_beta_n1[2] = {0.0,0.0};  
  double Lt_dN_alpha_n1[2] = {0.0,0.0};
  double Lt_dN_beta_n1[2] = {0.0,0.0};
#else  
  double d_phi_mT[9];
  double L[9];
  double E[9];
  double b_n[9];
  double dN_alpha_n1[3] = {0.0,0.0,0.0};
  double dN_beta_n1[3] = {0.0,0.0,0.0};
  double E_dN_alpha_n1[3] = {0.0,0.0,0.0};
  double E_dN_beta_n1[3] = {0.0,0.0,0.0};  
  double Lt_dN_alpha_n1[3] = {0.0,0.0,0.0};
  double Lt_dN_beta_n1[3] = {0.0,0.0,0.0};  
#endif

  // Compute the adjunt of the incremental deformation gradient
  STATUS = compute_adjunt__TensorLib__(d_phi_mT, d_phi);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,"" RED "Error in compute_adjunt__TensorLib__" RESET "\n");
    return EXIT_FAILURE;
  }

  // Compute rate tensors
  STATUS = spatial_velocity_gradient__Particles__(L, dFdt, D_phi_n1);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,"" RED "Error in spatial_velocity_gradient__Particles__" RESET "\n");
    return EXIT_FAILURE;
  }
  symmetrise__TensorLib__(E, L);


  // Do the projection of the shape function gradient to the n + 1 configuration
  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
      dN_alpha_n1[i] += d_phi_mT[i*Ndim + j]*dN_alpha_n[j];
      dN_beta_n1[i] += d_phi_mT[i*Ndim + j]*dN_beta_n[j];
    }
  }

  // 
  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
      E_dN_alpha_n1[i] += E[i*Ndim + j]*dN_alpha_n1[j];
      E_dN_beta_n1[i] += E[i*Ndim + j]*dN_beta_n1[j];

      Lt_dN_alpha_n1[i] += L[j*Ndim + i]*dN_alpha_n1[j];
      Lt_dN_beta_n1[i] += L[j*Ndim + i]*dN_beta_n1[j];

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
      Stiffness_Density[i*Ndim + j] = - c1 * dN_alpha_n1[i]*dN_beta_n1[j] +
                  c2 * dN_alpha_n1[j]*dN_beta_n1[i] +
                  2.0 * c0 * E_dN_alpha_n1[i]*dN_beta_n1[j] -
                  2.0 * c0 * E_dN_beta_n1[i]*dN_alpha_n1[j] +
                  alpha4 * c0 * (i == j) * lenght_0 -
                  c0 * lenght_0 * L[i*Ndim + j] -
                  c0 * dN_beta_n1[i]*Lt_dN_alpha_n1[j] +
                  (2.0 / 3.0) * c0 * dN_alpha_n1[i]*Lt_dN_beta_n1[j];
    }
  }

  return EXIT_SUCCESS;
}

/**************************************************************/
