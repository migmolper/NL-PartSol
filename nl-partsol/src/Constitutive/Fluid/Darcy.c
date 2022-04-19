#include <math.h>
#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/**************************************************************/




/**************************************************************/

int compute_stiffness_density_Darcy(
  double * Stiffness_Density,
  const double * dN_alpha_n,
  const double * dN_beta_n,
  double N_alpha_n,
  double N_beta_n,  
  State_Parameters IO_State_p,
  Mixture MixProp,
  Material MatProp)
{
 
    int STATUS = EXIT_SUCCESS;
    unsigned Ndim = NumberDimensions;
 
    // Material parameters
    double kappa_f = MatProp.Compressibility;

    // Mixture parameters
    double phi_f_0 = MixProp.phi_f_0;
    double phi_s_0 = MixProp.phi_s_0;    
    double * K = MixProp.Permeability;

    // State parameters
    double * D_phi_n1 = IO_State.D_phi_n1;
    double * D_phi_n = IO_State.D_phi_n;
    double * d_phi = IO_State.d_phi;
    double * dFdt = IO_State.dFdt;
    double J = IO_State.J;

    double * dpw_n = IO_State.D_pw_n;
    double pw_n1 = IO_State.pw;  
    double dpw_dt = IO_State.dpw_dt;    

    double phi_f = IO_State.Fluid_VolumeFraction;
    double intrinsic_rho_f = ;

    double * dyn_n1 = IO_State.dyn;
    double alpha4 = IO_State.alpha_4;
    double alpha1 = IO_State.alpha_1;

  // Compute tensorial variables
#if NumberDimensions == 2
    double d_phi_mT[4];
    double L[4];
    double dN_alpha_n1[2] = {0.0,0.0};
    double dN_beta_n1[2] = {0.0,0.0};
    double Lt_dN_beta_n1[2] = {0.0,0.0};
    double dpw_n1[2] = {0.0,0.0};
    double K_dN_alpha_n1[2] = {0.0,0.0};
    double K_dN_beta_n1[2] = {0.0,0.0};
    double K_dpw_n1[2] = {0.0,0.0};      
    double K_dyn_n1[2] = {0.0,0.0};
    double dN_alpha_o_dN_beta_n1[4] = {0.0,0.0,0.0,0.0};
    double skew_dN_alpha_o_dN_beta_n1[4] = {0.0,0.0,0.0,0.0};
    double dN_alpha_o_dN_beta_K_dpw_n1[2] = {0.0,0.0};
    double dN_beta_o_dN_alpha_K_dyn_n1[2] = {0.0,0.0};
    double skew_dN_alpha_o_dN_beta_K_dyn_n1[2] = {0.0,0.0};
#else  
    double d_phi_mT[9];
    double L[9];
    double dN_alpha_n1[3] = {0.0,0.0,0.0};  
    double dN_beta_n1[3] = {0.0,0.0,0.0};
    double Lt_dN_beta_n1[3] = {0.0,0.0,0.0};
    double dpw_n1[3] = {0.0,0.0,0.0};    
    double K_dN_alpha_n1[3] = {0.0,0.0,0.0};
    double K_dN_beta_n1[3] = {0.0,0.0,0.0};     
    double K_dpw_n1[3] = {0.0,0.0,0.0};   
    double K_dyn_n1[3] = {0.0,0.0,0.0};
    double dN_alpha_o_dN_beta_n1[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double skew_dN_alpha_o_dN_beta_n1[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double dN_alpha_o_dN_beta_K_dpw_n1[3] = {0.0,0.0,0.0};
    double dN_beta_o_dN_alpha_K_dyn_n1[3] = {0.0,0.0,0.0};
    double skew_dN_alpha_o_dN_beta_K_dyn_n1[3] = {0.0,0.0,0.0};
#endif

    STATUS = compute_adjunt__TensorLib__(d_phi_mT, d_phi);
    if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,"" RED "Error in compute_adjunt__TensorLib__" RESET "\n");
        return EXIT_FAILURE;
    }

    STATUS = spatial_velocity_gradient__Particles__(L, dFdt, D_phi_n1);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,"" RED "Error in spatial_velocity_gradient__Particles__" RESET "\n");
      return EXIT_FAILURE;
    }

#if NumberDimensions == 2
    double trL = L[0] + L[3];
#else  
    double trL = L[0] + L[4] + L[8];
#endif

  // Do the projection of the shape function gradient to the n + 1 configuration
  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
      dN_alpha_n1[i] += d_phi_mT[i*Ndim + j]*dN_alpha_n[j];
      dN_beta_n1[i] += d_phi_mT[i*Ndim + j]*dN_beta_n[j];
      dpw_n1[i] += d_phi_mT[i*Ndim + j]*dpw_n[j];
    }
  }

  // 
  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
        Lt_dN_beta_n1[i] += L[j*Ndim + i]*dN_beta_n1[j];
        K_dN_alpha_n1[i] += K[i*Ndim + j]*dN_alpha_n1[j];
        K_dN_beta_n1[i] += K[i*Ndim + j]*dN_beta_n1[j];
        K_dpw_n1[i] += K[i*Ndim + j]*dpw_n1[j];
        K_dyn_n1[i] += K[i*Ndim + j]*dyn_n1[j];
        dN_alpha_o_dN_beta_n1[i*Ndim + j] = dN_alpha_n1[i]*dN_beta_n1[j];
        skew_dN_alpha_o_dN_beta_n1[i*Ndim + j] = dN_alpha_n1[i]*dN_beta_n1[j] - dN_alpha_n1[j]*dN_beta_n1[i];        
    }
  }

  double K__x__dN_alpha_o_dN_beta_n1 = 0.0;
  double dN_alpha__x__K_dyn_n1 = 0.0;

  for (unsigned i = 0; i < Ndim; i++)
  {
    dN_alpha__x__K_dyn_n1 += dN_alpha_n1[i]*K_dyn_n1[i];
    for (unsigned j = 0; j < Ndim; j++)
    {
       dN_alpha_o_dN_beta_K_dpw_n1[i] += dN_alpha_o_dN_beta_n1[i*Ndim + j]*K_dpw_n1[j];
       dN_beta_o_dN_alpha_K_dyn_n1[i] += dN_alpha_o_dN_beta_n1[j*Ndim + i]*K_dyn_n1[j];
       skew_dN_alpha_o_dN_beta_K_dyn_n1[i] += skew_dN_alpha_o_dN_beta_n1[i*Ndim + j]*K_dyn_n1[j];
       K__x__dN_alpha_o_dN_beta_n1 += K[i*Ndim + j]*dN_alpha_o_dN_beta_n1[i*Ndim + j];
    }
  }

    

    // Compute scalar variables
    double c0 = dpw_dt * (1.0 - phi_f_0) * intrinsic_rho_f / (kappa_f * J);
    double c1 = dpw_dt * phi_f * intrinsic_rho_f * pw_n1 / (kappa_f * kappa_f * J);
    double c2 = intrinsic_rho_f * trL * pw_n1 / (kappa_f);
    double c3 = intrinsic_rho_f * J * trL;
    double c4 = intrinsic_rho_f * alpha4 * J;
    double c5 = intrinsic_rho_f * J;
    double c6 = dpw_dt * phi_f * intrinsic_rho_f / (J * kappa_f * kappa_f);
    double c7 = phi_f * intrinsic_rho_f * alpha4 / kappa_f;
    double c8 = intrinsic_rho_f * trL / kappa_f;

    double d0 = J * intrinsic_rho_f * alpha1 / g;
    double d1 = 1.0 / g;
    double d2 = J_p * intrinsic_rho_f_p * d1;
    double d3 = d1 * intrinsic_rho_f_p * pw_n1 / kappa_f;
    double d4 = d1 * intrinsic_rho_f_p / kappa_f;


  for (unsigned i = 0; i < Ndim; i++) {
    Stiffness_Density[i] =
        + c0 * N_alpha_n * dN_beta_n1[i] 
        - c1 * N_alpha_n * dN_beta_n1[i] 
        - c2 * N_alpha_n * dN_beta_n1[i] 
        + c3 * N_alpha_n * dN_beta_n1[i] 
        + c4 * N_alpha_n * dN_beta_n1[i] 
        - c5 * N_alpha_n * Lt_dN_beta_n1[i]
        - d0 * K_dN_alpha_n1[i] * N_beta_n
        + d1*dN_alpha_o_dN_beta_K_dpw_n1[i]
        + d1*K__x__dN_alpha_o_dN_beta_n1*dpw_n1[i]
        + d2*skew_dN_alpha_o_dN_beta_K_dyn_n1[i]
        + d3*dN_beta_o_dN_alpha_K_dyn_n1[i];
  }

  Stiffness_Density[Ndim] =
      + c6 * N_alpha_n * N_beta_n 
      + c7 * N_alpha_n * N_beta_n 
      + c8 * N_alpha_n * N_beta_n
      - d1*K__x__dN_alpha_o_dN_beta_n1
      - d4*dN_alpha__x__K_dyn_n1*N_beta_n;


    return STATUS;
}


/**************************************************************/