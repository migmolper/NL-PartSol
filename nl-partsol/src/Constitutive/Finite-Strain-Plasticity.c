#include "nl-partsol.h"

/**************************************************************/

static Tensor compute_kirchhoff_isotropic_linear_elasticity(Tensor, Tensor, Material);

/**************************************************************/

State_Parameters finite_strains_plasticity(
  Tensor P_p,
  State_Parameters Inputs_SP, 
  Material MatProp,
  State_Parameters(* infinitesimal_plasticity)(State_Parameters,Material))
/*
  Finite strains plasticity following the apporach of Ortiz and Camacho
*/
{
  int Ndim = NumberDimensions;

  /* Define auxiliar variables */
  State_Parameters Output_SP;
  Tensor F_m1_plastic = Inputs_SP.F_m1_plastic_p;
  Tensor F_total = memory_to_tensor__TensorLib__(Inputs_SP.F_n1_p,2);
  Tensor F_m1_total;
  Tensor F_trial_elastic;
  Tensor C_trial_elastic;
  Tensor E_trial_elastic;
  Tensor D_F_plastic;
  Tensor Fm1_plastic;
  Tensor T_p = alloc__TensorLib__(2);

  /* Compute the elastic right Cauchy-Green tensor using the intermediate configuration. */ 
  F_trial_elastic = matrix_product__TensorLib__(F_total,F_m1_plastic);

  C_trial_elastic = right_Cauchy_Green__Particles__(F_trial_elastic);

  /* Calculation of the small strain tensor */
  E_trial_elastic = logarithmic_strains__Particles__(C_trial_elastic);

  /* Calculation of the trial stress tensor using the trial small strain tensor */
  T_p = compute_kirchhoff_isotropic_linear_elasticity(T_p, E_trial_elastic, MatProp);

  /* Start plastic corrector algorithm in infinitesimal strains */
  Output_SP = infinitesimal_plasticity(Inputs_SP, MatProp);

  /* Use the Cuitiño & Ortiz exponential maping to compute the increment of plastic finite strains */
  update_plastic_deformation_gradient__Particles__(Output_SP.Increment_E_plastic, F_m1_plastic);

  /* Get the First Piola-Kirchhoff stress tensor (P_p) */
  F_m1_total = Inverse__TensorLib__(F_total);

  for(int i = 0 ; i < Ndim  ; i++)
  {
    for(int j = 0 ; j < Ndim  ; j++)    
    {

     P_p.N[i][j] = 0.0;

     for(int k = 0 ; k < Ndim  ; k++)
     {
        P_p.N[i][j] += T_p.N[i][k]*F_m1_total.N[k][j];
      }
    }
  }


  /* Free memory */
  free__TensorLib__(F_trial_elastic);
  free__TensorLib__(C_trial_elastic);
  free__TensorLib__(E_trial_elastic);
  free__TensorLib__(T_p);
  free__TensorLib__(F_m1_total);

  return Output_SP;
}

/**************************************************************/

static Tensor compute_kirchhoff_isotropic_linear_elasticity(
  Tensor Stress,
  Tensor Strain,
  Material Mat)
{

  int Ndim = NumberDimensions;
  double nu = Mat.nu; 
  double E = Mat.E;
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double traceStrain = I1__TensorLib__(Strain);

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      Stress.N[i][j] = 2*G*Strain.N[i][j] + (K - 2*G/3.0)*(i==j)*traceStrain;
    }
  }

  return Stress;  
}

/**************************************************************/
