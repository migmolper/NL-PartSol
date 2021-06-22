#include "nl-partsol.h"

/**************************************************************/

State_Parameters finite_strain_plasticity(
  State_Parameters Inputs_SP_finite, 
  Material MatProp,
  State_Parameters(* infinitesimal_plasticity)(State_Parameters,Material))
/*
  Finite strains plasticity following the apporach of Ortiz and Camacho
*/
{
  int Ndim = NumberDimensions;

  /* Define auxiliar variables */
  State_Parameters Output_SP;
  Tensor P_p = memory_to_tensor__TensorLib__(Inputs_SP_finite.Stress,2);
  Tensor F_m1_plastic = memory_to_tensor__TensorLib__(Inputs_SP_finite.F_m1_plastic_p,2);
  Tensor F_total = memory_to_tensor__TensorLib__(Inputs_SP_finite.F_n1_p,2);
  Tensor F_m1_total;
  Tensor F_trial_elastic;
  Tensor C_trial_elastic;
  Tensor E_trial_elastic;
  Tensor D_F_plastic;
  Tensor Fm1_plastic;
  Tensor kirchhoff_p;

  State_Parameters Input_SP_infinitesimal;
  State_Parameters Output_SP_infinitesimal;

  /* Defin the input parameters for the infinitesimal elasticity */
  Input_SP_infinitesimal.EPS = Inputs_SP_finite.EPS;
  Input_SP_infinitesimal.Back_stress = Inputs_SP_finite.Back_stress;
  Input_SP_infinitesimal.Stress = (double *)calloc(Ndim*Ndim,sizeof(double));
  Input_SP_infinitesimal.Strain = (double *)calloc(Ndim*Ndim,sizeof(double));

  /* Compute the elastic right Cauchy-Green tensor using the intermediate configuration. */ 
  F_trial_elastic = matrix_product__TensorLib__(F_total,F_m1_plastic);

  C_trial_elastic = right_Cauchy_Green__Particles__(F_trial_elastic);

  /* Calculation of the small strain tensor */
  E_trial_elastic = logarithmic_strains__Particles__(C_trial_elastic);

  /* Fill memory */
  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      Input_SP_infinitesimal.Strain[i*Ndim + j] = E_trial_elastic.N[i][j];
    }
  }
  
  /* Calculation of the trial stress tensor using the trial small strain tensor */
  Output_SP_infinitesimal = compute_kirchhoff_isotropic_linear_elasticity(Input_SP_infinitesimal, MatProp);

  /* Start plastic corrector algorithm in infinitesimal strains */
  Output_SP_infinitesimal = infinitesimal_plasticity(Input_SP_infinitesimal, MatProp);

  Tensor Increment_E_plastic = memory_to_tensor__TensorLib__(Output_SP_infinitesimal.Increment_E_plastic,2);

  /* Use the Cuitiño & Ortiz exponential maping to compute the increment of plastic finite strains */
  update_plastic_deformation_gradient__Particles__(Increment_E_plastic, F_m1_plastic);

  /* Get the First Piola-Kirchhoff stress tensor (P_p) */
  F_m1_total = Inverse__TensorLib__(F_total);
  kirchhoff_p = memory_to_tensor__TensorLib__(Input_SP_infinitesimal.Stress,2);

  for(int i = 0 ; i < Ndim  ; i++)
  {
    for(int j = 0 ; j < Ndim  ; j++)    
    {
      
      P_p.N[i][j] = 0.0;

     for(int k = 0 ; k < Ndim  ; k++)
     {
        P_p.N[i][j] += kirchhoff_p.N[i][k]*F_m1_total.N[k][j];
      }
    }
  }

  /* Free memory */
  free__TensorLib__(F_trial_elastic);
  free__TensorLib__(C_trial_elastic);
  free__TensorLib__(E_trial_elastic);
  free__TensorLib__(F_m1_total);
  free(Output_SP_infinitesimal.Increment_E_plastic);
  free(Input_SP_infinitesimal.Stress);
  free(Input_SP_infinitesimal.Strain);

  /*
    Update state paramers
  */
  Output_SP.EPS = Output_SP_infinitesimal.EPS;
  Output_SP.Back_stress = Output_SP_infinitesimal.Back_stress;


  return Output_SP;
}

/**************************************************************/
