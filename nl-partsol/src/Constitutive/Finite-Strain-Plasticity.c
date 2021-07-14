#include "nl-partsol.h"

static void elastic_trial(State_Parameters,Material);

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

  EigenTensor Eigen_C_trial_elastic;

  Tensor D_F_plastic_spectral = alloc__TensorLib__(2);
  Tensor D_F_plastic;

  Tensor F_m1_plastic_new;

  Tensor kirchhoff_spectral_p = alloc__TensorLib__(2);
  Tensor kirchhoff_p;

  State_Parameters Input_SP_infinitesimal;
  State_Parameters Output_SP_infinitesimal;

  /* Defin the input parameters for the infinitesimal elasticity */
  Input_SP_infinitesimal.EPS = Inputs_SP_finite.EPS;
  Input_SP_infinitesimal.Back_stress = Inputs_SP_finite.Back_stress;
  Input_SP_infinitesimal.Stress = (double *)calloc(3,sizeof(double));
  Input_SP_infinitesimal.Strain = (double *)calloc(3,sizeof(double));
  Input_SP_infinitesimal.Increment_E_plastic = (double *)calloc(3,sizeof(double));
  

  /*
    Compute the trial right Cauchy-Green tensor
  */
  if(MatProp.Locking_Control_Fbar)
  {    
    Tensor Fbar = memory_to_tensor__TensorLib__(Inputs_SP_finite.Fbar,2);
    F_trial_elastic = matrix_product__TensorLib__(Fbar,F_m1_plastic);
  }
  else
  {
    F_trial_elastic = matrix_product__TensorLib__(F_total,F_m1_plastic); 
  }

  C_trial_elastic = right_Cauchy_Green__Particles__(F_trial_elastic);

  /*
    Perform the spectral decomposition of the right Cauchy-Green tensor
  */
  Eigen_C_trial_elastic = Eigen_analysis__TensorLib__(C_trial_elastic);

  /*
    Compute the infinitesimal strain magnitude using the logarithmic mapping
  */ 
  for(int i = 0 ; i<Ndim ; i++)
  {
    Input_SP_infinitesimal.Strain[i] = 0.5*log(Eigen_C_trial_elastic.Value.n[i]);
  }

  /*
    Calculation of the trial stress tensor using the trial small strain tensor
    with a linear elastic material
  */
  elastic_trial(Input_SP_infinitesimal, MatProp);

  /*
    Start plastic corrector algorithm in infinitesimal strains
  */
  Output_SP_infinitesimal = infinitesimal_plasticity(Input_SP_infinitesimal, MatProp);

  /*
    Get the increment of plastic finite strains following the Cuitiño & Ortiz exponential maping
  */
  for(int i = 0 ; i<Ndim ; i++)
  {
    D_F_plastic_spectral.N[i][i] = exp(- Output_SP_infinitesimal.Increment_E_plastic[i]);
  }

  /*
    Rotate the increment of the plastic deformation gradient in the spectral representation
    to obtain its value in the cartesian representation. Note that we are employing the same
    directions to rotate D_F_platic as those employed for C_elastic
  */
  D_F_plastic = rotate__TensorLib__(D_F_plastic_spectral, Eigen_C_trial_elastic.Vector);

  /*
    Compute the new value of the plastic deformation gradient and update it
  */
  F_m1_plastic_new = matrix_product__TensorLib__(F_m1_plastic,D_F_plastic);

  for(int i = 0 ; i < Ndim  ; i++)
  {
    for(int j = 0 ; j < Ndim  ; j++)
    {
      F_m1_plastic.N[i][j] = F_m1_plastic_new.N[i][j];
    }
  }

  /*
    Rotate the kirchhoff stress tensor in the spectral representation
    to obtain its value in the cartesian representation. Note that we are employing the same
    directions to rotate D_F_platic as those employed for C_elastic
  */
  for(int i = 0 ; i<Ndim ; i++)
  {
    kirchhoff_spectral_p.N[i][i] = Input_SP_infinitesimal.Stress[i];
  }

  kirchhoff_p = rotate__TensorLib__(kirchhoff_spectral_p, Eigen_C_trial_elastic.Vector);

  /*
    Get the First Piola-Kirchhoff stress tensor from the Kirchhoff stress.
  */
  F_m1_total = Inverse__TensorLib__(F_total);

  for(int i = 0 ; i < Ndim  ; i++)
  {
    for(int j = 0 ; j < Ndim  ; j++)    
    {
      
      P_p.N[i][j] = 0.0;

      for(int k = 0 ; k < Ndim  ; k++)
      {
        P_p.N[i][j] += kirchhoff_p.N[i][k]*F_m1_total.N[j][k];
      }
    }
  }


  /* Free memory */
  free(Input_SP_infinitesimal.Stress);
  free(Input_SP_infinitesimal.Strain);
  free(Input_SP_infinitesimal.Increment_E_plastic);
  free__TensorLib__(F_trial_elastic);
  free__TensorLib__(C_trial_elastic);
  free__TensorLib__(Eigen_C_trial_elastic.Value);
  free__TensorLib__(Eigen_C_trial_elastic.Vector);
  free__TensorLib__(D_F_plastic_spectral);
  free__TensorLib__(D_F_plastic);
  free__TensorLib__(kirchhoff_spectral_p);
  free__TensorLib__(kirchhoff_p);
  free__TensorLib__(F_m1_plastic_new);
  free__TensorLib__(F_m1_total);

  /*
    Update state paramers
  */
  Output_SP.EPS = Output_SP_infinitesimal.EPS;
  Output_SP.Back_stress = Output_SP_infinitesimal.Back_stress;


  return Output_SP;
}

/**************************************************************/

static void elastic_trial(
  State_Parameters Intput_SP,
  Material MatProp_p)
{

  int Ndim = NumberDimensions;

  /* Get information from the state parameter */
  Tensor Strain = memory_to_tensor__TensorLib__(Intput_SP.Strain,2);
  Tensor Stress = memory_to_tensor__TensorLib__(Intput_SP.Stress,2);

  double nu = MatProp_p.nu; 
  double E = MatProp_p.E;
  double G = E/(2*(1+nu));
  double Lambda = nu*E/((1+nu)*(1-2*nu));
  double traceStrain = Intput_SP.Strain[0] + Intput_SP.Strain[1] + Intput_SP.Strain[2];

  Intput_SP.Stress[0] = 2*G*Intput_SP.Strain[0] - Lambda*traceStrain;
  Intput_SP.Stress[1] = 2*G*Intput_SP.Strain[1] - Lambda*traceStrain;
  Intput_SP.Stress[2] = 2*G*Intput_SP.Strain[2] - Lambda*traceStrain;
 
}

/**************************************************************/
