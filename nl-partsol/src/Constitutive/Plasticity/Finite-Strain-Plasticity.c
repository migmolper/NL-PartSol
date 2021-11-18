#include "nl-partsol.h"

static void elastic_trial(double *,double *,Material);

/**************************************************************/

int finite_strain_plasticity__Constitutive__(
  State_Parameters * ptr_SP_p, 
  Material MatProp,
  int (* infinitesimal_plasticity)(State_Parameters *,Material))
/*
  Finite strains plasticity following the apporach of Ortiz and Camacho
*/
{
  int Ndim = NumberDimensions;
  int status = 0;

  /* Define auxiliar variables */
  Tensor P_p = memory_to_tensor__TensorLib__(ptr_SP_p->Stress,2);
  Tensor F_m1_plastic = memory_to_tensor__TensorLib__(ptr_SP_p->F_m1_plastic_p,2);
  Tensor F_total = memory_to_tensor__TensorLib__(ptr_SP_p->F_n1_p,2);
  Tensor F_m1_total;
  Tensor F_trial_elastic;
  Tensor C_trial_elastic;

  EigenTensor Eigen_C_trial_elastic;

  Tensor D_F_plastic_spectral = alloc__TensorLib__(2);
  Tensor D_F_plastic;

  Tensor F_m1_plastic_new;

  Tensor kirchhoff_spectral_p = alloc__TensorLib__(2);
  Tensor kirchhoff_p;

  /* Set to zero memory */
  #if NumberDimensions == 3
  memset(ptr_SP_p->Stress, 0, 9 * sizeof(double));
  #endif
  #if NumberDimensions == 2
  memset(ptr_SP_p->Stress, 0, 5 * sizeof(double));
  #endif

  /* Define auxiliar memory */
  ptr_SP_p->Strain = (double *)calloc(3,sizeof(double));
  ptr_SP_p->Increment_E_plastic = (double *)calloc(3,sizeof(double));

  /*
    Compute the trial right Cauchy-Green tensor
  */
  if(MatProp.Locking_Control_Fbar)
  {
    Tensor Fbar = memory_to_tensor__TensorLib__(ptr_SP_p->Fbar,2);
    F_trial_elastic = matrix_product__TensorLib__(Fbar,F_m1_plastic);
  }
  else
  {
    F_total = memory_to_tensor__TensorLib__(ptr_SP_p->F_n1_p,2);
    F_trial_elastic = matrix_product__TensorLib__(F_total,F_m1_plastic); 
  }

  C_trial_elastic = right_Cauchy_Green__Particles__(F_trial_elastic);

  /*
    Perform the spectral decomposition of the right Cauchy-Green tensor
  */
  status = Eigen_analysis__TensorLib__(&Eigen_C_trial_elastic,C_trial_elastic);
  if(status)
  {
    fprintf(stderr,"%s %s \n%s %s\n",
    "Error in the function",__func__,
    "File",__FILE__);
    return EXIT_FAILURE;
  }
  
  /*
    Compute the infinitesimal strain magnitude using the logarithmic mapping
  */ 
  for(int i = 0 ; i<Ndim ; i++)
  {
    ptr_SP_p->Strain[i] = 0.5*log(Eigen_C_trial_elastic.Value.n[i]);
  }

  /*
    Calculation of the trial stress tensor using the trial small strain tensor
    with a linear elastic material
  */
  elastic_trial(ptr_SP_p->Stress,ptr_SP_p->Strain, MatProp);

  /*
    Start plastic corrector algorithm in infinitesimal strains
  */
  status = infinitesimal_plasticity(ptr_SP_p, MatProp);
  if(status)
  {
    fprintf(stderr,"%s %s \n%s %s\n",
    "Error in the function",__func__,
    "File",__FILE__);
    return EXIT_FAILURE;
  }

  /*
    Get the increment of plastic finite strains following the Cuitiño & Ortiz exponential maping
  */
  for(int i = 0 ; i<Ndim ; i++)
  {
    D_F_plastic_spectral.N[i][i] = exp( - ptr_SP_p->Increment_E_plastic[i]);
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

  if(I3__TensorLib__(F_m1_plastic) < 0)
  {
    fprintf(stderr,"%s %s: %s %i \n%s %s\n",
    "Error in the function",__func__,
    "negative Jacobian in the particle",ptr_SP_p->Particle_Idx,
    "File",__FILE__);
    return EXIT_FAILURE;
  }

  /*
    Rotate the kirchhoff stress tensor in the spectral representation
    to obtain its value in the cartesian representation. Note that we are employing the same
    directions to rotate D_F_platic as those employed for C_elastic
  */
  for(int i = 0 ; i<Ndim ; i++)
  {
    kirchhoff_spectral_p.N[i][i] = ptr_SP_p->Stress[i];
  }

  kirchhoff_p = rotate__TensorLib__(kirchhoff_spectral_p, Eigen_C_trial_elastic.Vector);

  /*
    Get the First Piola-Kirchhoff stress tensor from the Kirchhoff stress.
  */
  F_m1_total = Inverse__TensorLib__(F_total);


#if NumberDimensions == 2
  ptr_SP_p->Stress[4] = ptr_SP_p->Stress[2];
#endif

  for(int i = 0 ; i < Ndim  ; i++)
  {
    for(int j = 0 ; j < Ndim  ; j++)    
    {
      
      ptr_SP_p->Stress[i*Ndim + j] = 0.0;

      for(int k = 0 ; k < Ndim  ; k++)
      {
        ptr_SP_p->Stress[i*Ndim + j] += kirchhoff_p.N[i][k]*F_m1_total.N[j][k];
      }
    }
  }



  /* Free memory */
  free(ptr_SP_p->Strain);
  free(ptr_SP_p->Increment_E_plastic);
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


  return EXIT_SUCCESS;
}

/**************************************************************/

static void elastic_trial(
  double * Stress,
  double * Strain,
  Material MatProp_p)
{

  double nu = MatProp_p.nu; 
  double E = MatProp_p.E;
  double Lame_param = E*nu/((1 + nu)*(1 - 2*nu));
  double Shear_modulus = E/(2*(1 + nu));

  Stress[0] = (Lame_param + 2*Shear_modulus)*Strain[0] + Lame_param*Strain[1] + Lame_param*Strain[2];

  Stress[1] = Lame_param*Strain[0] + (Lame_param + 2*Shear_modulus)*Strain[1] + Lame_param*Strain[2];

  Stress[2] = Lame_param*Strain[0] + Lame_param*Strain[1] + (Lame_param + 2*Shear_modulus)*Strain[2];
}

/**************************************************************/
