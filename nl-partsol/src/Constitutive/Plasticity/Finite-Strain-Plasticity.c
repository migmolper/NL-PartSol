#include "nl-partsol.h"

static void elastic_trial(State_Parameters *,Material);

/**************************************************************/

int finite_strain_plasticity(
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
  memset(ptr_SP_p->Stress, 0, 12 * sizeof(double));

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
  Eigen_C_trial_elastic = Eigen_analysis__TensorLib__(C_trial_elastic);

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
  elastic_trial(ptr_SP_p, MatProp);

  /*
    Start plastic corrector algorithm in infinitesimal strains
  */
  status = infinitesimal_plasticity(ptr_SP_p, MatProp);

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
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in finite_strain_plasticity()",
	    "The Jacobian of the resulting plastic deformation gradient is less than 0");
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

#if NumberDimensions == 3

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

#endif

#if NumberDimensions == 2

  ptr_SP_p->Stress[0] = ;
  ptr_SP_p->Stress[1] = ;

  ptr_SP_p->Stress[2] = ;
  ptr_SP_p->Stress[3] = ;

  ptr_SP_p->Stress[4] = ;

#endif

  /*
    Plane strain conditions
  */
  if(Ndim == 2)
  {
    ptr_SP_p->Stress[4] = ptr_SP_p->Stress[2];
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


  return EXIT_FAILURE;
}

/**************************************************************/

static void elastic_trial(
  State_Parameters * ptr_SP_p,
  Material MatProp_p)
{

  double nu = MatProp_p.nu; 
  double E = MatProp_p.E;
  double Lame_param = E*nu/((1 + nu)*(1 - 2*nu));
  double Shear_modulus = E/(2*(1 + nu));

#if NumberDimensions == 3 

  ptr_SP_p->Stress[0] = (Lame_param + 2*Shear_modulus)*ptr_SP_p->Strain[0] + Lame_param*ptr_SP_p->Strain[1] + Lame_param*ptr_SP_p->Strain[2];

  ptr_SP_p->Stress[1] = Lame_param*ptr_SP_p->Strain[0] + (Lame_param + 2*Shear_modulus)*ptr_SP_p->Strain[1] + Lame_param*ptr_SP_p->Strain[2];

  ptr_SP_p->Stress[2] = Lame_param*ptr_SP_p->Strain[0] + Lame_param*ptr_SP_p->Strain[1] + (Lame_param + 2*Shear_modulus)*ptr_SP_p->Strain[2];

#endif

#if NumberDimensions == 2

  ptr_SP_p->Stress[0] = (Lame_param + 2*Shear_modulus)*ptr_SP_p->Strain[0] + Lame_param*ptr_SP_p->Strain[1];

  ptr_SP_p->Stress[1] = Lame_param*ptr_SP_p->Strain[0] + (Lame_param + 2*Shear_modulus)*ptr_SP_p->Strain[1];
  
  ptr_SP_p->Stress[2] = nu*(ptr_SP_p->Stress[0] + ptr_SP_p->Stress[1]);

#endif
}

/**************************************************************/
