#include "nl-partsol.h"

/*
  Call global variables
*/
double TOL_Radial_Returning;
int Max_Iterations_Radial_Returning;

/*
  Define local global variable for the relative error
*/
double Error0;

/*
  Auxiliar functions 
*/
static void   standard_error(char *);
static bool   check_convergence(double,double,int,int);

static double compute_yield_surface(double,double,double,Material);
static double compute_derivative_yield_surface(double, Material);

static double compute_K(double,Material);
static double compute_DK(double,Material);

static double update_increment_plastic_strain(double, double, double);
static double update_equivalent_plastic_strain(double, double, Material);

static Tensor compute_increment_plastic_strain_tensor(Tensor, double, double, Material);
static void   compute_finite_stress_tensor_plastic_region(Tensor, Tensor, Tensor, double, double, Material);
static void   compute_finite_stress_tensor_elastic_region(Tensor, Tensor, Tensor, Material);

/**************************************************************/

Plastic_status finite_strains_plasticity_Von_Mises(
  Tensor P_p,
  Tensor F_m1_plastic,
  Tensor F_total,
  Plastic_status Inputs_VarCons, 
  Material MatProp,
  double J)
/*
  Finite strains plasticity following the apporach of Ortiz and Camacho
*/
{
  int Ndim = NumberDimensions;

  /* Define auxiliar variables */
  Plastic_status Outputs_VarCons;
  Tensor F_trial_elastic;
  Tensor C_trial_elastic;
  Tensor E_trial_elastic;
  Tensor F_elastic;
  Tensor C_elastic;
  Tensor C_m1_elastic;
  Tensor Infinitesimal_Stress = alloc__TensorLib__(2);
  Tensor Increment_E_plastic;
  Tensor D_F_plastic;
  Tensor Fm1_plastic;
  Tensor M_p = alloc__TensorLib__(2);
  Tensor S_p = alloc__TensorLib__(2);

  /* Compute the elastic right Cauchy-Green tensor using the intermediate configuration. */ 
  F_trial_elastic = matrix_product__TensorLib__(F_total,F_m1_plastic);

  C_trial_elastic = right_Cauchy_Green__Particles__(F_trial_elastic);

  /* Calculation of the small strain tensor */
  E_trial_elastic = logarithmic_strains__Particles__(C_trial_elastic);

  /* Start plastic algorithm in infinitesimal strains */
  Outputs_VarCons = infinitesimal_strains_plasticity_Von_Mises(Infinitesimal_Stress, E_trial_elastic, Inputs_VarCons, MatProp);

  /* Update the logarithmic strain tensor */
  for(int i = 0 ; i < Ndim  ; i++)
  {
   for(int j = 0 ; j < Ndim  ; j++)
    {
      E_trial_elastic.N[i][j] -= Outputs_VarCons.Increment_E_plastic.N[i][j];
    }
  }

  /* Use the Cuitiño & Ortiz exponential maping to compute the increment of plastic finite strains */
  update_plastic_deformation_gradient__Particles__(Outputs_VarCons.Increment_E_plastic, F_m1_plastic);

  /* Compute the elastic stress tensor */
  F_elastic = matrix_product__TensorLib__(F_total,F_m1_plastic);

  /* Compute the inverse of the elastic right Cauchy-Green tensor */
  C_elastic = right_Cauchy_Green__Particles__(F_elastic);
  C_m1_elastic = Inverse__TensorLib__(C_elastic);

  /* Compute the Mandel stress tensor */
  for(int i = 0 ; i < Ndim  ; i++)
  {
   for(int j = 0 ; j < Ndim  ; j++)
     {
      /* Symmetric part */
      M_p.N[i][j] += Infinitesimal_Stress.N[i][j];

      /* Kew symetric part */
      for(int k = 0 ; k < Ndim  ; k++)
      {
        M_p.N[i][j] += E_trial_elastic.N[i][k]*Infinitesimal_Stress.N[k][j] - 
                       Infinitesimal_Stress.N[i][k]*E_trial_elastic.N[k][j];
      }
    }
  }

  /* Get the Second Piola-Kirchhoff stress tensor (S_p) */
  for(int i = 0 ; i < Ndim  ; i++)
  {
   for(int j = 0 ; j < Ndim  ; j++)
     {
      for(int k = 0 ; k < Ndim  ; k++)
      {
        S_p.N[i][j] += 0.5*(C_m1_elastic.N[i][k]*M_p.N[k][j] + M_p.N[k][i]*C_m1_elastic.N[k][j]);
      }
    }
  }

  /* Get the First Piola-Kirchhoff stress tensor (P_p) */
  for(int i = 0 ; i < Ndim  ; i++)
  {
   for(int j = 0 ; j < Ndim  ; j++)
     {
      P_p.N[i][j] = 0.0;

      for(int k = 0 ; k < Ndim  ; k++)
      {
        P_p.N[i][j] += F_total.N[i][k]*S_p.N[k][j];
      }
    }
  }

  /* Free memory */
  free__TensorLib__(F_trial_elastic);
  free__TensorLib__(C_trial_elastic);
  free__TensorLib__(E_trial_elastic);
  free__TensorLib__(F_elastic);
  free__TensorLib__(C_elastic);
  free__TensorLib__(C_m1_elastic);
  free__TensorLib__(Infinitesimal_Stress);
  free__TensorLib__(M_p);
  free__TensorLib__(S_p);
  free__TensorLib__(Outputs_VarCons.Increment_E_plastic);

  return Outputs_VarCons;
}

/**************************************************************/

Plastic_status infinitesimal_strains_plasticity_Von_Mises(
  Tensor sigma_k1,
  Tensor E_elastic,
  Plastic_status Inputs_VarCons,
  Material MatProp)
/*	
	Radial returning algorithm for the Von-Mises plastic criterium
*/
{

  int Ndim = NumberDimensions;

  double E_elastic_vol;
  Tensor E_elastic_dev;
  Tensor p_trial;
  Tensor s_trial;
  double s_trial_norm;
  double delta_Gamma;
  double Phi;
  double d_Phi;	
  Tensor Increment_E_plastic = alloc__TensorLib__(2);

  /* Load material marameters */
  double EPS = Inputs_VarCons.EPS;
  double EPS_k;
  double yield_stress_k = Inputs_VarCons.Yield_stress;

  /*
    Initialise convergence parameters for the solver
  */
  double TOL = TOL_Radial_Returning;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;
  bool Convergence = false;
  
  /*
    Elastic predictor : Volumetric and deviatoric stress measurements. Compute also
    the norm of the deviatoric tensor
  */
  E_elastic_vol = volumetric_component__TensorLib__(E_elastic);
  E_elastic_dev = deviatoric_component__TensorLib__(E_elastic,E_elastic_vol);

  p_trial = volumetric_stress__LinearElastic__(E_elastic_vol, MatProp);
  s_trial = deviatoric_stress__LinearElastic__(E_elastic_dev, MatProp);  

  s_trial_norm = EuclideanNorm__TensorLib__(s_trial); 

  /*
    Yield condition : Starting from incremental plastic strain equal to zero
  */
  EPS_k = EPS;
  delta_Gamma = 0;
  Phi = compute_yield_surface(s_trial_norm, delta_Gamma, yield_stress_k, MatProp);

  if(Phi < TOL)
  {
    compute_finite_stress_tensor_elastic_region(sigma_k1, s_trial, p_trial, MatProp);
  }
  else
  {     
    /*
      Newton-Rapson solver
    */
    while(Convergence == false)
    {
          
      Phi = compute_yield_surface(s_trial_norm, delta_Gamma, EPS_k, MatProp);

      d_Phi = compute_derivative_yield_surface(delta_Gamma, MatProp);

      delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

      EPS_k = update_equivalent_plastic_strain(EPS, delta_Gamma, MatProp);

      Convergence = check_convergence(Phi,TOL,Iter,MaxIter);
	     
      if(Convergence == false)
      {
        Iter++;
      }
    }

    /*
      Update plastic deformation gradient
    */
    Increment_E_plastic = compute_increment_plastic_strain_tensor(s_trial, s_trial_norm, delta_Gamma, MatProp);

    /*
      Update stress tensor in the deformed configuration
    */
    compute_finite_stress_tensor_plastic_region(sigma_k1, s_trial, p_trial, s_trial_norm, delta_Gamma, MatProp);

  }

  /*
    Free memory
  */
  free__TensorLib__(E_elastic_dev);
  free__TensorLib__(s_trial);
  free__TensorLib__(p_trial);

  /*
    Define output varible
  */
  Plastic_status Outputs_VarCons;
  Outputs_VarCons.EPS = EPS_k;
  Outputs_VarCons.Increment_E_plastic = Increment_E_plastic;
  
return Outputs_VarCons;
}

/***************************************************************************/

static void standard_error(char * Error_message)
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in plasticity_Von_Mises",Error_message);
    exit(EXIT_FAILURE);
}

/**************************************************************/

static bool check_convergence(double Error, double TOL, int Iter, int MaxIter)
{
  bool convergence = false;
  double Error_relative;
  char Error_message[MAXW];

  if(Iter > MaxIter)
    {
      sprintf(Error_message,"%s","Convergence not reached in the maximum number of iterations");
      printf("%e\n",Error_relative);
      standard_error(Error_message); 
    }

  /*
    Compute relative error
  */
  if(Iter == 0)
    {
      Error0 = Error;
      Error_relative = Error/Error0;
    }
    else
    {
      Error_relative = Error/Error0;
    }
      
    /*
      Check convergence using the relative error
    */
  if(Error_relative < TOL)
    {
      convergence = true;
    }

  return convergence;
}

/**************************************************************/

static double compute_yield_surface(
  double s_trial_norm,
  double delta_Gamma,
  double EPS_k,
  Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double K = compute_K(EPS_k,MatProp);
  double phi = - sqrt(2./3.)*K + s_trial_norm - 2*G*delta_Gamma;

  return phi;
}

/**************************************************************/


static double compute_derivative_yield_surface(
  double delta_Gamma,
  Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double DK = compute_DK(delta_Gamma,MatProp);

  return - 2*G - sqrt(2./3.)*DK;
}

/**************************************************************/

static double compute_K(
  double EPS_k1,
  Material MatProp)
/*
  Isotropic hardening function
*/
{
  double H = MatProp.hardening_modulus;
  double Sigma_y = MatProp.yield_stress_0;

  double K = Sigma_y + H*EPS_k1;

  return K;
}

/**************************************************************/

static double compute_DK(
  double delta_Gamma,
  Material MatProp)
/*
  derivative of the isotropic hardening function
*/
{
  double H = MatProp.hardening_modulus;

  double DK = sqrt(2./3.)*H*delta_Gamma;

  return DK;
}

/**************************************************************/

static double update_increment_plastic_strain(
  double delta_Gamma_k,
  double Phi,
  double d_Phi)
{
  double delta_Gamma_k1 = delta_Gamma_k - Phi/d_Phi;
  return delta_Gamma_k1;
}

/**************************************************************/

static double update_equivalent_plastic_strain(
  double EPS,
  double delta_Gamma,
  Material MatProp)
{
  double EPS_k = EPS + sqrt(2./3.)*delta_Gamma;
  return EPS_k;
}

/**************************************************************/

static Tensor compute_increment_plastic_strain_tensor(
  Tensor s_trial, 
  double s_trial_norm, 
  double delta_Gamma, 
  Material MatProp)
{
  int Ndim = NumberDimensions;
  double aux = delta_Gamma/s_trial_norm;
  Tensor D_E_plastic = alloc__TensorLib__(2); 

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        D_E_plastic.N[i][j] = aux*s_trial.N[i][j];
      }
    }


  return D_E_plastic;
}

/**************************************************************/

static void compute_finite_stress_tensor_plastic_region(Tensor sigma_k1, Tensor s_trial, Tensor p_trial, double s_trial_norm, double delta_Gamma, Material MatProp)
{
  int Ndim = NumberDimensions;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double aux = 1 - 2*G*delta_Gamma/s_trial_norm;

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        sigma_k1.N[i][j] = p_trial.N[i][j] + aux*s_trial.N[i][j];
      }
    }
}

/**************************************************************/

static void compute_finite_stress_tensor_elastic_region(Tensor sigma_k1, Tensor s_trial, Tensor p_trial, Material MatProp)
{
  int Ndim = NumberDimensions;

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        sigma_k1.N[i][j] = p_trial.N[i][j] + s_trial.N[i][j];
      }
    }

}

/**************************************************************/
