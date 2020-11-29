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
static double compute_yield_surface(double, double, double, Material);
static double compute_derivative_yield_surface(double, Material);
static double update_increment_plastic_strain(double, double, double);
static double update_equivalent_plastic_strain(double, double, Material);
static double update_yield_stress(double, Material);
static Tensor compute_increment_plastic_strain_tensor(Tensor, double, double, Material);
static Tensor compute_finite_stress_tensor_plastic_region(Tensor, Tensor, double, double, Material);
static Tensor compute_finite_stress_tensor_elastic_region(Tensor, Tensor, Material);

/**************************************************************/

Variables_Constitutive plasticity_Von_Mises(Tensor S_p, Tensor C_total, Tensor F_plastic, Tensor F_total, 
                					                  double J, Variables_Constitutive Inputs_VarCons, Material MatProp)
/*	
	Radial returning algorithm for the Von-Mises plastic criterium
*/
{

  int Ndim = NumberDimensions;

  double EPS_k;
  double yield_stress_k;

  Tensor C_elastic;
  Tensor E_elastic;
  double E_elastic_vol;
  Tensor E_elastic_dev;
  Tensor p_trial;
  Tensor s_trial;
  double s_trial_norm;
  double delta_Gamma;
  double Phi;
  double d_Phi;	
  double H;
  Tensor sigma_k1;
  Tensor Increment_E_plastic;

  /*
    Initialise convergence parameters for the solver
  */
  double TOL = TOL_Radial_Returning;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;
  bool Convergence = false;

  /*
    Define output varible
  */
  Variables_Constitutive Outputs_VarCons;

  /*  
    Get the value of the equivalent plastic stress and yield stress
  */
  EPS_k = Inputs_VarCons.EPS;
  yield_stress_k  = Inputs_VarCons.Yield_stress;

  

  /*
    Compute the trial elastic right Cauchy-Green tensor using the intermediate configuration.
  */
  C_elastic = alloc__TensorLib__(2); 
  covariant_push_forward_tensor__TensorLib__(C_elastic, C_total, F_plastic);

  /*  
    Calculation of the small strain tensor (the approach of Ortiz and Camacho)
  */
  E_elastic = logarithmic_strains__Particles__(C_elastic);
  
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
  delta_Gamma = 0;
  Phi = compute_yield_surface(s_trial_norm, delta_Gamma, yield_stress_k, MatProp);

  if(Phi >= TOL)
    {

      H = MatProp.hardening_modulus;
          
      /*
        Newton-Rapson solver
      */
      while(Convergence == false)
        {

          d_Phi = compute_derivative_yield_surface(H, MatProp);

          delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

          EPS_k = update_equivalent_plastic_strain(EPS_k, delta_Gamma, MatProp);

          yield_stress_k = update_yield_stress(EPS_k, MatProp);

          Phi = compute_yield_surface(s_trial_norm, delta_Gamma, yield_stress_k, MatProp);

          /*
		        Check convergence
	        */
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

      update_plastic_deformation_gradient__Particles__(Increment_E_plastic,F_plastic);

      /*
        Free increment of E_plastic
      */
      free__TensorLib__(Increment_E_plastic);

      /*
        Update stress tensor in the deformed configuration
      */
      sigma_k1 = compute_finite_stress_tensor_plastic_region(s_trial, p_trial, s_trial_norm, delta_Gamma, MatProp);

    }
    else
    {
      sigma_k1 = compute_finite_stress_tensor_elastic_region(s_trial, p_trial, MatProp);
    }

  /*
    Get the stress tensor in the reference configuration
  */
  contravariant_pull_back_tensor__TensorLib__(S_p, sigma_k1, F_total);
	
  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        S_p.N[i][j] *= J;
      }
    }

  /*
    Free memory
  */
  free__TensorLib__(C_elastic);
  free__TensorLib__(E_elastic);
  free__TensorLib__(E_elastic_dev);
  free__TensorLib__(s_trial);
  free__TensorLib__(p_trial);
  free__TensorLib__(sigma_k1);

  /*
    Return yield stress and EPS updated
  */
  Outputs_VarCons.EPS = EPS_k;
  Outputs_VarCons.Yield_stress = yield_stress_k;
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
  bool convergence;
  double Error_relative;
  char Error_message[MAXW];

  if(Iter > MaxIter)
    {
      sprintf(Error_message,"%s","Convergence not reached in the maximum number of iterations");
      standard_error(Error_message); 
    }
  else
    {

      /*
        Compute relative error
      */
      if(Iter == 0)
      {
        Error0 = Error;
        Error_relative = Error/Error0;
//        printf("Error iter %i : %1.4e ; %1.4e \n",Iter,Error,Error_relative);
      }
      else
      {
        Error_relative = Error/Error0;
//        printf("Error iter %i : %1.4e ; %1.4e \n",Iter,Error,Error_relative);
      }
      
      /*
        Check convergence using the relative error
      */
      if(Error_relative > TOL)
      {
        return false;
      }
      else
      {
        return true;
      }
    }
}

/**************************************************************/


static double compute_derivative_yield_surface(double H, Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));

  return - 2*G - H*sqrt(2/3.);
}

/**************************************************************/

static double update_increment_plastic_strain(double delta_Gamma_k, double Phi, double d_Phi)
{
  return delta_Gamma_k - Phi/d_Phi;
}

/**************************************************************/

static double update_equivalent_plastic_strain(double EPS_k, double delta_Gamma, Material MatProp)
{
  return EPS_k + delta_Gamma;
}

/**************************************************************/

static double update_yield_stress(double EPS_k1, Material MatProp)
{
  double yield_stress_0 = MatProp.yield_stress_0;
  double E_p0  = MatProp.E_plastic_reference;
  double aux = 1 + EPS_k1/E_p0;

  return yield_stress_0*aux;  
}

/**************************************************************/

static double compute_yield_surface(double s_trial_norm, double delta_Gamma, double yield_stress_k1, Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));

  return s_trial_norm - 2*G*delta_Gamma - sqrt(2/3.)*yield_stress_k1;
}

/**************************************************************/


static Tensor compute_increment_plastic_strain_tensor(Tensor s_trial, double s_trial_norm, double delta_Gamma, Material MatProp)
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

static Tensor compute_finite_stress_tensor_plastic_region(Tensor s_trial, Tensor p_trial, double s_trial_norm, double delta_Gamma, Material MatProp)
{
  int Ndim = NumberDimensions;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double aux = 1 - 2*G*delta_Gamma/s_trial_norm;
  Tensor sigma_k1 = alloc__TensorLib__(2);

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        sigma_k1.N[i][j] = p_trial.N[i][j] + aux*s_trial.N[i][j];
      }
    }

  return sigma_k1;
}

/**************************************************************/

static Tensor compute_finite_stress_tensor_elastic_region(Tensor s_trial, Tensor p_trial, Material MatProp)
{
  int Ndim = NumberDimensions;

  Tensor sigma_k1 = alloc__TensorLib__(2);

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        sigma_k1.N[i][j] = p_trial.N[i][j] + s_trial.N[i][j];
      }
    }

  return sigma_k1;
}

/**************************************************************/
