#include "nl-partsol.h"

/*
  Call global variables
*/
double TOL_Von_Mises;
int Max_Iterations_Von_Mises;

/*
  Auxiliar functions 
*/
static Tensor compute_small_strain_tensor(Tensor, Tensor);
static void   compute_volumetric_deviatoric_stress_tensor(double *, Tensor, Tensor, Material);
static double compute_yield_surface(double, double, double, Material);
static double compute_derivative_yield_surface(double, Material);
static double update_increment_plastic_strain(double, double, double);
static double update_equivalent_plastic_strain(double, double, Material);
static double update_yield_stress(double, Material);
static Tensor compute_increment_plastic_strain_tensor(Tensor, double, double, Material);
static Tensor compute_finite_stress_tensor(Tensor, double, double, double, Material);
extern void   update_plastic_deformation_gradient(Tensor, Tensor);

/**************************************************************/

Tensor plasticity_Von_Mises(Tensor grad_e, Tensor C, Tensor F_plastic, Tensor F, 
                					  double * ptr_EPS_k, double * ptr_yield_stress_k, double J, Material MatProp)
/*	
	Radial returning algorithm for the Von-Mises plastic criterium
*/
{

  int Ndim = NumberDimensions;
  int iter_NR = 0;

  Tensor E_elastic;
  double p_trial;
  Tensor s_trial;
  double s_trial_norm;
  double delta_Gamma;
  double Phi;
  double d_Phi;	
  double H;
  double EPS_k;
  Tensor sigma_k1;
  Tensor D_E_plastic;

  /*
    Initialise convergence parameters
  */
  double TOL_VM = TOL_Von_Mises;
  int iter_max_VM = Max_Iterations_Von_Mises;

  /*	
	  Calculation of the small strain tensor
  */
  E_elastic = compute_small_strain_tensor(Tensor C, Tensor F_plastic);

  /*
    Elastic predictor : Volumetric and deviatoric stress measurements. Compute also
    the norm of the deviatoric tensor
  */
  compute_volumetric_deviatoric_stress_tensor(&p_trial, s_trial, E_elastic, MatProp);

  s_trial_norm = EuclideanNorm__TensorLib__(s_trial);

  /*
    Yield condition : Starting from incremental plastic strain equal to zero
  */
  delta_Gamma = 0;
  Phi = compute_yield_surface(s_trial_norm, delta_Gamma, *ptr_yield_stress_k, MatProp);

  if(Phi >= TOL_VM)
    {

      H = MatProp.hardening_modulus;
          
      /*
        Newton-Rapson solver
      */
      while(iter_NR < iter_max_VM)
        {

          d_Phi = compute_derivative_yield_surface(H, MatProp);

          delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

          *ptr_EPS_k = update_equivalent_plastic_strain(*ptr_EPS_k, delta_Gamma, MatProp);

          *ptr_c_k = update_yield_stress(*ptr_EPS_k, MatProp);

          Phi = compute_yield_surface(s_trial_norm, delta_Gamma, *ptr_yield_stress_k, MatProp);

          /*
		        Check convergence
	        */
	        if(Phi < TOL_VM)
            {
              iter_NR++;
            }
          else
            {
              break;
            }
        }

      /*
        Update plastic deformation gradient
      */
	    D_E_plastic = compute_increment_plastic_strain_tensor(s_trial, s_trial_norm, delta_Gamma, MatProp);

      update_plastic_deformation_gradient(D_E_plastic,F_plastic);

    }

  /*
    Update stress tensor in the deformed configuration
  */
  sigma_k1 = compute_finite_stress_tensor(s_trial, p_trial, s_trial_norm, delta_Gamma, MatProp);

  /*
    Get the stress tensor in the reference configuration
  */
  push_backward_tensor__Particles__(grad_e, sigma_k1, F);
	
  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        grad_e.N[i][j] *= J;
      }
    }

  /*
    Free memory
  */
  free__TensorLib__(E_elastic);
  free__TensorLib__(s_trial);
  free__TensorLib__(sigma_k1);
  free__TensorLib__(Increment_E_plastic);

  /*
    Return stress tensor
  */
  return grad_e;
}

/**************************************************************/

static Tensor compute_small_strain_tensor(Tensor C, Tensor F_plastic)
{
  Tensor C_elastic;
  Tensor E_elastic;

  /*
    Compute the trial elastic right Cauchy-Green tensor using the intermediate configuration.
  */
  push_backward_tensor__Particles__(C_elastic, C, F_plastic);

  /*
    Use the approach of Ortiz and Camacho to compute the elastic infinitesimal strain tensor.
  */
  Eps_elastic = logarithmic_strains__Particles__(C_elastic);

  /*  
  Free memory
  */
  free__TensorLib__(C_elastic);


  return E_elastic;
}

/**************************************************************/

static void compute_volumetric_deviatoric_stress_tensor(double * p_trial, Tensor s_trial, Tensor E_elastic, Material MatProp)
{

  int Ndim = NumberDimensions;
  double nu = MatProp.nu; 
  double E = MatProp.E;
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double E_elastic_vol;
  Tensor E_elastic_dev;

  /*
    Split the deformation tensor into deviatoric and volumetric
  */
  volumetric_desviatoric_decomposition__TensorLib__(E_elastic,E_elastic_vol,E_elastic_dev);

  /*
    Compute deviatoric and volumetric stress tensor
  */
  p_trial = &(K*E_elastic_vol);

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        s_trial.N[i][j] = 2*G*E_elastic_dev.N[i][j];
      }
    }

  free__TensorLib__(E_elastic_dev);
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
  double N_exp = MatProp.hardening_exp;
  double basis = 1 + EPS_k1/E_p0;
  double exp   = 1/N_exp;

  return yield_stress_0*pow(basis,exp);  
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

static Tensor compute_finite_stress_tensor(Tensor s_trial, double p_trial, double s_trial_norm, double delta_Gamma, Material MatProp)
{
  int Ndim = NumberDimensions;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double aux = 1 - 2*G*delta_Gamma/s_trial_norm;
  Tensor sigma_k1 = alloc__TensorLib__(2);
  Tensor I = Identity__TensorLib__();

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        sigma_k1.N[i][j] = p_trial*I.N[i][j] + aux*s_trial.N[i][j];
      }
    }

  return sigma_k1;
}

/**************************************************************/

extern void update_plastic_deformation_gradient(Tensor D_E_plastic, Tensor F_plastic)
{
  int Ndim = NumberDimensions;

  /*
    Use the Cuitiño & Ortiz exponential maping
  */
  Tensor D_F_plastic = increment_Deformation_Gradient_exponential_strains__Particles__(D_E_plastic);

  /*
    Compute the new value of the plastic deformation gradient 
  */
  Tensor Aux_tensor = vector_linear_mapping__TensorLib__(D_F_plastic, F_plastic);

  /*
    Update value of the plastic deformation gradient
  */
  for(int i = 0 ; i < Ndim  ; i++)
    {
      for(int j = 0 ; j < Ndim  ; j++)
      {
        F_plastic.N[i][j] = Aux_tensor.N[i][j];
      }
    }

  /*
    Free memory 
  */
  free__TensorLib__(D_F_plastic);
  free__TensorLib__(Aux_tensor);

}

/**************************************************************/
