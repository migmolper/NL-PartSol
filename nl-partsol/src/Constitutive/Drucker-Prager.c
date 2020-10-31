#include "nl-partsol.h"


/*
Call global variables
*/
double TOL_Drucker_Prager;
int Max_Iterations_Drucker_Prager;

/*
  Auxiliar functions 
*/
static Tensor compute_small_strain_tensor(Tensor, Tensor);
static Tensor compute_volumetric_stress_tensor(double, Material);
static Tensor compute_deviatoric_stress_tensor(Tensor, Material);
static double compute_yield_surface(double, double, double, Material);
static double compute_hardening_modulus(double, Material);
static double compute_limit_between_classic_apex_algorithm(double, double, double, Material);
static double compute_derivative_yield_surface_classical(double, Material);
static double update_increment_plastic_strain(double, double, double);
static double update_equivalent_plastic_strain_classical(double, double, Material);
static double update_cohesion_modulus(double, Material);
static double compute_yield_surface_classical(double, double, double, double, Material);
static Tensor compute_increment_plastic_strain_classical(Tensor, double, double, Material);
static Tensor compute_finite_stress_tensor_classical(Tensor, Tensor, double, double, Material);
static double compute_derivative_yield_surface_apex(double, double, double, Material);
static double compute_yield_surface_apex(double, double, double, double, double, Material);
static double update_equivalent_plastic_strain_apex(double, double, double, Material);
static Tensor compute_increment_plastic_strain_apex(Tensor, double, double, Material);
static Tensor compute_finite_stress_tensor_apex(Tensor, double, Material);
static void   update_plastic_deformation_gradient(Tensor, Tensor);

/**************************************************************/

Tensor plasticity_Drucker_Prager_Sanavia(Tensor grad_e, Tensor C_total, Tensor F_plastic, Tensor F, 
                            					   double * ptr_EPS, double * ptr_c, double J, Material MatProp)
/*	
	Radial returning algorithm
*/
{

  int Ndim = NumberDimensions;
  int iter_NR = 0;

  Tensor E_elastic;
  double E_elastic_vol;
  Tensor E_elastic_dev;
  Tensor p_trial;
  Tensor s_trial;
  double delta_Gamma;
  double s_trial_norm;
  double p_trial_norm;
  double Phi;
  double d_Phi;	
  double H;
  double EPS_k;
  double c_k;
  double p_lim;
  Tensor sigma_k1;
  Tensor D_E_plastic;

  double TOL_DP = TOL_Drucker_Prager;
  int iter_max_DP = Max_Iterations_Drucker_Prager;

  /*  
    Get the value of the equivalent plastic stress and cohesion
  */
  EPS_k = *ptr_EPS;
  c_k  = *ptr_c;

  /*	
	calculation of the small strain tensor
  */
  E_elastic = compute_small_strain_tensor(C_total, F_plastic);

  /*
    Elastic predictor : Volumetric and deviatoric stress measurements. Compute also
    the norm of the deviatoric tensor
  */
  E_elastic_vol = volumetric_component__TensorLib__(E_elastic);
  E_elastic_dev = deviatoric_component__TensorLib__(E_elastic,E_elastic_vol);

  p_trial = compute_volumetric_stress_tensor(E_elastic_vol, MatProp);
  s_trial = compute_deviatoric_stress_tensor(E_elastic_dev, MatProp);  

  s_trial_norm = EuclideanNorm__TensorLib__(s_trial); 
  p_trial_norm = p_trial.N[0][0];

  /*
    Yield condition : Starting from incremental plastic strain equal to zero
  */
  delta_Gamma = 0;

  Phi = compute_yield_surface(p_trial_norm, s_trial_norm, c_k, MatProp);

  if(Phi >= 0)
    {

      H = compute_hardening_modulus(EPS_k, MatProp);

      p_lim = compute_limit_between_classic_apex_algorithm(s_trial_norm, c_k, H, MatProp);

      /*
	Classical plastic iterator
      */
      if(p_trial_norm <= p_lim)
	{
	  while(iter_NR < iter_max_DP)
	    {

	      d_Phi = compute_derivative_yield_surface_classical(H, MatProp);

	      delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

	      EPS_k = update_equivalent_plastic_strain_classical(EPS_k, delta_Gamma, MatProp);

        c_k = update_cohesion_modulus(EPS_k, MatProp);

        H = compute_hardening_modulus(EPS_k, MatProp);

	      Phi = compute_yield_surface_classical(s_trial_norm, p_trial_norm, delta_Gamma, c_k, MatProp);

	      /*
		Check convergence
	      */
	      if(Phi < TOL_DP)
		{
		  iter_NR++;
		}
	      else
		{
		  break;
		}

	    }

	  /*
	    update
	  */
	  D_E_plastic = compute_increment_plastic_strain_classical(s_trial, s_trial_norm, delta_Gamma, MatProp);

	  sigma_k1 = compute_finite_stress_tensor_classical(s_trial, p_trial, s_trial_norm, delta_Gamma, MatProp);

	}
      /*
	Apex plastic iterator
      */
      else
	{

	  while(iter_NR < iter_max_DP)
	    {

	      d_Phi = compute_derivative_yield_surface_apex(H, s_trial_norm, delta_Gamma, MatProp);

	      delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

	      Phi = compute_yield_surface_apex(p_trial_norm, s_trial_norm, H, delta_Gamma, c_k, MatProp);

	      /*
		Check convergence
	      */
	      if(Phi < TOL_DP)
		{
		  iter_NR++;
		}
	      else
		{
		  break;
		}

	    }

	  /*
	    update
	  */
	  EPS_k = update_equivalent_plastic_strain_apex(EPS_k, s_trial_norm, delta_Gamma, MatProp);

	  D_E_plastic = compute_increment_plastic_strain_apex(s_trial, s_trial_norm, delta_Gamma,MatProp);

	  sigma_k1 = compute_finite_stress_tensor_apex(p_trial, delta_Gamma, MatProp);

	}

    }

  /*
    Update plastic deformation gradient
  */
  update_plastic_deformation_gradient(D_E_plastic,F_plastic);

  /*
    Get the stress tensor in the deformed configuration
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
    Update equivalent plastic stress and cohesion
  */
  *ptr_EPS = EPS_k;
  *ptr_c = c_k;

  /*
    Free memory
  */
  free__TensorLib__(E_elastic);
  free__TensorLib__(E_elastic_dev);
  free__TensorLib__(s_trial);
  free__TensorLib__(p_trial);
  free__TensorLib__(sigma_k1);
  free__TensorLib__(D_E_plastic);

  /*
    Return stress tensor
  */
  return grad_e;
}

/**************************************************************/

static Tensor compute_small_strain_tensor(Tensor C_total, Tensor F_plastic)
{
  Tensor C_elastic;
  Tensor E_elastic;

  /*
    Compute the trial elastic right Cauchy-Green tensor using the intermediate configuration.
  */
  push_backward_tensor__Particles__(C_elastic, C_total, F_plastic);

  /*
    Use the approach of Ortiz and Camacho to compute the elastic infinitesimal strain tensor.
  */
  E_elastic = logarithmic_strains__Particles__(C_elastic);

  /*	
	Free memory
  */
  free__TensorLib__(C_elastic);


  return E_elastic;
}

/**************************************************************/

static Tensor compute_volumetric_stress_tensor(double E_elastic_vol, Material MatProp)
{

  int Ndim = NumberDimensions;
  double nu = MatProp.nu; 
  double E = MatProp.E;
  double K = E/(3*(1-2*nu));
  double aux = K*E_elastic_vol;
  Tensor p_trial = alloc__TensorLib__(2);

  /*
    Compute deviatoric stress tensor
  */
  for(int i = 0 ; i<Ndim ; i++)
  {
    p_trial.N[i][i] = aux;
  }

  return p_trial;
}

/**************************************************************/

static Tensor compute_deviatoric_stress_tensor(Tensor E_elastic_dev, Material MatProp)
{

  int Ndim = NumberDimensions;
  double nu = MatProp.nu; 
  double E = MatProp.E;
  double G = E/(2*(1+nu));
  Tensor s_trial = alloc__TensorLib__(2);
  /*
    Compute deviatoric stress tensor
  */

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        s_trial.N[i][j] = 2*G*E_elastic_dev.N[i][j];
      }
    }

  return s_trial;
}

/**************************************************************/

static double compute_yield_surface(double p_trial_norm, double s_trial_norm, double c_k, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;

  double Phi = 3*alpha_F*p_trial_norm + s_trial_norm - beta*c_k;

  return Phi;
}

/**************************************************************/

static double compute_hardening_modulus(double EPS_k, Material MatProp)
{
  double E0_p  = MatProp.E_plastic_reference;
  double N_exp = MatProp.hardening_exp;
  double c0    = MatProp.cohesion_reference;

  double H = (c0/(N_exp*E0_p))*pow((1 + EPS_k/E0_p),(1/N_exp - 1));

  return H;
}

/**************************************************************/

static double compute_limit_between_classic_apex_algorithm(double s_trial_norm, double c_k, double H, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; 
  double E = MatProp.E;
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));

  double p_lim = 3*alpha_Q*K*s_trial_norm/(2*G) + beta/(3*alpha_F)*(s_trial_norm*H/(2*G)*sqrt(1+3*DSQR(alpha_Q)) + c_k);

  return p_lim;

}

/**************************************************************/

static double compute_derivative_yield_surface_classical(double H, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));

  double d_Phi = - 9*K*alpha_F*alpha_Q - 2*G - H*beta*sqrt(3*DSQR(alpha_Q) + 1);

  return d_Phi;
}

/**************************************************************/

static double update_increment_plastic_strain(double delta_Gamma_k, double Phi, double d_Phi)
{
  double delta_Gamma_k1 = delta_Gamma_k - Phi/d_Phi;

  return delta_Gamma_k1;
}

/**************************************************************/

static double update_equivalent_plastic_strain_classical(double EPS_k, double delta_Gamma, Material MatProp)
{
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;

  double EPS_k1 = EPS_k + delta_Gamma*sqrt(3*DSQR(alpha_Q) + 1);

  return EPS_k1;
}

/**************************************************************/

static double update_cohesion_modulus(double EPS_k1, Material MatProp)
{
  double c_0   = MatProp.cohesion_reference;
  double E0_p  = MatProp.E_plastic_reference;
  double N_exp = MatProp.hardening_exp;
  double exp   = 1/N_exp - 1;
  double basis = 1 + EPS_k1/E0_p;

  double c_k1 = (c_0/(N_exp*E0_p))*pow(basis,exp);

  return c_k1;	
}


/**************************************************************/

static double compute_yield_surface_classical(double s_trial_norm, double p_trial_norm, double delta_Gamma, double c_k1, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));

  double Phi = s_trial_norm - 2*G*delta_Gamma + 3*alpha_F*(p_trial_norm - 3*K*alpha_Q*delta_Gamma) - beta*c_k1;

  return Phi;
}

/**************************************************************/

static Tensor compute_increment_plastic_strain_classical(Tensor s_trial, double s_trial_norm, double delta_Gamma, Material MatProp)
{

  int Ndim = NumberDimensions;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;

  double aux1 = alpha_Q*delta_Gamma;
  double aux2 = delta_Gamma/s_trial_norm;

  Tensor D_E_plastic = alloc__TensorLib__(2);	
  Tensor I = Identity__TensorLib__();

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
	{
	  D_E_plastic.N[i][j] = aux1*I.N[i][j] + aux2*s_trial.N[i][j];
	}
    }

  free__TensorLib__(I);

  return D_E_plastic;
}

/**************************************************************/

static Tensor compute_finite_stress_tensor_classical(Tensor s_trial, Tensor p_trial, double s_trial_norm, double delta_Gamma, Material MatProp)
{
  int Ndim = NumberDimensions;

  Tensor sigma_k1 = alloc__TensorLib__(2);
  Tensor I = Identity__TensorLib__();

  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));

  double aux1 = 3*K*alpha_Q*delta_Gamma;
  double aux2 = 1 - 2*G*delta_Gamma/s_trial_norm;

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
	{
	  sigma_k1.N[i][j] = p_trial.N[i][j] - aux1*I.N[i][j] + aux2*s_trial.N[i][j];
	}
    }

  return sigma_k1;
}

/**************************************************************/

static double compute_derivative_yield_surface_apex(double H, double s_trial_norm, double delta_Gamma2, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/G;

  double aux1 = 3*alpha_Q*K;
  double aux2 = 3*H*beta*DSQR(alpha_Q)*(delta_Gamma1 + delta_Gamma2);
  double aux3 = 3*DSQR(alpha_Q);
  double aux4 = DSQR(delta_Gamma1 + delta_Gamma2);
  double aux5 = 3*alpha_F*sqrt(DSQR(delta_Gamma1) + aux3*aux4);

  double d_Phi = aux1 + aux2/aux5;

  return d_Phi;
}

/**************************************************************/

static double compute_yield_surface_apex(double p_trial_norm, double s_trial_norm, double H, double delta_Gamma2, double c_k, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/G;

  double aux1 = beta/(3*alpha_F);
  double aux2 = 3*DSQR(alpha_Q);
  double aux3 = DSQR(delta_Gamma1 + delta_Gamma2);
  double aux4 = H*sqrt(DSQR(delta_Gamma1) + aux2*aux3); 
  double aux5 = 3*K*alpha_Q*(delta_Gamma1 + delta_Gamma2);

  double Phi = aux1*(c_k + aux4 - p_trial_norm + aux5);

  return Phi;
}

/**************************************************************/

static double update_equivalent_plastic_strain_apex(double EPS_k, double s_trial_norm, double delta_Gamma2, Material MatProp)
{
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/G;

  double aux1 = DSQR(delta_Gamma1);
  double aux2 = 3*DSQR(alpha_Q);
  double aux3 = DSQR(delta_Gamma1 + delta_Gamma2);

  double EPS_k1 = EPS_k + sqrt(aux1 + aux2*aux3);

  return EPS_k1;
}

/**************************************************************/


static Tensor compute_increment_plastic_strain_apex(Tensor s_trial, double s_trial_norm, double delta_Gamma2, Material MatProp)
{

  int Ndim = NumberDimensions;

  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/G;

  double aux1 = alpha_Q*(delta_Gamma1 + delta_Gamma2);
  double aux2 = delta_Gamma1/s_trial_norm;

  Tensor D_E_plastic = alloc__TensorLib__(2);	
  Tensor I = Identity__TensorLib__();

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
	{
	  D_E_plastic.N[i][j] = aux1*I.N[i][j] + aux2*s_trial.N[i][j];
	}
    }

  free__TensorLib__(I);

  return D_E_plastic;
}

/**************************************************************/

static Tensor compute_finite_stress_tensor_apex(Tensor p_trial, double delta_Gamma2, Material MatProp)
{
  int Ndim = NumberDimensions;

  Tensor sigma_k1 = alloc__TensorLib__(2);
  Tensor I = Identity__TensorLib__();

  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));

  double aux1 = 3*K*alpha_Q*delta_Gamma2;

  for(int i = 0 ; i<Ndim ; i++)
    {
	   sigma_k1.N[i][i] = p_trial.N[i][i] - aux1*I.N[i][i];
    }

  return sigma_k1;
}

/**************************************************************/

static void update_plastic_deformation_gradient(Tensor D_E_plastic, Tensor F_plastic)
{
  int Ndim = NumberDimensions;
  double aux1;

  Tensor D_F_plastic;
  Tensor Aux_tensor = alloc__TensorLib__(2);

  D_F_plastic = increment_Deformation_Gradient_exponential_strains__Particles__(D_E_plastic);

    
  /*
    Compute the new value of the plastic deformation gradient 
  */
  for(int i = 0 ; i < Ndim  ; i++)
    {
      for(int j = 0 ; j < Ndim  ; j++)
	{

	  /*
	    Compute row-column multiplication
	  */
	  aux1 = 0;
	  for(int k = 0 ; k < Ndim  ; k++)
	    {
	      aux1 += D_F_plastic.N[i][k]*F_plastic.N[k][j];
	    }

	  /*
	    New value
	  */
	  Aux_tensor.N[i][j] = aux1;
	}
    }

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



