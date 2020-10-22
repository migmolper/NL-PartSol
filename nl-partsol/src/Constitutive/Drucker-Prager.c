#include "nl-partsol.h"

/*
  Auxiliar functions 
*/
static Tensor compute_small_strain_tensor(Tensor, Tensor);
static void   compute_volumetric_deviatoric_stress_tensor(double *, Tensor, Tensor, Material);
static double compute_yield_surface(double, double, double, Material);
static double compute_hardening_modulus(double, Material);
static double compute_limit_between_classic_apex_algorithm(double, double, double, Material);
static double compute_derivative_yield_surface_classical(double, Material);
static double compute_increment_plastic_strain(double, double, double);
static double update_equivalent_plastic_strain_classical(double, double, Material);
static double update_cohesion_modulus(double, Material);
static double compute_yield_surface_classical(double, double, double, double, Material);
static Tensor compute_increment_plastic_strain_classical(Tensor, double, double, Material);
static Tensor compute_finite_stress_tensor_classical(Tensor, double, double, double, Material);
static double compute_derivative_yield_surface_apex(double, double, double, Material);
static double compute_yield_surface_apex(double, double, double, double, Material);
static double update_equivalent_plastic_strain_apex(double, double, double, Material);
static Tensor compute_increment_plastic_strain_apex(Tensor, double, double, Material);
static Tensor compute_finite_stress_tensor_apex(double, double, Material);
extern void   update_plastic_deformation_gradient(Tensor, Tensor);

/**************************************************************/

Tensor plasticity_Drucker_Prager_Sanavia(Tensor grad_e, Tensor C, Tensor F_plastic, Tensor F, 
                            					   double * ptr_EPS_k, double * ptr_c_k, double J, Material MatProp)
/*	
	Radial returning algorithm
*/
{

  int Ndim = NumberDimensions;

  Tensor E_elastic;
  double p_trial;
  Tensor s_trial;
  double delta_Gamma;
  double s_trial_norm;
  double Phi;
  double d_Phi;	
  double H;
  double EPS_k;
  double c_k;
  double p_lim;
  Tensor sigma_k1;
  Tensor D_E_plastic;

  double TOL_DP = MatProp.TOL_Drucker_Prager;
  int iter_max_DP = MatProp.Max_Iterations_Drucker_Prager;

  /*	
	calculation of the small strain tensor
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

  Phi = compute_yield_surface(p_trial, s_trial_norm, c_k, MatProp);

  if(Phi >= 0)
    {
      if(strcmp(MatProp_p.Type,"Von-Mises") == 0)
        {
          H = MatProp.hardening_modulus;
        }
        else
        {
          H = compute_hardening_modulus(*ptr_EPS_k, MatProp);
        }

      p_lim = compute_limit_between_classic_apex_algorithm(s_trial_norm, *ptr_c_k, H, MatProp);

      /*
	Classical plastic iterator
      */
      if(p_trial <= p_lim)
	{
	  while(iter < iter_max_DP)
	    {

	      d_Phi = compute_derivative_yield_surface_classical(H, MatProp);

	      delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

	      ptr_EPS_k = &(update_equivalent_plastic_strain_classical(*ptr_EPS_k, delta_Gamma, MatProp));

        if(strcmp(MatProp_p.Type,"Von-Mises") == 0)
        {
          H = MatProp.hardening_modulus;

          c_k = MatProp.yield_stress;
        }
        else
        {
          ptr_c_k = &(update_cohesion_modulus(*ptr_EPS_k, MatProp));

          H = compute_hardening_modulus(*ptr_EPS_k, MatProp);
        }

	      Phi = compute_yield_surface_classical(s_trial_norm, p_trial, delta_Gamma, *ptr_c_k, MatProp);

	      /*
		Check convergence
	      */
	      if(Phi < TOL_DP)
		{
		  iter++;
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

	  while(iter < iter_max_DP)
	    {

	      d_Phi = compute_derivative_yield_surface_apex(H, s_trial_norm, delta_Gamma, MatProp);

	      delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

	      Phi = compute_yield_surface_apex(p_trial, H, delta_Gamma, *ptr_c_k, MatProp);

	      /*
		Check convergence
	      */
	      if(Phi < TOL_DP)
		{
		  iter++;
		}
	      else
		{
		  break;
		}

	    }

	  /*
	    update
	  */
	  ptr_EPS_k = &(update_equivalent_plastic_strain_apex(*ptr_EPS_k, s_trial_norm, delta_Gamma, MatProp));

	  D_E_plastic = compute_increment_plastic_strain_apex(s_trial, s_trial_norm, delta_Gamma,MatProp);

	  sigma_k1 = compute_finite_stress_tensor_apex(p_trial, delta_Gamma2, MatProp);

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
  double nu = Mat.nu; 
  double E = Mat.E;
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

static double compute_yield_surface(double p_trial, double s_trial_norm, double c_k, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;

  double Phi = 3*alpha_F*p_trial + s_trial_norm - beta*c_k;

  return Phi;
}

/**************************************************************/

static double compute_hardening_modulus(double EPS_k1, Material MatProp)
{
  double E0_p  = MatProp.E_plastic_reference_Drucker_Prager;
  double N_exp = MatProp.hardening_exp_Drucker_Prager;
  double c0    = MatProp.cohesion_reference;

  double H = (c0/(N_exp*E0))*pow((1 + EPS_k/E0_p),(1/N_exp - 1));

  return H;
}

/**************************************************************/

static double compute_limit_between_classic_apex_algorithm(double s_trial_norm, double c_k, double H, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = Mat.nu; 
  double E = Mat.E;
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
  double nu = Mat.nu; /* Poisson modulus */
  double E = Mat.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));

  double d_Phi = - 9*K*alpha_F*alpha_Q - 2*G - H*beta*sqrt(3*DSQR(alpha_Q) + 1);

  return d_Phi;
}

/**************************************************************/

static double compute_increment_plastic_strain(double delta_Gamma_k, double Phi, double d_Phi)
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
  double c0    = MatProp.cohesion_reference;
  double E0_p  = MatProp.E_plastic_reference_Drucker_Prager;
  double N_exp = MatProp.hardening_exp_Drucker_Prager;
  double exp   = 1/N_exp - 1;
  double basis = 1 + EPS_k1/E0_p;

  double c_k1 = (c_0/(N_exp*E0_p))*pow(basis,exp);

  return c_k1;	
}


/**************************************************************/

static double compute_yield_surface_classical(double s_trial_norm, double p_trial, double delta_Gamma, double c_k1, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = Mat.nu; /* Poisson modulus */
  double E = Mat.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));

  double Phi = s_trial_norm - 2*G*delta_Gamma + 3*alpha_F*(p_trial - 3*K*alpha_Q*delta_Gamma) - beta*c_k1;

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

static Tensor compute_finite_stress_tensor_classical(Tensor s_trial, double p_trial, double s_trial_norm, double delta_Gamma, Material MatProp)
{
  int Ndim = NumberDimensions;

  Tensor sigma_k1 = alloc__TensorLib__(2);
  Tensor I = Identity__TensorLib__();

  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = Mat.nu; /* Poisson modulus */
  double E = Mat.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));

  double aux1 = p_trial - 3*K*alpha_Q*delta_Gamma;
  double aux2 = 1 - 2*G*delta_Gamma/s_trial_norm;

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
	{
	  sigma_k1.N[i][j] = aux1*I.N[i][j] + aux2*s_trial.N[i][j];
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
  double nu = Mat.nu; /* Poisson modulus */
  double E = Mat.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/G

    double aux1 = 3*alpha_Q*K;
  double aux2 = 3*H*beta*DSQR(alpha_Q)*(delta_Gamma1 + delta_Gamma2);
  double aux3 = 3*alpha_F*sqrt(DSQR(delta_Gamma1) + 3*DSQR(alpha_Q)*DSQR(delta_Gamma1 + delta_Gamma2));

  double d_Phi = aux1 + aux2/aux3;

  return d_Phi;
}

/**************************************************************/

static double compute_yield_surface_apex(double p_trial, double H, double delta_Gamma2, double c_k, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = Mat.nu; /* Poisson modulus */
  double E = Mat.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/G;

  double aux1 = beta/(3*alpha_F);
  double aux2 = H*sqrt(DSQR(delta_Gamma1) + 3*DSQR(alpha_Q)*DSQR(delta_Gamma1 + delta_Gamma2));
  double aux3 = 3*K*alpha_Q*(delta_Gamma1 + delta_Gamma2);

  double Phi = aux1*(c_k + aux2 - p_trial + aux3);

  return Phi;
}

/**************************************************************/

static double update_equivalent_plastic_strain_apex(double EPS_k, double s_trial_norm, double delta_Gamma2, Material MatProp)
{
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double nu = Mat.nu; /* Poisson modulus */
  double E = Mat.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/G;

  double aux1 = DSQR(delta_Gamma1) + 3*DSQR(alpha_Q)*DSQR(delta_Gamma1 + delta_Gamma2);

  double EPS_k1 = EPS_k + sqrt(aux1);

  return EPS_k1;
}

/**************************************************************/


static Tensor compute_increment_plastic_strain_apex(Tensor s_trial, double s_trial_norm, double delta_Gamma2, Material MatProp)
{

  int Ndim = NumberDimensions;

  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double nu = Mat.nu; /* Poisson modulus */
  double E = Mat.E; /* Elastic modulus */
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

static Tensor compute_finite_stress_tensor_apex(double p_trial, double delta_Gamma2, Material MatProp)
{
  int Ndim = NumberDimensions;

  Tensor sigma_k1 = alloc__TensorLib__(2);
  Tensor I = Identity__TensorLib__();

  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double nu = Mat.nu; /* Poisson modulus */
  double E = Mat.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));

  double aux1 = p_trial - 3*K*alpha_Q*delta_Gamma2;

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
	{
	  sigma_k1.N[i][j] = aux1*I.N[i][j];
	}
    }

  return sigma_k1;
}

/**************************************************************/

extern void update_plastic_deformation_gradient(Tensor D_E_plastic, Tensor F_plastic)
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



