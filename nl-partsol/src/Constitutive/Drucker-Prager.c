#include "nl-partsol.h"


/*
Call global variables
*/
double TOL_Radial_Returning;
int Max_Iterations_Radial_Returning;

/*
  Auxiliar functions 
*/
static void   standard_error(char *);
static bool   check_convergence(double,double,int,int);
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
static Tensor compute_finite_stress_tensor_elastic_region(Tensor, Tensor, Material);

/**************************************************************/

Plastic_status plasticity_Drucker_Prager_Sanavia(Tensor S_p, Tensor C_total, Tensor F_total, Tensor F_plastic, 
                                                         double J, Plastic_status Inputs_VarCons, Material MatProp)
/*	
	Radial returning algorithm
*/
{

  int Ndim = NumberDimensions;

  Tensor C_elastic;
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
  double p_lim;
  Tensor sigma_k1;
  Tensor Increment_E_plastic;

  /* Load material parameters */
  double EPS_k = Inputs_VarCons.EPS;
  double cohesion_k  = Inputs_VarCons.Cohesion;
  double H;

  /* Initialise convergence prameters for the solver */
  double TOL = TOL_Radial_Returning;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;
  bool Convergence = false;

  /* Define output varible */
  Plastic_status Outputs_VarCons;

  /* Compute the trial elastic right Cauchy-Green tensor using the intermediate configuration. */
  C_elastic = alloc__TensorLib__(2); 
  covariant_push_forward_tensor__TensorLib__(C_elastic, C_total, F_plastic);

  /* Calculation of the small strain tensor (the approach of Ortiz and Camacho) */
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
  p_trial_norm = p_trial.N[0][0];

  /*
    Yield condition : Starting from incremental plastic strain equal to zero
  */
  delta_Gamma = 0;
  Phi = compute_yield_surface(p_trial_norm, s_trial_norm, cohesion_k, MatProp);

  if(Phi < TOL)
  {
    sigma_k1 = compute_finite_stress_tensor_elastic_region(s_trial, p_trial, MatProp);
  }
  else
  {

    H = compute_hardening_modulus(EPS_k, MatProp);

    p_lim = compute_limit_between_classic_apex_algorithm(s_trial_norm, cohesion_k, H, MatProp);

      /*
	Classical plastic iterator
      */
    if(p_trial_norm <= p_lim)
    {
     while(Convergence == false)
     {

       d_Phi = compute_derivative_yield_surface_classical(H, MatProp);

       delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

       EPS_k = update_equivalent_plastic_strain_classical(EPS_k, delta_Gamma, MatProp);

       cohesion_k = update_cohesion_modulus(EPS_k, MatProp);

       H = compute_hardening_modulus(EPS_k, MatProp);

       Phi = compute_yield_surface_classical(s_trial_norm, p_trial_norm, delta_Gamma, cohesion_k, MatProp);

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
	    update
	  */
    Increment_E_plastic = compute_increment_plastic_strain_classical(s_trial, s_trial_norm, delta_Gamma, MatProp);

    update_plastic_deformation_gradient__Particles__(Increment_E_plastic,F_plastic);

    free__TensorLib__(Increment_E_plastic);
    
    sigma_k1 = compute_finite_stress_tensor_classical(s_trial, p_trial, s_trial_norm, delta_Gamma, MatProp);


  }
      /*
	Apex plastic iterator
      */
  else
  {

   while(Convergence == false)
   {

    d_Phi = compute_derivative_yield_surface_apex(H, s_trial_norm, delta_Gamma, MatProp);

    delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

    Phi = compute_yield_surface_apex(p_trial_norm, s_trial_norm, H, delta_Gamma, cohesion_k, MatProp);

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
	    update
	  */
  EPS_k = update_equivalent_plastic_strain_apex(EPS_k, s_trial_norm, delta_Gamma, MatProp);

  Increment_E_plastic = compute_increment_plastic_strain_apex(s_trial, s_trial_norm, delta_Gamma,MatProp);

  update_plastic_deformation_gradient__Particles__(Increment_E_plastic,F_plastic);

  free__TensorLib__(Increment_E_plastic);

  sigma_k1 = compute_finite_stress_tensor_apex(p_trial, delta_Gamma, MatProp);

}

}


/*
  Get the stress tensor in the reference configuration (S_p) using the Piola transformation
*/
compute_Piola_transformation(S_p, sigma_k1, F_total, J);



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
    Return cohesion and EPS updated
  */
  Outputs_VarCons.EPS = EPS_k;
  Outputs_VarCons.Cohesion = cohesion_k;
return Outputs_VarCons;
}
/***************************************************************************/

static void standard_error(char * Error_message)
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in plasticity_Drucker_Prager_Sanavia",Error_message);
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
      standard_error(Error_message); 
    }

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
  if(Error_relative < TOL)
    {
      convergence = true;
    }

  return convergence;
}

/**************************************************************/

static double compute_yield_surface(double p_trial_norm, double s_trial_norm, double cohesion_k, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;

  double Phi = 3*alpha_F*p_trial_norm + s_trial_norm - beta*cohesion_k;

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

static double compute_limit_between_classic_apex_algorithm(double s_trial_norm, double cohesion_k, double H, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; 
  double E = MatProp.E;
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));

  double p_lim = 3*alpha_Q*K*s_trial_norm/(2*G) + beta/(3*alpha_F)*(s_trial_norm*H/(2*G)*sqrt(1+3*DSQR(alpha_Q)) + cohesion_k);

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

  double cohesion_k1 = (c_0/(N_exp*E0_p))*pow(basis,exp);

  return cohesion_k1;	
}


/**************************************************************/

static double compute_yield_surface_classical(double s_trial_norm, double p_trial_norm,
                                              double delta_Gamma, double cohesion_k1, Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));

  double Phi = s_trial_norm - 2*G*delta_Gamma + 3*alpha_F*(p_trial_norm - 3*K*alpha_Q*delta_Gamma) - beta*cohesion_k1;

  return Phi;
}

/**************************************************************/

static Tensor compute_increment_plastic_strain_classical(Tensor s_trial, double s_trial_norm, double delta_Gamma, Material MatProp)
{

  int Ndim = NumberDimensions;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;

  double aux1 = alpha_Q*delta_Gamma;
  double aux2 = delta_Gamma/s_trial_norm;

  Tensor Increment_E_plastic = alloc__TensorLib__(2);	
  Tensor I = Identity__TensorLib__();

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
     Increment_E_plastic.N[i][j] = aux1*I.N[i][j] + aux2*s_trial.N[i][j];
   }
 }

 free__TensorLib__(I);

 return Increment_E_plastic;
}

/**************************************************************/

static Tensor compute_finite_stress_tensor_classical(Tensor s_trial, Tensor p_trial, double s_trial_norm,
                                                     double delta_Gamma, Material MatProp)
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

static double compute_derivative_yield_surface_apex(double H, double s_trial_norm, 
                                                    double delta_Gamma2, Material MatProp)
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

static double compute_yield_surface_apex(double p_trial_norm, double s_trial_norm, double H, 
                                         double delta_Gamma2, double cohesion_k, Material MatProp)
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

  double Phi = aux1*(cohesion_k + aux4 - p_trial_norm + aux5);

  return Phi;
}

/**************************************************************/

static double update_equivalent_plastic_strain_apex(double EPS_k, double s_trial_norm,
                                                    double delta_Gamma2, Material MatProp)
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


static Tensor compute_increment_plastic_strain_apex(Tensor s_trial, double s_trial_norm, 
                                                    double delta_Gamma2, Material MatProp)
{

  int Ndim = NumberDimensions;

  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/G;

  double aux1 = alpha_Q*(delta_Gamma1 + delta_Gamma2);
  double aux2 = delta_Gamma1/s_trial_norm;

  Tensor Increment_E_plastic = alloc__TensorLib__(2);	
  Tensor I = Identity__TensorLib__();

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
     Increment_E_plastic.N[i][j] = aux1*I.N[i][j] + aux2*s_trial.N[i][j];
   }
 }

 free__TensorLib__(I);

 return Increment_E_plastic;
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





