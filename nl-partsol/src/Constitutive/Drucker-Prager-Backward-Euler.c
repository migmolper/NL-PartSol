#include "nl-partsol.h"


/*
Call global variables
*/
double TOL_Radial_Returning;
int Max_Iterations_Radial_Returning;

double Error0;

/*
  Auxiliar functions 
*/
static void   standard_error(char *);
static bool   check_convergence(double,double,int,int);

static double compute_yield_surface(double, double, double, Material);
static double compute_D_K(double, Material);
static double compute_limit_between_classic_apex_algorithm(double, double, double, Material);
static double compute_derivative_yield_surface_classical(double, Material);
static double update_increment_plastic_strain(double, double, double);
static double update_equivalent_plastic_strain_classical(double, double, Material);
static double compute_K(double, double, Material);
static double compute_yield_surface_classical(double, double, double, double, Material);
static Tensor compute_increment_plastic_strain_classical(Tensor, double, double, Material);
static void   compute_finite_stress_tensor_classical(Tensor, Tensor, Tensor, double, double, Material);
static double compute_derivative_yield_surface_apex(double, double, double, Material);
static double compute_yield_surface_apex(double, double, double, double, double, Material);
static double update_equivalent_plastic_strain_apex(double, double, double, Material);
static Tensor compute_increment_plastic_strain_apex(Tensor, double, double, Material);
static void   compute_finite_stress_tensor_apex(Tensor, Tensor, double, double, Material);
static void   compute_finite_stress_tensor_elastic_region(Tensor, Tensor, Tensor, Material);


/**************************************************************/

State_Parameters Drucker_Prager_backward_euler(
  State_Parameters Inputs_SP,
  Material MatProp)
/*	
	Radial returning algorithm for the Drucker-Prager de Sanavia algorithm
*/
{

  int Ndim = NumberDimensions;

  Tensor p_trial;
  Tensor s_trial;
  double delta_Gamma;
  double s_trial_norm;
  double p_trial_norm;
  double Phi;
  double d_Phi;	
  double p_lim;
  Tensor Increment_E_plastic = alloc__TensorLib__(2);

  /*
    Load state parameters
  */
  Tensor sigma_k1 = memory_to_tensor__TensorLib__(Inputs_SP.Stress,2);
  double EPS_k = Inputs_SP.EPS;

  /* Initialise convergence prameters for the solver */
  double TOL = TOL_Radial_Returning;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;
  bool Convergence = false;

  /*
    Decompose the elastic trial in volumetric and deviatoric components
  */
  p_trial = volumetric_component__TensorLib__(sigma_k1);
  s_trial = deviatoric_component__TensorLib__(sigma_k1,p_trial);

  s_trial_norm = EuclideanNorm__TensorLib__(s_trial);
  p_trial_norm = EuclideanNorm__TensorLib__(p_trial); 

  /*
    Yield condition : Starting from incremental plastic strain equal to zero
  */
  delta_Gamma = 0.0;
  Phi = compute_yield_surface(p_trial_norm, s_trial_norm, EPS, MatProp);

  if(Phi < TOL)
  {
    compute_finite_stress_tensor_elastic_region(sigma_k1, s_trial, p_trial, MatProp);
  }
  else
  {
    p_lim = compute_limit_between_classic_apex_algorithm(s_trial_norm, cohesion_k, EPS_k, MatProp);

    /*
	    Classical plastic iterator
    */
    if(p_trial_norm <= p_lim)
    {

     while(Convergence == false)
     {

      Phi = compute_yield_surface_classical(s_trial_norm, p_trial_norm, delta_Gamma, cohesion_k, MatProp);
      
      d_Phi = compute_derivative_yield_surface_classical(hardening_k, MatProp);

      delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

      EPS_k = update_equivalent_plastic_strain_classical(EPS_k, delta_Gamma, MatProp);

      hardening_k = compute_D_K(EPS_k, MatProp);

      Convergence = check_convergence(Phi,TOL,Iter,MaxIter);

      if(Convergence == false)
      {
        Iter++;
      }

    }


    Increment_E_plastic = compute_increment_plastic_strain_classical(s_trial, s_trial_norm, delta_Gamma, MatProp);
    
    compute_finite_stress_tensor_classical(sigma_k1, s_trial, p_trial, s_trial_norm, delta_Gamma, MatProp);


  }
  /*
    Apex plastic iterator
  */
  else
  {

   while(Convergence == false)
   {

    Phi = compute_yield_surface_apex(p_trial_norm, s_trial_norm, hardening_k, delta_Gamma, cohesion_k, MatProp);

    d_Phi = compute_derivative_yield_surface_apex(hardening_k, s_trial_norm, delta_Gamma, MatProp);

    delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

    Convergence = check_convergence(Phi,TOL,Iter,MaxIter);

    if(Convergence == false)
    {
      Iter++;
    }

  }

  EPS_k = update_equivalent_plastic_strain_apex(EPS_k, s_trial_norm, delta_Gamma, MatProp);

  Increment_E_plastic = compute_increment_plastic_strain_apex(s_trial, s_trial_norm, delta_Gamma,MatProp);

  compute_finite_stress_tensor_apex(sigma_k1, p_trial, s_trial_norm, delta_Gamma, MatProp);

}

}

  /*
    Free memory
  */
free__TensorLib__(s_trial);
free__TensorLib__(p_trial);

  /*
    Define output varible
  */
  State_Parameters Outputs_VarCons;
  Outputs_VarCons.EPS = EPS_k;
  Outputs_VarCons.Increment_E_plastic = Increment_E_plastic.N[0];

return Outputs_VarCons;

}
/***************************************************************************/

static void standard_error(
  char * Error_message)
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in plasticity_Drucker_Prager_Sanavia",Error_message);
    exit(EXIT_FAILURE);
}

/**************************************************************/

static bool check_convergence(
  double Error,
  double TOL,
  int Iter, 
  int MaxIter)
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
      Error0 = fabs(Error);
      Error_relative = fabs(Error/Error0);
    }
    else
    {
      Error_relative = fabs(Error/Error0);
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
  double p_trial_norm,
  double s_trial_norm,
  double EPS,
  Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double K = compute_K(EPS, MatProp); // Isotropic hardening

  double Phi = 3*alpha_F*p_trial_norm + s_trial_norm - beta*K;

  return Phi;
}


/**************************************************************/

static double compute_limit_between_classic_apex_algorithm(
  double s_trial_norm,
  double cohesion_k,
  double EPS_k,
  Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; 
  double E = MatProp.E;
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double H = compute_D_K(EPS_k, MatProp);

  double p_lim = 3*alpha_Q*K*s_trial_norm/(2*G) + (beta/(3*alpha_F))*((s_trial_norm/(2*G))*H*sqrt(1+3*DSQR(alpha_Q)) + cohesion_k);

  return p_lim;

}

/**************************************************************/

static double compute_derivative_yield_surface_classical(
  double hardening_k,
  Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));

  double d_Phi = - 2*G - 9*K*alpha_F*alpha_Q - beta*hardening_k*sqrt(3*DSQR(alpha_Q) + 1);

  return d_Phi;
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

static double update_equivalent_plastic_strain_classical(
  double EPS_k,
  double delta_Gamma, 
  Material MatProp)
{
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;

  double EPS_k1 = EPS_k + delta_Gamma*sqrt(3*DSQR(alpha_Q) + 1);

  return EPS_k1;
}


/**************************************************************/

static double compute_yield_surface_classical(
  double s_trial_norm,
  double p_trial_norm,
  double delta_Gamma,
  double cohesion_k1,
  Material MatProp)
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

static Tensor compute_increment_plastic_strain_classical(
  Tensor s_trial,
  double s_trial_norm,
  double delta_Gamma, 
  Material MatProp)
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

static double compute_derivative_yield_surface_apex(
  double hardening_k,
  double s_trial_norm,
  double delta_Gamma2,
  Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/(2*G);

  double aux1 = 3*alpha_Q*K;
  double aux2 = 3*beta*DSQR(alpha_Q)*(delta_Gamma1 + delta_Gamma2);
  double aux3 = 3*DSQR(alpha_Q);
  double aux4 = DSQR(delta_Gamma1 + delta_Gamma2);
  double aux5 = 3*alpha_F*sqrt(DSQR(delta_Gamma1) + aux3*aux4);

  double d_Phi = aux1 + hardening_k*aux2/aux5;

  return d_Phi;
}

/**************************************************************/

static double compute_yield_surface_apex(
  double p_trial_norm,
  double s_trial_norm,
  double hardening_k, 
  double delta_Gamma2,
  double cohesion_k,
  Material MatProp)
{
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta = MatProp.beta_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/(2*G);

  double aux1 = beta/(3*alpha_F);
  double aux2 = 3*DSQR(alpha_Q);
  double aux3 = DSQR(delta_Gamma1 + delta_Gamma2);
  double aux4 = sqrt(DSQR(delta_Gamma1) + aux2*aux3); 
  double aux5 = 3*K*alpha_Q*(delta_Gamma1 + delta_Gamma2);

  double Phi = aux1*(cohesion_k + hardening_k*aux4) - p_trial_norm + aux5;

  return Phi;
}

/**************************************************************/

static double update_equivalent_plastic_strain_apex(
  double EPS_k,
  double s_trial_norm,
  double delta_Gamma2,
  Material MatProp)
{
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/(2*G);

  double aux1 = DSQR(delta_Gamma1);
  double aux2 = 3*DSQR(alpha_Q);
  double aux3 = DSQR(delta_Gamma1 + delta_Gamma2);

  double EPS_k1 = EPS_k + sqrt(aux1 + aux2*aux3);

  return EPS_k1;
}


/**************************************************************/

static double compute_K(
  double EPS,
  Material MatProp)
{
  double Sigma_y = MatProp.yield_stress_0;

  if(MatProp.Hardening_Ortiz)
  {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double Exponent = MatProp.Exponent_Hardening_Ortiz;
    double Ref_PS = MatProp.Reference_Plastic_Strain_Ortiz;

    return Sigma_y*pow(1 + EPS/Ref_PS,1.0/Exponent);
  }
  else
  {
    return Sigma_y;
  }
}

/**************************************************************/

static double compute_D_K(
  double EPS,
  Material MatProp)
{
  if(MatProp.Hardening_Ortiz)
  {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double Exponent = MatProp.Exponent_Hardening_Ortiz;
    double Ref_PS = MatProp.Reference_Plastic_Strain_Ortiz;
    double Sigma_y = MatProp.yield_stress_0;

    return (Sigma_y/(Ref_PS*Exponent))*pow(1 + EPS/Ref_PS,1.0/Exponent - 1);
  }
  else
  {
    return 0.0;
  }
}


/**************************************************************/


static Tensor compute_increment_plastic_strain_apex(
  Tensor s_trial,
  double s_trial_norm, 
  double delta_Gamma2,
  Material MatProp)
{

  int Ndim = NumberDimensions;

  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/(2*G);

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

static void compute_finite_stress_tensor_classical(
  Tensor sigma_k1,
  Tensor s_trial,
  Tensor p_trial,
  double s_trial_norm,
  double delta_Gamma, Material MatProp)
{
  int Ndim = NumberDimensions;

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

  free__TensorLib__(I);
}

/**************************************************************/

static void compute_finite_stress_tensor_apex(
  Tensor sigma_k1,
  Tensor p_trial,
  double s_trial_norm,
  double delta_Gamma2,
  Material MatProp)
{
  int Ndim = NumberDimensions;

  Tensor I = Identity__TensorLib__();

  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double delta_Gamma1 = s_trial_norm/(2*G);

  double aux1 = 3*K*alpha_Q*(delta_Gamma1 + delta_Gamma2);

  for(int i = 0 ; i<Ndim ; i++)
  {
    sigma_k1.N[i][i] = p_trial.N[i][i] - aux1*I.N[i][i];
  }

  free__TensorLib__(I);
}

/**************************************************************/

static void compute_finite_stress_tensor_elastic_region(
  Tensor sigma_k1,
  Tensor s_trial,
  Tensor p_trial,
  Material MatProp)
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





