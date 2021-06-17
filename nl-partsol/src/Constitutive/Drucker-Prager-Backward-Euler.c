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

static Tensor compute_plastic_flow_direction(Tensor *,double);

static double compute_yield_surface(double,double,double,double,double,Material);
static double compute_limit_between_classic_apex_algorithm(double,double,double,double,double,Material);

static double compute_objective_function(double,double,double,double,double,double,double,double,Material);
static double compute_derivative_objective_function(double,double,double,double,double,double,Material);

static double compute_Yield(double, Material);
static double compute_D_Yield(Material);

static double update_increment_plastic_strain(double, double, double);
static double update_equivalent_plastic_strain(double, double, double, Material);

static void compute_increment_plastic_strain_tensor(Tensor *,Tensor *, double);
static void apply_plastic_corrector_stress_tensor(Tensor *,Tensor *, double,Material);


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
  Tensor plastic_flow_direction;
  double delta_Gamma;
  double s_trial_norm;
  double p_trial_norm;
  double Phi;
  double G;
  double d_G;	
  double p_lim;

  /*
    Load state parameters
  */
  Tensor sigma_k1 = memory_to_tensor__TensorLib__(Inputs_SP.Stress,2);

  /*
    Compute material parameters for the D-P
  */
  double rad_friction_angle  = (PI__MatrixLib__/180)*MatProp.friction_angle;
  double rad_dilatancy_angle = (PI__MatrixLib__/180)*MatProp.dilatancy_angle;
  
  double alpha_F, alpha_Q, beta;

  if(strcmp(MatProp.Type,"Drucker-Prager-Plane-Strain") == 0)
  {
    alpha_F = sqrt(2/3.)*tan(rad_friction_angle)/sqrt(3+4*DSQR(tan(rad_friction_angle)));
    alpha_Q = sqrt(2/3.)*tan(rad_dilatancy_angle)/sqrt(3+4*DSQR(tan(rad_dilatancy_angle)));
    beta = sqrt(2/3.)*3/sqrt(3+4*DSQR(tan(rad_friction_angle)));
  }
  else if(strcmp(MatProp.Type,"Drucker-Prager-Outer-Cone") == 0)
  {
    alpha_F = sqrt(2/3.)*2*sin(rad_friction_angle)/(3-sin(rad_friction_angle));
    alpha_Q = sqrt(2/3.)*2*sin(rad_dilatancy_angle)/(3-sin(rad_dilatancy_angle));
    beta = sqrt(2/3.)*6*cos(rad_friction_angle)/(3-sin(rad_friction_angle));
  }

  /*
    Define output state parameters
  */
  State_Parameters Outputs_VarCons;
  Outputs_VarCons.Increment_E_plastic = (double *)calloc(Ndim*Ndim,sizeof(double));

  Tensor Increment_E_plastic = memory_to_tensor__TensorLib__(Outputs_VarCons.Increment_E_plastic,2);

  /*
    Initialise solver parameters
  */
  double EPS_k = Inputs_SP.EPS;
  double delta_Gamma_k = 0;
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
    Compute plastic flow
  */
  plastic_flow_direction = compute_plastic_flow_direction(&s_trial,s_trial_norm);

  /*
    Check yield condition
  */
  Phi = compute_yield_surface(p_trial_norm, s_trial_norm, EPS_k, alpha_F, beta, MatProp);

  if(Phi > TOL)
  {

    /*
      Compute parameter to determine apex region
    */
    p_lim = compute_limit_between_classic_apex_algorithm(s_trial_norm, EPS_k, alpha_F, alpha_Q, beta, MatProp);

    /*
      Newton-Rapson solver
    */
    while(Convergence == false)
    {
      Iter++;

      G = compute_objective_function(s_trial_norm, p_trial_norm, p_lim, delta_Gamma, EPS_k, alpha_F, alpha_Q, beta, MatProp);

      Convergence = check_convergence(Phi,TOL,Iter,MaxIter);

      d_G = compute_derivative_objective_function(p_trial_norm,p_lim,EPS_k,alpha_F,alpha_Q,beta,MatProp);

      delta_Gamma = update_increment_plastic_strain(delta_Gamma, G, d_G);

      EPS_k = update_equivalent_plastic_strain(alpha_Q, EPS_k, delta_Gamma, MatProp);

    }

    /*
      Update plastic deformation gradient
    */
    compute_increment_plastic_strain_tensor(&Increment_E_plastic, &plastic_flow_direction, delta_Gamma_k);

    /*
      Update stress tensor with the plastic corrector
    */
    apply_plastic_corrector_stress_tensor(&sigma_k1,&plastic_flow_direction,delta_Gamma_k,MatProp);

  }
  else
  {
    Outputs_VarCons.EPS = Inputs_SP.EPS;
  }

  /*
    Free memory
  */
  free__TensorLib__(s_trial);
  free__TensorLib__(p_trial);
  free__TensorLib__(plastic_flow_direction);


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

static Tensor compute_plastic_flow_direction(
  Tensor * deviatoric_stress,
  double deviatoric_stress_norm)
{
  int Ndim = NumberDimensions;
  Tensor plastic_flow_direction = alloc__TensorLib__(2); 

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      plastic_flow_direction.N[i][j] = deviatoric_stress->N[i][j]/deviatoric_stress_norm;
    }
  }

  return plastic_flow_direction;
}

/**************************************************************/

static double compute_yield_surface(
  double p_trial_norm,
  double s_trial_norm,
  double EPS_k,
  double alpha_F,
  double beta,
  Material MatProp)
{

  double Yield = compute_Yield(EPS_k,MatProp);

  double Phi = 3*alpha_F*p_trial_norm + s_trial_norm - beta*Yield;

  return Phi;
}

/**************************************************************/

static double compute_limit_between_classic_apex_algorithm(

  double s_trial_norm,
  double EPS_k,
  double alpha_F,
  double alpha_Q,
  double beta,
  Material MatProp)
{
  double nu = MatProp.nu; 
  double E = MatProp.E;
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double Hardening_modulus = MatProp.Hardening_modulus;
  double D_EPS = sqrt(1+3*DSQR(alpha_Q));

  double Yield = compute_Yield(EPS_k,MatProp);

  double constant_1 = 3*alpha_Q*K/(2*G);
  double constant_2 = beta*sqrt(2./3.)/(3*alpha_F);
  double e_trial_norm = s_trial_norm/(2*G);

  return constant_1*s_trial_norm + constant_2*(e_trial_norm*Hardening_modulus*D_EPS + Yield);
}

/**************************************************************/

static double compute_Yield(
  double EPS_k,
  Material MatProp)
{
  double Sigma_y = MatProp.yield_stress_0;
  double Hardening_modulus = MatProp.Hardening_modulus;

  return Sigma_y + Hardening_modulus*EPS_k;
}

/**************************************************************/

static double compute_D_Yield(
  Material MatProp)
{
  double Hardening_modulus = MatProp.Hardening_modulus;

  return Hardening_modulus;
}

/**************************************************************/

static double compute_objective_function(
  double s_trial_norm,
  double p_trial_norm,
  double p_lim,
  double delta_Gamma,
  double EPS_k,
  double alpha_F,
  double alpha_Q,
  double beta,
  Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double Yield = compute_Yield(EPS_k,MatProp);

  if(p_trial_norm <= p_lim)
  {
    double p_k = p_trial_norm - 3*K*alpha_Q*delta_Gamma;
    
    return s_trial_norm - 2*G*delta_Gamma + 3*alpha_F*p_k - beta*Yield;
  }
  else
  {
    double e_trial_norm = s_trial_norm/(2*G);
    double p_k = p_trial_norm - 3*K*alpha_Q*(e_trial_norm + delta_Gamma);

    return 3*alpha_F*p_k - beta*Yield;
  }

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

static double compute_derivative_objective_function(
  double p_trial_norm,
  double p_lim,
  double EPS_k,
  double alpha_F,
  double alpha_Q,
  double beta,
  Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  
  if(p_trial_norm <= p_lim)
  {
    double D_Hardening = compute_D_Yield(MatProp);
    double D_EPS = sqrt(3*DSQR(alpha_Q) + 1);
    double d_p_k = - 3*K*alpha_Q;
    
    return - 2*G + 3*alpha_F*d_p_k - beta*D_Hardening*D_EPS;
  }
  else
  {
    return 0.0;
  }
  
}

/**************************************************************/

static double update_equivalent_plastic_strain(
  double alpha_Q,
  double EPS_k,
  double delta_Gamma, 
  Material MatProp)
{
  double EPS_k1 = EPS_k + delta_Gamma*sqrt(3*DSQR(alpha_Q) + 1);

  return EPS_k1;
}

/**************************************************************/

static void compute_increment_plastic_strain_tensor(
  Tensor * D_E_plastic,
  Tensor * plastic_flow_direction, 
  double delta_Gamma)
{
  int Ndim = NumberDimensions;
 
  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      D_E_plastic->N[i][j] = delta_Gamma*plastic_flow_direction->N[i][j];
    }
  }

}

/**************************************************************/

static void apply_plastic_corrector_stress_tensor(
  Tensor * sigma_k1, 
  Tensor * plastic_flow_direction, 
  double delta_Gamma, 
  Material MatProp)
{
  int Ndim = NumberDimensions;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      sigma_k1->N[i][j] -= 2*G*delta_Gamma*plastic_flow_direction->N[i][j];
    }
  }
}

/**************************************************************/




