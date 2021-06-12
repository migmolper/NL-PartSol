#include "nl-partsol.h"

/*
  Call global variables
*/
double TOL_Radial_Returning;
int Max_Iterations_Radial_Returning;
double DeltaTimeStep;


/*
  Define local global variable for the relative error
*/
double Error0;

/*
  Auxiliar functions 
*/
static void   standard_error(char *);
static bool   check_convergence(double,double,int,int);

static Tensor compute_plastic_flow_direction(Tensor *,double);

static double compute_yield_surface(double,double,Material);

static double compute_objective_function(double,double,double,Material);
static double compute_derivative_objective_function(double,Material);

static double compute_K(double,Material);
static double compute_D_K(double,Material);

static double compute_DeltaH(double,Material);
static double compute_D_DeltaH(Material);

static double update_increment_plastic_strain(double, double, double);
static double update_equivalent_plastic_strain(double, double);
static  void  update_back_stress(Tensor *,Tensor *,double, Material);

static void compute_increment_plastic_strain_tensor(Tensor *,Tensor *,double);
static void   apply_plastic_corrector_stress_tensor(Tensor *,Tensor *,double,Material);

/**************************************************************/

State_Parameters Von_Mises_backward_euler(
  State_Parameters Inputs_SP,
  Material MatProp)
/*	
	Radial returning algorithm for the Von-Mises plastic criterium
*/
{

  int Ndim = NumberDimensions;

  Tensor p_trial;
  Tensor s_trial;
  Tensor relative_stress;
  Tensor plastic_flow_direction;
  double s_trial_norm;
  double relative_stress_norm;
  double Phi;
  double G;
  double d_G;
  
  /*
    Load state parameters
  */
  Tensor sigma_k1 = memory_to_tensor__TensorLib__(Inputs_SP.Stress,2);
  Tensor Back_stress = memory_to_tensor__TensorLib__(Inputs_SP.Back_stress,2);

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

  /*
    Compute the relative stress
  */
  relative_stress = subtraction__TensorLib__(s_trial,Back_stress);
  relative_stress_norm = EuclideanNorm__TensorLib__(relative_stress);

  /*
    Compute plastic flow
  */
  plastic_flow_direction = compute_plastic_flow_direction(&relative_stress,relative_stress_norm);

  /*
    Yield condition : Starting from incremental plastic strain equal to zero
  */
  Phi = compute_yield_surface(relative_stress_norm, Inputs_SP.EPS, MatProp);

  if(Phi > TOL)
  {     
    /*
      Newton-Rapson solver
    */
    while(Convergence == false)
    {
          
      G = compute_objective_function(relative_stress_norm, delta_Gamma_k, EPS_k, MatProp);

      Convergence = check_convergence(G,TOL,Iter,MaxIter);
    
      d_G = compute_derivative_objective_function(EPS_k,MatProp);

      delta_Gamma_k = update_increment_plastic_strain(delta_Gamma_k, G, d_G);

      EPS_k = update_equivalent_plastic_strain(Inputs_SP.EPS, delta_Gamma_k);

      Iter++;

    }

    /*
      Update plastic deformation gradient
    */
    compute_increment_plastic_strain_tensor(&Increment_E_plastic,&plastic_flow_direction, delta_Gamma_k);

    /*
      Update back stress
    */
    update_back_stress(&Back_stress,&plastic_flow_direction,delta_Gamma_k,MatProp);

    /*
      Update stress tensor in the deformed configuration
    */
    apply_plastic_corrector_stress_tensor(&sigma_k1,&plastic_flow_direction,delta_Gamma_k,MatProp);

    /*
      Update equivalent plastic strain
    */
    Outputs_VarCons.EPS = EPS_k;

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
  free__TensorLib__(relative_stress);
  free__TensorLib__(plastic_flow_direction);

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
  Tensor * relative_stress,
  double relative_stress_norm)
{
  int Ndim = NumberDimensions;
  Tensor plastic_flow_direction = alloc__TensorLib__(2); 

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      plastic_flow_direction.N[i][j] = relative_stress->N[i][j]/relative_stress_norm;
    }
  }

  return plastic_flow_direction;
}

/**************************************************************/

static double compute_yield_surface(
  double relative_stress_norm,
  double EPS,
  Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double K = compute_K(EPS,MatProp); /* Isotropic hardening */

  return - sqrt(2./3.)*K + relative_stress_norm;
}

/**************************************************************/

static double compute_objective_function(
  double relative_stress_norm,
  double delta_Gamma,
  double EPS_k,
  Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double K = compute_K(EPS_k,MatProp); /* Isotropic hardening */
  double DeltaH = compute_DeltaH(delta_Gamma,MatProp); /* Kinematic hardening */

  if(MatProp.Viscous_regularization)
  {
    double fluidity_param = MatProp.fluidity_param;
    double DT = 1.0;

    return - sqrt(2./3.)*(K + DeltaH) + relative_stress_norm - DMAX(0,delta_Gamma*fluidity_param/DT) - 2*G*delta_Gamma;
  }
  else
  {
      return - sqrt(2./3.)*(K + DeltaH) + relative_stress_norm - 2*G*delta_Gamma;
  }
}

/**************************************************************/


static double compute_derivative_objective_function(
  double EPS_k,
  Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double D_K = compute_D_K(EPS_k,MatProp); /* Derivative of the isotropic hardening */
  double D_DeltaH = compute_D_DeltaH(MatProp); /* Derivative of the kinematic hardening */

  if(MatProp.Viscous_regularization)
  {
    double fluidity_param = MatProp.fluidity_param;
    double DT = 1.0;

    return - sqrt(2./3.)*(D_K + D_DeltaH) - 2*G - fluidity_param/DT;
  }
  else
  {
    return - sqrt(2./3.)*(D_K + D_DeltaH) - 2*G;
  }

}

/**************************************************************/

static double compute_K(
  double EPS,
  Material MatProp)
/*
  Isotropic hardening function
*/
{
  double Sigma_y = MatProp.yield_stress_0;

  if(MatProp.Hardening_Hughes)
  {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.Parameter_Hardening_Hughes;

    return Sigma_y + theta*Hardening_modulus*EPS;
  }
  else if(MatProp.Hardening_Cervera)
  {
    double Hardening_modulus = MatProp.Hardening_modulus;

    return Sigma_y*exp(-2*Hardening_modulus*EPS/Sigma_y);
  }
  else if(MatProp.Hardening_Ortiz)
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
/*
  Derivative of the isotropic hardening function
*/
{
  
  if(MatProp.Hardening_Hughes)
  {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.Parameter_Hardening_Hughes;

    return sqrt(2./3.)*theta*Hardening_modulus;
  }
  else if(MatProp.Hardening_Cervera)
  {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double Sigma_y = MatProp.yield_stress_0;

    return -2.0*sqrt(2./3.)*Hardening_modulus*exp(-2.0*Hardening_modulus*EPS/Sigma_y);
  }
  else if(MatProp.Hardening_Ortiz)
  {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double Exponent = MatProp.Exponent_Hardening_Ortiz;
    double Ref_PS = MatProp.Reference_Plastic_Strain_Ortiz;
    double Sigma_y = MatProp.yield_stress_0;

    return sqrt(2./3.)*(Sigma_y/(Ref_PS*Exponent))*pow(1 + EPS/Ref_PS,1.0/Exponent - 1);
  }
  else
  {
    return 0.0;
  }

}

/**************************************************************/

static double compute_DeltaH(
  double delta_Gamma,
  Material MatProp)
/*
  Kinematic hardening function
*/
{

  if(MatProp.Hardening_Hughes)
  {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.Parameter_Hardening_Hughes;

    return sqrt(2./3.)*(1-theta)*Hardening_modulus*delta_Gamma;
  }
  else
  {
    return 0.0;
  }

}

/**************************************************************/

static double compute_D_DeltaH(
  Material MatProp)
/*
  Derivative of the kinematic hardening function
*/
{

  if(MatProp.Hardening_Hughes)
  {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.Parameter_Hardening_Hughes;

    return sqrt(2./3.)*(1-theta)*Hardening_modulus;
  }
  else
  {
    return 0.0;
  }

}

/**************************************************************/

static double update_increment_plastic_strain(
  double delta_Gamma_k,
  double G,
  double d_G)
{
  return delta_Gamma_k - G/d_G;
}

/**************************************************************/

static double update_equivalent_plastic_strain(
  double EPS,
  double delta_Gamma)
{
  return EPS + sqrt(2./3.)*delta_Gamma;
}

/**************************************************************/

static void update_back_stress(
  Tensor * Back_stress,
  Tensor * plastic_flow_direction,
  double delta_Gamma,
  Material MatProp)
{
  int Ndim = NumberDimensions;
  double DeltaH = compute_DeltaH(delta_Gamma,MatProp);

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      Back_stress->N[i][j] += sqrt(2./3.)*DeltaH*plastic_flow_direction->N[i][j];
    }
  }

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
