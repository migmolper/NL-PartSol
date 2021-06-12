#include "nl-partsol.h"

/*
  Call global variables
*/
double DeltaTimeStep;

/*
  Auxiliar functions 
*/
static Tensor compute_plastic_flow_direction(Tensor *,double);
static double compute_yield_surface(double,double,Material);
static double compute_DeltaH(double,Material);
static double compute_increment_flow_rule(double, Material);
static double update_equivalent_plastic_strain(double, double);
static  void  update_back_stress(Tensor *,Tensor *,double, Material);
static void compute_increment_plastic_strain_tensor(Tensor *,Tensor *,double);
static void   apply_plastic_corrector_stress_tensor(Tensor *,Tensor *,double,Material);

/**************************************************************/

State_Parameters Von_Mises_forward_euler(
  State_Parameters Inputs_SP,
  Material MatProp)
/*	
	Radial returning algorithm for the Von-Mises plastic criterium
*/
{

  int Ndim = NumberDimensions;

  Tensor p_trial;
  Tensor s_trial;
  Tensor relative_stress_tr;
  Tensor plastic_flow_direction;
  double s_trial_norm;
  double norm_relative_stress_tr;
  double Phi_tr;

  /* Load input state marameters */
  Tensor sigma_k1 = memory_to_tensor__TensorLib__(Inputs_SP.Stress,2);
  Tensor Back_stress = memory_to_tensor__TensorLib__(Inputs_SP.Back_stress,2);

  /*
    Define output state parameters
  */
  State_Parameters Outputs_VarCons;
  Outputs_VarCons.Increment_E_plastic = (double *)calloc(Ndim*Ndim,sizeof(double));

  Tensor Increment_E_plastic = memory_to_tensor__TensorLib__(Outputs_VarCons.Increment_E_plastic,2);

  /* Auxiliar variables during iterations */
  double delta_Gamma;

  /*
    Initialise parameters for the solver
  */
  double TOL = TOL_Radial_Returning;
  
  /*
    Decompose the elastic trial in volumetric and deviatoric components
  */
  p_trial = volumetric_component__TensorLib__(sigma_k1);
  s_trial = deviatoric_component__TensorLib__(sigma_k1,p_trial);

  s_trial_norm = EuclideanNorm__TensorLib__(s_trial); 

  /*
    Compute the trial relative stress
  */
  relative_stress_tr = subtraction__TensorLib__(s_trial,Back_stress);
  norm_relative_stress_tr = EuclideanNorm__TensorLib__(relative_stress_tr);

  /*
    Compute plastic flow
  */
  plastic_flow_direction = compute_plastic_flow_direction(&relative_stress_tr,norm_relative_stress_tr);

  /*
    Yield condition : Starting from incremental plastic strain equal to zero
  */
  Phi_tr = compute_yield_surface(norm_relative_stress_tr, Inputs_SP.EPS, MatProp);

  if(Phi_tr > TOL)
  {

    delta_Gamma = compute_increment_flow_rule(Phi_tr,MatProp);

    Outputs_VarCons.EPS = update_equivalent_plastic_strain(Inputs_SP.EPS, delta_Gamma);

    compute_increment_plastic_strain_tensor(&Increment_E_plastic,&plastic_flow_direction, delta_Gamma);

    update_back_stress(&Back_stress,&plastic_flow_direction,delta_Gamma,MatProp);

    apply_plastic_corrector_stress_tensor(&sigma_k1,&plastic_flow_direction,delta_Gamma,MatProp);

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
  free__TensorLib__(relative_stress_tr);
  free__TensorLib__(plastic_flow_direction);
  
  return Outputs_VarCons;
}

/**************************************************************/

static Tensor compute_plastic_flow_direction(
  Tensor * relative_stress_tr,
  double norm_relative_stress_tr)
{
  int Ndim = NumberDimensions;
  Tensor plastic_flow_direction = alloc__TensorLib__(2); 

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      plastic_flow_direction.N[i][j] = relative_stress_tr->N[i][j]/norm_relative_stress_tr;
    }
  }

  return plastic_flow_direction;
}

/**************************************************************/

static double compute_yield_surface(
  double norm_relative_stress_tr,
  double EPS_n,
  Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double Sigma_y = MatProp.yield_stress_0;

  if(MatProp.Hardening_Hughes)
  {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.Parameter_Hardening_Hughes;
    
    return norm_relative_stress_tr - sqrt(2./3.)*(Sigma_y + Hardening_modulus*theta*EPS_n);
  }
  else
  {
    return norm_relative_stress_tr - sqrt(2./3.)*Sigma_y;
  }
  
}


/**************************************************************/

static double compute_increment_flow_rule(
  double Phi_tr,
  Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));

  if(MatProp.Hardening_Hughes)
  {
    double Hardening_modulus = MatProp.Hardening_modulus;

    if(MatProp.Viscous_regularization)
    {
      double fluidity_param = MatProp.fluidity_param;
      double DT = 1.0;

      return (Phi_tr/(2*G))/(fluidity_param/(DT*2*G) + 1 + Hardening_modulus/(3*G));
    }
    else
    {
      return (Phi_tr/(2*G))/(1 + Hardening_modulus/(3*G));
    }
    
  }
  else
  {
    
    if(MatProp.Viscous_regularization)
    {
      double fluidity_param = MatProp.fluidity_param; 
      double DT = 1.0;

      return (Phi_tr/(2*G))/(fluidity_param/(DT*2*G) + 1);
    }
    else
    {
      return (Phi_tr/(2*G));
    }

  }

}


/**************************************************************/

static double update_equivalent_plastic_strain(
  double EPS_n,
  double delta_Gamma)
{
  return EPS_n + sqrt(2./3.)*delta_Gamma;
}

/**************************************************************/

static void update_back_stress(
  Tensor * Back_stress,
  Tensor * plastic_flow_direction,
  double delta_Gamma,
  Material MatProp)
{
  int Ndim = NumberDimensions;
  double DeltaH;

  if(MatProp.Hardening_Hughes)
  { 
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.Parameter_Hardening_Hughes;
    
    DeltaH = sqrt(2./3.)*(1-theta)*Hardening_modulus*delta_Gamma;
  }
  else
  {
    DeltaH = 0.0;
  }


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
  Tensor * Increment_E_plastic,
  Tensor * plastic_flow_direction, 
  double delta_Gamma)
{
  int Ndim = NumberDimensions;

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      Increment_E_plastic->N[i][j] = delta_Gamma*plastic_flow_direction->N[i][j];
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
