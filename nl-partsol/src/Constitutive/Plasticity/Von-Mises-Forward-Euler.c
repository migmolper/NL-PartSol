#include "nl-partsol.h"

/*
  Call global variables
*/
double DeltaTimeStep;

/*
  Auxiliar functions 
*/

static void compute_relative_stress(double *,double *, double *);
static void compute_plastic_flow_direction(double *, double *,double);
static double compute_yield_surface(double,double,Material);
static double compute_DeltaH(double,Material);
static double compute_increment_flow_rule(double, Material);
static double update_equivalent_plastic_strain(double, double);
static  void  update_back_stress(double *,double *,double, Material);
static void compute_increment_plastic_strain_tensor(double *, double *,double);
static void   apply_plastic_corrector_stress_tensor(double *, double *,double,Material);

/**************************************************************/

State_Parameters Von_Mises_forward_euler(
  State_Parameters Inputs_SP,
  Material MatProp)
/*	
	Radial returning algorithm for the Von-Mises plastic criterium
*/
{

  int Ndim = NumberDimensions;

  double relative_stress[3];
  double plastic_flow_direction[3];
  double relative_stress_norm;
  double Phi_tr;

  /*
    Define output state parameters
  */
  State_Parameters Outputs_VarCons;
  
  /* Auxiliar variables during iterations */
  double delta_Gamma;

  /*
    Initialise parameters for the solver
  */
  double TOL = TOL_Radial_Returning;
  
  /*
    Compute the relative stress
  */
  compute_relative_stress(relative_stress,Inputs_SP.Stress,Inputs_SP.Back_stress);
  relative_stress_norm = sqrt(relative_stress[0]*relative_stress[0] + relative_stress[1]*relative_stress[1] + relative_stress[2]*relative_stress[2]);

  /*
    Compute plastic flow
  */
  compute_plastic_flow_direction(plastic_flow_direction,relative_stress,relative_stress_norm);

  /*
    Yield condition : Starting from incremental plastic strain equal to zero
  */
  Phi_tr = compute_yield_surface(relative_stress_norm, Inputs_SP.Equiv_Plast_Str, MatProp);

  if(Phi_tr > TOL)
  {

    delta_Gamma = compute_increment_flow_rule(Phi_tr,MatProp);

    Outputs_VarCons.Equiv_Plast_Str = update_equivalent_plastic_strain(Inputs_SP.Equiv_Plast_Str, delta_Gamma);

    compute_increment_plastic_strain_tensor(Inputs_SP.Increment_E_plastic,plastic_flow_direction, delta_Gamma);

    update_back_stress(Inputs_SP.Back_stress,plastic_flow_direction,delta_Gamma,MatProp);

    apply_plastic_corrector_stress_tensor(Inputs_SP.Stress,plastic_flow_direction,delta_Gamma,MatProp);

    Outputs_VarCons.Increment_E_plastic = Inputs_SP.Increment_E_plastic;
  }
  else
  {
    Outputs_VarCons.Equiv_Plast_Str = Inputs_SP.Equiv_Plast_Str;
    Outputs_VarCons.Increment_E_plastic = Inputs_SP.Increment_E_plastic;
  }
  
  return Outputs_VarCons;
}

/**************************************************************/

static void compute_relative_stress(
  double * Relative_Stress,
  double * Stress,
  double * Back_stress)
{
  double vol_Stress = (Stress[0] + Stress[1] + Stress[2])/3.0;

  Relative_Stress[0] = Stress[0] - vol_Stress - Back_stress[0];
  Relative_Stress[1] = Stress[1] - vol_Stress - Back_stress[1];
  Relative_Stress[2] = Stress[2] - vol_Stress - Back_stress[2];
}

/**************************************************************/

static void compute_plastic_flow_direction(
  double * plastic_flow_direction,
  double * Relative_Stress,
  double Relative_Stress_norm)
{
  plastic_flow_direction[0] = Relative_Stress[0]/Relative_Stress_norm;
  plastic_flow_direction[1] = Relative_Stress[1]/Relative_Stress_norm;
  plastic_flow_direction[2] = Relative_Stress[2]/Relative_Stress_norm;
}

/**************************************************************/

static double compute_yield_surface(
  double relative_stress_norm,
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
    
    return relative_stress_norm - sqrt(2./3.)*(Sigma_y + Hardening_modulus*theta*EPS_n);
  }
  else
  {
    return relative_stress_norm - sqrt(2./3.)*Sigma_y;
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
      double DT = DeltaTimeStep;

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
      double DT = DeltaTimeStep;

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
  double * Back_Stress,
  double * Plastic_Flow_Direction,
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

  Back_Stress[0] += sqrt(2./3.)*DeltaH*Plastic_Flow_Direction[0];
  Back_Stress[1] += sqrt(2./3.)*DeltaH*Plastic_Flow_Direction[1];
  Back_Stress[2] += sqrt(2./3.)*DeltaH*Plastic_Flow_Direction[2];

}

/**************************************************************/

static void compute_increment_plastic_strain_tensor(
  double * Increment_E_plastic,
  double * plastic_flow_direction, 
  double delta_Gamma)
{

  Increment_E_plastic[0] = delta_Gamma*plastic_flow_direction[0];
  Increment_E_plastic[1] = delta_Gamma*plastic_flow_direction[1];
  Increment_E_plastic[2] = delta_Gamma*plastic_flow_direction[2];
}

/**************************************************************/

static void apply_plastic_corrector_stress_tensor(
  double * Stress, 
  double * Plastic_Flow_Direction, 
  double delta_Gamma, 
  Material MatProp)
{
  int Ndim = NumberDimensions;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));

  Stress[0] -= 2*G*delta_Gamma*Plastic_Flow_Direction[0];
  Stress[1] -= 2*G*delta_Gamma*Plastic_Flow_Direction[1];
  Stress[2] -= 2*G*delta_Gamma*Plastic_Flow_Direction[2];
}

/**************************************************************/
