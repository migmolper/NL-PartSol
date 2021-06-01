#include "nl-partsol.h"

/*
  Call global variables
*/
double DeltaTimeStep;
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

static Tensor compute_plastic_flow_direction(Tensor, double);

static double compute_yield_surface(double,double,double,Material);
static double compute_derivative_yield_surface(Material);

static double compute_K(double,Material);
static double compute_D_K(Material);

static double compute_DeltaH(double,Material);
static double compute_D_DeltaH(Material);

static double update_increment_plastic_strain(double, double, double);
static double update_equivalent_plastic_strain(double, double);
static  void  update_back_stress(Tensor,Tensor,double,Material);

static Tensor compute_increment_viscoplastic_plastic_strain_tensor(Tensor,double);
static void   apply_plastic_corrector_stress_tensor(Tensor,Tensor,double,Material);

/**************************************************************/

Plastic_status finite_strains_viscoplasticity_Von_Mises_Perzyna(
  Tensor P_p,
  Plastic_status Inputs_VarCons, 
  Material MatProp)
/*
  Finite strains plasticity following the apporach of Ortiz and Camacho
*/
{
  int Ndim = NumberDimensions;

  /* Define auxiliar variables */
  Plastic_status Outputs_VarCons;
  Tensor F_m1_plastic = Inputs_VarCons.F_m1_plastic_p;
  Tensor F_total = Inputs_VarCons.F_n1_p;
  Tensor F_m1_total;
  Tensor F_trial_elastic;
  Tensor C_trial_elastic;
  Tensor E_trial_elastic;
  Tensor D_F_plastic;
  Tensor Fm1_plastic;
  Tensor T_p = alloc__TensorLib__(2);

  /* Compute the elastic right Cauchy-Green tensor using the intermediate configuration. */ 
  F_trial_elastic = matrix_product__TensorLib__(F_total,F_m1_plastic);

  C_trial_elastic = right_Cauchy_Green__Particles__(F_trial_elastic);

  /* Calculation of the small strain tensor */
  E_trial_elastic = logarithmic_strains__Particles__(C_trial_elastic);

  /* Calculation of the trial stress tensor using the trial small strain tensor */
  T_p = LinearElastic(T_p, E_trial_elastic, MatProp);

  /* Start plastic algorithm in infinitesimal strains */
  Outputs_VarCons = infinitesimal_strains_viscoplasticity_Von_Mises_Perzyna(T_p, Inputs_VarCons, MatProp);

  /* Use the Cuitiño & Ortiz exponential maping to compute the increment of plastic finite strains */
  update_plastic_deformation_gradient__Particles__(Outputs_VarCons.Increment_E_plastic, F_m1_plastic);

  /* Get the First Piola-Kirchhoff stress tensor (P_p) */
  F_m1_total = Inverse__TensorLib__(F_total);

  for(int i = 0 ; i < Ndim  ; i++)
  {
    for(int j = 0 ; j < Ndim  ; j++)    
    {

      P_p.N[i][j] = 0.0;

     for(int k = 0 ; k < Ndim  ; k++)
     {
        P_p.N[i][j] += T_p.N[i][k]*F_m1_total.N[k][j];
      }
    }
  }

  /* Free memory */
  free__TensorLib__(F_trial_elastic);
  free__TensorLib__(C_trial_elastic);
  free__TensorLib__(E_trial_elastic);
  free__TensorLib__(T_p);
  free__TensorLib__(F_m1_total);

  return Outputs_VarCons;
}

/**************************************************************/

Plastic_status infinitesimal_strains_viscoplasticity_Von_Mises_Perzyna(
  Tensor sigma_k1,
  Plastic_status Inputs_VarCons,
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
  double d_Phi;
  Tensor Increment_E_plastic;

  /* Load material marameters */
  double EPS = Inputs_VarCons.EPS;
  Tensor Back_stress = Inputs_VarCons.Back_stress;

  /* Auxiliar variables during iterations */
  double EPS_k;
  double delta_Gamma_k;
 
  /*
    Initialise convergence parameters for the solver
  */
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
  plastic_flow_direction = compute_plastic_flow_direction(relative_stress,relative_stress_norm);

  /*
    Yield condition : Starting from incremental plastic strain equal to zero
  */
  EPS_k = EPS;
  delta_Gamma_k = 0;
  Phi = compute_yield_surface(relative_stress_norm, delta_Gamma_k, EPS_k, MatProp);

  if(Phi > TOL)
  {     
    /*
      Newton-Rapson solver
    */
    while(Convergence == false)
    {
          
      Phi = compute_yield_surface(relative_stress_norm, delta_Gamma_k, EPS_k, MatProp);

      d_Phi = compute_derivative_yield_surface(MatProp);

      delta_Gamma_k = update_increment_plastic_strain(delta_Gamma_k, Phi, d_Phi);

      EPS_k = update_equivalent_plastic_strain(EPS, delta_Gamma_k);

      Convergence = check_convergence(Phi,TOL,Iter,MaxIter);
	     
      if(Convergence == false)
      {
        Iter++;
      }
    }

    /*
      Update plastic deformation gradient
    */
    Increment_E_plastic = compute_increment_viscoplastic_plastic_strain_tensor(plastic_flow_direction, delta_Gamma_k);

    /*
      Update back stress
    */
    update_back_stress(Back_stress,plastic_flow_direction,delta_Gamma_k,MatProp);

    /*
      Update stress tensor in the deformed configuration
    */
    apply_plastic_corrector_stress_tensor(sigma_k1, plastic_flow_direction, delta_Gamma_k, MatProp);

  }
  else
  {
    Increment_E_plastic = alloc__TensorLib__(2);
  }

  /*
    Free memory
  */
  free__TensorLib__(s_trial);
  free__TensorLib__(p_trial);
  free__TensorLib__(relative_stress);
  free__TensorLib__(plastic_flow_direction);

  /*
    Define output varible
  */
  Plastic_status Outputs_VarCons;
  Outputs_VarCons.EPS = EPS_k;
  Outputs_VarCons.Increment_E_plastic = Increment_E_plastic;
  
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
      Error0 = Error;
      Error_relative = Error/Error0;
    }
    else
    {
      Error_relative = Error/Error0;
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
  Tensor relative_stress,
  double relative_stress_norm)
{
  int Ndim = NumberDimensions;
  Tensor plastic_flow_direction = alloc__TensorLib__(2); 

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      plastic_flow_direction.N[i][j] = relative_stress.N[i][j]/relative_stress_norm;
    }
  }

  return plastic_flow_direction;
}

/**************************************************************/

static double compute_yield_surface(
  double relative_stress_norm,
  double delta_Gamma,
  double EPS_k,
  Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double fluidity_param = MatProp.fluidity_param;
  double K = compute_K(EPS_k,MatProp);
  double DeltaH = compute_DeltaH(delta_Gamma,MatProp);

  return - sqrt(2./3.)*(K + DeltaH) + relative_stress_norm - (2*G + fluidity_param/DeltaTimeStep)*delta_Gamma;
}

/**************************************************************/


static double compute_derivative_yield_surface(
  Material MatProp)
{
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E; /* Elastic modulus */
  double G = E/(2*(1+nu));
  double fluidity_param = MatProp.fluidity_param;
  double D_K = compute_D_K(MatProp);
  double D_DeltaH = compute_D_DeltaH(MatProp);

  return - sqrt(2./3.)*(D_K + D_DeltaH) - (2*G + fluidity_param/DeltaTimeStep);
}

/**************************************************************/

static double compute_K(
  double EPS_k1,
  Material MatProp)
/*
  Isotropic hardening function
*/
{
  double H = MatProp.isotropic_hardening_modulus;
  double theta = MatProp.isotropic_hardening_theta;
  double Sigma_y = MatProp.yield_stress_0;

  double K = Sigma_y + theta*H*EPS_k1;

  return K;
}

/**************************************************************/

static double compute_D_K(
  Material MatProp)
/*
  derivative of the isotropic hardening function
*/
{
  double H = MatProp.isotropic_hardening_modulus;
  double theta = MatProp.isotropic_hardening_theta;

  double DK = sqrt(2./3.)*theta*H;

  return DK;
}

/**************************************************************/

static double compute_DeltaH(
  double delta_Gamma,
  Material MatProp)
/*
  Kinematic hardening function
*/
{
  double H = MatProp.kinematic_hardening_modulus;
  double beta = MatProp.kinematic_hardening_beta;

  double DeltaH = sqrt(2./3.)*(1-beta)*H*delta_Gamma;

  return DeltaH;
}

/**************************************************************/

static double compute_D_DeltaH(
  Material MatProp)
/*
  derivative of the kinematic hardening function
*/
{
  double H = MatProp.kinematic_hardening_modulus;
  double beta = MatProp.kinematic_hardening_beta;

  double D_DeltaH = sqrt(2./3.)*(1-beta)*H;

  return D_DeltaH;
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

static double update_equivalent_plastic_strain(
  double EPS,
  double delta_Gamma)
{
  double EPS_k = EPS + sqrt(2./3.)*delta_Gamma;
  return EPS_k;
}

/**************************************************************/

static void update_back_stress(
  Tensor Back_stress,
  Tensor plastic_flow_direction,
  double delta_Gamma,
  Material MatProp)
{
  int Ndim = NumberDimensions;
  double DeltaH = compute_DeltaH(delta_Gamma,MatProp);

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      Back_stress.N[i][j] += sqrt(2./3.)*DeltaH*plastic_flow_direction.N[i][j];
    }
  }

}

/**************************************************************/

static Tensor compute_increment_viscoplastic_plastic_strain_tensor(
  Tensor plastic_flow_direction, 
  double delta_Gamma)
{
  int Ndim = NumberDimensions;
  Tensor D_E_plastic = alloc__TensorLib__(2); 

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      D_E_plastic.N[i][j] = delta_Gamma*plastic_flow_direction.N[i][j];
    }
  }


  return D_E_plastic;
}

/**************************************************************/

static void apply_plastic_corrector_stress_tensor(
  Tensor sigma_k1,
  Tensor plastic_flow_direction, 
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
      sigma_k1.N[i][j] -= 2*G*delta_Gamma*plastic_flow_direction.N[i][j];
    }
  }
}

/**************************************************************/
