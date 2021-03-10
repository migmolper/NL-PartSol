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
static double compute_yield_surface(double, double, double, double, double);
static double compute_hardening_modulus(double,double,double,double);
static double compute_limit_between_classic_apex_algorithm(double,double,double,double,double,double,double,double);
static double update_increment_plastic_strain(double, double, double);
static double update_cohesion_modulus(double,double,double,double);
static double compute_yield_surface_classical(double,double,double,double,double,double,double,double,double);
static double compute_yield_surface_apex(double,double,double,double,double,double,double,double,double,double);
static double compute_derivative_yield_surface_classical(double,double,double,double,double,double);
static double compute_derivative_yield_surface_apex(double,double,double,double,double,double,double,double);
static double update_equivalent_plastic_strain_classical(double,double,double);
static double update_equivalent_plastic_strain_apex(double,double,double,double,double);
static Tensor compute_increment_plastic_strain_classical(Tensor,double,double,double);
static Tensor compute_increment_plastic_strain_apex(Tensor,double,double,double);
static void   update_pq_classical(double,Tensor,double,double,double,double,double,double);
static void   update_pq_apex(double,Tensor,double,double,double);

/**************************************************************/

Plastic_status finite_strains_plasticity_Drucker_Prager_Sanavia(
  Tensor SPK,
  Tensor C_total,
  Tensor F_plastic,
  Plastic_status Inputs_VarCons,
  Material MatProp)
/*
  Finite strains plasticity following the apporach of Ortiz and Camacho
*/
{

  /* Define output */
  Plastic_status Outputs_VarCons;

  double p;
  Tensor q;

  Tensor C_elastic = alloc__TensorLib__(2);
  Tensor log_Strain;
  Tensor log_Stress = alloc__TensorLib__(2);
  Tensor M_Mandel;
  Tensor D_plastic;
  Tensor D_F_plastic;

  /* Compute the elastic right Cauchy-Green tensorn. */ 
  covariant_push_forward_tensor__TensorLib__(C_elastic, C_total, F_plastic);

  /* Calculation of the logarithmic/Henky strain tensor */
  log_Strain = logarithmic_strains__Particles__(C_elastic);

  /* Compute the logaritmic stress meassure using Hencky model */
  log_Stress = compute_log_stress_tensor__Hencky__(log_Stress, log_Strain, MatProp);

  /* Compute volumetric and deviatoric decomposition of the "trial elastic Kirchhoff stress tensor" */
  p = volumetric_component__TensorLib__(log_Stress);
  q = deviatoric_component__TensorLib__(log_Stress,p);

  /* Start plastic algorithm in infinitesimal strains */
  Outputs_VarCons = infinitesimal_strains_plasticity_Drucker_Prager_Sanavia(p, q, Inputs_VarCons, MatProp);

  /* Update the "trial elastic Kirchhoff stress tensor" */
  assemble_volumetric_deviatoric__TensorLib__(log_Stress,p,q);

  /* In the special case of elastic isotropy, M == T, since skew-symmetric part of M vanishes */
  M_Mandel = log_Stress;

  /* Use the Cuitiño & Ortiz exponential maping to compute the increment of plastic finite strains */
  D_plastic = Outputs_VarCons.Increment_E_plastic;
  D_F_plastic = increment_Deformation_Gradient_exponential_strains__Particles__(D_plastic);

  /* Compute the new plastic deformation gradient */
  update_plastic_deformation_gradient__Particles__(D_F_plastic, F_plastic);

  /* Compute the updated elastic right Cauchy-Green tensor. */ 
  covariant_push_forward_tensor__TensorLib__(C_elastic, C_total, F_plastic);

  /* Get second piola-kirchhoff stress from the Mandel stress tensor */
  Mandel_to_second_Piola_Kirchhoff__Particles__(SPK,M_Mandel,C_elastic);

  /* Free memory */
  free__TensorLib__(C_elastic);
  free__TensorLib__(log_Strain);
  free__TensorLib__(log_Stress);
  free__TensorLib__(q);
  free__TensorLib__(D_plastic);
  free__TensorLib__(D_F_plastic);

  return Outputs_VarCons;
}

/**************************************************************/

Plastic_status infinitesimal_strains_plasticity_Drucker_Prager_Sanavia(
  double p,
  Tensor q,
  Plastic_status Inputs_VarCons,
  Material MatProp)
/*	
	Radial returning algorithm for the Drucker-Prager de Sanavia algorithm
*/
{

  int Ndim = NumberDimensions;

  /* Auxiliar parameters for the computations */
  double delta_Gamma;
  double delta_Gamma_1;
  double delta_Gamma_2;
  double q_norm;
  double p_norm;
  double Phi;
  double d_Phi;	
  double p_lim;
  Tensor D_plastic = alloc__TensorLib__(2);

  /* Load material parameters */
  double alpha_F = MatProp.alpha_F_Drucker_Prager;
  double alpha_Q = MatProp.alpha_Q_Drucker_Prager;
  double beta    = MatProp.beta_Drucker_Prager;
  double c_0     = MatProp.cohesion_reference;
  double E_0p    = MatProp.E_plastic_reference;
  double N_exp   = MatProp.hardening_exp;
  double nu = MatProp.nu;
  double E = MatProp.E;
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double EPS_k = Inputs_VarCons.EPS;
  double cohesion_k = Inputs_VarCons.Cohesion;
  double H;
  if(MatProp.Hardening_Ortiz)
  {
    H = 0;
  }
  else
  {
    H = MatProp.hardening_modulus;
  }

  /*
    Initialise convergence prameters for the solver
  */
  double TOL = TOL_Radial_Returning;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;
  bool Convergence = false;

  /*
    Compute the norm of the deviatoric tensor
  */ 
  p_norm = p;
  q_norm = EuclideanNorm__TensorLib__(q); 

  /*
    Yield condition : Starting from incremental plastic strain equal to zero
  */
  delta_Gamma = 0.0;
  Phi = compute_yield_surface(p_norm, q_norm, cohesion_k, alpha_F, beta);

  if(Phi > TOL)
  {

    if(MatProp.Hardening_Ortiz)
    {
      H = compute_hardening_modulus(EPS_k, c_0, E_0p, N_exp);
    }

    p_lim = compute_limit_between_classic_apex_algorithm(q_norm, cohesion_k, H, K, G, alpha_F, alpha_Q, beta);

    /*
	    Classical plastic iterator
    */
    if(p_norm <= p_lim)
    {

     while(Convergence == false)
     {

       d_Phi = compute_derivative_yield_surface_classical(H, K, G, alpha_F, alpha_Q, beta);

       delta_Gamma = update_increment_plastic_strain(delta_Gamma, Phi, d_Phi);

       EPS_k = update_equivalent_plastic_strain_classical(EPS_k, delta_Gamma, alpha_Q);

       cohesion_k = update_cohesion_modulus(EPS_k, c_0, E_0p, N_exp);

       if(MatProp.Hardening_Ortiz)
       {
          H = compute_hardening_modulus(EPS_k, c_0, E_0p, N_exp);
       }

       Phi = compute_yield_surface_classical(q_norm, p_norm, delta_Gamma, cohesion_k, K, G, alpha_F, alpha_Q, beta);

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
    D_plastic = compute_increment_plastic_strain_classical(q,q_norm,delta_Gamma,alpha_Q);
    update_pq_classical(p,q,q_norm,delta_Gamma,alpha_Q,beta,K,G);

  }
  /*
    Apex plastic iterator
  */
  else
  {
    delta_Gamma_1 = q_norm/(2*G);
    delta_Gamma_2 = 0.0;
    delta_Gamma = delta_Gamma_1 + delta_Gamma_2;

    while(Convergence == false)
    {

      d_Phi = compute_derivative_yield_surface_apex(H, q_norm, delta_Gamma, K, G, alpha_F, alpha_Q, beta);

      delta_Gamma_2 = update_increment_plastic_strain(delta_Gamma_2, Phi, d_Phi);

      delta_Gamma = delta_Gamma_1 + delta_Gamma_2;

      Phi = compute_yield_surface_apex(p_norm, q_norm, H, delta_Gamma, cohesion_k, K, G, alpha_F, alpha_Q, beta);

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
  EPS_k = update_equivalent_plastic_strain_apex(EPS_k, q_norm, delta_Gamma, G, alpha_Q);
  D_plastic = compute_increment_plastic_strain_apex(q,delta_Gamma,G,alpha_Q);


  update_pq_apex(p, q, delta_Gamma, alpha_Q, K);

}

}

  /*
    Define output varible
  */
  Plastic_status Outputs_VarCons;
  Outputs_VarCons.EPS = EPS_k;
  Outputs_VarCons.Cohesion = cohesion_k;
  Outputs_VarCons.Increment_E_plastic = D_plastic;

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

static double compute_yield_surface(
  double p_norm,
  double q_norm,
  double cohesion_k,
  double alpha_F, 
  double beta)
{

  double Phi = 3*alpha_F*p_norm + q_norm - beta*cohesion_k;

  return Phi;
}

/**************************************************************/

static double compute_hardening_modulus(
  double EPS_k,  
  double c_0,
  double E_0p,
  double N_exp)
{
  double basis = 1 + EPS_k/E_0p;
  double exp = 1/N_exp - 1;

  double H = (c_0/(N_exp*E_0p))*pow(basis,exp);

  return H;
}

/**************************************************************/

static double compute_limit_between_classic_apex_algorithm(
  double q_norm,
  double cohesion_k,
  double H,
  double K,
  double G,
  double alpha_F,
  double alpha_Q,
  double beta)
{

  double aux1 = 3*alpha_Q*K*q_norm/(2*G);
  double aux2 = (beta/(3*alpha_F))*((q_norm*H/(2*G))*sqrt(1+3*DSQR(alpha_Q)) + cohesion_k);

  double p_lim = - aux1 - aux2;

  return p_lim;

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

static double update_cohesion_modulus(
  double EPS_k1,
  double c_0,
  double E_0p,
  double N_exp)
{
  double exp   = 1/N_exp;
  double basis = 1 + EPS_k1/E_0p;

  double cohesion_k1 = c_0*pow(basis,exp);

  return cohesion_k1;	
}


/**************************************************************/

static double compute_yield_surface_classical(
  double q_norm,
  double p_norm,
  double delta_Gamma,
  double cohesion_k1,
  double K,
  double G,
  double alpha_F,
  double alpha_Q,
  double beta)
{

  double Phi = q_norm - 2*G*delta_Gamma - 3*alpha_F*(p_norm + 3*K*alpha_Q*delta_Gamma) - beta*cohesion_k1;

  return Phi;
}

/**************************************************************/

static double compute_yield_surface_apex(
  double p_norm,
  double q_norm,
  double H, 
  double delta_Gamma,
  double cohesion_k,
  double K,
  double G,
  double alpha_F,
  double alpha_Q,
  double beta)
{

  double aux = sqrt(DSQR(q_norm/(2*G)) + 3*DSQR(alpha_Q)*DSQR(delta_Gamma)); 

  double Phi = (beta/(3*alpha_F))*(cohesion_k + H*aux) - p_norm + 3*K*alpha_Q*delta_Gamma;

  return Phi;
}

/**************************************************************/

static double compute_derivative_yield_surface_classical(
  double H,
  double K,
  double G,
  double alpha_F,
  double alpha_Q,
  double beta)
{

  double d_Phi = - 9*K*alpha_F*alpha_Q - 2*G - H*beta*sqrt(3*DSQR(alpha_Q) + 1);

  return d_Phi;
}

/**************************************************************/

static double compute_derivative_yield_surface_apex(
  double H,
  double q_norm,
  double delta_Gamma,
  double K,
  double G,
  double alpha_F,
  double alpha_Q,
  double beta)
{
  double aux1 = 3*H*beta*DSQR(alpha_Q)*delta_Gamma;
  double aux2 = 3*alpha_F*sqrt(DSQR(q_norm/(2*G)) + 3*DSQR(alpha_Q)*DSQR(delta_Gamma));

  double d_Phi = 3*alpha_Q*K + aux1/aux2;

  return d_Phi;
}

/**************************************************************/

static double update_equivalent_plastic_strain_classical(
  double EPS_k,
  double delta_Gamma,
  double alpha_Q)
{

  double EPS_k1 = EPS_k + delta_Gamma*sqrt(3*DSQR(alpha_Q) + 1);

  return EPS_k1;
}

/**************************************************************/

static double update_equivalent_plastic_strain_apex(
  double EPS_k,
  double q_norm,
  double delta_Gamma,
  double G,
  double alpha_Q)
{

  double EPS_k1 = EPS_k + sqrt(DSQR(q_norm/(2*G)) + 3*DSQR(alpha_Q)*DSQR(delta_Gamma));

  return EPS_k1;
}

/**************************************************************/

static Tensor compute_increment_plastic_strain_classical(
  Tensor q,
  double q_norm,
  double delta_Gamma,
  double alpha_Q)
{

  int Ndim = NumberDimensions;

  Tensor D_plastic = alloc__TensorLib__(2); 

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
     D_plastic.N[i][j] = alpha_Q*delta_Gamma*(i==j) + (delta_Gamma/q_norm)*q.N[i][j];
   }
 }

 return D_plastic;
}

/**************************************************************/

static Tensor compute_increment_plastic_strain_apex(
  Tensor q,
  double delta_Gamma,
  double G,
  double alpha_Q)
{

  int Ndim = NumberDimensions;

  double aux1 = alpha_Q*delta_Gamma;
  double aux2 = 1/(2*G);

  Tensor D_plastic = alloc__TensorLib__(2);	

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
     D_plastic.N[i][j] = aux1*(i==j) + aux2*q.N[i][j];
   }
 }

 return D_plastic;
}


/**************************************************************/

static void update_pq_classical(
  double p,
  Tensor q, 
  double q_norm,
  double delta_Gamma,
  double alpha_Q,
  double beta,
  double K,
  double G)
{
  int Ndim = NumberDimensions;

  for(int i = 0 ; i<Ndim ; i++)
  {
    p = p - 3*K*alpha_Q*delta_Gamma;

    for(int j = 0 ; j<Ndim ; j++)
    {
      q.N[i][j] = (1 - 2*G*delta_Gamma/q_norm)*q.N[i][j];
    }
 }

}

/**************************************************************/

static void update_pq_apex(
  double p,
  Tensor q,
  double delta_Gamma,
  double alpha_Q,
  double K)
{
  int Ndim = NumberDimensions;

  for(int i = 0 ; i<Ndim ; i++)
  {
    p = p - 3*K*alpha_Q*delta_Gamma;
  
    for(int j = 0 ; j<Ndim ; j++)
    {
      q.N[i][j] = 0.0;          
    } 
  }

}

/**************************************************************/





