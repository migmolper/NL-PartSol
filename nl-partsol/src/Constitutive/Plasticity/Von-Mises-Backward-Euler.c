#include <math.h>
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
static void standard_error(char *);
static bool check_convergence(double, double, int, int);

static void compute_relative_stress(double *, double *, double *);

static void compute_plastic_flow_direction(double *, double *, double);

static double compute_yield_surface(double, double, Material);

static double compute_objective_function(double, double, double, Material);
static double compute_derivative_objective_function(double, Material);

static double compute_K(double, Material);
static double compute_D_K(double, Material);

static double compute_DeltaH(double, Material);
static double compute_D_DeltaH(Material);

static double update_increment_plastic_strain(double, double, double);
static double update_equivalent_plastic_strain(double, double);
static void update_back_stress(double *, double *, double, Material);

static void update_increment_plastic_strain_tensor(double *, double *, double);
static void apply_plastic_corrector_stress_tensor(double *, double *, double,
                                                  Material);

/**************************************************************/

State_Parameters Von_Mises_backward_euler(State_Parameters Inputs_SP,
                                          Material MatProp)
/*
        Radial returning algorithm for the Von-Mises plastic criterium
*/
{

  int Ndim = NumberDimensions;

  double relative_stress[3];
  double plastic_flow_direction[3];
  double relative_stress_norm;
  double Phi;
  double G;
  double d_G;

  /*
    Define output state parameters
  */
  State_Parameters Outputs_VarCons;

  /*
    Initialise solver parameters
  */
  double alpha_n = *Inputs_SP.EPS;
  double alpha_k = alpha_n;
  double delta_Gamma_k = 0;
  double TOL = TOL_Radial_Returning;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;
  bool Convergence = false;

  /*
    Compute the relative stress
  */
  compute_relative_stress(relative_stress, Inputs_SP.Stress,
                          Inputs_SP.Back_stress);
  relative_stress_norm = sqrt(relative_stress[0] * relative_stress[0] +
                              relative_stress[1] * relative_stress[1] +
                              relative_stress[2] * relative_stress[2]);

  /*
    Compute plastic flow
  */
  compute_plastic_flow_direction(plastic_flow_direction, relative_stress,
                                 relative_stress_norm);

  /*
    Check yield condition
  */
  Phi = compute_yield_surface(relative_stress_norm, alpha_n, MatProp);

  if (Phi > TOL) {

    /*
      Newton-Rapson solver
    */
    while (Convergence == false) {

      G = compute_objective_function(relative_stress_norm, delta_Gamma_k,
                                     alpha_k, MatProp);

      Convergence = check_convergence(G, TOL, Iter, MaxIter);

      d_G = compute_derivative_objective_function(alpha_k, MatProp);

      delta_Gamma_k = update_increment_plastic_strain(delta_Gamma_k, G, d_G);

      alpha_k = update_equivalent_plastic_strain(alpha_n, delta_Gamma_k);

      Iter++;
    }

    /*
      Update plastic deformation gradient
    */
    update_increment_plastic_strain_tensor(
        Inputs_SP.Increment_E_plastic, plastic_flow_direction, delta_Gamma_k);

    /*
      Update back stress
    */
    update_back_stress(Inputs_SP.Back_stress, plastic_flow_direction,
                       delta_Gamma_k, MatProp);

    /*
      Update stress tensor with the plastic corrector
    */
    apply_plastic_corrector_stress_tensor(
        Inputs_SP.Stress, plastic_flow_direction, delta_Gamma_k, MatProp);

    /*
      Update equivalent plastic strain and increment of plastic deformation
    */
    *Outputs_VarCons.EPS = alpha_k;
    Outputs_VarCons.Stress = Inputs_SP.Stress;
    Outputs_VarCons.Increment_E_plastic = Inputs_SP.Increment_E_plastic;
  } else {
    Outputs_VarCons.EPS = Inputs_SP.EPS;
    Outputs_VarCons.Stress = Inputs_SP.Stress;
    Outputs_VarCons.Increment_E_plastic = Inputs_SP.Increment_E_plastic;
  }

  return Outputs_VarCons;
}

/***************************************************************************/

static void standard_error(char *Error_message) {
  fprintf(stderr, "%s : %s !!! \n", "Error in plasticity_Von_Mises",
          Error_message);
  exit(EXIT_FAILURE);
}

/**************************************************************/

static bool check_convergence(double Error, double TOL, int Iter, int MaxIter) {
  bool convergence = false;
  double Error_relative;
  char Error_message[MAXW];

  if (Iter > MaxIter) {
    sprintf(Error_message, "%s",
            "Convergence not reached in the maximum number of iterations");
    printf("%e\n", Error_relative);
    standard_error(Error_message);
  }

  /*
    Compute relative error
  */
  if (Iter == 0) {

    Error0 = fabs(Error);
    Error_relative = fabs(Error / Error0);

    if (Error0 < TOL) {
      convergence = true;
    }

  } else {
    Error_relative = fabs(Error / Error0);

    if (Error_relative < TOL) {
      convergence = true;
    }
  }

  return convergence;
}

/**************************************************************/

static void compute_relative_stress(double *Relative_Stress, double *Stress,
                                    double *Back_stress) {
  double vol_Stress = (Stress[0] + Stress[1] + Stress[2]) / 3.0;

  Relative_Stress[0] = Stress[0] - vol_Stress - Back_stress[0];
  Relative_Stress[1] = Stress[1] - vol_Stress - Back_stress[1];
  Relative_Stress[2] = Stress[2] - vol_Stress - Back_stress[2];
}

/**************************************************************/

static void compute_plastic_flow_direction(double *plastic_flow_direction,
                                           double *Relative_Stress,
                                           double Relative_Stress_norm) {
  plastic_flow_direction[0] = Relative_Stress[0] / Relative_Stress_norm;
  plastic_flow_direction[1] = Relative_Stress[1] / Relative_Stress_norm;
  plastic_flow_direction[2] = Relative_Stress[2] / Relative_Stress_norm;
}

/**************************************************************/

static double compute_yield_surface(double relative_stress_norm, double alpha,
                                    Material MatProp) {
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E;   /* Elastic modulus */
  double G = E / (2 * (1 + nu));
  double K = compute_K(alpha, MatProp); /* Isotropic hardening */

  return relative_stress_norm - sqrt(2. / 3.) * K;
}

/**************************************************************/

static double compute_objective_function(double relative_stress_norm,
                                         double delta_Gamma, double alpha_k,
                                         Material MatProp) {
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E;   /* Elastic modulus */
  double G = E / (2 * (1 + nu));
  double K = compute_K(alpha_k, MatProp); /* Isotropic hardening */
  double DeltaH =
      compute_DeltaH(delta_Gamma, MatProp); /* Kinematic hardening */

  if (MatProp.Viscous_regularization) {
    double fluidity_param = MatProp.fluidity_param;
    double DT = 1.0;

    return relative_stress_norm - sqrt(2. / 3.) * (K + DeltaH) -
           DMAX(0, delta_Gamma * fluidity_param / DT) - 2 * G * delta_Gamma;
  } else {
    return relative_stress_norm - sqrt(2. / 3.) * (K + DeltaH) -
           2 * G * delta_Gamma;
  }
}

/**************************************************************/

static double compute_derivative_objective_function(double alpha_k,
                                                    Material MatProp) {
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E;   /* Elastic modulus */
  double G = E / (2 * (1 + nu));
  double D_K =
      compute_D_K(alpha_k, MatProp); /* Derivative of the isotropic hardening */
  double D_DeltaH =
      compute_D_DeltaH(MatProp); /* Derivative of the kinematic hardening */

  if (MatProp.Viscous_regularization) {
    double fluidity_param = MatProp.fluidity_param;
    double DT = 1.0;

    return -sqrt(2. / 3.) * (D_K + D_DeltaH) - 2 * G - fluidity_param / DT;
  } else {
    return -sqrt(2. / 3.) * (D_K + D_DeltaH) - 2 * G;
  }
}

/**************************************************************/

static double compute_K(double alpha, Material MatProp)
/*
  Isotropic hardening function
*/
{
  double Sigma_y = MatProp.kappa_0;

  if (MatProp.Hardening_Hughes) {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.Parameter_Hardening_Hughes;

    return Sigma_y + theta * Hardening_modulus * alpha;
  } else if (MatProp.Hardening_Cervera) {
    double Hardening_modulus = MatProp.Hardening_modulus;

    return Sigma_y * exp(-2 * Hardening_modulus * alpha / Sigma_y);
  } else if (MatProp.Hardening_Ortiz) {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double Exponent = MatProp.Exponent_Hardening_Ortiz;
    double Ref_PS = MatProp.Plastic_Strain_0;

    return Sigma_y * pow(1 + alpha / Ref_PS, 1.0 / Exponent);
  } else if (MatProp.Hardening_Voce) {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.theta_Hardening_Voce;
    double K_0 = MatProp.K_0_Hardening_Voce;
    double K_inf = MatProp.K_inf_Hardening_Voce;
    double delta = MatProp.delta_Hardening_Voce;

    return Sigma_y + theta * Hardening_modulus * alpha +
           (K_inf - K_0) * (1 - exp(-delta * alpha));
  } else {
    return Sigma_y;
  }
}

/**************************************************************/

static double compute_D_K(double alpha, Material MatProp)
/*
  Derivative of the isotropic hardening function
*/
{

  if (MatProp.Hardening_Hughes) {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.Parameter_Hardening_Hughes;

    return sqrt(2. / 3.) * theta * Hardening_modulus;
  } else if (MatProp.Hardening_Cervera) {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double Sigma_y = MatProp.kappa_0;

    return -2.0 * sqrt(2. / 3.) * Hardening_modulus *
           exp(-2.0 * Hardening_modulus * alpha / Sigma_y);
  } else if (MatProp.Hardening_Ortiz) {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double Exponent = MatProp.Exponent_Hardening_Ortiz;
    double Ref_PS = MatProp.Plastic_Strain_0;
    double Sigma_y = MatProp.kappa_0;

    return sqrt(2. / 3.) * (Sigma_y / (Ref_PS * Exponent)) *
           pow(1 + alpha / Ref_PS, 1.0 / Exponent - 1);
  } else if (MatProp.Hardening_Voce) {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.theta_Hardening_Voce;
    double K_0 = MatProp.K_0_Hardening_Voce;
    double K_inf = MatProp.K_inf_Hardening_Voce;
    double delta = MatProp.delta_Hardening_Voce;

    return sqrt(2. / 3.) * (theta * Hardening_modulus +
                            delta * (K_inf - K_0) * exp(-delta * alpha));
  } else {
    return 0.0;
  }
}

/**************************************************************/

static double compute_DeltaH(double delta_Gamma, Material MatProp)
/*
  Kinematic hardening function
*/
{

  if (MatProp.Hardening_Hughes) {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.Parameter_Hardening_Hughes;

    return sqrt(2. / 3.) * (1 - theta) * Hardening_modulus * delta_Gamma;
  } else if (MatProp.Hardening_Voce) {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.theta_Hardening_Voce;

    return sqrt(2. / 3.) * (1 - theta) * Hardening_modulus * delta_Gamma;
  } else {
    return 0.0;
  }
}

/**************************************************************/

static double compute_D_DeltaH(Material MatProp)
/*
  Derivative of the kinematic hardening function
*/
{

  if (MatProp.Hardening_Hughes) {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.Parameter_Hardening_Hughes;

    return sqrt(2. / 3.) * (1 - theta) * Hardening_modulus;
  } else if (MatProp.Hardening_Voce) {
    double Hardening_modulus = MatProp.Hardening_modulus;
    double theta = MatProp.theta_Hardening_Voce;

    return sqrt(2. / 3.) * (1 - theta) * Hardening_modulus;
  } else {
    return 0.0;
  }
}

/**************************************************************/

static double update_increment_plastic_strain(double delta_Gamma_k, double G,
                                              double d_G) {
  return delta_Gamma_k - G / d_G;
}

/**************************************************************/

static double update_equivalent_plastic_strain(double alpha,
                                               double delta_Gamma) {
  return alpha + sqrt(2. / 3.) * delta_Gamma;
}

/**************************************************************/

static void update_back_stress(double *Back_Stress,
                               double *Plastic_Flow_Direction,
                               double delta_Gamma, Material MatProp) {
  double DeltaH = compute_DeltaH(delta_Gamma, MatProp);

  Back_Stress[0] += sqrt(2. / 3.) * DeltaH * Plastic_Flow_Direction[0];
  Back_Stress[1] += sqrt(2. / 3.) * DeltaH * Plastic_Flow_Direction[1];
  Back_Stress[2] += sqrt(2. / 3.) * DeltaH * Plastic_Flow_Direction[2];
}

/**************************************************************/

static void
update_increment_plastic_strain_tensor(double *Increment_E_plastic,
                                       double *Plastic_Flow_Direction,
                                       double delta_Gamma) {
  Increment_E_plastic[0] = delta_Gamma * Plastic_Flow_Direction[0];
  Increment_E_plastic[1] = delta_Gamma * Plastic_Flow_Direction[1];
  Increment_E_plastic[2] = delta_Gamma * Plastic_Flow_Direction[2];
}

/**************************************************************/

static void
apply_plastic_corrector_stress_tensor(double *Stress,
                                      double *Plastic_Flow_Direction,
                                      double delta_Gamma, Material MatProp) {
  int Ndim = NumberDimensions;
  double nu = MatProp.nu; /* Poisson modulus */
  double E = MatProp.E;   /* Elastic modulus */
  double G = E / (2 * (1 + nu));

  Stress[0] -= 2 * G * delta_Gamma * Plastic_Flow_Direction[0];
  Stress[1] -= 2 * G * delta_Gamma * Plastic_Flow_Direction[1];
  Stress[2] -= 2 * G * delta_Gamma * Plastic_Flow_Direction[2];
}

/**************************************************************/
