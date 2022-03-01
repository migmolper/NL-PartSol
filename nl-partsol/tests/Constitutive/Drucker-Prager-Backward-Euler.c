#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/*! \struct Material
 * Properties of a general material model
 */
typedef struct {

  /*!
   * Index of the material
   */
  int Id;

  /*!
   * Name of the material
   */
  char Type[100];

  /*!
   * Material celerity
   */
  double Cel;

  /*!
   * Initial density (mixture/fluid/solid)
   */
  double rho;

  /*!
   * Elastic modulus
   */
  double E;

  /*!
   * Poisson ratio
   */
  double nu;

  /*!
   * Compressibility
   */
  double Compressibility;

  /*!
   * Fluid parameters
   */
  double ReferencePressure;
  double Viscosity;
  double n_Macdonald_model;

  /*!
   * Activate eigenerosion-fracture modulus (Eigenerosion/Eigensoftening)
   */
  bool Eigenerosion;
  bool Eigensoftening;
  double Ceps; /*! Normalizing constant (Eigenerosion/Eigensoftening) */
  double Gf;   /*! Failure energy (Eigenerosion) */
  double ft;   /*! Tensile strengt of the material (Eigensoftening) */
  double heps; /*! Bandwidth of the cohesive fracture (Eigensoftening) */
  double Wc;   /*! Critical opening displacement (Eigensoftening) */

  /*!
   * Integration algorithm for plasticity
   */
  char Plastic_Solver[100];

  /*!
   * General plastic parameters
   */
  double yield_stress_0;
  double Hardening_modulus;
  double atmospheric_pressure;

  /*!
   * Frictional material (Borja et al. 2003)
   * */
  char Yield_Function_Frictional[100];
  double m_Frictional;
  double c0_Frictional;
  double phi_Frictional;
  double psi_Frictional;

  /*!
   * Hardening Hughes
   */
  bool Hardening_Hughes;
  double Parameter_Hardening_Hughes;

  /*!
   * Hardening Cervera
   * */
  bool Hardening_Cervera;

  /*!
   * Hardening Ortiz
   * */
  bool Hardening_Ortiz;
  double Exponent_Hardening_Ortiz;
  double Reference_Plastic_Strain_Ortiz;

  /*!
   * Hardening Voce
   * */
  bool Hardening_Voce;
  double K_0_Hardening_Voce;
  double K_inf_Hardening_Voce;
  double delta_Hardening_Voce;
  double theta_Hardening_Voce;

  /*!
   * Hardening Borja et al. 2003
   * */
  bool Hardening_Borja;
  double a_Hardening_Borja[3];
  double alpha_Hardening_Borja;

  /*!
   * Viscoplasticity parameters
   * */
  bool Viscous_regularization;
  double fluidity_param;

  /*!
   * Activate auxiliar techniques
   * */
  bool Locking_Control_Fbar;
  double alpha_Fbar;

} Material;

/*!
 * \struct State_Parameters
 */
typedef struct {
  /*!
   * Particle identifier
   * */
  int Particle_Idx;

  /*!
   * Stress/strain parameters
   * */
  double *Stress;
  double *Strain;
  double Pressure;

  /*!
   * Finite strain kinematic parameters
   * */
  double *d_phi;
  double *D_phi_n1;
  double *Fbar;
  double *rate_D_phi_n1;
  double J;

  /*!
   * Plasticity parameters
   * */
  double *Back_stress;
  double *D_phi_e;
  double *Increment_E_plastic;

  double Cohesion;
  double Yield_stress;

  // Internal variables
  double *Equiv_Plast_Str; // Equivalent plastic strain
  double *Kappa;           // Hardening Parameter

} State_Parameters;

static int __compute_trial_D_phi_e(
    double *e_hencky /**< [out] Trial elastic strain tensor (henky). */,
    double *EE /**< [out] Principal directions (right). */,
    double *EE_m1 /**< [out] Principal directions (left). */,
    double *D_phi_e_tr /**< [out] (trial) Elastic deformation gradient.*/,
    const double *d_phi /**< [in] Incremental deformation gradient. */,
    const double *D_phi_e /**< [in] Elastic deformation gradient. */);

static int __compute_D_phi_e(
    double *D_phi_e /**< [out] Elastic deformation gradient. */,
    const double *D_phi_e_tr /**< [out] (trial) Elastic deformation gradient.*/,
    const double *EE /**< [out] Principal directions (right). */,
    const double *EE_m1 /**< [out] Principal directions (left). */,
    const double *Increment_E_plastic /**< [in] Increment plastic strain */);

static int __trial_elastic(
    double *kirchhoff_tr_vol /**< [in/out] Volumetric elastic stress tensor. */,
    double *kirchhoff_tr_dev /**< [in/out] Deviatoric elastic stress tensor. */,
    double *pressure /**< [out] First invariant of the stress tensor */,
    double *J2 /**< [out] Second invariant of the deviatoric stress tensor */,
    const double *E_hencky_trial, /**< [in] Trial elastic strain tensor. */
    double K /**< [in] First Lamé invariant. */,
    double G /**< [in] Second Lamé invariant. */,
    double p_ref /**< [in] Reference pressure. */);

static int __update_internal_variables_elastic(
    State_Parameters IO_State /**< [in/out] List of internal variables */,
    const double
        *kirchhoff_tr_vol /**< [in] Volumetric elastic stress tensor */,
    const double
        *kirchhoff_tr_dev /**< [in] Deviatoric elastic stress tensor */);

static int __compute_plastic_flow_direction(
    double *n /**< [out] Plastic flow direction */,
    const double
        *kirchhoff_tr_dev /**< [in] Deviatoric elastic stress tensor */,
    double J2 /**< [in] Second invariant of the deviatoric stress tensor */);

static int __compute_eps(
    double *eps_k /**< [out] Equivalent plastic strain*/,
    double d_gamma_k /**< [in] Discrete plastic multiplier */,
    double eps_n /**< [in] Equivalent plastic strain in the last step */,
    double alpha_Q /**< [in] Plastic potential parameter */);

static int __compute_kappa(double *kappa_k /**< [out] Hardening function. */,
                           double kappa_0 /**< [in] Reference hardening */,
                           double exp_param /**< [in] Hardening exponential*/,
                           double eps_k /**< [in] Equivalent plastic strain*/,
                           double eps_0 /**< [in] Reference plastic strain */);

static int __compute_d_kappa(
    double *d_kappa /**< [out] Derivative of the hardening function. */,
    double kappa_0 /**< [in] Reference hardening */,
    double eps_k /**< [in] Equivalent plastic strain*/,
    double eps_0 /**< [in] Reference plastic strain */,
    double exp_param /**< [in] Hardening exponential*/);

static int __compute_pressure_limit(
    double *pressure_limit /**< [out] Limit for the apex region. */,
    double J2 /**< [in] Second invariant of the deviatoric stress tensor */,
    double kappa_n /**< [in] Hardening function. */,
    double d_kappa /**< [in] Derivative of the hardening function. */,
    double K /**< [in] First Lamé invariant. */,
    double G /**< [in] Second Lamé invariant. */,
    double alpha_F /**< [in] Yield surface parameter I. */,
    double alpha_Q /**< [in] Plastic potential parameter. */,
    double beta /**< [in] Yield surface parameter II. */);

static double __yield_function_classical(
    double pressure /**< [in] First invariant of the stress tensor */,
    double J2 /**< [in] Second invariant of the deviatoric stress tensor */,
    double d_gamma_k /**< [in] Discrete plastic multiplier */,
    double kappa_k /**< [in] Hardening function. */,
    double alpha_F /**< [in] Yield surface parameter I. */,
    double alpha_Q /**< [in] Plastic potential parameter. */,
    double beta /**< [in] Yield surface parameter II. */,
    double K /**< [in] First Lamé invariant. */,
    double G /**< [in] Second Lamé invariant. */);

static double __d_yield_function_classical(
    double d_kappa_k /**< [in] Derivative of the hardening function. */,
    double K /**< [in] First Lamé invariant. */,
    double G /**< [in] Second Lamé invariant. */,
    double alpha_F /**< [in] Yield surface parameter I. */,
    double alpha_Q /**< [in] Plastic potential parameter. */,
    double beta /**< [in] Yield surface parameter II. */);

static int __update_internal_variables_classical(
    State_Parameters IO_State /**< [in/out] List of internal variables. */,
    const double
        *kirchhoff_tr_vol /**< [in] Volumetric elastic stress tensor. */,
    const double
        *kirchhoff_tr_dev /**< [in] Deviatoric elastic stress tensor. */,
    const double *n /**< [out] Plastic flow direction. */,
    double d_gamma_k /**< [in] Discrete plastic multiplier */,
    double alpha_Q /**< [in] Plastic potential parameter. */,
    double K /**< [in] First Lamé invariant. */,
    double G /**< [in] Second Lamé invariant. */,
    double eps_k /**< [in] Equivalent plastic strain*/,
    double kappa_k /**< [in] Hardening function. */);

static double __yield_function_apex(
    double pressure /**< [in] First invariant of the stress tensor */,
    double d_gamma_k /**< [in] Discrete plastic multiplier */,
    double d_gamma_1 /**< [in] Discrete plastic multiplier I */,
    double kappa_k /**< [in] Hardening function. */,
    double d_kappa_k /**< [in] Discrete plastic multiplier */,
    double K /**< [in] First Lamé invariant. */,
    double alpha_F /**< [in] Yield surface parameter I. */,
    double alpha_Q /**< [in] Plastic potential parameter. */,
    double beta /**< [in] Yield surface parameter II. */);

static double __d_yield_function_apex(
    double d_gamma_k /**< [in] Discrete plastic multiplier */,
    double d_gamma_1 /**< [in] Discrete plastic multiplier I */,
    double d_kappa_k /**< [in] Discrete plastic multiplier */,
    double K /**< [in] First Lamé invariant. */,
    double alpha_F /**< [in] Yield surface parameter I. */,
    double alpha_Q /**< [in] Plastic potential parameter. */,
    double beta /**< [in] Yield surface parameter II. */);

static int __update_internal_variables_apex(
    State_Parameters IO_State /**< [in/out] List of internal variables. */,
    const double
        *kirchhoff_tr_vol /**< [in] Volumetric elastic stress tensor. */,
    const double *n /**< [out] Plastic flow direction. */,
    double d_gamma_k /**< [in] Discrete plastic multiplier */,
    double d_gamma_1 /**< [in] Discrete plastic multiplier I */,
    double alpha_Q /**< [in] Plastic potential parameter. */,
    double K /**< [in] First Lamé invariant. */,
    double eps_k /**< [in] Equivalent plastic strain*/,
    double kappa_k /**< [in] Hardening function. */);

void ERROR() { fprintf(stderr, "\033[1;31m"); }
void SUCCES() { fprintf(stdout, "\033[1;33m"); }
void RESET_ERROR() { fprintf(stderr, "\033[0m"); }

int Drucker_Prager_backward_euler(State_Parameters IO_State, Material MatProp);

/**************************************************************/

/**************************************************************/
/******************** Solver Parameters ***********************/
/**************************************************************/
#define TOL_Radial_Returning 10E-14
#define Max_Iterations_Radial_Returning 10
#define PI__MatrixLib__ 3.14159265358979323846

static float sqr_arg;
#define SQR(a) ((sqr_arg = (a)) == 0.0 ? 0.0 : sqr_arg * sqr_arg)
static double dsqr_arg;
#define DSQR(a) ((dsqr_arg = (a)) == 0.0 ? 0.0 : dsqr_arg * dsqr_arg)
static double dmax_arg1, dmax_arg2;
#define DMAX(a, b)                                                             \
  (dmax_arg1 = (a), dmax_arg2 = (b),                                           \
   (dmax_arg1) > (dmax_arg2) ? (dmax_arg1) : (dmax_arg2))
static double dmin_arg1, dmin_arg2;
#define DMIN(a, b)                                                             \
  (dmin_arg1 = (a), dmin_arg2 = (b),                                           \
   (dmin_arg1) < (dmin_arg2) ? (dmin_arg1) : (dmin_arg2))
static float max_arg1, max_arg2;
#define FMAX(a, b)                                                             \
  (max_arg1 = (a), max_arg2 = (b),                                             \
   (max_arg1) > (max_arg2) ? (max_arg1) : (max_arg2))
static float min_arg1, min_arg2;
#define FMIN(a, b)                                                             \
  (min_arg1 = (a), min_arg2 = (b),                                             \
   (min_arg1) < (min_arg2) ? (min_arg1) : (min_arg2))
static long lmax_arg1, lmax_arg2;
#define LMAX(a, b)                                                             \
  (lmax_arg1 = (a), lmax_arg2 = (b),                                           \
   (lmax_arg1) > (lmax_arg2) ? (lmax_arg1) : (lmax_arg2))
static long lmin_arg1, lmin_arg2;
#define LMIN(a, b)                                                             \
  (lmin_arg1 = (a), lmin_arg2 = (b),                                           \
   (lmin_arg1) < (lmin_arg2) ? (lmin_arg1) : (lmin_arg2))
static int imax_arg1, imax_arg2;
#define IMAX(a, b)                                                             \
  (imax_arg1 = (a), imax_arg2 = (b),                                           \
   (imax_arg1) > (imax_arg2) ? (imax_arg1) : (imax_arg2))
static int imin_arg1, imin_arg2;
#define IMIN(a, b)                                                             \
  (imin_arg1 = (a), imin_arg2 = (b),                                           \
   (imin_arg1) < (imin_arg2) ? (imin_arg1) : (imin_arg2))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) s

/**************************************************************/
/******************* Material Parameters **********************/
/**************************************************************/
#define NumberDimensions 2
#define YoungMouduls 100E3
#define PoissonRatio 0.3
#define param_kappa_0 10.0
#define phi 30.0
#define psi 30.0
#define hardening_modulus -10000.0
#define hardening_exp -4.0
#define NumberSteps 50 // 50 // 3500 // 2500 //
#define Delta_strain_II -0.000001
#define dF_yy 0.99999

int main() {
  int STATUS = EXIT_SUCCESS;
  State_Parameters IO_State;
  Material MatProp;

  // Define material variables
  MatProp.E = YoungMouduls;
  MatProp.nu = PoissonRatio;
  MatProp.yield_stress_0 = param_kappa_0;
  MatProp.phi_Frictional = phi;
  MatProp.psi_Frictional = psi;
  MatProp.Exponent_Hardening_Ortiz = hardening_exp;
  MatProp.Reference_Plastic_Strain_Ortiz =
      (param_kappa_0 / (hardening_exp * hardening_modulus)) *
      pow(1, (1.0 / hardening_exp - 1.0));
  MatProp.ReferencePressure = 0.0;

  // Initialize state variables
  double *stress = (double *)calloc(3 * NumberSteps, sizeof(double));
  double d_phi[5] = {1.0, 0.0, 0.0, dF_yy, 1.0};
  double *D_phi = (double *)calloc(5 * NumberSteps, sizeof(double));
  double *D_phi_e = (double *)calloc(5 * NumberSteps, sizeof(double));
  double *kappa1 = (double *)calloc(NumberSteps, sizeof(double));
  double *Equiv_Plast_Str = (double *)calloc(NumberSteps, sizeof(double));
  double *Increment_E_plastic =
      (double *)calloc(3 * NumberSteps, sizeof(double));

  // Set initial values
  double stress_tr[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; i++) {
    stress[i] = MatProp.ReferencePressure;
  }

  kappa1[0] = param_kappa_0;

  D_phi[0] = 1.0;
  D_phi[3] = 1.0;
  D_phi[4] = 1.0;

  D_phi_e[0] = 1.0;
  D_phi_e[3] = 1.0;
  D_phi_e[4] = 1.0;

  // Start time integration
  for (int i = 1; i < NumberSteps; i++) {

    printf("Step: %i \n", i);

    // Asign variables to the solver
    IO_State.Particle_Idx = 0;
    IO_State.Stress = &stress[i * 3];
    IO_State.D_phi_e = &D_phi_e[i * 5];
    IO_State.Increment_E_plastic = &Increment_E_plastic[i * 3];
    IO_State.Equiv_Plast_Str = &Equiv_Plast_Str[i];
    IO_State.Kappa = &kappa1[i];
    IO_State.d_phi = d_phi;

    // Initialize
    IO_State.D_phi_e[0] = D_phi_e[(i - 1) * 5 + 0];
    IO_State.D_phi_e[1] = D_phi_e[(i - 1) * 5 + 1];
    IO_State.D_phi_e[2] = D_phi_e[(i - 1) * 5 + 2];
    IO_State.D_phi_e[3] = D_phi_e[(i - 1) * 5 + 3];
    IO_State.D_phi_e[4] = D_phi_e[(i - 1) * 5 + 4];

    D_phi[i * 5 + 3] = D_phi[(i - 1) * 5 + 3] * d_phi[3];

    *IO_State.Equiv_Plast_Str = Equiv_Plast_Str[i - 1];
    *IO_State.Kappa = kappa1[i - 1];

    // Run solver
    STATUS = Drucker_Prager_backward_euler(IO_State, MatProp);

    if (STATUS) {
      ERROR();
      fprintf(stderr, "Failure\n");
      RESET_ERROR();
      free(stress);
      free(D_phi);
      free(D_phi_e);
      free(kappa1);
      free(Equiv_Plast_Str);
      return STATUS;
    }
  }

  // Save data in a csv file
  FILE *DP_output = fopen("DP_output.csv", "w");
  for (int i = 0; i < NumberSteps; i++) {
    fprintf(DP_output, "%e, %e\n",
            -0.5 * (D_phi[i * 5 + 3] * D_phi[i * 5 + 3] - 1.0),
            sqrt(0.5 * (pow((stress[i * 3 + 0] - stress[i * 3 + 1]), 2.0) +
                        pow((stress[i * 3 + 1] - stress[i * 3 + 2]), 2.0) +
                        pow((stress[i * 3 + 2] - stress[i * 3 + 0]), 2.0))));
  }
  fclose(DP_output);

  // Print data with gnuplot
  FILE *gnuplot = popen("gnuplot -persistent", "w");
  fprintf(gnuplot, "set termoption enhanced \n");
  fprintf(gnuplot, "set datafile separator ',' \n");
  fprintf(gnuplot,
          "set xlabel '- {/Symbol e}_{II}' font 'Times,20' enhanced \n");
  fprintf(gnuplot, "set ylabel '{/Symbol s}_{VM}' font "
                   "'Times,20' enhanced \n");
  fprintf(gnuplot, "plot %s \n",
          "'DP_output.csv' title 'Us' with line lt rgb 'red' lw 2");
  fflush(gnuplot);

  // Free memory
  free(stress);
  free(D_phi);
  free(D_phi_e);
  free(kappa1);
  free(Equiv_Plast_Str);

  return STATUS;
}

/**************************************************************/

int Drucker_Prager_backward_euler(State_Parameters IO_State, Material MatProp)
/*
  Backward Euler algorithm for the Drucker-Prager (Lorenzo Sanavia)
*/
{
  int STATUS = EXIT_SUCCESS;

  // Read input/output parameters
  double E_hencky_trial[3] = {0.0, 0.0, 0.0};
  double EE[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double EE_m1[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double D_phi_e_tr[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  double kirchhoff_tr_vol[3] = {0.0, 0.0, 0.0};
  double kirchhoff_tr_dev[3] = {0.0, 0.0, 0.0};

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  puts("t = n elastic deformation gradient");
  printf("%f %f %f \n", IO_State.D_phi_e[0], IO_State.D_phi_e[1], 0.0);
  printf("%f %f %f \n", IO_State.D_phi_e[2], IO_State.D_phi_e[3], 0.0);
  printf("%f %f %f \n", 0.0, 0.0, IO_State.D_phi_e[4]);

#endif
#endif

  STATUS = __compute_trial_D_phi_e(E_hencky_trial, EE, EE_m1, D_phi_e_tr,
                                   IO_State.d_phi, IO_State.D_phi_e);
  if (STATUS) {
    ERROR();
    fprintf(stderr, "__compute_trial_D_phi_e\n");
    RESET_ERROR();
    return STATUS;
  }

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  puts("trial elastic deformation gradient");
  printf("%f %f %f \n", D_phi_e_tr[0], D_phi_e_tr[1], D_phi_e_tr[2]);
  printf("%f %f %f \n", D_phi_e_tr[3], D_phi_e_tr[4], D_phi_e_tr[5]);
  printf("%f %f %f \n", D_phi_e_tr[6], D_phi_e_tr[7], D_phi_e_tr[8]);

  printf("E_hencky_trial: [%f, %f, %f] \n", E_hencky_trial[0],
         E_hencky_trial[1], E_hencky_trial[2]);

  puts("Principal directions");
  printf("%f %f %f \n", EE[0], EE[1], EE[2]);
  printf("%f %f %f \n", EE[3], EE[4], EE[5]);
  printf("%f %f %f \n", EE[6], EE[7], EE[8]);

#endif
#endif

  // Material parameters
  double K = MatProp.E / (3.0 * (1.0 - 2.0 * MatProp.nu));
  double G = MatProp.E / (2.0 * (1.0 + MatProp.nu));
  double p_ref = MatProp.ReferencePressure;
  double rad_friction_angle =
      (PI__MatrixLib__ / 180.0) * MatProp.phi_Frictional;
  double rad_dilatancy_angle =
      (PI__MatrixLib__ / 180.0) * MatProp.psi_Frictional;
  double exp_param = MatProp.Exponent_Hardening_Ortiz;
  double kappa_0 = MatProp.yield_stress_0;
  double eps_0 = MatProp.Reference_Plastic_Strain_Ortiz;

#if NumberDimensions == 2
  double alpha_F = sqrt(2. / 3.) * tan(rad_friction_angle) /
                   sqrt(3. + 4. * DSQR(tan(rad_friction_angle)));
  double alpha_Q = sqrt(2. / 3.) * tan(rad_dilatancy_angle) /
                   sqrt(3. + 4. * DSQR(tan(rad_dilatancy_angle)));
  double beta =
      sqrt(2. / 3.) * 3. / sqrt(3. + 4. * DSQR(tan(rad_friction_angle)));
#else
  double alpha_F = sqrt(2 / 3.) * 2 * sin(rad_friction_angle) /
                   (3 - sin(rad_friction_angle));
  double alpha_Q = sqrt(2 / 3.) * 2 * sin(rad_dilatancy_angle) /
                   (3 - sin(rad_dilatancy_angle));
  double beta = sqrt(2 / 3.) * 6 * cos(rad_friction_angle) /
                (3 - sin(rad_friction_angle));
#endif

  // Define tensorial internal variables
  double n[3] = {0.0, 0.0, 0.0};

  // Define scalar internal variables
  double PHI, PHI_0 = 0.0;
  double d_PHI = 0.0;
  double J2 = 0.0;
  double pressure = 0.0;
  double pressure_limit = 0.0;
  double d_gamma_k = 0;
  double eps_n = *IO_State.Equiv_Plast_Str;
  double eps_k = eps_n;
  double kappa_k = *IO_State.Kappa;
  double d_kappa_k = 0.0;

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  printf("Equivalent Plastic Strain: %f \n", *IO_State.Equiv_Plast_Str);
  printf("Kappa: %f \n", *IO_State.Kappa);
#endif
#endif

  // Initialise solver parameters
  double TOL = TOL_Radial_Returning;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;
  bool Convergence = false;

  STATUS = __trial_elastic(kirchhoff_tr_vol, kirchhoff_tr_dev, &pressure, &J2,
                           E_hencky_trial, K, G, p_ref);
  if (STATUS) {
    ERROR();
    fprintf(stderr, "__trial_elastic\n");
    RESET_ERROR();
    return STATUS;
  }

  PHI_0 = __yield_function_classical(pressure, J2, d_gamma_k, kappa_k, alpha_F,
                                     alpha_Q, beta, K, G);

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  printf("Initial value of the yield function: %f \n", PHI_0);
#endif
#endif

  // Elastic
  if (PHI_0 <= 0.0) {

    STATUS = __update_internal_variables_elastic(IO_State, kirchhoff_tr_vol,
                                                 kirchhoff_tr_dev);
    if (STATUS) {
      ERROR();
      fprintf(stderr, "__update_internal_variables_elastic\n");
      RESET_ERROR();
      return STATUS;
    }

  }
  // Plastic (check yield condition)
  else {

    PHI = PHI_0;

    STATUS = __compute_plastic_flow_direction(n, kirchhoff_tr_dev, J2);
    if (STATUS) {
      ERROR();
      fprintf(stderr, "__compute_plastic_flow_direction\n");
      RESET_ERROR();
      return STATUS;
    }

    STATUS = __compute_d_kappa(&d_kappa_k, kappa_0, eps_n, eps_0, exp_param);
    if (STATUS) {
      ERROR();
      fprintf(stderr, "__compute_d_kappa\n");
      RESET_ERROR();
      return STATUS;
    }

    STATUS = __compute_pressure_limit(&pressure_limit, J2, kappa_k, d_kappa_k,
                                      K, G, alpha_F, alpha_Q, beta);
    if (STATUS) {
      ERROR();
      fprintf(stderr, "__compute_pressure_limit\n");
      RESET_ERROR();
      return STATUS;
    }

    // Classic radial returning algorithm
    if (-pressure < pressure_limit) {

      while (fabs(PHI / PHI_0) >= TOL) {

        Iter++;

        if (Iter == MaxIter)
          break;

        d_PHI = __d_yield_function_classical(d_kappa_k, K, G, alpha_F, alpha_Q,
                                             beta);

        d_gamma_k += -PHI / d_PHI;

        STATUS = __compute_eps(&eps_k, d_gamma_k, eps_n, alpha_Q);
        if (STATUS) {
          ERROR();
          fprintf(stderr, "__compute_eps\n");
          RESET_ERROR();
          return STATUS;
        }

        STATUS = __compute_kappa(&kappa_k, kappa_0, exp_param, eps_k, eps_0);
        if (STATUS) {
          ERROR();
          fprintf(stderr, "__compute_kappa\n");
          RESET_ERROR();
          return STATUS;
        }

        if (kappa_k > 0) {
          STATUS =
              __compute_d_kappa(&d_kappa_k, kappa_0, eps_k, eps_0, exp_param);
          if (STATUS) {
            ERROR();
            fprintf(stderr, "__compute_d_kappa\n");
            RESET_ERROR();
            return STATUS;
          }

        } else {
          kappa_k = 0.0;
          d_kappa_k = 0.0;
        }

        PHI = __yield_function_classical(pressure, J2, d_gamma_k, kappa_k,
                                         alpha_F, alpha_Q, beta, K, G);
      }

      STATUS = __update_internal_variables_classical(
          IO_State, kirchhoff_tr_vol, kirchhoff_tr_dev, n, d_gamma_k, alpha_Q,
          K, G, eps_k, kappa_k);
      if (STATUS) {
        ERROR();
        fprintf(stderr, "__update_internal_variables_classical\n");
        RESET_ERROR();
        return STATUS;
      }

    }
    // Apex radial returning algorithm
    else {

      double d_gamma_1 = J2 / (2.0 * G);
      double d_gamma_2_k = 0.0;
      d_gamma_k = d_gamma_1 + d_gamma_2_k;

      while (fabs(PHI / PHI_0) >= TOL) {

        Iter++;

        if (Iter == MaxIter)
          break;

        d_PHI = __d_yield_function_apex(d_gamma_k, d_gamma_1, d_kappa_k, K,
                                        alpha_F, alpha_Q, beta);

        d_gamma_2_k += -PHI / d_PHI;

        d_gamma_k = d_gamma_1 + d_gamma_2_k;

        PHI = __yield_function_apex(pressure, d_gamma_k, d_gamma_1, kappa_k,
                                    d_kappa_k, K, alpha_F, alpha_Q, beta);
      }

      STATUS = __compute_eps(&eps_k, d_gamma_k, eps_n, alpha_Q);
      if (STATUS) {
        ERROR();
        fprintf(stderr, "__compute_eps\n");
        RESET_ERROR();
        return STATUS;
      }

      STATUS = __update_internal_variables_apex(IO_State, kirchhoff_tr_vol, n,
                                                d_gamma_k, d_gamma_1, alpha_Q,
                                                K, eps_k, kappa_k);
      if (STATUS) {
        ERROR();
        fprintf(stderr, "__update_internal_variables_apex\n");
        RESET_ERROR();
        return STATUS;
      }
    }
  }

  STATUS = __compute_D_phi_e(IO_State.D_phi_e, D_phi_e_tr, EE, EE_m1,
                             IO_State.Increment_E_plastic);
  if (STATUS) {
    ERROR();
    fprintf(stderr, "__compute_D_phi_e\n");
    RESET_ERROR();
    return STATUS;
  }

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  printf("Numer of iterations: %i, error: %e \n", Iter, PHI);
  if (Iter == MaxIter) {
    ERROR();
    printf("Number of iterations: %i \n", Iter);
    printf("Final value of the yield function: %e \n", PHI);
    RESET_ERROR();
  }

  printf("Out stress tensor: [%f, %f, %f] \n", IO_State.Stress[0],
         IO_State.Stress[1], IO_State.Stress[2]);
  printf("Increment of the plastic tensor: [%e, %e, %e] \n",
         IO_State.Increment_E_plastic[0], IO_State.Increment_E_plastic[1],
         IO_State.Increment_E_plastic[2]);
  puts("Out elastic deformation gradient: ");
  printf("%e, %e, %e \n", IO_State.D_phi_e[0], IO_State.D_phi_e[1], 0.0);
  printf("%e, %e, %e \n", IO_State.D_phi_e[2], IO_State.D_phi_e[3], 0.0);
  printf("%e, %e, %e \n", 0.0, 0.0, IO_State.D_phi_e[4]);
#endif
#endif

  return STATUS;
}

/**************************************************************/

static int __compute_trial_D_phi_e(double *e_hencky, double *EE, double *EE_m1,
                                   double *D_phi_e_tr, const double *d_phi,
                                   const double *D_phi_e) {

  double C_e_tr[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

#if NumberDimensions == 2

  D_phi_e_tr[0] = d_phi[0] * D_phi_e[0] + d_phi[1] * D_phi_e[2];
  D_phi_e_tr[1] = d_phi[0] * D_phi_e[1] + d_phi[1] * D_phi_e[3];
  D_phi_e_tr[3] = d_phi[2] * D_phi_e[0] + d_phi[3] * D_phi_e[2];
  D_phi_e_tr[4] = d_phi[2] * D_phi_e[1] + d_phi[3] * D_phi_e[3];
  D_phi_e_tr[8] = d_phi[4] * D_phi_e[4];

  C_e_tr[0] = D_phi_e_tr[0] * D_phi_e_tr[0] + D_phi_e_tr[3] * D_phi_e_tr[3];
  C_e_tr[1] = D_phi_e_tr[0] * D_phi_e_tr[1] + D_phi_e_tr[3] * D_phi_e_tr[4];
  C_e_tr[3] = D_phi_e_tr[1] * D_phi_e_tr[0] + D_phi_e_tr[4] * D_phi_e_tr[3];
  C_e_tr[4] = D_phi_e_tr[1] * D_phi_e_tr[1] + D_phi_e_tr[4] * D_phi_e_tr[4];
  C_e_tr[8] = D_phi_e_tr[8] * D_phi_e_tr[8];

#else
  No esta implementado
#endif

  /* Locals */
  int n = 3;
  int lda = 3;
  int ldvl = 3;
  int ldvr = 3;
  int info;
  int lwork;
  double wkopt;
  double *work;

  /* Local arrays */
  double wr[3];
  double wi[3];
  double vl[9];
  double vr[9];

  /*
    Query and allocate the optimal workspace
  */
  lwork = -1;
  dgeev_("V", "V", &n, C_e_tr, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt,
         &lwork, &info);
  lwork = (int)wkopt;
  work = (double *)malloc(lwork * sizeof(double));

  /* Check for convergence */
  if (info > 0) {
    free(work);
    printf("Error in Eigen_analysis__TensorLib__() : The algorithm failed to "
           "compute eigenvalues.\n");
    return EXIT_FAILURE;
  }

  /* Solve eigenproblem */
  dgeev_("V", "V", &n, C_e_tr, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork,
         &info);

  /* Check for convergence */
  if (info > 0) {
    free(work);
    printf("Error in Eigen_analysis__TensorLib__() : The algorithm failed to "
           "compute eigenvalues.\n");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    free(work);
    printf("Error in Eigen_analysis__TensorLib__() : the %i-th argument had an "
           "illegal value.\n",
           abs(info));
    return EXIT_FAILURE;
  }

  free(work);

  e_hencky[0] = 0.5 * log(wr[0]);
  e_hencky[1] = 0.5 * log(wr[1]);
  e_hencky[2] = 0.5 * log(wr[2]);

  EE[0] = vr[0];
  EE[1] = vr[3];
  EE[2] = vr[6];
  EE[3] = vr[1];
  EE[4] = vr[4];
  EE[5] = vr[7];
  EE[6] = vr[2];
  EE[7] = vr[5];
  EE[8] = vr[8];

  EE_m1[0] = vl[0];
  EE_m1[1] = vl[1];
  EE_m1[2] = vl[2];
  EE_m1[3] = vl[3];
  EE_m1[4] = vl[4];
  EE_m1[5] = vl[5];
  EE_m1[6] = vl[6];
  EE_m1[7] = vl[7];
  EE_m1[8] = vl[8];

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __compute_D_phi_e(double *D_phi_e, const double *D_phi_e_tr,
                             const double *EE, const double *EE_m1,
                             const double *Increment_E_plastic) {

  double exp_Increment_E_plastic[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                       0.0, 0.0, 0.0, 0.0};
  double d_phi_e_corr[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  exp_Increment_E_plastic[0] = exp(-Increment_E_plastic[0]);
  exp_Increment_E_plastic[4] = exp(-Increment_E_plastic[1]);
  exp_Increment_E_plastic[8] = exp(-Increment_E_plastic[2]);

#if NumberDimensions == 2
  d_phi_e_corr[0] = EE[0] * exp_Increment_E_plastic[0] * EE_m1[0] +
                    EE[1] * exp_Increment_E_plastic[4] * EE_m1[3];
  d_phi_e_corr[1] = EE[0] * exp_Increment_E_plastic[0] * EE_m1[1] +
                    EE[1] * exp_Increment_E_plastic[4] * EE_m1[4];

  d_phi_e_corr[3] = EE[3] * exp_Increment_E_plastic[0] * EE_m1[0] +
                    EE[4] * exp_Increment_E_plastic[4] * EE_m1[3];
  d_phi_e_corr[4] = EE[3] * exp_Increment_E_plastic[0] * EE_m1[1] +
                    EE[4] * exp_Increment_E_plastic[4] * EE_m1[4];

  d_phi_e_corr[8] = exp_Increment_E_plastic[8];

  D_phi_e[0] = D_phi_e_tr[0];
  D_phi_e[1] = D_phi_e_tr[1];
  D_phi_e[2] = D_phi_e_tr[3];
  D_phi_e[3] = D_phi_e_tr[4];
  D_phi_e[4] = D_phi_e_tr[8];
  /*
    D_phi_e[0] =
        D_phi_e_tr[0] * d_phi_e_corr[0] + D_phi_e_tr[1] * d_phi_e_corr[3];

    D_phi_e[1] =
        D_phi_e_tr[0] * d_phi_e_corr[1] + D_phi_e_tr[1] * d_phi_e_corr[4];

    D_phi_e[2] =
        D_phi_e_tr[3] * d_phi_e_corr[0] + D_phi_e_tr[4] * d_phi_e_corr[3];

    D_phi_e[3] =
        D_phi_e_tr[3] * d_phi_e_corr[1] + D_phi_e_tr[4] * d_phi_e_corr[4];

    D_phi_e[4] = D_phi_e_tr[8] * d_phi_e_corr[8];
  */

#else
  No esta implementado
#endif

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  puts("Corrector elastic deformation gradient");
  printf("%f %f %f \n", d_phi_e_corr[0], d_phi_e_corr[1], d_phi_e_corr[2]);
  printf("%f %f %f \n", d_phi_e_corr[3], d_phi_e_corr[4], d_phi_e_corr[5]);
  printf("%f %f %f \n", d_phi_e_corr[6], d_phi_e_corr[7], d_phi_e_corr[8]);
#endif
#endif

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __trial_elastic(double *kirchhoff_tr_vol, double *kirchhoff_tr_dev,
                           double *pressure, double *J2,
                           const double *E_hencky_trial, double K, double G,
                           double p_ref) {

  double E_hencky_trial_vol[3];
  double tr_E_hencky_trial =
      E_hencky_trial[0] + E_hencky_trial[1] + E_hencky_trial[2];

  E_hencky_trial_vol[0] = (1.0 / 3.0) * tr_E_hencky_trial;
  E_hencky_trial_vol[1] = (1.0 / 3.0) * tr_E_hencky_trial;
  E_hencky_trial_vol[2] = (1.0 / 3.0) * tr_E_hencky_trial;

  kirchhoff_tr_vol[0] = -K * E_hencky_trial_vol[0] - p_ref;
  kirchhoff_tr_vol[1] = -K * E_hencky_trial_vol[1] - p_ref;
  kirchhoff_tr_vol[2] = -K * E_hencky_trial_vol[2] - p_ref;

  kirchhoff_tr_dev[0] = 2 * G * (E_hencky_trial[0] - E_hencky_trial_vol[0]);
  kirchhoff_tr_dev[1] = 2 * G * (E_hencky_trial[1] - E_hencky_trial_vol[1]);
  kirchhoff_tr_dev[2] = 2 * G * (E_hencky_trial[2] - E_hencky_trial_vol[2]);

  *pressure =
      (kirchhoff_tr_vol[0] + kirchhoff_tr_vol[1] + kirchhoff_tr_vol[2]) / 3.0;
  *J2 = sqrt(kirchhoff_tr_dev[0] * kirchhoff_tr_dev[0] +
             kirchhoff_tr_dev[1] * kirchhoff_tr_dev[1] +
             kirchhoff_tr_dev[2] * kirchhoff_tr_dev[2]);

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __update_internal_variables_elastic(State_Parameters IO_State,
                                               const double *kirchhoff_tr_vol,
                                               const double *kirchhoff_tr_dev) {

  IO_State.Increment_E_plastic[0] = 0.0;
  IO_State.Increment_E_plastic[1] = 0.0;
  IO_State.Increment_E_plastic[2] = 0.0;

  IO_State.Stress[0] = -kirchhoff_tr_vol[0] + kirchhoff_tr_dev[0];
  IO_State.Stress[1] = -kirchhoff_tr_vol[1] + kirchhoff_tr_dev[1];
  IO_State.Stress[2] = -kirchhoff_tr_vol[2] + kirchhoff_tr_dev[2];

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __compute_plastic_flow_direction(double *n,
                                            const double *kirchhoff_tr_dev,
                                            double J2) {

  n[0] = kirchhoff_tr_dev[0] / J2;
  n[1] = kirchhoff_tr_dev[1] / J2;
  n[2] = kirchhoff_tr_dev[2] / J2;

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __compute_eps(double *eps_k, double d_gamma_k, double eps_n,
                         double alpha_Q) {

  *eps_k = eps_n + d_gamma_k * sqrt(3 * alpha_Q * alpha_Q + 1);

  if (*eps_k < 0.0)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __compute_kappa(double *kappa_k, double kappa_0, double exp_param,
                           double eps_k, double eps_0) {
  double base = 1.0 + eps_k / eps_0;
  double exp = 1.0 / exp_param;

  if (base < 0.0)
    return EXIT_FAILURE;

  *kappa_k = kappa_0 * pow(base, exp);

  if (*kappa_k < 0.0)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __compute_d_kappa(double *d_kappa, double kappa_0, double eps_k,
                             double eps_0, double exp_param) {

  double base = 1.0 + eps_k / eps_0;
  double exp = 1.0 / exp_param - 1.0;

  if (base < 0.0)
    return EXIT_FAILURE;

  *d_kappa = (kappa_0 / (exp_param * eps_0)) * pow(base, exp);

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __compute_pressure_limit(double *pressure_limit, double J2,
                                    double kappa_n, double d_kappa, double K,
                                    double G, double alpha_F, double alpha_Q,
                                    double beta) {

  double ads = sqrt(1.0 + 3.0 * alpha_Q * alpha_Q);

  if (alpha_F == 0.0)
    return EXIT_FAILURE;

  *pressure_limit =
      3.0 * alpha_Q * K / (2.0 * G) * J2 +
      beta / (3.0 * alpha_F) * ((J2 / (2.0 * G)) * d_kappa * ads + kappa_n);

  return EXIT_SUCCESS;
}

/**************************************************************/

static double __yield_function_classical(double pressure, double J2,
                                         double d_gamma_k, double kappa_k,
                                         double alpha_F, double alpha_Q,
                                         double beta, double K, double G) {

  double PHI = (J2 - 2.0 * G * d_gamma_k -
                3.0 * alpha_F * (pressure - 3.0 * K * alpha_Q * d_gamma_k) -
                beta * kappa_k);

  return PHI;
}

/***************************************************************************/

static double __d_yield_function_classical(double d_kappa_k, double K, double G,
                                           double alpha_F, double alpha_Q,
                                           double beta) {

  double ads = sqrt(1.0 + 3.0 * alpha_Q * alpha_Q);

  double d_PHI =
      +9.0 * K * alpha_F * alpha_Q - 2.0 * G - beta * d_kappa_k * ads;

  return d_PHI;
}

/***************************************************************************/

static int __update_internal_variables_classical(
    State_Parameters IO_State, const double *kirchhoff_tr_vol,
    const double *kirchhoff_tr_dev, const double *n, double d_gamma_k,
    double alpha_Q, double K, double G, double eps_k, double kappa_k) {

  IO_State.Increment_E_plastic[0] = d_gamma_k * (alpha_Q + n[0]);
  IO_State.Increment_E_plastic[1] = d_gamma_k * (alpha_Q + n[1]);
  IO_State.Increment_E_plastic[2] = d_gamma_k * (alpha_Q + n[2]);

  IO_State.Stress[0] = -kirchhoff_tr_vol[0] + kirchhoff_tr_dev[0] -
                       d_gamma_k * (3 * K * alpha_Q + 2 * G * n[0]);
  IO_State.Stress[1] = -kirchhoff_tr_vol[1] + kirchhoff_tr_dev[1] -
                       d_gamma_k * (3 * K * alpha_Q + 2 * G * n[1]);
  IO_State.Stress[2] = -kirchhoff_tr_vol[2] + kirchhoff_tr_dev[2] -
                       d_gamma_k * (3 * K * alpha_Q + 2 * G * n[2]);

  *IO_State.Equiv_Plast_Str = eps_k;
  *IO_State.Kappa = kappa_k;

  return EXIT_SUCCESS;
}

/***************************************************************************/

static double __yield_function_apex(double pressure, double d_gamma_k,
                                    double d_gamma_1, double kappa_k,
                                    double d_kappa_k, double K, double alpha_F,
                                    double alpha_Q, double beta) {

  double PHI = (beta / (3.0 * alpha_F) *
                    (kappa_k + d_kappa_k * sqrt((d_gamma_1 * d_gamma_1) +
                                                3.0 * (alpha_Q * alpha_Q) *
                                                    (d_gamma_k * d_gamma_k))) -
                pressure + 3.0 * K * alpha_Q * d_gamma_k);

  return PHI;
}

/***************************************************************************/

static double __d_yield_function_apex(double d_gamma_k, double d_gamma_1,
                                      double d_kappa_k, double K,
                                      double alpha_F, double alpha_Q,
                                      double beta) {

  double d_PHI =
      3.0 * alpha_Q * K +
      3.0 * d_kappa_k * beta * (alpha_Q * alpha_Q) * d_gamma_k /
          (3.0 * alpha_F *
           sqrt((d_gamma_1 * d_gamma_1) +
                3.0 * (alpha_Q * alpha_Q) * (d_gamma_k * d_gamma_k)));

  return d_PHI;
}

/***************************************************************************/

static int __update_internal_variables_apex(State_Parameters IO_State,
                                            const double *kirchhoff_tr_vol,
                                            const double *n, double d_gamma_k,
                                            double d_gamma_1, double alpha_Q,
                                            double K, double eps_k,
                                            double kappa_k) {

  IO_State.Increment_E_plastic[0] = d_gamma_k * alpha_Q + d_gamma_1 * n[0];
  IO_State.Increment_E_plastic[1] = d_gamma_k * alpha_Q + d_gamma_1 * n[1];
  IO_State.Increment_E_plastic[2] = d_gamma_k * alpha_Q + d_gamma_1 * n[2];

  IO_State.Stress[0] = -kirchhoff_tr_vol[0] - d_gamma_k * 3 * K * alpha_Q;
  IO_State.Stress[1] = -kirchhoff_tr_vol[1] - d_gamma_k * 3 * K * alpha_Q;
  IO_State.Stress[2] = -kirchhoff_tr_vol[2] - d_gamma_k * 3 * K * alpha_Q;

  *IO_State.Equiv_Plast_Str = eps_k;
  *IO_State.Kappa = kappa_k;

  return EXIT_SUCCESS;
}

/***************************************************************************/
