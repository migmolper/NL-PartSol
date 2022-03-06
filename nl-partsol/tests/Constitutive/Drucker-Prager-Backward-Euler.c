#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

#define RESET "\033[0m"
#define RED "\033[31m" /* Red */

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
  double J2_degradated;

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
  double *D_phi;
  double *Fbar;
  double *rate_D_phi;
  double J;

  /*!
   * Plasticity parameters
   * */
  double *Back_stress;
  double *b_e;

  double Cohesion;
  double Yield_stress;

  // Internal hardening variables
  double *Equiv_Plast_Str; // Equivalent plastic strain
  double *Kappa;           // Hardening Parameter

  //
  bool *Failure;

} State_Parameters;

static int __compute_trial_b_e(
    double *eigval_b_e_tr /**< [out] Eigenvalues of b elastic trial. */,
    double *eigvec_b_e_tr /**< [out] Eigenvector of b elastic trial. */,
    const double *b_e /**< [in] (n) Elastic left Cauchy-Green.*/,
    const double *d_phi /**< [in] Incremental deformation gradient. */);

static int __corrector_b_e(
    double *b_e /**< [out] (n+1) Elastic deformation gradient. */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */,
    const double *E_hencky_trial /**< [in] Corrected Henky strain */);

static int __trial_elastic(
    double *T_tr_vol /**< [in/out] Volumetric elastic stress tensor. */,
    double *T_tr_dev /**< [in/out] Deviatoric elastic stress tensor. */,
    double *pressure /**< [out] First invariant of the stress tensor */,
    double *J2 /**< [out] Second invariant of the deviatoric stress tensor */,
    const double *E_hencky_trial, /**< [in] Trial elastic strain tensor. */
    double K /**< [in] First Lamé invariant. */,
    double G /**< [in] Second Lamé invariant. */,
    double p_ref /**< [in] Reference pressure. */);

static int __update_internal_variables_elastic(
    double *Stress /**< [in/out] Nominal stress tensor */,
    const double *D_phi /**< [in] Total deformation gradient. */,
    const double *T_tr_vol /**< [in] Volumetric elastic stress tensor */,
    const double *T_tr_dev /**< [in] Deviatoric elastic stress tensor */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */);

static int __compute_plastic_flow_direction(
    double *n /**< [out] Plastic flow direction */,
    const double *T_tr_dev /**< [in] Deviatoric elastic stress tensor */,
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
    double *Increment_E_plastic /**< [in/out] Increment plastic strain */,
    double *Stress /**< [in/out] Nominal stress tensor */,
    double *eps_n1 /**< [in/out] Equivalent plastic strain*/,
    double *kappa_n1 /**< [in/out] Hardening function. */,
    const double *D_phi /**< [in] Total deformation gradient. */,
    const double *T_tr_vol /**< [in] Volumetric elastic stress tensor. */,
    const double *T_tr_dev /**< [in] Deviatoric elastic stress tensor. */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */,
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
    double *Increment_E_plastic /**< [in/out] Increment plastic strain */,
    double *Stress /**< [in/out] Nominal stress tensor */,
    double *eps_n1 /**< [in/out] Equivalent plastic strain*/,
    double *kappa_n1 /**< [in/out] Hardening function. */,
    const double *D_phi /**< [in] Total deformation gradient. */,
    const double *T_tr_vol /**< [in] Volumetric elastic stress tensor. */,
    const double *T_tr_dev /**< [in] Deviatoric elastic stress tensor. */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */,
    const double *n /**< [out] Plastic flow direction. */,
    double d_gamma_k /**< [in] Discrete plastic multiplier */,
    double d_gamma_1 /**< [in] Discrete plastic multiplier I */,
    double alpha_Q /**< [in] Plastic potential parameter. */,
    double K /**< [in] First Lamé invariant. */,
    double G /**< [in] Second Lamé invariant. */,
    double eps_k /**< [in] Equivalent plastic strain*/,
    double kappa_k /**< [in] Hardening function. */);

static bool __degradated(
    const double *Stress /**< [in/out] Nominal stress tensor */,
    double J2_degradated /**< [in] Critical value of the J2 invariant */);

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
#define YoungMouduls 10.0E3
#define PoissonRatio 0.2
#define kappa0 40.0
#define phi 39.0
#define psi -30.0
#define H -200.0 //-300.0
#define m -0.2
#define NumberSteps 30000 // 50 // 3500 // 2500 //
#define dF_yy 0.99999
#define J2_degradated_value 5.0
#define Confining_pressure -20

int main() {
  int STATUS = EXIT_SUCCESS;
  State_Parameters IO_State;
  Material MatProp;

  // Define material variables
  MatProp.E = YoungMouduls;
  MatProp.nu = PoissonRatio;
  MatProp.yield_stress_0 = kappa0;
  MatProp.phi_Frictional = phi;
  MatProp.psi_Frictional = psi;
  MatProp.Hardening_modulus = H;
  MatProp.Exponent_Hardening_Ortiz = m;
  MatProp.Reference_Plastic_Strain_Ortiz =
      (kappa0 / (m * H)) * pow(1.0, (1.0 / m - 1.0));
  MatProp.ReferencePressure = Confining_pressure;
  MatProp.J2_degradated = J2_degradated_value;

  // Initialize state variables
  bool Status_particle = false;
  double *stress = (double *)calloc(5 * NumberSteps, sizeof(double));
  double d_phi[5] = {1.0, 0.0, 0.0, dF_yy, 1.0};
  double *D_phi = (double *)calloc(5 * NumberSteps, sizeof(double));
  double *b_e = (double *)calloc(5 * NumberSteps, sizeof(double));
  double *kappa1 = (double *)calloc(NumberSteps, sizeof(double));
  double *Equiv_Plast_Str = (double *)calloc(NumberSteps, sizeof(double));

  // Set initial values
  stress[0] = MatProp.ReferencePressure;
  stress[3] = MatProp.ReferencePressure;
  stress[4] = MatProp.ReferencePressure;

  kappa1[0] = kappa0;

  D_phi[0] = 1.0;
  D_phi[3] = 1.0;
  D_phi[4] = 1.0;

  b_e[0] = 1.0;
  b_e[3] = 1.0;
  b_e[4] = 1.0;

  // Start time integration
  for (int i = 1; i < NumberSteps; i++) {

    printf("Step: %i \n", i);

    // Asign variables to the solver
    IO_State.Particle_Idx = 0;
    IO_State.Stress = &stress[i * 5];
    IO_State.b_e = &b_e[i * 5];
    IO_State.Equiv_Plast_Str = &Equiv_Plast_Str[i];
    IO_State.Kappa = &kappa1[i];
    IO_State.d_phi = d_phi;
    IO_State.D_phi = &D_phi[i * 5];
    IO_State.Failure = &Status_particle;

    IO_State.b_e[0] = b_e[(i - 1) * 5 + 0];
    IO_State.b_e[1] = b_e[(i - 1) * 5 + 1];
    IO_State.b_e[2] = b_e[(i - 1) * 5 + 2];
    IO_State.b_e[3] = b_e[(i - 1) * 5 + 3];
    IO_State.b_e[4] = b_e[(i - 1) * 5 + 4];

    IO_State.D_phi[0] = D_phi[(i - 1) * 5 + 0] * d_phi[0];
    IO_State.D_phi[3] = D_phi[(i - 1) * 5 + 3] * d_phi[3];
    IO_State.D_phi[4] = D_phi[(i - 1) * 5 + 4] * d_phi[4];

    *IO_State.Equiv_Plast_Str = Equiv_Plast_Str[i - 1];
    *IO_State.Kappa = kappa1[i - 1];

    // Run solver
    STATUS = Drucker_Prager_backward_euler(IO_State, MatProp);

    if (STATUS) {
      fprintf(stderr, "Failure\n");
      free(stress);
      free(D_phi);
      free(b_e);
      free(kappa1);
      free(Equiv_Plast_Str);
      return STATUS;
    }
  }

  // Save data in a csv file
  FILE *E2_vs_VM = fopen("E2_vs_VM.csv", "w");
  for (int i = 0; i < NumberSteps; i++) {
    fprintf(E2_vs_VM, "%e, %e\n",
            -100 * 0.5 * (D_phi[i * 5 + 3] * D_phi[i * 5 + 3] - 1.0),
            sqrt(0.5 * (pow((stress[i * 5 + 0] - stress[i * 5 + 3]), 2.0) +
                        pow((stress[i * 5 + 3] - stress[i * 5 + 4]), 2.0) +
                        pow((stress[i * 5 + 4] - stress[i * 5 + 0]), 2.0))));
  }
  fclose(E2_vs_VM);

  FILE *p_vs_q = fopen("p_vs_q.csv.csv", "w");
  for (int i = 0; i < NumberSteps; i++) {
    fprintf(p_vs_q, "%e, %e\n", -0.5 * (stress[i * 5 + 0] + stress[i * 5 + 3]),
            0.5 * ((stress[i * 5 + 0] - stress[i * 5 + 3])));
  }
  fclose(p_vs_q);

  FILE *kappa_vs_eps = fopen("kappa_vs_eps.csv", "w");
  for (int i = 0; i < NumberSteps; i++) {
    fprintf(kappa_vs_eps, "%e, %e\n", Equiv_Plast_Str[i], kappa1[i]);
  }
  fclose(kappa_vs_eps);

  // Print data with gnuplot
  FILE *gnuplot;

  gnuplot = popen("gnuplot -persistent", "w");
  fprintf(gnuplot, "set termoption enhanced \n");
  fprintf(gnuplot, "set datafile separator ',' \n");
  fprintf(gnuplot, "set xlabel '- p' font 'Times,20' enhanced \n");
  fprintf(gnuplot, "set ylabel 'q' font "
                   "'Times,20' enhanced \n");
  fprintf(gnuplot, "plot %s \n",
          "'p_vs_q.csv.csv' title 'Us' with line lt rgb 'red' lw 2");
  fflush(gnuplot);

  // Print data with gnuplot
  gnuplot = popen("gnuplot -persistent", "w");
  fprintf(gnuplot, "set termoption enhanced \n");
  fprintf(gnuplot, "set datafile separator ',' \n");
  fprintf(gnuplot,
          "set xlabel '- {/Symbol e}_{II}' font 'Times,20' enhanced \n");
  fprintf(gnuplot, "set ylabel '{/Symbol s}_{VM}' font "
                   "'Times,20' enhanced \n");
  //  fprintf(gnuplot, "plot %s \n",
  //          "'E2_vs_VM.csv' title 'Us' with line lt rgb 'red' lw 2");
  fprintf(gnuplot, "plot %s \n", "'E2_vs_VM.csv' title 'Us'  with points pt 7");
  fflush(gnuplot);

  // Print data with gnuplot
  gnuplot = popen("gnuplot -persistent", "w");
  fprintf(gnuplot, "set termoption enhanced \n");
  fprintf(gnuplot, "set datafile separator ',' \n");
  fprintf(gnuplot, "set xlabel 'eps' font 'Times,20' enhanced \n");
  fprintf(gnuplot, "set ylabel 'kappa' font "
                   "'Times,20' enhanced \n");
  fprintf(gnuplot, "plot %s \n",
          "'kappa_vs_eps.csv' title 'Us' with line lt rgb 'red' lw 2");
  fflush(gnuplot);

  // Free memory
  free(stress);
  free(D_phi);
  free(b_e);
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
  double eigval_b_e_tr[3] = {0.0, 0.0, 0.0};
  double eigvec_b_e_tr[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double E_hencky_trial[3] = {0.0, 0.0, 0.0};
  double T_tr_vol[3] = {0.0, 0.0, 0.0};
  double T_tr_dev[3] = {0.0, 0.0, 0.0};

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  puts("t = n elastic left Cauchy-Green tensor");
  printf("%f %f %f \n", IO_State.b_e[0], IO_State.b_e[1], 0.0);
  printf("%f %f %f \n", IO_State.b_e[2], IO_State.b_e[3], 0.0);
  printf("%f %f %f \n", 0.0, 0.0, IO_State.b_e[4]);

#endif
#endif

  STATUS = __compute_trial_b_e(eigval_b_e_tr, eigvec_b_e_tr, IO_State.b_e,
                               IO_State.d_phi);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __compute_trial_b_e" RESET "\n");
    return EXIT_FAILURE;
  }

  E_hencky_trial[0] = 0.5 * log(eigval_b_e_tr[0]);
  E_hencky_trial[1] = 0.5 * log(eigval_b_e_tr[1]);
  E_hencky_trial[2] = 0.5 * log(eigval_b_e_tr[2]);

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  printf("Eigenvalues left Cauchy-Green tensor: [%f, %f, %f] \n",
         eigval_b_e_tr[0], eigval_b_e_tr[1], eigval_b_e_tr[2]);

  puts("Eigenvectors (files) left Cauchy-Green tensor");
  printf("%f %f %f \n", eigvec_b_e_tr[0], eigvec_b_e_tr[1], eigvec_b_e_tr[2]);
  printf("%f %f %f \n", eigvec_b_e_tr[3], eigvec_b_e_tr[4], eigvec_b_e_tr[5]);
  printf("%f %f %f \n", eigvec_b_e_tr[6], eigvec_b_e_tr[7], eigvec_b_e_tr[8]);

  printf("E_hencky_trial: [%f, %f, %f] \n", E_hencky_trial[0],
         E_hencky_trial[1], E_hencky_trial[2]);

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
  double Increment_E_plastic[3] = {0.0, 0.0, 0.0};

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

  STATUS = __trial_elastic(T_tr_vol, T_tr_dev, &pressure, &J2, E_hencky_trial,
                           K, G, p_ref);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __trial_elastic" RESET "\n");
    return EXIT_FAILURE;
  }

  PHI = PHI_0 = __yield_function_classical(pressure, J2, d_gamma_k, kappa_k,
                                           alpha_F, alpha_Q, beta, K, G);

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  printf("Initial value of the yield function: %f \n", PHI_0);
  printf("T trial: [%f, %f, %f] \n", -T_tr_vol[0] + T_tr_dev[0],
         -T_tr_vol[1] + T_tr_dev[1], -T_tr_vol[2] + T_tr_dev[2]);

#endif
#endif

  // Elastic
  if (PHI_0 <= 0.0) {

    STATUS = __update_internal_variables_elastic(
        IO_State.Stress, IO_State.D_phi, T_tr_vol, T_tr_dev, eigvec_b_e_tr);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in __update_internal_variables_elastic()" RESET
              "\n");
      return EXIT_FAILURE;
    }
  }
  // Plastic (check yield condition)
  else {

    STATUS = __compute_plastic_flow_direction(n, T_tr_dev, J2);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED " Error in __compute_plastic_flow_direction" RESET "\n");
      return EXIT_FAILURE;
    }

    STATUS = __compute_d_kappa(&d_kappa_k, kappa_0, eps_n, eps_0, exp_param);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __compute_d_kappa" RESET "\n");
      return EXIT_FAILURE;
    }

    STATUS = __compute_pressure_limit(&pressure_limit, J2, kappa_k, d_kappa_k,
                                      K, G, alpha_F, alpha_Q, beta);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __compute_pressure_limit" RESET "\n");
      return EXIT_FAILURE;
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
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr,
                  "" RED "Error in __compute_eps (apex loop)" RESET "\n");
          return EXIT_FAILURE;
        }

        STATUS = __compute_kappa(&kappa_k, kappa_0, exp_param, eps_k, eps_0);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr, "" RED "Error in __compute_kappa" RESET "\n");
          return EXIT_FAILURE;
        }

        STATUS =
            __compute_d_kappa(&d_kappa_k, kappa_0, eps_k, eps_0, exp_param);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr, "" RED "Error in __compute_d_kappa" RESET "\n");
          return EXIT_FAILURE;
        }

        PHI = __yield_function_classical(pressure, J2, d_gamma_k, kappa_k,
                                         alpha_F, alpha_Q, beta, K, G);
      }

      STATUS = __update_internal_variables_classical(
          Increment_E_plastic, IO_State.Stress, IO_State.Equiv_Plast_Str,
          IO_State.Kappa, IO_State.D_phi, T_tr_vol, T_tr_dev, eigvec_b_e_tr, n,
          d_gamma_k, alpha_Q, K, G, eps_k, kappa_k);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __update_internal_variables_classical" RESET
                "\n");
        return EXIT_FAILURE;
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
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __compute_eps (apex loop)" RESET "\n");
        return EXIT_FAILURE;
      }

      STATUS = __update_internal_variables_apex(
          Increment_E_plastic, IO_State.Stress, IO_State.Equiv_Plast_Str,
          IO_State.Kappa, IO_State.D_phi, T_tr_vol, T_tr_dev, eigvec_b_e_tr, n,
          d_gamma_k, d_gamma_1, alpha_Q, K, G, eps_k, kappa_k);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in__update_internal_variables_apex" RESET "\n");
        return EXIT_FAILURE;
      }
    }
  }

  E_hencky_trial[0] -= Increment_E_plastic[0];
  E_hencky_trial[1] -= Increment_E_plastic[1];
  E_hencky_trial[2] -= Increment_E_plastic[2];

  // Update elastic left Cauchy-Green tensor
  STATUS = __corrector_b_e(IO_State.b_e, eigvec_b_e_tr, E_hencky_trial);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __corrector_b_e" RESET "\n");
    return EXIT_FAILURE;
  }

  printf("Numer of iterations: %i, error: %e \n", Iter, PHI);
  if (Iter == MaxIter) {
    printf("Number of iterations: %i \n", Iter);
    printf("Final value of the yield function: %e \n", PHI);
  }

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  printf("Increment of the plastic tensor: [%e, %e, %e] \n",
         Increment_E_plastic[0], Increment_E_plastic[1],
         Increment_E_plastic[2]);
  puts("t = n + 1 elastic left Cauchy-Green tensor");
  printf("%e, %e, %e \n", IO_State.b_e[0], IO_State.b_e[1], 0.0);
  printf("%e, %e, %e \n", IO_State.b_e[2], IO_State.b_e[3], 0.0);
  printf("%e, %e, %e \n", 0.0, 0.0, IO_State.b_e[4]);
#endif
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __compute_trial_b_e(double *eigval_b_e_tr, double *eigvec_b_e_tr,
                               const double *b_e, const double *d_phi) {

  double b_e_tr[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#if NumberDimensions == 2

  b_e_tr[0] = d_phi[0] * b_e[0] * d_phi[0] + d_phi[0] * b_e[1] * d_phi[1] +
              d_phi[1] * b_e[2] * d_phi[0] + d_phi[1] * b_e[3] * d_phi[1];

  b_e_tr[1] = d_phi[0] * b_e[0] * d_phi[2] + d_phi[0] * b_e[1] * d_phi[3] +
              d_phi[1] * b_e[2] * d_phi[2] + d_phi[1] * b_e[3] * d_phi[3];

  b_e_tr[3] = d_phi[2] * b_e[0] * d_phi[0] + d_phi[2] * b_e[1] * d_phi[1] +
              d_phi[3] * b_e[2] * d_phi[0] + d_phi[3] * b_e[3] * d_phi[1];

  b_e_tr[4] = d_phi[2] * b_e[0] * d_phi[2] + d_phi[2] * b_e[1] * d_phi[3] +
              d_phi[3] * b_e[2] * d_phi[2] + d_phi[3] * b_e[3] * d_phi[3];

  b_e_tr[8] = b_e[4];

#else
  No esta implementado
#endif

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  puts("Trial elastic left Cauchy-Gree");
  printf("%f %f %f \n", b_e_tr[0], b_e_tr[1], b_e_tr[2]);
  printf("%f %f %f \n", b_e_tr[3], b_e_tr[4], b_e_tr[5]);
  printf("%f %f %f \n", b_e_tr[6], b_e_tr[7], b_e_tr[8]);

#endif
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
  double wi[3];
  double vl[9];

  /*
    Query and allocate the optimal workspace
  */
  lwork = -1;
  dgeev_("N", "V", &n, b_e_tr, &lda, eigval_b_e_tr, wi, vl, &ldvl,
         eigvec_b_e_tr, &ldvr, &wkopt, &lwork, &info);
  lwork = (int)wkopt;
  work = (double *)malloc(lwork * sizeof(double));

  /* Check for convergence */
  if (info > 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dgeev_(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }

  /* Solve eigenproblem */
  dgeev_("N", "V", &n, b_e_tr, &lda, eigval_b_e_tr, wi, vl, &ldvl,
         eigvec_b_e_tr, &ldvr, work, &lwork, &info);

  /* Check for convergence */
  if (info > 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dgeev_(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dgeev_(): the %i-th argument had an "
            "illegal value." RESET "\n",
            abs(info));
    return EXIT_FAILURE;
  }

  free(work);

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __corrector_b_e(double *b_e, const double *eigvec_b_e_tr,
                           const double *E_hencky_trial) {

  double eigval_b_e[3] = {0.0, 0.0, 0.0};

  eigval_b_e[0] = exp(2 * E_hencky_trial[0]);
  eigval_b_e[1] = exp(2 * E_hencky_trial[1]);
  eigval_b_e[2] = exp(2 * E_hencky_trial[2]);

#if NumberDimensions == 2

  b_e[0] = 0.0;
  b_e[1] = 0.0;
  b_e[2] = 0.0;
  b_e[3] = 0.0;
  b_e[4] = 0.0;

  for (unsigned i = 0; i < 3; i++) {
    b_e[0] +=
        eigval_b_e[i] * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 0];
    b_e[1] +=
        eigval_b_e[i] * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 1];
    b_e[2] +=
        eigval_b_e[i] * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 0];
    b_e[3] +=
        eigval_b_e[i] * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 1];
    b_e[4] +=
        eigval_b_e[i] * eigvec_b_e_tr[i * 3 + 2] * eigvec_b_e_tr[i * 3 + 2];
  }

#else
  No esta implementado
#endif

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __trial_elastic(double *T_tr_vol, double *T_tr_dev, double *pressure,
                           double *J2, const double *E_hencky_trial, double K,
                           double G, double p_ref) {

  double E_hencky_trial_vol[3];
  double tr_E_hencky_trial =
      E_hencky_trial[0] + E_hencky_trial[1] + E_hencky_trial[2];

  E_hencky_trial_vol[0] = (1.0 / 3.0) * tr_E_hencky_trial;
  E_hencky_trial_vol[1] = (1.0 / 3.0) * tr_E_hencky_trial;
  E_hencky_trial_vol[2] = (1.0 / 3.0) * tr_E_hencky_trial;

  T_tr_vol[0] = -p_ref; //-K * E_hencky_trial_vol[0];
  T_tr_vol[1] = -p_ref; //-K * E_hencky_trial_vol[1];
  T_tr_vol[2] = -p_ref; //-K * E_hencky_trial_vol[2];

  T_tr_dev[0] = 0.0; // 2 * G * (E_hencky_trial[0] - E_hencky_trial_vol[0]);
  T_tr_dev[1] =
      YoungMouduls * (E_hencky_trial[1]); // 2 * G * (E_hencky_trial[1] -
                                          // E_hencky_trial_vol[1]);
  T_tr_dev[2] = 0.0; // 2 * G * (E_hencky_trial[2] - E_hencky_trial_vol[2]);

  *pressure = (T_tr_vol[0] + T_tr_vol[1] + T_tr_vol[2]) / 3.0;
  *J2 = sqrt(T_tr_dev[0] * T_tr_dev[0] + T_tr_dev[1] * T_tr_dev[1] +
             T_tr_dev[2] * T_tr_dev[2]);

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __update_internal_variables_elastic(double *Stress,
                                               const double *D_phi,
                                               const double *T_tr_vol,
                                               const double *T_tr_dev,
                                               const double *eigvec_b_e_tr) {

  // Compute the transpose of D_phi
  double D_phi_mT[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#if NumberDimensions == 2

  D_phi_mT[0] = D_phi[0];
  D_phi_mT[1] = D_phi[2];
  D_phi_mT[3] = D_phi[1];
  D_phi_mT[4] = D_phi[3];
  D_phi_mT[8] = D_phi[4];

#else
  No esta implementado
#endif

  // compute the inverse of D_phi
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N, &N, D_phi_mT, &LDA, IPIV, &INFO);
  // Check output of dgetrf
  if (INFO != 0) {
    if (INFO < 0) {
      printf(
          "" RED
          "Error in dgetrf_(): the %i-th argument had an illegal value " RESET
          "\n",
          abs(INFO));
    } else if (INFO > 0) {

      printf("" RED
             "Error in dgetrf_(): D_phi_mT(%i,%i) %s \n %s \n %s \n %s " RESET
             "\n",
             INFO, INFO, "is exactly zero. The factorization",
             "has been completed, but the factor D_phi_mT is exactly",
             "singular, and division by zero will occur if it is used",
             "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  dgetri_(&N, D_phi_mT, &LDA, IPIV, WORK, &LWORK, &INFO);
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(stderr, "" RED "%s: the %i-th argument %s" RESET "\n",
              "Error in dgetri_()", abs(INFO), "had an illegal value");
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetri_(): D_phi_mT(%i,%i) %s \n %s \n %s \n %s " RESET
              "\n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor D_phi_mT is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  puts("Adjunt of the deformation gradient");
  printf("%f %f %f \n", D_phi_mT[0], D_phi_mT[1], D_phi_mT[2]);
  printf("%f %f %f \n", D_phi_mT[3], D_phi_mT[4], D_phi_mT[5]);
  printf("%f %f %f \n", D_phi_mT[6], D_phi_mT[7], D_phi_mT[8]);

#endif
#endif

  double T[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for (unsigned i = 0; i < 3; i++) {

    double T_i = -T_tr_vol[i] + T_tr_dev[i];

    T[0] += T_i * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 0];
    T[1] += T_i * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 1];
    T[3] += T_i * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 0];
    T[4] += T_i * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 1];
    T[8] += T_i * eigvec_b_e_tr[i * 3 + 2] * eigvec_b_e_tr[i * 3 + 2];
  }

#if NumberDimensions == 2

  Stress[0] = T[0] * D_phi_mT[0] + T[1] * D_phi_mT[3];
  Stress[1] = T[0] * D_phi_mT[1] + T[1] * D_phi_mT[4];
  Stress[2] = T[3] * D_phi_mT[0] + T[4] * D_phi_mT[3];
  Stress[3] = T[3] * D_phi_mT[1] + T[4] * D_phi_mT[4];
  Stress[4] = T[8] * D_phi_mT[8];
#else
  No esta implementado
#endif

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
#if NumberDimensions == 2
  puts("Nominal stress tensor");
  printf("%f %f %f \n", Stress[0], Stress[1], 0.0);
  printf("%f %f %f \n", Stress[2], Stress[3], 0.0);
  printf("%f %f %f \n", 0.0, 0.0, Stress[4]);
#endif
#endif
#endif

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __compute_plastic_flow_direction(double *n, const double *T_tr_dev,
                                            double J2) {

  n[0] = T_tr_dev[0] / J2;
  n[1] = T_tr_dev[1] / J2;
  n[2] = T_tr_dev[2] / J2;

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __compute_eps(double *eps_k, double d_gamma_k, double eps_n,
                         double alpha_Q) {

  *eps_k = eps_n + d_gamma_k * sqrt(3.0 * alpha_Q * alpha_Q + 1.0);

  if (*eps_k < 0.0) {
    fprintf(stderr, "" RED "Negative value of the EPS_k: %f " RESET "\n",
            *eps_k);
    return EXIT_FAILURE;
  }

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
    double *Increment_E_plastic, double *Stress, double *eps_n1,
    double *kappa_n1, const double *D_phi, const double *T_tr_vol,
    const double *T_tr_dev, const double *eigvec_b_e_tr, const double *n,
    double d_gamma_k, double alpha_Q, double K, double G, double eps_k,
    double kappa_k) {

  // Update hardening parameters
  *eps_n1 = eps_k;
  *kappa_n1 = kappa_k;

  // Update the increment of the plastic strain tensor
  Increment_E_plastic[0] = d_gamma_k * (alpha_Q + n[0]);
  Increment_E_plastic[1] = d_gamma_k * (alpha_Q + n[1]);
  Increment_E_plastic[2] = d_gamma_k * (alpha_Q + n[2]);

  // Compute the transpose of D_phi
  double D_phi_mT[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#if NumberDimensions == 2

  D_phi_mT[0] = D_phi[0];
  D_phi_mT[1] = D_phi[2];
  D_phi_mT[3] = D_phi[1];
  D_phi_mT[4] = D_phi[3];
  D_phi_mT[8] = D_phi[4];

#else
  No esta implementado
#endif

  // compute the inverse of D_phi
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N, &N, D_phi_mT, &LDA, IPIV, &INFO);
  // Check output of dgetrf
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(
          stderr,
          "" RED
          "Error in dgetrf_(): the %i-th argument had an illegal value" RESET
          "",
          abs(INFO));
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetrf_(): D_phi_mT(%i,%i) %s \n %s \n %s \n %s" RESET
              " \n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor D_phi_mT is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  dgetri_(&N, D_phi_mT, &LDA, IPIV, WORK, &LWORK, &INFO);
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(stderr,
              "" RED "Error in dgetri_(): the %i-th argument of dgetrf_ had an "
              "illegal value" RESET "\n",
              abs(INFO));
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetri_(): D_phi_mT(%i,%i) %s \n %s \n %s \n %s " RESET
              "\n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor D_phi_mT is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  puts("Adjunt of the deformation gradient");
  printf("%f %f %f \n", D_phi_mT[0], D_phi_mT[1], D_phi_mT[2]);
  printf("%f %f %f \n", D_phi_mT[3], D_phi_mT[4], D_phi_mT[5]);
  printf("%f %f %f \n", D_phi_mT[6], D_phi_mT[7], D_phi_mT[8]);

#endif
#endif

  double T[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for (unsigned i = 0; i < 3; i++) {

    double T_i = -T_tr_vol[i] + T_tr_dev[i] -
                 d_gamma_k * (3 * K * alpha_Q + 2 * G * n[i]);

    T[0] += T_i * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 0];
    T[1] += T_i * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 1];
    T[3] += T_i * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 0];
    T[4] += T_i * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 1];
    T[8] += T_i * eigvec_b_e_tr[i * 3 + 2] * eigvec_b_e_tr[i * 3 + 2];
  }

#if NumberDimensions == 2

  Stress[0] = T[0] * D_phi_mT[0] + T[1] * D_phi_mT[3];
  Stress[1] = T[0] * D_phi_mT[1] + T[1] * D_phi_mT[4];
  Stress[2] = T[3] * D_phi_mT[0] + T[4] * D_phi_mT[3];
  Stress[3] = T[3] * D_phi_mT[1] + T[4] * D_phi_mT[4];
  Stress[4] = T[8] * D_phi_mT[8];

#else
  No esta implementado
#endif

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
#if NumberDimensions == 2
  puts("Nominal stress tensor");
  printf("%f %f %f \n", Stress[0], Stress[1], 0.0);
  printf("%f %f %f \n", Stress[2], Stress[3], 0.0);
  printf("%f %f %f \n", 0.0, 0.0, Stress[4]);
#endif
#endif
#endif

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

static int __update_internal_variables_apex(
    double *Increment_E_plastic, double *Stress, double *eps_n1,
    double *kappa_n1, const double *D_phi, const double *T_tr_vol,
    const double *T_tr_dev, const double *eigvec_b_e_tr, const double *n,
    double d_gamma_k, double d_gamma_1, double alpha_Q, double K, double G,
    double eps_k, double kappa_k) {

  // Update hardening parameters
  *eps_n1 = eps_k;
  *kappa_n1 = kappa_k;

  // Update the increment of the plastic strain tensor
  Increment_E_plastic[0] = d_gamma_k * alpha_Q + d_gamma_1 * n[0];
  Increment_E_plastic[1] = d_gamma_k * alpha_Q + d_gamma_1 * n[1];
  Increment_E_plastic[2] = d_gamma_k * alpha_Q + d_gamma_1 * n[2];

  // Compute the transpose of D_phi
  double D_phi_mT[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#if NumberDimensions == 2

  D_phi_mT[0] = D_phi[0];
  D_phi_mT[1] = D_phi[2];
  D_phi_mT[3] = D_phi[1];
  D_phi_mT[4] = D_phi[3];
  D_phi_mT[8] = D_phi[4];

#else
  No esta implementado
#endif

  // compute the inverse of D_phi
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N, &N, D_phi_mT, &LDA, IPIV, &INFO);
  // Check output of dgetrf
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(
          stderr,
          "" RED
          "Error in dgetrf_(): the %i-th argument had an illegal value" RESET
          "",
          abs(INFO));
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetrf_(): D_phi_mT(%i,%i) %s \n %s \n %s \n %s" RESET
              " \n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor D_phi_mT is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  dgetri_(&N, D_phi_mT, &LDA, IPIV, WORK, &LWORK, &INFO);
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(stderr,
              "" RED "Error in dgetri_(): the %i-th argument of dgetrf_ had an "
              "illegal value" RESET "\n",
              abs(INFO));
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetri_(): D_phi_mT(%i,%i) %s \n %s \n %s \n %s " RESET
              "\n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor D_phi_mT is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  puts("Adjunt of the deformation gradient");
  printf("%f %f %f \n", D_phi_mT[0], D_phi_mT[1], D_phi_mT[2]);
  printf("%f %f %f \n", D_phi_mT[3], D_phi_mT[4], D_phi_mT[5]);
  printf("%f %f %f \n", D_phi_mT[6], D_phi_mT[7], D_phi_mT[8]);

#endif
#endif

  double T[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for (unsigned i = 0; i < 3; i++) {

    double T_i = -T_tr_vol[i] - d_gamma_k * 3 * K * alpha_Q;

    T[0] += T_i * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 0];
    T[1] += T_i * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 1];
    T[3] += T_i * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 0];
    T[4] += T_i * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 1];
    T[8] += T_i * eigvec_b_e_tr[i * 3 + 2] * eigvec_b_e_tr[i * 3 + 2];
  }

#if NumberDimensions == 2

  Stress[0] = T[0] * D_phi_mT[0] + T[1] * D_phi_mT[3];
  Stress[1] = T[0] * D_phi_mT[1] + T[1] * D_phi_mT[4];
  Stress[2] = T[3] * D_phi_mT[0] + T[4] * D_phi_mT[3];
  Stress[3] = T[3] * D_phi_mT[1] + T[4] * D_phi_mT[4];
  Stress[4] = T[8] * D_phi_mT[8];

#else
  No esta implementado
#endif

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
#if NumberDimensions == 2
  puts("Nominal stress tensor");
  printf("%f %f %f \n", Stress[0], Stress[1], 0.0);
  printf("%f %f %f \n", Stress[2], Stress[3], 0.0);
  printf("%f %f %f \n", 0.0, 0.0, Stress[4]);
#endif
#endif
#endif

  return EXIT_SUCCESS;
}

/***************************************************************************/

static bool __degradated(const double *Stress, double J2_degradated) {

  double J2 = 0.0;

#if NumberDimensions == 2
  J2 = sqrt(0.5 * (pow((Stress[0] - Stress[3]), 2.0) +
                   pow((Stress[3] - Stress[4]), 2.0) +
                   pow((Stress[4] - Stress[0]), 2.0)));
#else
  No esta implementado
#endif

  if (J2 <= J2_degradated) {
    return true;
  } else {
    return false;
  }
}

/***************************************************************************/