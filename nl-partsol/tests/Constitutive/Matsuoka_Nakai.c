/**************************************************************/
/******************* Material Parameters **********************/
/**************************************************************/
// Kenichi
/*
#define NumberDimensions 2
#define YoungMouduls 10.0E3
#define PoissonRatio 0.2
#define AtmosphericPressure -100.0
#define a1_Parameter 370.0
#define a2_Parameter 0.01
#define a//3_Parameter 5.0
#define phi_Parameter 37.0
#define alpha_Parameter 0.10 // psi: 6.0
#define c_cotphi_value 8.0
#define NumberSteps 2000 // 1000 // 5000 // 3500 // 2500 //
#define Delta_strain_II -0.00001
#define Confining_pressure -20.0
*/
/*
// Borja
#define NumberDimensions 2
#define YoungMouduls 100.0E3
#define PoissonRatio 0.2
#define AtmosphericPressure -100.0
#define c_cotphi_value 0.0
#define a1_Parameter 20000.0
#define a2_Parameter 0.005
#define a3_Parameter 35.0
#define alpha_Parameter 0.5
#define NumberSteps 3000 // 1000 // 5000 // 3500 // 2500 //
#define Delta_strain_II -0.00001
#define Confining_pressure -200.0
#define EPS_0 0.0
#define kappa_0 0.0
*/

// Borja
#define NumberDimensions 2
#define YoungMouduls 10.0E3
#define PoissonRatio 0.2
#define AtmosphericPressure -100.0
#define c_cotphi_value 0.0
#define a1_Parameter 10.0
#define a2_Parameter 0.0
#define a3_Parameter 0.8
#define alpha_Parameter 0.162
#define NumberSteps 20000 // 1000 // 5000 // 3500 // 2500 //
#define Delta_strain_II -0.00001
#define Confining_pressure -20.0
#define EPS_0 1.065199
#define kappa_0 4.543

/*
Unitary test for the smooth Mohr-Coulomb model.
One single point is initially confined with an
hydrostatic stress state of -200 kPa. Before that,
the stress in the II direction is incresed using
strain control.
-------------------------------------------------------------------
Written by Miguel Molinos, 2021
Universidad Politécnica de Madrid
-------------------------------------------------------------------
Useful information is in Borja et al. (2003):
'On the numerical integration of three-invariant elastoplastic constitutive
models' https://doi.org/10.1016/S0045-7825(02)00620-5
-------------------------------------------------------------------
Requires:
  * accelerate library or LAPACK (solver)
  * gnuplot (plots)
-------------------------------------------------------------------
Compilation recipy (MacOSX):
gcc -framework Accelerate Frictional-Monolithic.c -o Frictional-Monolithic
-------------------------------------------------------------------
Compilation recipy (Linux)
gcc Frictional-Monolithic.c -o Frictional-Monolithic -llapack -lm
-------------------------------------------------------------------
*/

#ifdef __linux__
#include <lapacke.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
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
  double Cohesion;

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
  double *Strain_e;
  double Pressure;

  /*!
   * Finite strain kinematic parameters
   * */
  double *d_phi;
  double *D_phi_n1;
  double *D_phi_n;
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
  double *EPS;   // Equivalent plastic strain
  double *Kappa; // Hardening Parameter

  //
  bool *Failure;
  bool compute_C_ep;
  double *C_ep;

} State_Parameters;

/*
  Auxiliar functions
*/

static void __elastic_tangent(double *CC /**< [out] Elastic compliance */,
                              double *AA /**< [out] Elastic matrix */,
                              double E /**< [in] Young modulus */,
                              double nu /**< [in] Poisson ratio */,
                              double Lame /**< [in] Lamé parameter */,
                              double G /**< [in] Shear modulus */);

static void
__trial_elastic(double *T_tr /**< [out] Trial elastic stress tensor*/,
                const double *E_hencky_trial /**< [in] Henky strain (trial) */,
                const double *AA /**< [in] Elastic matrix */,
                double c_cotphi /**< [in] Cohesion parameter */);

static void
__E_hencky(double *E_hencky_k /**< [out] Henky strain (iter k) */,
           const double *T_k /**< [in] Local stress tensor (iter k) */,
           const double *CC /**< [in] Elastic compliance */,
           double c_cotphi /**< [in] Cohesion parameter */);

static void
__kappa(double *kappa /**< [out] Hardening vector */,
        const double *a /**< [in] Vector with fit parameters (hardening) */,
        double Lambda /**< [in] Total plastic multiplier */,
        double I1 /**< [in] First invariant of the stress tensor */,
        double alpha /**< [in] Dilatance parameter*/);

static void __d_kappa_phi_d_stress(
    double *d_kappa_phi_d_stress /**< [out] Stress derivative of kappa[0] */,
    const double *a /**< [in] Vector with fit parameters (hardening) */,
    double Lambda /**< [in] Total plastic multiplier */,
    double I1 /**< [in] First invariant of the stress tensor */);

static void __d_kappa_phi_d_lambda(
    double *d_kappa_phi_d_lambda /**< [out] Lambda derivative of kappa[0] */,
    const double *a /**< [in] Vector with fit parameters (hardening) */,
    double Lambda /**< [in] Total plastic multiplier */,
    double I1 /**< [in] First invariant of the stress tensor */);

static double __F(double *F /**< [out] Yield function evaluation */,
                  double kappa_phi /**< [in] Friction angle hardening */,
                  double I1 /**< [in] First invariant of the stress tensor */,
                  double I2 /**< [in] Second invariant of the stress tensor */,
                  double I3 /**< [in] Third invariant of the stress tensor */);

static void __d_F_d_stress(
    double *d_F_d_stress /**< [out] Yield function derivative (stress) */,
    const double *T_k /**< [in] Local stress tensor (iter k) */,
    double I1 /**< [in] First invariant of the stress tensor */,
    double I2 /**< [in] Second invariant of the stress tensor */,
    double I3 /**< [in] Third invariant of the stress tensor */,
    double kappa_phi /**< [in] Friction angle hardening */);

static void __d_F_d_kappa_phi(
    double *d_F_d_kappa_phi /**< [out] Yield function derivative (kappa[0]) */,
    double I1 /**< [in] First invariant of the stress tensor */,
    double I3 /**< [in] Third invariant of the stress tensor */,
    double kappa_phi /**< [in] Friction angle hardening */);

static void
__d_G_d_stress(double *d_G_d_stress /**< [out] Plastic potential function
                                       derivative (stress) */
               ,
               const double *T_k /**< [in] Local stress tensor (iter k) */,
               double I1 /**< [in] First invariant of the stress tensor */,
               double I2 /**< [in] Second invariant of the stress tensor */,
               double I3 /**< [in] Third invariant of the stress tensor */,
               double kappa_psi /**< [in] Dilatance angle hardening */);

static void __dd_G_dd_stress(
    double *dd_G_dd_stress /**< [out] Plastic potential hessian (stress) */,
    const double *T_k /**< [in] Local stress tensor (iter k) */,
    double kappa_psi /**< [in] Dilatance angle hardening */,
    double I1 /**< [in] First invariant of the stress tensor */,
    double I2 /**< [in] Second invariant of the stress tensor */,
    double I3 /**< [in] Third invariant of the stress tensor */);

static void __dd_G_d_stress_d_kappa_psi(
    double *dd_G_d_stress_d_kappa_psi /**< [out] Plastic potential deriv */,
    const double *T_k /**< [in] Local stress tensor (iter k) */,
    double I1 /**< [in] First invariant of the stress tensor */,
    double I3 /**< [in] Third invariant of the stress tensor */,
    double kappa_psi /**< [in] Dilatance angle hardening */);

static int __residual(
    double *Residual /**< [out] Residual of the problem */,
    double *Error_k /**< [out] Norm of the residual */,
    const double *E_hencky_trial /**< [in] Henky strain (trial) */,
    const double *E_hencky_k /**< [in] Henky strain (iter k) */,
    const double *
        d_G_d_stress /**< [in] Plastic potential function derivative (stress) */
    ,
    const double *kappa_k /**< [in] Hardening vector (iter k) */,
    const double *kappa_hat /**< [in] Hardening vector (eval) */,
    double F_k /**< [in] Yield function evaluation (iter k) */,
    double delta_lambda_k /**< [in] Discrete plastic multiplier (iter k) */);

static int
__solver(double *Tangent_Matrix /**< [in/out] Tangent matrix of the problem */,
         double *Residual /**< [in/out] Residual of the problem */);

static int __reciprocal_condition_number(
    double *RCOND /**< [out] Condition number of the tangent matrix */,
    double *Tangent_Matrix /**< [in/out] Tangent matrix of the problem */);

int compute_1PK_Matsuoka_Nakai__Constitutive__(State_Parameters IO_State, Material MatProp);

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

int main() {
  int STATUS = EXIT_SUCCESS;
  State_Parameters IO_State;
  Material MatProp;

  // Define material variables
  MatProp.E = YoungMouduls;
  MatProp.nu = PoissonRatio;
  MatProp.atmospheric_pressure = AtmosphericPressure;
  MatProp.a_Hardening_Borja[0] = a1_Parameter;
  MatProp.a_Hardening_Borja[1] = a2_Parameter;
  MatProp.a_Hardening_Borja[2] = a3_Parameter;
  MatProp.alpha_Hardening_Borja = alpha_Parameter;
  MatProp.Cohesion = 0.0;

  // Initialize state variables
  bool Status_particle = false;
  double *stress = (double *)calloc(3 * NumberSteps, sizeof(double));
  double *strain = (double *)calloc(3 * NumberSteps, sizeof(double));
  double *strain_e = (double *)calloc(3 * NumberSteps, sizeof(double));
  double *kappa1 = (double *)calloc(NumberSteps, sizeof(double));
  double *EPS = (double *)calloc(NumberSteps, sizeof(double));

  // Set initial value stress
  stress[0] = Confining_pressure;
  stress[1] = Confining_pressure;
  stress[2] = Confining_pressure;
  //
  //  double rad_friction_angle = (PI__MatrixLib__ / 180.0) * phi_Parameter;
  double a1 = MatProp.a_Hardening_Borja[0];
  double a2 = MatProp.a_Hardening_Borja[1];
  double a3 = MatProp.a_Hardening_Borja[2];
  double f, df;
  double I1 = 3.0 * Confining_pressure;
  int iter = 0;

  kappa1[0] = kappa_0;
  EPS[0] = EPS_0;

  // Start time integration
  for (int i = 1; i < NumberSteps; i++) {

    printf("*********************\n");
    printf("Step: %i \n", i);

    // Asign variables to the solver
    IO_State.Particle_Idx = 0;
    IO_State.Stress = &stress[i * 3];
    IO_State.Strain = &strain[i * 3];
    IO_State.Strain_e = &strain_e[i * 3];
    IO_State.EPS = &EPS[i];
    IO_State.Kappa = &kappa1[i];
    IO_State.Failure = &Status_particle;
    IO_State.compute_C_ep = false;

    // Trial strain
    IO_State.Strain[1] = strain[(i - 1) * 3 + 1] + Delta_strain_II;
    IO_State.Strain_e[1] = IO_State.Strain[1];

    // Trial stress
    IO_State.Stress[0] = Confining_pressure;
    IO_State.Stress[1] =
        stress[(i - 1) * 3 + 1] + YoungMouduls * Delta_strain_II;
    IO_State.Stress[2] = Confining_pressure;

    *IO_State.EPS = EPS[i - 1];
    *IO_State.Kappa = kappa1[i - 1];

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

    printf("E_hencky_trial: [%e, %e, %e] \n", IO_State.Strain[0],
           IO_State.Strain[1], IO_State.Strain[2]);

#endif
#endif

    // Run solver
    STATUS = compute_1PK_Matsuoka_Nakai__Constitutive__(IO_State, MatProp);

    if (STATUS) {
      fprintf(stderr, "Failure\n");
      free(stress);
      free(strain);
      free(kappa1);
      free(EPS);
      return STATUS;
    }
  }

  // Save data in a csv file
  FILE *E2_vs_VM = fopen("E2_vs_VM.csv", "w");
  for (int i = 0; i < NumberSteps; i++) {
    fprintf(E2_vs_VM, "%e, %e\n", -100 * strain[i * 3 + 1],
            fabs(stress[i * 3 + 1] - stress[i * 3 + 0]));
  }
  fclose(E2_vs_VM);

  FILE *p_vs_q = fopen("p_vs_q.csv.csv", "w");
  for (int i = 0; i < NumberSteps; i++) {
    fprintf(p_vs_q, "%e, %e\n", -0.5 * (stress[i * 3 + 0] + stress[i * 3 + 1]),
            0.5 * ((stress[i * 3 + 0] - stress[i * 3 + 1])));
  }
  fclose(p_vs_q);

  FILE *kappa_vs_eps = fopen("kappa_vs_eps.csv", "w");
  for (int i = 0; i < NumberSteps; i++) {
    fprintf(kappa_vs_eps, "%e, %e, %e\n", EPS[i], kappa1[i],
            (180.0 / PI__MatrixLib__) *
                asin(sqrt(kappa1[i] / (kappa1[i] + 8))));
  }
  fclose(kappa_vs_eps);

  // Print data with gnuplot
  FILE *gnuplot;

  gnuplot = popen("gnuplot -persistent", "w");
  fprintf(gnuplot, "set termoption enhanced \n");
  fprintf(gnuplot, "set datafile separator ',' \n");
  fprintf(gnuplot, "set xlabel '-E2' font 'Times,20' enhanced \n");
  fprintf(gnuplot, "set ylabel 'abs(s2 - s1)' font "
                   "'Times,20' enhanced \n");
  fprintf(gnuplot, "plot %s \n",
          "'E2_vs_VM.csv' title 'Us' with line lt rgb 'red' lw 2");
  fflush(gnuplot);

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
  fprintf(gnuplot, "set xlabel 'eps' font 'Times,20' enhanced \n");
  fprintf(gnuplot, "set ylabel 'kappa' font "
                   "'Times,20' enhanced \n");
  fprintf(
      gnuplot, "plot %s \n",
      "'kappa_vs_eps.csv' using 1:3 title 'Us' with line lt rgb 'red' lw 2");
  fflush(gnuplot);

  // Free memory
  free(stress);
  free(strain);
  free(kappa1);
  free(EPS);

  return STATUS;
}

/**************************************************************/

int compute_1PK_Matsuoka_Nakai__Constitutive__(State_Parameters IO_State, Material MatProp)
/*
  Monolithic algorithm for a smooth Mohr-Coulomb plastic criterium
*/
{

  int STATUS = EXIT_SUCCESS;

  // Read input/output parameters
  double eigval_b_e_tr[3] = {0.0, 0.0, 0.0};
#if NumberDimensions == 2
  double eigvec_b_e_tr[4] = {0.0, 0.0, 0.0, 0.0};
#else
  double eigvec_b_e_tr[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif
  double E_hencky_trial[3] = {0.0, 0.0, 0.0};
  double E_hencky_k1[3] = {0.0, 0.0, 0.0};
  double E_hencky_k2[3] = {0.0, 0.0, 0.0};

  E_hencky_trial[0] = IO_State.Strain[0]; // 0.5 * log(eigval_b_e_tr[0]);
  E_hencky_trial[1] = IO_State.Strain[1]; // 0.5 * log(eigval_b_e_tr[1]);
  E_hencky_trial[2] = IO_State.Strain[2]; // 0.5 * log(eigval_b_e_tr[2]);

  // Material parameters
  double E = MatProp.E;
  double nu = MatProp.nu;
  double Lame = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  double G = E / (2.0 * (1.0 + nu));
  double friction_angle = MatProp.phi_Frictional;
  //  double rad_friction_angle = (PI__MatrixLib__ / 180.0) * friction_angle;
  double c_cotphi = c_cotphi_value;
  double alpha = MatProp.alpha_Hardening_Borja;
  double a[3] = {0.0, 0.0, 0.0};
  a[0] = MatProp.a_Hardening_Borja[0];
  a[1] = MatProp.a_Hardening_Borja[1];
  a[2] = MatProp.a_Hardening_Borja[2];

  bool Activate_CutOff = false;
  double CutOff = 0.0;

  double CC[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double AA[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  __elastic_tangent(CC, AA, E, nu, Lame, G);

  // Define scalar variables internal variables
  double F_k1, F_k2, F_0;
  double I1, I2, I3;
  double Lambda_n = *IO_State.EPS;
  double Lambda_k1;
  double Lambda_k2;
  double delta_lambda_k0 = 0.0;
  double delta_lambda_k1;
  double delta_lambda_k2;

  // Define tensorial internal variables
  double T_tr[3] = {0.0, 0.0, 0.0}; // Trial elastic stress
  double T_k1[3];                   // Stress iteration k
  double T_k2[3];                   // Stress iteration k (line search)
  double kappa_n[2];
  kappa_n[0] = (*IO_State.Kappa);
  kappa_n[1] = alpha * (*IO_State.Kappa);
  double kappa_k1[2]; // Hardening iteration k
  double kappa_k2[2]; // Hardening iteration k (line search)
  double kappa_hat[3] = {0.0, 0.0, 0.0};
  double d_G_d_stress[3] = {0.0, 0.0, 0.0}; // Plastic flow
  double dd_G_dd_stress[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dd_G_d_stress_d_kappa_psi[3] = {0.0, 0.0, 0.0};
  double d_kappa_phi_d_stress[3] = {0.0, 0.0, 0.0};
  double d_kappa_phi_d_lambda = 0.0;
  double d_F_d_stress[3] = {0.0, 0.0, 0.0};
  double d_F_d_kappa_phi = 0;

  /*
    Initialize Newton-Raphson solver
  */
  double TOL = 1E-10; // TOL_Radial_Returning;
  double Residual_k1[5] = {0.0, 0.0, 0.0};
  double Residual_k2[5] = {0.0, 0.0, 0.0};
  double Tangent_Matrix[25];
  double rcond;
  double Norm_Residual_k0;
  double Norm_Residual_k1;
  double Norm_Residual_k2;
  double delta = 1;
  int MaxIter_k1 = Max_Iterations_Radial_Returning;
  int MaxIter_k2 = 10 * Max_Iterations_Radial_Returning;
  int Iter_k1 = 0;
  int Iter_k2 = 0;

  T_tr[0] = IO_State.Stress[0] - c_cotphi;
  T_tr[1] = IO_State.Stress[1] - c_cotphi;
  T_tr[2] = IO_State.Stress[2] - c_cotphi;

  I1 = T_tr[0] + T_tr[1] + T_tr[2];
  I2 = T_tr[0] * T_tr[1] + T_tr[1] * T_tr[2] + T_tr[0] * T_tr[2];
  I3 = T_tr[0] * T_tr[1] * T_tr[2];

  // Check yield
  __F(&F_0, kappa_n[0], I1, I2, I3);

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  printf("Initial value of the yield function: %e \n", F_0);
  printf("T trial: [%e, %e, %e] \n", T_tr[0], T_tr[1], T_tr[2]);
#endif
#endif

  // Assign trial (elastic) stress-strain values
  T_k1[0] = T_tr[0];
  T_k1[1] = T_tr[1];
  T_k1[2] = T_tr[2];

  // Elastic
  if (F_0 <= 0.0) {

    IO_State.Stress[0] = T_k1[0] + c_cotphi;
    IO_State.Stress[1] = T_k1[1] + c_cotphi;
    IO_State.Stress[2] = T_k1[2] + c_cotphi;
    IO_State.Strain_e[0] = IO_State.Strain[0];
    IO_State.Strain_e[1] = IO_State.Strain[1];
    IO_State.Strain_e[2] = IO_State.Strain[2];
    //    *IO_State.EPS = Lambda_n;

  }
  // Plastic (monolithic solver with line search)
  else {

    __E_hencky(E_hencky_k1, T_k1, CC, c_cotphi);

    E_hencky_trial[0] = E_hencky_k1[0];
    E_hencky_trial[1] = E_hencky_k1[1];
    E_hencky_trial[2] = E_hencky_k1[2];

    __kappa(kappa_hat, a, Lambda_n, I1, alpha);

    __d_G_d_stress(d_G_d_stress, T_tr, I1, I2, I3, kappa_n[1]);

    STATUS =
        __residual(Residual_k1, &Norm_Residual_k0, E_hencky_trial, E_hencky_k1,
                   d_G_d_stress, kappa_n, kappa_hat, F_0, delta_lambda_k0);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __residual (trial)" RESET "\n");
      return EXIT_FAILURE;
    }

    // Assign values to the k iteration
    kappa_k1[0] = kappa_n[0];
    kappa_k1[1] = kappa_n[1];
    F_k1 = F_0;
    delta_lambda_k1 = delta_lambda_k0;
    Lambda_k1 = Lambda_n;
    Norm_Residual_k1 = Norm_Residual_k0;
    Iter_k1 = 0;

    while (fabs(Norm_Residual_k1 / Norm_Residual_k0) >= TOL) {

      delta = 1.0;

      // Evaluate hardening derivatives
      __d_kappa_phi_d_stress(d_kappa_phi_d_stress, a, Lambda_k1, I1);
      __d_kappa_phi_d_lambda(&d_kappa_phi_d_lambda, a, Lambda_k1, I1);
      // Evaluate yield function derivatives
      __d_F_d_stress(d_F_d_stress, T_k1, I1, I2, I3, kappa_k1[0]);
      __d_F_d_kappa_phi(&d_F_d_kappa_phi, I1, I3, kappa_k1[0]);

      // Evaluate plastic flow rule derivatives
      __dd_G_dd_stress(dd_G_dd_stress, T_k1, kappa_k1[1], I1, I2, I3);
      __dd_G_d_stress_d_kappa_psi(dd_G_d_stress_d_kappa_psi, T_k1, I1, I3,
                                  kappa_k1[1]);

      // Assemble tangent matrix
      Tangent_Matrix[0] = CC[0] + delta_lambda_k1 * dd_G_dd_stress[0];
      Tangent_Matrix[1] = CC[1] + delta_lambda_k1 * dd_G_dd_stress[1];
      Tangent_Matrix[2] = CC[2] + delta_lambda_k1 * dd_G_dd_stress[2];
      Tangent_Matrix[3] =
          alpha * delta_lambda_k1 * dd_G_d_stress_d_kappa_psi[0];
      Tangent_Matrix[4] = d_G_d_stress[0];

      Tangent_Matrix[5] = CC[3] + delta_lambda_k1 * dd_G_dd_stress[3];
      Tangent_Matrix[6] = CC[4] + delta_lambda_k1 * dd_G_dd_stress[4];
      Tangent_Matrix[7] = CC[5] + delta_lambda_k1 * dd_G_dd_stress[5];
      Tangent_Matrix[8] =
          alpha * delta_lambda_k1 * dd_G_d_stress_d_kappa_psi[1];
      Tangent_Matrix[9] = d_G_d_stress[1];

      Tangent_Matrix[10] = CC[6] + delta_lambda_k1 * dd_G_dd_stress[6];
      Tangent_Matrix[11] = CC[7] + delta_lambda_k1 * dd_G_dd_stress[7];
      Tangent_Matrix[12] = CC[8] + delta_lambda_k1 * dd_G_dd_stress[8];
      Tangent_Matrix[13] =
          alpha * delta_lambda_k1 * dd_G_d_stress_d_kappa_psi[2];
      Tangent_Matrix[14] = d_G_d_stress[2];

      Tangent_Matrix[15] = -d_kappa_phi_d_stress[0];
      Tangent_Matrix[16] = -d_kappa_phi_d_stress[1];
      Tangent_Matrix[17] = -d_kappa_phi_d_stress[2];
      Tangent_Matrix[18] = 1.0;
      Tangent_Matrix[19] = -d_kappa_phi_d_lambda;

      Tangent_Matrix[20] = d_F_d_stress[0];
      Tangent_Matrix[21] = d_F_d_stress[1];
      Tangent_Matrix[22] = d_F_d_stress[2];
      Tangent_Matrix[23] = d_F_d_kappa_phi;
      Tangent_Matrix[24] = 0.0;

      // Compute increments and update variables
      STATUS = __solver(Tangent_Matrix, Residual_k1);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "__solver" RESET "\n");
        return EXIT_FAILURE;
      }

      // Update values for the next step (line search)
      T_k2[0] = T_k1[0] - delta * Residual_k1[0];
      T_k2[1] = T_k1[1] - delta * Residual_k1[1];
      T_k2[2] = T_k1[2] - delta * Residual_k1[2];
      kappa_k2[0] = kappa_k1[0] - delta * Residual_k1[3];
      kappa_k2[1] = alpha * kappa_k2[0];
      delta_lambda_k2 = delta_lambda_k1 - delta * Residual_k1[4];
      Lambda_k2 = Lambda_n + delta_lambda_k2;
      Iter_k2 = 0;

      if (Activate_CutOff) {
        if ((T_k2[0] > CutOff) && (T_k2[1] > CutOff) && (T_k2[2] > CutOff)) {
          T_k2[0] = T_k2[1] = T_k2[2] = CutOff;
          Lambda_k2 = Lambda_n;
          kappa_k2[0] = kappa_n[0];
          kappa_k2[1] = alpha * kappa_k2[0];
          delta_lambda_k2 = 0.0;
          break;
        }
      }

      I1 = T_k2[0] + T_k2[1] + T_k2[2];
      I2 = T_k2[0] * T_k2[1] + T_k2[1] * T_k2[2] + T_k2[0] * T_k2[2];
      I3 = T_k2[0] * T_k2[1] * T_k2[2];

      if (Lambda_k2 < 0.0) {
        fprintf(stderr, "" RED "Negative value of Lambda: %f " RESET "\n",
                Lambda_k2);
        return EXIT_FAILURE;
      }

      if (I1 > 0.0) {
        fprintf(stderr, "" RED "Positive value of I1: %f " RESET "\n", I1);
        return EXIT_FAILURE;
      }

      // Compute the residual of the next step
      __E_hencky(E_hencky_k2, T_k2, CC, c_cotphi);

      __kappa(kappa_hat, a, Lambda_k2, I1, alpha);

      __d_G_d_stress(d_G_d_stress, T_k2, I1, I2, I3, kappa_k2[1]);

      __F(&F_k2, kappa_k2[0], I1, I2, I3);

      STATUS = __residual(Residual_k2, &Norm_Residual_k2, E_hencky_trial,
                          E_hencky_k2, d_G_d_stress, kappa_k2, kappa_hat, F_k2,
                          delta_lambda_k2);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __residual" RESET "\n");
        return EXIT_FAILURE;
      }

      while ((Norm_Residual_k2 - Norm_Residual_k1) > TOL) {
        delta =
            pow(delta, 2.0) * 0.5 * Norm_Residual_k1 /
            (Norm_Residual_k2 - delta * Norm_Residual_k1 + Norm_Residual_k1);

        if (delta < TOL)
          break;

        T_k2[0] = T_k1[0] - delta * Residual_k2[0];
        T_k2[1] = T_k1[1] - delta * Residual_k2[1];
        T_k2[2] = T_k1[2] - delta * Residual_k2[2];
        kappa_k2[0] = kappa_k1[0] - delta * Residual_k2[3];
        kappa_k2[1] = alpha * kappa_k2[0];
        delta_lambda_k2 = delta_lambda_k1 - delta * Residual_k2[4];
        Lambda_k2 = Lambda_n + delta_lambda_k2;

        if (Activate_CutOff) {
          if ((T_k2[0] > CutOff) && (T_k2[1] > CutOff) && (T_k2[2] > CutOff)) {
            T_k2[0] = T_k2[1] = T_k2[2] = CutOff;
            Lambda_k2 = Lambda_n;
            kappa_k2[0] = kappa_n[0];
            kappa_k2[1] = alpha * kappa_k2[0];
            delta_lambda_k2 = 0.0;
            break;
          }
        }

        I1 = T_k2[0] + T_k2[1] + T_k2[2];
        I2 = T_k2[0] * T_k2[1] + T_k2[1] * T_k2[2] + T_k2[0] * T_k2[2];
        I3 = T_k2[0] * T_k2[1] * T_k2[2];

        if (Lambda_k2 < 0.0) {
          fprintf(stderr,
                  "" RED "Negative value of Lambda (line search): %f " RESET
                  "\n",
                  Lambda_k2);
          return EXIT_FAILURE;
        }

        __E_hencky(E_hencky_k2, T_k2, CC, c_cotphi);

        __kappa(kappa_hat, a, Lambda_k2, I1, alpha);

        __d_G_d_stress(d_G_d_stress, T_k2, I1, I2, I3, kappa_k2[1]);

        __F(&F_k2, kappa_k2[0], I1, I2, I3);

        STATUS = __residual(Residual_k2, &Norm_Residual_k2, E_hencky_trial,
                            E_hencky_k2, d_G_d_stress, kappa_k2, kappa_hat,
                            F_k2, delta_lambda_k2);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr,
                  "" RED "Error in __residual (line search loop)" RESET "\n");
          return EXIT_FAILURE;
        }

        Iter_k2++;

        if (Iter_k2 == MaxIter_k2) {
          break;
        }
      }

      // Update variables for the next evaluation of the residual
      T_k1[0] = T_k2[0];
      T_k1[1] = T_k2[1];
      T_k1[2] = T_k2[2];
      E_hencky_k1[0] = E_hencky_k2[0];
      E_hencky_k1[1] = E_hencky_k2[1];
      E_hencky_k1[2] = E_hencky_k2[2];
      kappa_k1[0] = kappa_k2[0];
      kappa_k1[1] = kappa_k2[1];
      Lambda_k1 = Lambda_k2;
      F_k1 = F_k2;
      delta_lambda_k1 = delta_lambda_k2;
      Residual_k1[0] = Residual_k2[0];
      Residual_k1[1] = Residual_k2[1];
      Residual_k1[2] = Residual_k2[2];
      Residual_k1[3] = Residual_k2[3];
      Residual_k1[4] = Residual_k2[4];
      Norm_Residual_k1 = Norm_Residual_k2;
      Iter_k1++;

      if (Iter_k1 == MaxIter_k1) {
        STATUS = __reciprocal_condition_number(&rcond, Tangent_Matrix);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr,
                  "" RED "Error in __reciprocal_condition_number" RESET "\n");
          return EXIT_FAILURE;
        }

        if (rcond < 1E-10) {
          fprintf(stderr,
                  "" RED "Reciprocal condition number below 1E-10: %e" RESET
                  "\n",
                  rcond);
          return EXIT_FAILURE;
        }
        break;
      }
    }

    /*
      Update equivalent plastic strain and increment of plastic deformation
    */
    IO_State.Stress[0] = T_k1[0] + c_cotphi;
    IO_State.Stress[1] = T_k1[1] + c_cotphi;
    IO_State.Stress[2] = T_k1[2] + c_cotphi;
    IO_State.Strain_e[0] = E_hencky_k1[0];
    IO_State.Strain_e[1] = E_hencky_k1[1];
    IO_State.Strain_e[2] = E_hencky_k1[2];
    *IO_State.EPS = Lambda_k1;
    *IO_State.Kappa = kappa_k1[0];

    if (Iter_k1 == MaxIter_k1) {
      *IO_State.EPS = Lambda_n;
      *IO_State.Kappa = kappa_n[0];
    }
  }

  printf("Iter: %i\n", Iter_k1);

  return EXIT_SUCCESS;
}

/**************************************************************/

static void __elastic_tangent(double *CC, double *AA, double E, double nu,
                              double Lame, double G) {
  CC[0] = 1.0 / E;
  CC[1] = -nu / E;
  CC[2] = -nu / E;

  CC[3] = -nu / E;
  CC[4] = 1.0 / E;
  CC[5] = -nu / E;

  CC[6] = -nu / E;
  CC[7] = -nu / E;
  CC[8] = 1.0 / E;

  AA[0] = Lame + 2 * G;
  AA[1] = Lame;
  AA[2] = Lame;

  AA[3] = Lame;
  AA[4] = Lame + 2 * G;
  AA[5] = Lame;

  AA[6] = Lame;
  AA[7] = Lame;
  AA[8] = Lame + 2 * G;
}

/**************************************************************/

static void __trial_elastic(double *T_tr, const double *E_hencky_trial,
                            const double *AA, double c_cotphi) {

  T_tr[0] = AA[0] * E_hencky_trial[0] + AA[1] * E_hencky_trial[1] +
            AA[2] * E_hencky_trial[2] - c_cotphi;
  T_tr[1] = AA[3] * E_hencky_trial[0] + AA[4] * E_hencky_trial[1] +
            AA[5] * E_hencky_trial[2] - c_cotphi;
  T_tr[2] = AA[6] * E_hencky_trial[0] + AA[7] * E_hencky_trial[1] +
            AA[8] * E_hencky_trial[2] - c_cotphi;
}

/**************************************************************/

static void __E_hencky(double *E_hencky_k, const double *T_k, const double *CC,
                       double c_cotphi) {

  E_hencky_k[0] = CC[0] * (T_k[0] + c_cotphi) + CC[1] * (T_k[1] + c_cotphi) +
                  CC[2] * (T_k[2] + c_cotphi);
  E_hencky_k[1] = CC[3] * (T_k[0] + c_cotphi) + CC[4] * (T_k[1] + c_cotphi) +
                  CC[5] * (T_k[2] + c_cotphi);
  E_hencky_k[2] = CC[6] * (T_k[0] + c_cotphi) + CC[7] * (T_k[1] + c_cotphi) +
                  CC[8] * (T_k[2] + c_cotphi);
}

/**************************************************************/

static void __kappa(double *kappa, const double *a, double Lambda, double I1,
                    double alpha) {

  kappa[0] = a[0] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  kappa[1] = alpha * kappa[0];
}

/**************************************************************/

static void __d_kappa_phi_d_stress(double *d_kappa_phi_d_stress,
                                   const double *a, double Lambda, double I1) {

  d_kappa_phi_d_stress[0] =
      a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  d_kappa_phi_d_stress[1] =
      a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  d_kappa_phi_d_stress[2] =
      a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
}

/**************************************************************/

static void __d_kappa_phi_d_lambda(double *d_kappa_phi_d_lambda,
                                   const double *a, double Lambda, double I1) {

  *d_kappa_phi_d_lambda =
      (1 - a[2] * Lambda) * a[0] * exp(a[1] * I1) * exp(-a[2] * Lambda);
}

/**************************************************************/

static double __F(double *F, double kappa_phi, double I1, double I2,
                  double I3) {

  double K1 = 9.0 + kappa_phi;

  *F = cbrt(K1 * I3) - cbrt(I1 * I2);
}

/**************************************************************/

static void __d_F_d_stress(double *d_F_d_stress, const double *T_k, double I1,
                           double I2, double I3, double kappa_phi) {
  double Grad_f[3];

  double K1 = 9.0 + kappa_phi;

  Grad_f[0] = (I1 * (I1 - T_k[0]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2.0));
  Grad_f[1] = (I1 * (I1 - T_k[1]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2.0));
  Grad_f[2] = (I1 * (I1 - T_k[2]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2.0));

  d_F_d_stress[0] = cbrt(K1 * I3) / (3.0 * T_k[0]) - Grad_f[0];
  d_F_d_stress[1] = cbrt(K1 * I3) / (3.0 * T_k[1]) - Grad_f[1];
  d_F_d_stress[2] = cbrt(K1 * I3) / (3.0 * T_k[2]) - Grad_f[2];
}

/**************************************************************/

static void __d_F_d_kappa_phi(double *d_F_d_kappa_phi, double I1, double I3,
                              double kappa_phi) {

  double K1 = 9.0 + kappa_phi;

  *d_F_d_kappa_phi = (1.0 / 3.0) * pow(cbrt(K1), -2.0) * cbrt(I3);
}

/**************************************************************/

static void __d_G_d_stress(double *d_G_d_stress, const double *T_k, double I1,
                           double I2, double I3, double kappa_psi) {

  double Grad_g[3] = {0.0, 0.0, 0.0};
  double K2 = 9.0 + kappa_psi;

  Grad_g[0] = (I1 * (I1 - T_k[0]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2.0));
  Grad_g[1] = (I1 * (I1 - T_k[1]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2.0));
  Grad_g[2] = (I1 * (I1 - T_k[2]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2.0));

  d_G_d_stress[0] = cbrt(K2 * I3) / (3.0 * T_k[0]) - Grad_g[0];
  d_G_d_stress[1] = cbrt(K2 * I3) / (3.0 * T_k[1]) - Grad_g[1];
  d_G_d_stress[2] = cbrt(K2 * I3) / (3.0 * T_k[2]) - Grad_g[2];
}

/**************************************************************/

static void __dd_G_dd_stress(double *dd_G_dd_stress, const double *T_k,
                             double kappa_psi, double I1, double I2,
                             double I3) {

  double K2 = 9.0 + kappa_psi;

  double d_g_d_stress[3];
  d_g_d_stress[0] = (I1 * (I1 - T_k[0]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2.0));
  d_g_d_stress[1] = (I1 * (I1 - T_k[1]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2.0));
  d_g_d_stress[2] = (I1 * (I1 - T_k[2]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2.0));

  double dd_g_dd_stress[9];
  for (unsigned A = 0; A < 3; A++) {
    for (unsigned B = 0; B < 3; B++) {
      dd_g_dd_stress[A * 3 + B] =
          pow(cbrt(I1 * I2), -2.0) / 3.0 *
              (3.0 * I1 - T_k[A] - T_k[B] - I1 * (A == B)) -
          (2.0 / cbrt(I1 * I2)) * d_g_d_stress[A] * d_g_d_stress[B];
    }
  }

  for (unsigned A = 0; A < 3; A++) {
    for (unsigned B = 0; B < 3; B++) {
      dd_G_dd_stress[A * 3 + B] = (1.0 / 3.0) * cbrt(K2 * I3) *
                                      (1.0 / (3.0 * T_k[A] * T_k[B]) -
                                       1.0 * (A == B) / pow(T_k[A], 2.0)) -
                                  dd_g_dd_stress[A * 3 + B];
    }
  }
}
/**************************************************************/

static void __dd_G_d_stress_d_kappa_psi(double *dd_G_d_stress_d_kappa_psi,
                                        const double *T_k, double I1, double I3,
                                        double kappa_psi) {

  double K2 = 9.0 + kappa_psi;

  dd_G_d_stress_d_kappa_psi[0] =
      (cbrt(I3) / (3.0 * T_k[0])) / (3.0 * pow(cbrt(K2), 2));
  dd_G_d_stress_d_kappa_psi[1] =
      (cbrt(I3) / (3.0 * T_k[1])) / (3.0 * pow(cbrt(K2), 2));
  dd_G_d_stress_d_kappa_psi[2] =
      (cbrt(I3) / (3.0 * T_k[2])) / (3.0 * pow(cbrt(K2), 2));
}

/**************************************************************/

static int __residual(double *Residual, double *Error_k,
                      const double *E_hencky_trial, const double *E_hencky_k,
                      const double *d_G_d_stress, const double *kappa_k,
                      const double *kappa_hat, double F_k,
                      double delta_lambda_k) {

  Residual[0] =
      E_hencky_k[0] - E_hencky_trial[0] + delta_lambda_k * d_G_d_stress[0];
  Residual[1] =
      E_hencky_k[1] - E_hencky_trial[1] + delta_lambda_k * d_G_d_stress[1];
  Residual[2] =
      E_hencky_k[2] - E_hencky_trial[2] + delta_lambda_k * d_G_d_stress[2];
  Residual[3] = kappa_k[0] - kappa_hat[0];
  Residual[4] = F_k;

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  printf("*********************************\n");
  printf("__residual(): \n");
  printf("d_lambda_k: %e \n", delta_lambda_k);
  printf("F_k: %e \n", F_k);
  printf("E_hencky_tr: [%e, %e, %e] \n", E_hencky_trial[0], E_hencky_trial[1],
         E_hencky_trial[2]);
  printf("E_hencky_k: [%e, %e, %e] \n", E_hencky_k[0], E_hencky_k[1],
         E_hencky_k[2]);
  printf("Kappa_phi: %e, kappa_phi_hat: %e \n", kappa_k[0], kappa_hat[0]);
  printf("Plastic flow: [%e, %e, %e] \n", d_G_d_stress[0], d_G_d_stress[1],
         d_G_d_stress[2]);
  printf("*********************************\n");
#endif
#endif

  /*
    Compute absolute error from the residual
  */
  double norm_R = 0.0;
  for (unsigned A = 0; A < 5; A++) {
    norm_R += DSQR(Residual[A]);
  }

  *Error_k = pow(norm_R, 0.5);

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __reciprocal_condition_number(double *RCOND, double *Tangent_Matrix)
/*
  C = rcond(Tangent_Matrix) returns an estimate for the reciprocal condition of
  Tangent_Matrix in 1-norm.
*/
{

  double ANORM;
  int INFO;
  int N_rows = 5;
  int N_cols = 5;
  int LDA = 5;
  double WORK_ANORM[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double WORK_RCOND[20];
  int IPIV[5] = {0, 0, 0, 0, 0};
  int IWORK_RCOND[5] = {0, 0, 0, 0, 0};

  // Compute 1-norm
  ANORM = dlange_("1", &N_rows, &N_cols, Tangent_Matrix, &LDA, WORK_ANORM);

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N_rows, &N_cols, Tangent_Matrix, &LDA, IPIV, &INFO);

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

  // Compute the Reciprocal condition number
  dgecon_("1", &N_rows, Tangent_Matrix, &LDA, &ANORM, RCOND, WORK_RCOND,
          IWORK_RCOND, &INFO);

  if (INFO < 0) {
    fprintf(
        stderr,
        "" RED
        "Error in dgecon_() : the %i-th argument had an illegal value " RESET
        "\n",
        abs(INFO));
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __solver(double *Tangent_Matrix, double *Residual) {
  int Order = 5;
  int LDA = 5;
  int LDB = 5;
  char TRANS = 'T'; /* (Transpose) */
  int INFO = 3;
  int IPIV[5] = {0, 0, 0, 0, 0};
  int NRHS = 1;

  /*
    Compute the LU factorization
  */
  dgetrf_(&Order, &Order, Tangent_Matrix, &LDA, IPIV, &INFO);
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

  /*
    Solve
  */
  dgetrs_(&TRANS, &Order, &NRHS, Tangent_Matrix, &LDA, IPIV, Residual, &LDB,
          &INFO);
  if (INFO) {
    fprintf(stderr, "" RED "Error in dgetrs_() " RESET "\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/*************************************************************/