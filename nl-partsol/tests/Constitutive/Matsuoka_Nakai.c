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

/*
  Auxiliar functions
*/
static int __compute_trial_b_e(
    double *eigval_b_e_tr /**< [out] Eigenvalues of b elastic trial. */,
    double *eigvec_b_e_tr /**< [out] Eigenvector of b elastic trial. */,
    const double *b_e /**< [in] (n) Elastic left Cauchy-Green.*/,
    const double *d_phi /**< [in] Incremental deformation gradient. */);

static int __corrector_b_e(
    double *b_e /**< [out] (n+1) Elastic deformation gradient. */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */,
    const double *E_hencky_trial /**< [in] Corrected Henky strain */);

static int __elastic_tangent(double *CC /**< [out] Elastic compliance */,
                             double *AA /**< [out] Elastic matrix */,
                             double E /**< [in] Young modulus */,
                             double nu /**< [in] Poisson ratio */,
                             double K /**< [in] Lamé parameter */,
                             double G /**< [in] Shear modulus */);

static int
__trial_elastic(double *T_tr /**< [out] Trial elastic stress tensor*/,
                const double *E_hencky_trial /**< [in] Henky strain (trial) */,
                const double *AA /**< [in] Elastic matrix */);

static int
__E_hencky(double *E_hencky_k /**< [out] Henky strain (iter k) */,
           const double *T_k /**< [in] Local stress tensor (iter k) */,
           const double *CC /**< [in] Elastic compliance */);

static int __update_internal_variables_elastic(
    double *Stress /**< [in/out] Nominal stress tensor */,
    const double *D_phi /**< [in] Total deformation gradient. */,
    const double *T_tr /**< [in] Elastic stress tensor */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */);

static int
__kappa(double *kappa /**< [out] Hardening vector */,
        const double *a /**< [in] Vector with fit parameters (hardening) */,
        double Lambda /**< [in] Total plastic multiplier */,
        double I1 /**< [in] First invariant of the stress tensor */,
        double alpha /**< [in] Dilatance parameter*/);

static int __d_kappa_phi_d_stress(
    double *d_kappa_phi_d_stress /**< [out] Stress derivative of kappa[0] */,
    const double *a /**< [in] Vector with fit parameters (hardening) */,
    double Lambda /**< [in] Total plastic multiplier */,
    double I1 /**< [in] First invariant of the stress tensor */);

static int __d_kappa_phi_d_lambda(
    double *d_kappa_phi_d_lambda /**< [out] Lambda derivative of kappa[0] */,
    const double *a /**< [in] Vector with fit parameters (hardening) */,
    double Lambda /**< [in] Total plastic multiplier */,
    double I1 /**< [in] First invariant of the stress tensor */);

static double __F(double c0 /**< [in] Yield function fit parameter */,
                  double kappa_phi /**< [in] Friction angle hardening */,
                  double pa /**< [in] Atmospheric pressure */,
                  double I1 /**< [in] First invariant of the stress tensor */,
                  double I2 /**< [in] Second invariant of the stress tensor */,
                  double I3 /**< [in] Third invariant of the stress tensor */,
                  double m /**< [in] Yield function fit parameter */);

static int __d_F_d_stress(
    double *d_F_d_stress /**< [out] Yield function derivative (stress) */,
    const double *T_k /**< [in] Local stress tensor (iter k) */,
    double I1 /**< [in] First invariant of the stress tensor */,
    double I2 /**< [in] Second invariant of the stress tensor */,
    double I3 /**< [in] Third invariant of the stress tensor */,
    double c0 /**< [in] Yield function fit parameter */,
    double kappa_phi /**< [in] Friction angle hardening */,
    double pa /**< [in] Atmospheric pressure */,
    double m /**< [in] Yield function fit parameter */);

static int __d_F_d_kappa_phi(
    double *d_F_d_kappa_phi /**< [out] Yield function derivative (kappa[0]) */,
    double I1 /**< [in] First invariant of the stress tensor */,
    double I3 /**< [in] Third invariant of the stress tensor */,
    double c0 /**< [in] Yield function fit parameter */,
    double m /**< [in] Yield function fit parameter */,
    double pa /**< [in] Atmospheric pressure */,
    double kappa_phi /**< [in] Friction angle hardening */);

static double __G(double c0 /**< [in] Yield function fit parameter */,
                  double kappa_psi /**< [in] Dilatance angle hardening */,
                  double pa /**< [in] Atmospheric pressure */,
                  double I1 /**< [in] First invariant of the stress tensor */,
                  double I2 /**< [in] Second invariant of the stress tensor */,
                  double I3 /**< [in] Third invariant of the stress tensor */,
                  double m /**< [in] Yield function fit parameter */);

static int
__d_G_d_T(double *d_G_d_T /**< [out] Plastic potential function
                                       derivative (stress) */
          ,
          const double *T_k /**< [in] Local stress tensor (iter k) */,
          double I1 /**< [in] First invariant of the stress tensor */,
          double I2 /**< [in] Second invariant of the stress tensor */,
          double I3 /**< [in] Third invariant of the stress tensor */,
          double c0 /**< [in] Yield function fit parameter */,
          double kappa_psi /**< [in] Dilatance angle hardening */,
          double pa /**< [in] Atmospheric pressure */,
          double m /**< [in] Yield function fit parameter */);

static int
__dd_G_dd_T(double *dd_G_dd_T /**< [out] Plastic potential hessian (stress) */,
            const double *T_k /**< [in] Local stress tensor (iter k) */,
            double kappa_psi /**< [in] Dilatance angle hardening */,
            double I1 /**< [in] First invariant of the stress tensor */,
            double I2 /**< [in] Second invariant of the stress tensor */,
            double I3 /**< [in] Third invariant of the stress tensor */,
            double m /**< [in] Yield function fit parameter */,
            double pa /**< [in] Atmospheric pressure */,
            double c0 /**< [in] Yield function fit parameter */);

static int __dd_G_d_T_d_kappa_psi(
    double *dd_G_d_T_d_kappa_psi /**< [out] Plastic potential deriv */,
    const double *T_k /**< [in] Local stress tensor (iter k) */,
    double I1 /**< [in] First invariant of the stress tensor */,
    double I3 /**< [in] Third invariant of the stress tensor */,
    double m /**< [in] Yield function fit parameter */,
    double pa /**< [in] Atmospheric pressure */,
    double c0 /**< [in] Yield function fit parameter */,
    double kappa_psi /**< [in] Dilatance angle hardening */);

static int __residual(
    double *Residual /**< [out] Residual of the problem */,
    double *Error_k /**< [out] Norm of the residual */,
    const double *E_hencky_trial /**< [in] Henky strain (trial) */,
    const double *E_hencky_k /**< [in] Henky strain (iter k) */,
    const double
        *d_G_d_T /**< [in] Plastic potential function derivative (stress) */
    ,
    const double *kappa_k /**< [in] Hardening vector (iter k) */,
    const double *kappa_hat /**< [in] Hardening vector (eval) */,
    double F_k /**< [in] Yield function evaluation (iter k) */,
    double delta_lambda_k /**< [in] Discrete plastic multiplier (iter k) */);

static int __tangent_matrix(
    double *Tangent_Matrix /**< [out] Tangent matrix of the problem */,
    const double *CC /**< [in] Elastic compliance */,
    const double *d_F_d_stress /**< [in] Yield gradient (stress) */,
    double d_F_d_kappa_phi /**< [in] Yield gradient (kappa-phi) */,
    const double *d_G_d_T /**< [in] Plastic potential gradient (stress) */,
    const double *dd_G_dd_T /**< [in] Plastic potential hessian (stress) */,
    const double *dd_G_d_T_d_kappa2 /**< [in] Plastic potential hessian
                                            (stress-kappa) */
    ,
    const double
        *d_kappa_phi_d_stress /**< [in] Hardening friction gradient (stress) */,
    double
        d_kappa_phi_d_lambda /**< [in] Hardening friction gradient (lambda) */,
    double alpha /**< [in] Dilatance parameter*/,
    double delta_lambda_k /**< [in] Discrete plastic multiplier (iter k) */);

static int __update_internal_variables_plastic(
    double *Increment_E_plastic /**< [out] Plastic corrector */,
    double *Stress /**< [out] Nominal stress tensor */,
    double *eps_n1 /**< [out] Equivalent plastic strain */,
    double *kappa_n1 /**< [out] Friction angle hardening */,
    const double *D_phi /**< [in] Total deformation gradient. */,
    const double *T_tr_k /**< [in] Stress tensor (iter k). */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */,
    const double
        *d_G_d_T /**< [in] Plastic potential function derivative (stress) */
    ,
    double Lambda_k /**< [in] Total plastic multiplier (iter k) */,
    double delta_lambda_k /**< [in] Discrete plastic multiplier (iter k) */,
    double kappa_phi_k /**< [out] Friction angle hardening (iter k)*/);

static int
__solver(double *Tangent_Matrix /**< [in/out] Tangent matrix of the problem */,
         double *Residual /**< [in/out] Residual of the problem */);

static int __condition_number(
    double *RCOND /**< [out] Condition number of the tangent matrix */,
    double *Tangent_Matrix /**< [in/out] Tangent matrix of the problem */);

int compute_1PK_Matsuoka_Nakai(State_Parameters IO_State, Material MatProp);

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
#define PoissonRatio 0.2
#define AtmosphericPressure -100
#define m_Parameter 0.0
#define c0_Parameter 9.0
#define a1_Parameter 20000
#define a2_Parameter 0.005
#define a3_Parameter 35
#define alpha_Parameter 0.5
#define NumberSteps 30 // 5000 // 3500 // 2500 //
#define Delta_strain_II -0.00001
#define dF_yy 0.9999
#define Confining_pressure -200

int main() {
  int STATUS = EXIT_SUCCESS;
  State_Parameters IO_State;
  Material MatProp;

  // Define material variables
  MatProp.E = YoungMouduls;
  MatProp.nu = PoissonRatio;
  MatProp.atmospheric_pressure = AtmosphericPressure;
  MatProp.m_Frictional = m_Parameter;
  MatProp.c0_Frictional = c0_Parameter;
  MatProp.a_Hardening_Borja[0] = a1_Parameter;
  MatProp.a_Hardening_Borja[1] = a2_Parameter;
  MatProp.a_Hardening_Borja[2] = a3_Parameter;
  MatProp.alpha_Hardening_Borja = alpha_Parameter;

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

  D_phi[0] = 1.0;
  D_phi[3] = 1.0;
  D_phi[4] = 1.0;

  b_e[0] = 1.0;
  b_e[3] = 1.0;
  b_e[4] = 1.0;

  // Start time integration
  for (int i = 1; i < NumberSteps; i++) {

    printf("*********************\n");
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
    STATUS = compute_1PK_Matsuoka_Nakai(IO_State, MatProp);

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

int compute_1PK_Matsuoka_Nakai(State_Parameters IO_State, Material MatProp)
/*
  Monolithic algorithm for a smooth Mohr-Coulomb plastic criterium
*/
{

  int STATUS = EXIT_SUCCESS;

  // Read input/output parameters
  double eigval_b_e_tr[3] = {0.0, 0.0, 0.0};
  double eigvec_b_e_tr[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double E_hencky_trial[3] = {0.0, 0.0, 0.0};
  double E_hencky_k1[3] = {0.0, 0.0, 0.0};
  double E_hencky_k2[3] = {0.0, 0.0, 0.0};

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
  double E = MatProp.E;
  double nu = MatProp.nu;
  double K = E / (3.0 * (1.0 - 2.0 * nu));
  double G = E / (2.0 * (1.0 + nu));
  double alpha = MatProp.alpha_Hardening_Borja;
  double pa = MatProp.atmospheric_pressure;
  double rad_friction_angle =
      (PI__MatrixLib__ / 180.0) * MatProp.phi_Frictional;
  double m = 0.0;
  double c0 = 9.0;
  double a[3] = {0.0, 0.0, 0.0};
  a[0] = MatProp.a_Hardening_Borja[0];
  a[1] = MatProp.a_Hardening_Borja[1];
  a[2] = MatProp.a_Hardening_Borja[2];

  double CC[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double AA[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  __elastic_tangent(CC, AA, E, nu, K, G);

  // Define scalar variables internal variables
  double F_k1, F_k2, F_0 = 0.0;
  double I1, I2, I3 = 0.0;
  double Lambda_n = *IO_State.Equiv_Plast_Str;
  double Lambda_k1;
  double Lambda_k2;
  double delta_lambda_k0 = 0.0;
  double delta_lambda_k1;
  double delta_lambda_k2;

  // Define tensorial internal variables
  double Increment_E_plastic[3] = {0.0, 0.0, 0.0};
  double T_tr[3] = {0.0, 0.0, 0.0}; // Trial elastic stress
  double T_k1[3];                   // Stress iteration k
  double T_k2[3];                   // Stress iteration k (line search)
  double kappa_n[2];
  kappa_n[0] = (*IO_State.Kappa);
  kappa_n[1] = alpha * (*IO_State.Kappa);
  double kappa_k1[2]; // Hardening iteration k
  double kappa_k2[2]; // Hardening iteration k (line search)
  double kappa_hat[3] = {0.0, 0.0, 0.0};
  double d_G_d_T[3] = {0.0, 0.0, 0.0}; // Plastic flow
  double dd_G_dd_T[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dd_G_d_T_d_kappa_psi[3] = {0.0, 0.0, 0.0};
  double d_kappa_phi_d_stress[3] = {0.0, 0.0, 0.0};
  double d_kappa_phi_d_lambda = 0.0;
  double d_F_d_stress[3] = {0.0, 0.0, 0.0};
  double d_F_d_kappa_phi = 0;

  /*
    Initialize Newton-Raphson solver
  */
  double TOL = TOL_Radial_Returning;
  double Residual_k1[5] = {0.0, 0.0, 0.0};
  double Residual_k2[5] = {0.0, 0.0, 0.0};
  double Tangent_Matrix[25];
  double rcond;
  double Norm_Residual_1;
  double Norm_Residual_2;
  double delta = 1;
  int MaxIter_k1 = Max_Iterations_Radial_Returning;
  int MaxIter_k2 = 10 * Max_Iterations_Radial_Returning;
  int Iter_k1;
  int Iter_k2;

  __trial_elastic(T_tr, E_hencky_trial, AA);

  I1 = T_tr[0] + T_tr[1] + T_tr[2];
  I2 = T_tr[0] * T_tr[1] + T_tr[1] * T_tr[2] + T_tr[0] * T_tr[2];
  I3 = T_tr[0] * T_tr[1] * T_tr[2];

  F_k1 = F_0 = __F(c0, kappa_k1[0], pa, I1, I2, I3, m);

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
  printf("Initial value of the yield function: %f \n", F_0);
  printf("T trial: [%f, %f, %f] \n", T_tr[0], T_tr[1], T_tr[2]);
#endif
#endif

  // Elastic
  if (F_0 <= 0.0) {

    STATUS = __update_internal_variables_elastic(
        IO_State.Stress, IO_State.D_phi, T_tr, eigvec_b_e_tr);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in __update_internal_variables_elastic()" RESET
              "\n");
      return EXIT_FAILURE;
    }

  }
  // Plastic (monolithic solver with line search)
  else {

    // Compute trial plastic flow rule
    STATUS = __d_G_d_T(d_G_d_T, T_tr, I1, I2, I3, c0, kappa_n[1], pa, m);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __d_G_d_T" RESET "\n");
      return EXIT_FAILURE;
    }

    // Compute trial residual
    STATUS = __residual(Residual_k1, &Norm_Residual_1, E_hencky_trial,
                        E_hencky_trial, d_G_d_T, kappa_n, kappa_n, F_0,
                        delta_lambda_k0);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __residual" RESET "\n");
      return EXIT_FAILURE;
    }

    // Assign values to the k iteration
    T_k1[0] = T_tr[0];
    T_k1[1] = T_tr[1];
    T_k1[2] = T_tr[2];
    kappa_k1[0] = kappa_n[0];
    kappa_k1[1] = kappa_n[1];
    delta_lambda_k1 = 0.0;
    Lambda_k1 = Lambda_n;
    Iter_k1 = 0;

    while (fabs(F_k1 / F_0) >= TOL) {

      delta = 1;

      // Evaluate hardening derivatives
      STATUS = __d_kappa_phi_d_stress(d_kappa_phi_d_stress, a, Lambda_k1, I1);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __d_kappa_phi_d_stress" RESET "\n");
        return EXIT_FAILURE;
      }
      STATUS = __d_kappa_phi_d_lambda(&d_kappa_phi_d_lambda, a, Lambda_k1, I1);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __d_kappa_phi_d_lambda" RESET "\n");
        return EXIT_FAILURE;
      }

      // Evaluate yield function derivatives
      STATUS = __d_F_d_stress(d_F_d_stress, T_k1, I1, I2, I3, c0, kappa_k1[0],
                              pa, m);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __d_F_d_stress" RESET "\n");
        return EXIT_FAILURE;
      }
      STATUS =
          __d_F_d_kappa_phi(&d_F_d_kappa_phi, I1, I3, c0, m, pa, kappa_k1[0]);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __d_F_d_kappa_phi" RESET "\n");
        return EXIT_FAILURE;
      }

      // Evaluate plastic flow rule derivatives
      STATUS = __dd_G_dd_T(dd_G_dd_T, T_k1, kappa_k1[1], I1, I2, I3, m, pa, c0);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __dd_G_dd_T" RESET "\n");
        return EXIT_FAILURE;
      }
      STATUS = __dd_G_d_T_d_kappa_psi(dd_G_d_T_d_kappa_psi, T_k1, I1, I3, m, pa,
                                      c0, kappa_k1[1]);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __dd_G_d_T_d_kappa_psi" RESET "\n");
        return EXIT_FAILURE;
      }

      // Assemble tangent matrix
      STATUS = __tangent_matrix(Tangent_Matrix, CC, d_F_d_stress,
                                d_F_d_kappa_phi, d_G_d_T, dd_G_dd_T,
                                dd_G_d_T_d_kappa_psi, d_kappa_phi_d_stress,
                                d_kappa_phi_d_lambda, alpha, delta_lambda_k1);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __tangent_matrix" RESET "\n");
        return EXIT_FAILURE;
      }

      // Compute increments and update variables
      __solver(Tangent_Matrix, Residual_k1);

      // Update values for the next step (line search)
      T_k2[0] = T_k1[0] - delta * Residual_k1[0];
      T_k2[1] = T_k1[1] - delta * Residual_k1[1];
      T_k2[2] = T_k1[2] - delta * Residual_k1[2];
      kappa_k2[0] = kappa_k1[0] - delta * Residual_k1[3];
      kappa_k2[1] = alpha * kappa_k2[0];
      delta_lambda_k2 = delta_lambda_k1 - delta * Residual_k1[4];
      Lambda_k2 = Lambda_n + delta_lambda_k2;
      Iter_k2 = 0;

      I1 = T_k2[0] + T_k2[1] + T_k2[2];
      I2 = T_k2[0] * T_k2[1] + T_k2[1] * T_k2[2] + T_k2[0] * T_k2[2];
      I3 = T_k2[0] * T_k2[1] * T_k2[2];

      // Compute the residual of the next step (line search)
      STATUS = __E_hencky(E_hencky_k2, T_k2, CC);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __E_hencky (line search loop)" RESET "\n");
        return EXIT_FAILURE;
      }

      STATUS = __kappa(kappa_hat, a, Lambda_k2, I1, alpha);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __kappa (line search loop)" RESET "\n");
        return EXIT_FAILURE;
      }

      STATUS = __d_G_d_T(d_G_d_T, T_k2, I1, I2, I3, c0, kappa_k2[1], pa, m);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __d_G_d_T (line search loop)" RESET "\n");
        return EXIT_FAILURE;
      }

      F_k2 = __F(c0, kappa_k2[0], pa, I1, I2, I3, m);

      STATUS =
          __residual(Residual_k2, &Norm_Residual_2, E_hencky_trial, E_hencky_k2,
                     d_G_d_T, kappa_k2, kappa_hat, F_k2, delta_lambda_k2);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __residual (line search loop)" RESET "\n");
        return EXIT_FAILURE;
      }

      while ((Norm_Residual_2 - Norm_Residual_1) > TOL_Radial_Returning) {

        delta = pow(delta, 2) * 0.5 * Norm_Residual_1 /
                (Norm_Residual_2 - delta * Norm_Residual_1 + Norm_Residual_1);

        if (delta < TOL_Radial_Returning)
          break;

        T_k2[0] = T_k1[0] - delta * Residual_k1[0];
        T_k2[1] = T_k1[1] - delta * Residual_k1[1];
        T_k2[2] = T_k1[2] - delta * Residual_k1[2];
        kappa_k2[0] = kappa_k1[0] - delta * Residual_k1[3];
        kappa_k2[1] = alpha * kappa_k2[0];
        delta_lambda_k2 = delta_lambda_k1 - delta * Residual_k1[4];
        Lambda_k2 = Lambda_n + delta_lambda_k2;

        Iter_k2++;

        if (Iter_k2 == MaxIter_k2) {
          break;
        }

        STATUS = __residual(Residual_k2, &Norm_Residual_2, E_hencky_trial,
                            E_hencky_k2, d_G_d_T, kappa_k2, kappa_hat, F_k2,
                            delta_lambda_k2);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr,
                  "" RED "Error in __residual (line search loop)" RESET "\n");
          return EXIT_FAILURE;
        }
      }

      // Update variables for the next evaluation of the residual
      T_k1[0] = T_k2[0];
      T_k1[1] = T_k2[1];
      T_k1[2] = T_k2[2];
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
      Norm_Residual_1 = Norm_Residual_2;
      Iter_k1++;

      if (Iter_k1 == MaxIter_k1) {
        STATUS = __condition_number(&rcond, Tangent_Matrix);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr, "" RED "Error in __condition_number" RESET "\n");
          return EXIT_FAILURE;
        }

        if (rcond < 1E-10) {
          fprintf(stderr, "%s: %s, %s: %e\n", "Error in Frictional_Monolithic",
                  "Tangent_Matrix is near to singular", "rcond", rcond);
        }
        break;
      }
    }

    /*
      Update equivalent plastic strain and increment of plastic deformation
    */
    __update_internal_variables_plastic(
        Increment_E_plastic, IO_State.Stress, IO_State.Equiv_Plast_Str,
        IO_State.Kappa, IO_State.D_phi, T_k1, eigvec_b_e_tr, d_G_d_T, Lambda_k1,
        delta_lambda_k1, kappa_k1[0]);
  }

  // Plastic corrector step for the left Cauchy-Green tensor
  E_hencky_trial[0] -= Increment_E_plastic[0];
  E_hencky_trial[1] -= Increment_E_plastic[1];
  E_hencky_trial[2] -= Increment_E_plastic[2];

  // Update elastic left Cauchy-Green tensor
  STATUS = __corrector_b_e(IO_State.b_e, eigvec_b_e_tr, E_hencky_trial);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __corrector_b_e" RESET "\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __compute_trial_b_e(double *eigval_b_e_tr, double *eigvec_b_e_tr,
                               const double *b_e, const double *d_phi) {

#if NumberDimensions == 2

  eigvec_b_e_tr[0] =
      d_phi[0] * b_e[0] * d_phi[0] + d_phi[0] * b_e[1] * d_phi[1] +
      d_phi[1] * b_e[2] * d_phi[0] + d_phi[1] * b_e[3] * d_phi[1];

  eigvec_b_e_tr[1] =
      d_phi[0] * b_e[0] * d_phi[2] + d_phi[0] * b_e[1] * d_phi[3] +
      d_phi[1] * b_e[2] * d_phi[2] + d_phi[1] * b_e[3] * d_phi[3];

  eigvec_b_e_tr[3] =
      d_phi[2] * b_e[0] * d_phi[0] + d_phi[2] * b_e[1] * d_phi[1] +
      d_phi[3] * b_e[2] * d_phi[0] + d_phi[3] * b_e[3] * d_phi[1];

  eigvec_b_e_tr[4] =
      d_phi[2] * b_e[0] * d_phi[2] + d_phi[2] * b_e[1] * d_phi[3] +
      d_phi[3] * b_e[2] * d_phi[2] + d_phi[3] * b_e[3] * d_phi[3];

  eigvec_b_e_tr[8] = b_e[4];

#else
  No esta implementado
#endif

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  puts("Trial elastic left Cauchy-Green");
  printf("%e %e %e \n", eigvec_b_e_tr[0], eigvec_b_e_tr[1], eigvec_b_e_tr[2]);
  printf("%e %e %e \n", eigvec_b_e_tr[3], eigvec_b_e_tr[4], eigvec_b_e_tr[5]);
  printf("%e %e %e \n", eigvec_b_e_tr[6], eigvec_b_e_tr[7], eigvec_b_e_tr[8]);

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
  int IPIV[3] = {0, 0, 0};
  double wi[3];
  double vl[9];

  /*
    Query and allocate the optimal workspace
  */
  lwork = -1;
  dsyev_("V", "L", &n, eigvec_b_e_tr, &lda, eigval_b_e_tr, &wkopt, &lwork,
         &info);
  lwork = (int)wkopt;
  work = (double *)malloc(lwork * sizeof(double));

  /* Check for convergence */
  if (info > 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): the %i-th argument had an "
            "illegal value." RESET "\n",
            abs(info));
    return EXIT_FAILURE;
  }

  dsyev_("V", "L", &n, eigvec_b_e_tr, &lda, eigval_b_e_tr, work, &lwork, &info);
  /* Check for convergence */
  if (info > 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    free(work);
    fprintf(stderr,
            "" RED "Error in dsyev_(): the %i-th argument had an "
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

/**************************************************************/

static int __elastic_tangent(double *CC, double *AA, double E, double nu,
                             double K, double G) {
  CC[0] = 1.0 / E;
  CC[1] = -nu / E;
  CC[2] = -nu / E;

  CC[3] = -nu / E;
  CC[4] = 1.0 / E;
  CC[5] = -nu / E;

  CC[6] = -nu / E;
  CC[7] = -nu / E;
  CC[8] = 1.0 / E;

  AA[0] = K + 2 * G;
  AA[1] = K;
  AA[2] = K;

  AA[3] = K;
  AA[4] = K + 2 * G;
  AA[5] = K;

  AA[6] = K;
  AA[7] = K;
  AA[8] = K + 2 * G;

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __trial_elastic(double *T_tr, const double *E_hencky_trial,
                           const double *AA) {

#if NumberDimensions == 2
  //  T_tr[0] = AA[0] * E_hencky_trial[0] + AA[1] * E_hencky_trial[1] +
  //            AA[2] * E_hencky_trial[2];
  //  T_tr[1] = AA[3] * E_hencky_trial[0] + AA[4] * E_hencky_trial[1] +
  //            AA[5] * E_hencky_trial[2];
  //  T_tr[2] = AA[6] * E_hencky_trial[0] + AA[7] * E_hencky_trial[1] +
  //            AA[8] * E_hencky_trial[2];

  T_tr[0] = Confining_pressure;
  T_tr[1] = Confining_pressure + AA[4] * E_hencky_trial[1];
  T_tr[2] = Confining_pressure;

#else
  No esta implementado
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __E_hencky(double *E_hencky_k, const double *T_k, const double *CC) {
#if NumberDimensions == 2

  E_hencky_k[0] = CC[0] * T_k[0] + CC[1] * T_k[1] + CC[2] * T_k[2];
  E_hencky_k[1] = CC[3] * T_k[0] + CC[4] * T_k[1] + CC[5] * T_k[2];
  E_hencky_k[2] = CC[6] * T_k[0] + CC[7] * T_k[1] + CC[8] * T_k[2];

#else
  No esta implementado
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __update_internal_variables_elastic(double *Stress,
                                               const double *D_phi,
                                               const double *T_tr,
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

    T[0] += T_tr[i] * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 0];
    T[1] += T_tr[i] * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 1];
    T[3] += T_tr[i] * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 0];
    T[4] += T_tr[i] * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 1];
    T[8] += T_tr[i] * eigvec_b_e_tr[i * 3 + 2] * eigvec_b_e_tr[i * 3 + 2];
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

/**************************************************************/

static int __kappa(double *kappa, const double *a, double Lambda, double I1,
                   double alpha) {

  kappa[0] = a[0] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  kappa[1] = alpha * kappa[0];

  if (kappa[0] < 0.0) {
    fprintf(stderr, "" RED "Negative value of kappa: %f " RESET "\n", kappa[0]);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __d_kappa_phi_d_stress(double *d_kappa_phi_d_stress, const double *a,
                                  double Lambda, double I1) {

#if NumberDimensions == 2

  d_kappa_phi_d_stress[0] =
      a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  d_kappa_phi_d_stress[1] =
      a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  d_kappa_phi_d_stress[2] =
      a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);

#else
  No esta implementado
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __d_kappa_phi_d_lambda(double *d_kappa_phi_d_lambda, const double *a,
                                  double Lambda, double I1) {

  *d_kappa_phi_d_lambda =
      (1 - a[2] * Lambda) * a[0] * exp(a[1] * I1) * exp(-a[2] * Lambda);

  return EXIT_SUCCESS;
}

/**************************************************************/

static double __F(double c0, double kappa_phi, double pa, double I1, double I2,
                  double I3, double m) {

  double K1 = c0 + kappa_phi * pow(pa / I1, m);

  double PHI = cbrt(K1 * I3) - cbrt(I1 * I2);

  return PHI;
}

/**************************************************************/

static int __d_F_d_stress(double *d_F_d_stress, const double *T_k, double I1,
                          double I2, double I3, double c0, double kappa_phi,
                          double pa, double m) {
  double Grad_f[3];

  double K1 = c0 + kappa_phi * pow(pa / I1, m);
  double b1 = m * kappa_phi * (pow(pa / I1, m)) * (cbrt(I3) / I1);

  Grad_f[0] = (I1 * (I1 - T_k[0]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2));
  Grad_f[1] = (I1 * (I1 - T_k[1]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2));
  Grad_f[2] = (I1 * (I1 - T_k[2]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2));

  d_F_d_stress[0] = cbrt(K1 * I3) / (3.0 * T_k[0]) -
                    b1 / (3.0 * pow(cbrt(K1), 2.0)) - Grad_f[0];
  d_F_d_stress[1] = cbrt(K1 * I3) / (3.0 * T_k[1]) -
                    b1 / (3.0 * pow(cbrt(K1), 2.0)) - Grad_f[1];
  d_F_d_stress[2] = cbrt(K1 * I3) / (3.0 * T_k[2]) -
                    b1 / (3.0 * pow(cbrt(K1), 2.0)) - Grad_f[2];

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __d_F_d_kappa_phi(double *d_F_d_kappa_phi, double I1, double I3,
                             double c0, double m, double pa, double kappa_phi) {

  double K1 = c0 + kappa_phi * pow(pa / I1, m);

  *d_F_d_kappa_phi = cbrt(I3) / (3.0 * pow(cbrt(K1), 2.0)) * pow(pa / I1, m);

  return EXIT_SUCCESS;
}

/**************************************************************/

static double __G(double c0, double kappa2, double pa, double I1, double I2,
                  double I3, double m) {

  double K2 = c0 + kappa2 * pow(pa / I1, m);

  double G = cbrt(K2 * I3) - cbrt(I1 * I2);

  return G;
}

/**************************************************************/

static int __d_G_d_T(double *d_G_d_T, const double *T_k, double I1, double I2,
                     double I3, double c0, double kappa_psi, double pa,
                     double m) {

  double Grad_g[3] = {0.0, 0.0, 0.0};
  double K2 = c0 + kappa_psi * pow(pa / I1, m);
  double b2 = m * kappa_psi * (pow(pa / I1, m)) * (cbrt(I3) / I1);

  Grad_g[0] = (I1 * (I1 - T_k[0]) + I2) / (3 * pow(cbrt(I1 * I2), 2));
  Grad_g[1] = (I1 * (I1 - T_k[1]) + I2) / (3 * pow(cbrt(I1 * I2), 2));
  Grad_g[2] = (I1 * (I1 - T_k[2]) + I2) / (3 * pow(cbrt(I1 * I2), 2));

  d_G_d_T[0] =
      cbrt(K2 * I3) / (3 * T_k[0]) - b2 / (3 * pow(cbrt(K2), 2)) - Grad_g[0];
  d_G_d_T[1] =
      cbrt(K2 * I3) / (3 * T_k[1]) - b2 / (3 * pow(cbrt(K2), 2)) - Grad_g[1];
  d_G_d_T[2] =
      cbrt(K2 * I3) / (3 * T_k[2]) - b2 / (3 * pow(cbrt(K2), 2)) - Grad_g[2];

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __dd_G_dd_T(double *dd_G_dd_T, const double *T_k, double kappa_psi,
                       double I1, double I2, double I3, double m, double pa,
                       double c0) {

  double K2 = c0 + kappa_psi * pow(pa / I1, m);

  double b2 = m * kappa_psi * (pow(pa / I1, m)) * (cbrt(I3) / I1);

  double Grad_K2[3] = {0.0, 0.0, 0.0};
  Grad_K2[0] = -(m * kappa_psi / I1) * pow(pa / I1, m);
  Grad_K2[1] = -(m * kappa_psi / I1) * pow(pa / I1, m);
  Grad_K2[2] = -(m * kappa_psi / I1) * pow(pa / I1, m);

  double Grad_b2[3] = {0.0, 0.0, 0.0};
  Grad_b2[0] = (b2 / I1) * (I1 / (3 * T_k[0]) - m - 1.0);
  Grad_b2[1] = (b2 / I1) * (I1 / (3 * T_k[1]) - m - 1.0);
  Grad_b2[2] = (b2 / I1) * (I1 / (3 * T_k[2]) - m - 1.0);

  double Hess_g[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double d_g_d_stress[3];

  d_g_d_stress[0] = (I1 * (I1 - T_k[0]) + I2) / (3 * pow(cbrt(I1 * I2), 2));
  d_g_d_stress[1] = (I1 * (I1 - T_k[1]) + I2) / (3 * pow(cbrt(I1 * I2), 2));
  d_g_d_stress[2] = (I1 * (I1 - T_k[2]) + I2) / (3 * pow(cbrt(I1 * I2), 2));

  for (unsigned A = 0; A < 3; A++) {
    for (unsigned B = 0; B < 3; B++) {
      Hess_g[A * 3 + B] =
          (3.0 * I1 - T_k[A] - T_k[B] - I1 * (A == B)) /
              (3.0 * pow(cbrt(I1 * I2), 2)) -
          (2.0 / cbrt(I1 * I2)) * d_g_d_stress[A] * d_g_d_stress[B];
    }
  }

  for (unsigned A = 0; A < 3; A++) {
    for (unsigned B = 0; B < 3; B++) {
      dd_G_dd_T[A * 3 + B] =
          (1.0 / 3.0) * cbrt(K2 * I3) *
              (1.0 / (3 * T_k[A] * T_k[B]) - 1.0 * (A == B) / pow(T_k[A], 2)) +
          (cbrt(I3) / T_k[A] + 2.0 * b2 / K2) * Grad_K2[B] /
              (9.0 * pow(cbrt(K2), 2)) -
          Grad_b2[B] / (3.0 * pow(cbrt(K2), 2)) - Hess_g[A * 3 + B];
    }
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __dd_G_d_T_d_kappa_psi(double *dd_G_d_T_d_kappa_psi,
                                  const double *T_k, double I1, double I3,
                                  double m, double pa, double c0,
                                  double kappa_psi) {

  double K2 = c0 + kappa_psi * pow(pa / I1, m);
  double b2 = m * kappa_psi * (pow(pa / I1, m)) * (cbrt(I3) / I1);

  dd_G_d_T_d_kappa_psi[0] =
      pow((pa / I1), m) *
      (cbrt(I3) / (3.0 * T_k[0]) + 2.0 * b2 / (3.0 * K2) - m * cbrt(I3) / I1) /
      (3.0 * pow(cbrt(K2), 2));
  dd_G_d_T_d_kappa_psi[1] =
      pow((pa / I1), m) *
      (cbrt(I3) / (3.0 * T_k[1]) + 2.0 * b2 / (3.0 * K2) - m * cbrt(I3) / I1) /
      (3.0 * pow(cbrt(K2), 2));
  dd_G_d_T_d_kappa_psi[2] =
      pow((pa / I1), m) *
      (cbrt(I3) / (3.0 * T_k[2]) + 2.0 * b2 / (3.0 * K2) - m * cbrt(I3) / I1) /
      (3.0 * pow(cbrt(K2), 2));
  return EXIT_SUCCESS;
}

/**************************************************************/

static int __residual(double *Residual, double *Error_k,
                      const double *E_hencky_trial, const double *E_hencky_k,
                      const double *d_G_d_T, const double *kappa_k,
                      const double *kappa_hat, double F_k,
                      double delta_lambda_k) {

  Residual[0] = E_hencky_k[0] - E_hencky_trial[0] + delta_lambda_k * d_G_d_T[0];
  Residual[1] = E_hencky_k[1] - E_hencky_trial[1] + delta_lambda_k * d_G_d_T[1];
  Residual[2] = E_hencky_k[2] - E_hencky_trial[2] + delta_lambda_k * d_G_d_T[2];
  Residual[3] = kappa_k[0] - kappa_hat[0];
  Residual[4] = F_k;

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

static int __tangent_matrix(double *Tangent_Matrix, const double *CC,
                            const double *d_F_d_stress, double d_F_d_kappa_phi,
                            const double *d_G_d_T, const double *dd_G_dd_T,
                            const double *dd_G_d_T_d_kappa2,
                            const double *d_kappa_phi_d_stress,
                            double d_kappa_phi_d_lambda, double alpha,
                            double delta_lambda_k) {

  /* First row */
  Tangent_Matrix[0] = CC[0] + delta_lambda_k * dd_G_dd_T[0];
  Tangent_Matrix[1] = CC[1] + delta_lambda_k * dd_G_dd_T[1];
  Tangent_Matrix[2] = CC[2] + delta_lambda_k * dd_G_dd_T[2];
  Tangent_Matrix[3] = alpha * delta_lambda_k * dd_G_d_T_d_kappa2[0];
  Tangent_Matrix[4] = d_G_d_T[0];

  /* Second row */
  Tangent_Matrix[5] = CC[3] + delta_lambda_k * dd_G_dd_T[3];
  Tangent_Matrix[6] = CC[4] + delta_lambda_k * dd_G_dd_T[4];
  Tangent_Matrix[7] = CC[5] + delta_lambda_k * dd_G_dd_T[5];
  Tangent_Matrix[8] = alpha * delta_lambda_k * dd_G_d_T_d_kappa2[1];
  Tangent_Matrix[9] = d_G_d_T[1];

  /* Third row */
  Tangent_Matrix[10] = CC[6] + delta_lambda_k * dd_G_dd_T[6];
  Tangent_Matrix[11] = CC[7] + delta_lambda_k * dd_G_dd_T[7];
  Tangent_Matrix[12] = CC[8] + delta_lambda_k * dd_G_dd_T[8];
  Tangent_Matrix[13] = alpha * delta_lambda_k * dd_G_d_T_d_kappa2[2];
  Tangent_Matrix[14] = d_G_d_T[2];

  /* Four row */
  Tangent_Matrix[15] = -d_kappa_phi_d_stress[0];
  Tangent_Matrix[16] = -d_kappa_phi_d_stress[1];
  Tangent_Matrix[17] = -d_kappa_phi_d_stress[2];
  Tangent_Matrix[18] = 1.0;
  Tangent_Matrix[19] = -d_kappa_phi_d_lambda;

  /* Five row */
  Tangent_Matrix[20] = d_F_d_stress[0];
  Tangent_Matrix[21] = d_F_d_stress[1];
  Tangent_Matrix[22] = d_F_d_stress[2];
  Tangent_Matrix[23] = d_F_d_kappa_phi;
  Tangent_Matrix[24] = 0.0;

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __condition_number(double *RCOND, double *Tangent_Matrix)
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
  double WORK_RCOND[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  int IPIV[5] = {0, 0, 0, 0, 0};
  int IWORK_RCOND[5] = {0, 0, 0, 0, 0};

  // Compute 1-norm
  ANORM = dlange_("1", &N_rows, &N_cols, Tangent_Matrix, &LDA, WORK_ANORM);

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N_rows, &N_cols, Tangent_Matrix, &LDA, IPIV, &INFO);

  // Check output of dgetrf
  if (INFO != 0) {
    if (INFO < 0) {
      printf("%s : \n", "Error in Frictional_Monolithic()");
      printf("the %i-th argument had an illegal value", abs(INFO));
    } else if (INFO > 0) {
      printf("%s :\n", "Error in Frictional_Monolithic()");
      printf(" M(%i,%i) %s \n %s \n %s \n %s \n", INFO, INFO,
             "is exactly zero. The factorization",
             "has been completed, but the factor M is exactly",
             "singular, and division by zero will occur if it is used",
             "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  // Compute the Reciprocal condition number
  dgecon_("1", &N_rows, Tangent_Matrix, &LDA, &ANORM, RCOND, WORK_RCOND,
          IWORK_RCOND, &INFO);

  if (INFO < 0) {
    printf("Error in Frictional_Monolithic() : the %i-th argument of dgecon_ "
           "had an illegal value\n",
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

  /*
    Check error messages in the LAPACK LU descomposition
  */
  if (INFO) {
    fprintf(stderr, "%s : %s %s %s \n", "Error in Frictional_Monolithic",
            "The function", "dgetrf_", "returned an error message !!!");
    return EXIT_FAILURE;
  }

  /*
    Solve
  */
  dgetrs_(&TRANS, &Order, &NRHS, Tangent_Matrix, &LDA, IPIV, Residual, &LDB,
          &INFO);
  if (INFO) {
    fprintf(stderr, "%s : %s %s %s \n", "Error in Frictional_Monolithic",
            "The function", "dgetrs_", "returned an error message !!!");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __update_internal_variables_plastic(
    double *Increment_E_plastic, double *Stress, double *eps_n1,
    double *kappa_n1, const double *D_phi, const double *T_tr_k,
    const double *eigvec_b_e_tr, const double *d_G_d_T, double Lambda_k,
    double delta_lambda_k, double kappa_phi_k) {

  // Update hardening parameters
  *eps_n1 = Lambda_k;
  *kappa_n1 = kappa_phi_k;

  // Update the increment of the plastic strain tensor
  Increment_E_plastic[0] = delta_lambda_k * d_G_d_T[0];
  Increment_E_plastic[1] = delta_lambda_k * d_G_d_T[1];
  Increment_E_plastic[2] = delta_lambda_k * d_G_d_T[2];

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

    T[0] += T_tr_k[i] * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 0];
    T[1] += T_tr_k[i] * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 1];
    T[3] += T_tr_k[i] * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 0];
    T[4] += T_tr_k[i] * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 1];
    T[8] += T_tr_k[i] * eigvec_b_e_tr[i * 3 + 2] * eigvec_b_e_tr[i * 3 + 2];
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