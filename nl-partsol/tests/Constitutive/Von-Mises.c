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
  double kappa_0;
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
  double *EPS;   // Equivalent plastic strain
  double *Kappa; // Hardening Parameter

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
    const double *T_back /**< [in] Back-stress (kinematic hardening) */,
    double *J2 /**< [out] Second invariant of the deviatoric stress tensor */,
    const double *E_hencky_trial, /**< [in] Trial elastic strain tensor. */
    double K /**< [in] First Lamé invariant. */,
    double G /**< [in] Second Lamé invariant. */);

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

static int __kappa(
    double *kappa_k /**< [out] Hardening function (isotropic,kinematic). */,
    double sigma_y /**< [in] Reference hardening */,
    double eps_k /**< [in] Equivalent plastic strain */,
    double H /**< [in] Hardening modulus */,
    double theta /**< [in] Ratio isotropic/kinematic hardening */,
    double
        K_0 /**< [in] Reference non-linear hardening (saturation parameter) */,
    double K_inf /**< [in] Saturation non-linear hardening (saturation
                    parameter) */
    ,
    double delta /**< [in] Saturation eps (saturation parameter) */);

static int __d_kappa(
    double *d_kappa /**< [out] Derivative of the hardening function
                       (isotropic,kinematic) */
    ,
    double eps_k /**< [in] Equivalent plastic strain*/,
    double H /**< [in] Hardening modulus */,
    double theta /**< [in] Ratio isotropic/kinematic hardening */,
    double
        K_0 /**< [in] Reference non-linear hardening (saturation parameter) */,
    double K_inf /**< [in] Saturation non-linear hardening (saturation
                    parameter) */
    ,
    double delta /**< [in] Saturation eps (saturation parameter) */);

static double __yield_function(
    const double
        *kappa_k /**< [in] Hardening function (isotropic,kinematic), t = k */,
    const double
        *kappa_n /**< [in] Hardening function (isotropic,kinematic), t = n */,
    double J2 /**< [in] Second invariant of the deviatoric stress tensor */,
    double d_gamma_k /**< [in] Increment of the discrete plastic multiplier */,
    double eps_k /**< [in] Equivalent plastic strain*/,
    double G /**< [in] Second Lamé invariant. */);

static double __d_yield_function(
    const double *d_kappa_k /**< [in] Hardening function derivative
                               (isotropic,kinematic), t = k */
    ,
    double G /**< [in] Second Lamé invariant. */);

static int __update_internal_variables_plastic(
    double *Increment_E_plastic /**< [in/out] Increment plastic strain */,
    double *Stress /**< [in/out] Nominal stress tensor */,
    double *T_back /**< [in] Back-stress (kinematic hardening) */,
    double *eps_n1 /**< [in/out] Equivalent plastic strain*/,
    const double *D_phi /**< [in] Total deformation gradient. */,
    const double *T_tr_vol /**< [in] Volumetric elastic stress tensor. */,
    const double *T_tr_dev /**< [in] Deviatoric elastic stress tensor. */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */,
    const double *n /**< [out] Plastic flow direction. */,
    double d_gamma_k /**< [in] Discrete plastic multiplier */,
    double G /**< [in] Second Lamé invariant. */,
    double eps_k /**< [in] Equivalent plastic strain */,
    double d_K_kin /**< [in] Increment of the kinematic hardening */);

int compute_Kirchhoff_Stress_Von_Mises__Constitutive__(State_Parameters IO_State, Material MatProp);

/**************************************************************/
/******************** Solver Parameters ***********************/
/**************************************************************/
#define TOL_Radial_Returning 10E-14
#define Max_Iterations_Radial_Returning 10
#define PI__MatrixLib__ 3.14159265358979323846

/**************************************************************/
/******************* Material Parameters **********************/
/**************************************************************/
#define NumberDimensions 2
#define YoungMouduls 10.0E3
#define PoissonRatio 0.2
#define yield_value 10
#define H_value 0.0
#define theta_value 0.0
#define K_0_value 0.0
#define K_inf_value 0.0
#define delta_value 0.0
#define NumberSteps 100 // 50 // 3500 // 2500 //
#define dF_yy 0.9999
#define Confining_pressure 0.0

int main() {
  int STATUS = EXIT_SUCCESS;
  State_Parameters IO_State;
  Material MatProp;

  // Define material variables
  MatProp.E = YoungMouduls;
  MatProp.nu = PoissonRatio;
  MatProp.kappa_0 = yield_value;
  MatProp.Hardening_modulus = H_value;
  MatProp.theta_Hardening_Voce = theta_value;
  MatProp.K_0_Hardening_Voce = K_0_value;
  MatProp.K_inf_Hardening_Voce = K_inf_value;
  MatProp.delta_Hardening_Voce = delta_value;
  MatProp.ReferencePressure = Confining_pressure;

  // Initialize state variables
  bool Status_particle = false;
  double *stress = (double *)calloc(5 * NumberSteps, sizeof(double));
  double d_phi[5] = {1.0, 0.0, 0.0, dF_yy, 1.0};
  double *D_phi = (double *)calloc(5 * NumberSteps, sizeof(double));
  double *b_e = (double *)calloc(5 * NumberSteps, sizeof(double));
  double *beta = (double *)calloc(3 * NumberSteps, sizeof(double));
  double *EPS = (double *)calloc(NumberSteps, sizeof(double));

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

    printf("Step: %i \n", i);

    // Asign variables to the solver
    IO_State.Particle_Idx = 0;
    IO_State.Stress = &stress[i * 5];
    IO_State.b_e = &b_e[i * 5];
    IO_State.EPS = &EPS[i];
    IO_State.Back_stress = &beta[i];
    IO_State.d_phi = d_phi;
    IO_State.D_phi = &D_phi[i * 5];

    IO_State.b_e[0] = b_e[(i - 1) * 5 + 0];
    IO_State.b_e[1] = b_e[(i - 1) * 5 + 1];
    IO_State.b_e[2] = b_e[(i - 1) * 5 + 2];
    IO_State.b_e[3] = b_e[(i - 1) * 5 + 3];
    IO_State.b_e[4] = b_e[(i - 1) * 5 + 4];

    IO_State.D_phi[0] = D_phi[(i - 1) * 5 + 0] * d_phi[0];
    IO_State.D_phi[3] = D_phi[(i - 1) * 5 + 3] * d_phi[3];
    IO_State.D_phi[4] = D_phi[(i - 1) * 5 + 4] * d_phi[4];

    *IO_State.EPS = EPS[i - 1];
    *IO_State.Back_stress = beta[i - 1];

    // Run solver
    STATUS = compute_Kirchhoff_Stress_Von_Mises__Constitutive__(IO_State, MatProp);

    if (STATUS) {
      fprintf(stderr, "Failure\n");
      free(stress);
      free(D_phi);
      free(b_e);
      free(beta);
      free(EPS);
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

  // Free memory
  free(stress);
  free(D_phi);
  free(b_e);
  free(beta);
  free(EPS);

  return STATUS;
}

/**************************************************************/
int compute_Kirchhoff_Stress_Von_Mises__Constitutive__(State_Parameters IO_State, Material MatProp)
/*
  Radial returning algorithm for the Von-Mises (Simo and Hughes)
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
#if NumberDimensions == 2
  printf("%f %f %f \n", eigvec_b_e_tr[0], eigvec_b_e_tr[1], 0.0);
  printf("%f %f %f \n", eigvec_b_e_tr[2], eigvec_b_e_tr[3], 0.0);
  printf("%f %f %f \n", 0.0, 0.0, 1.0);
#else
  printf("%f %f %f \n", eigvec_b_e_tr[0], eigvec_b_e_tr[1], eigvec_b_e_tr[2]);
  printf("%f %f %f \n", eigvec_b_e_tr[3], eigvec_b_e_tr[4], eigvec_b_e_tr[5]);
  printf("%f %f %f \n", eigvec_b_e_tr[6], eigvec_b_e_tr[7], eigvec_b_e_tr[8]);
#endif

  printf("E_hencky_trial: [%f, %f, %f] \n", E_hencky_trial[0],
         E_hencky_trial[1], E_hencky_trial[2]);

#endif
#endif

  // Material parameters
  double K = MatProp.E / (3.0 * (1.0 - 2.0 * MatProp.nu));
  double G = MatProp.E / (2.0 * (1.0 + MatProp.nu));
  double sigma_y = MatProp.kappa_0;
  double H = MatProp.Hardening_modulus;
  double theta = MatProp.theta_Hardening_Voce;
  double K_0 = MatProp.K_0_Hardening_Voce;
  double K_inf = MatProp.K_inf_Hardening_Voce;
  double delta = MatProp.delta_Hardening_Voce;

  // Define tensorial internal variables
  double n[3] = {0.0, 0.0, 0.0};
  double Increment_E_plastic[3] = {0.0, 0.0, 0.0};
  double T_back[3];
  T_back[0] = IO_State.Back_stress[0];
  T_back[1] = IO_State.Back_stress[1];
  T_back[2] = IO_State.Back_stress[2];
  double kappa_n[2];
  double kappa_k[2];
  double d_kappa_k[2];

  // Define scalar internal variables
  double PHI, PHI_0 = 0.0;
  double d_PHI = 0.0;
  double J2 = 0.0;
  double eps_n = *IO_State.EPS;
  double eps_k = eps_n;
  double d_gamma_k = 0.0;

  // Initialise solver parameters
  double TOL = TOL_Radial_Returning;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;

  STATUS =
      __trial_elastic(T_tr_vol, T_tr_dev, T_back, &J2, E_hencky_trial, K, G);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __trial_elastic" RESET "\n");
    return EXIT_FAILURE;
  }

  STATUS = __kappa(kappa_n, sigma_y, eps_n, H, theta, K_0, K_inf, delta);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __kappa" RESET "\n");
    return EXIT_FAILURE;
  }

  //  Check yield condition
  PHI_0 = __yield_function(kappa_n, kappa_n, J2, d_gamma_k, eps_n, G);

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

    kappa_k[0] = kappa_n[0];
    kappa_k[1] = kappa_n[1];

    PHI = PHI_0;

    while (fabs(PHI / PHI_0) >= TOL) {

      Iter++;

      if (Iter == MaxIter)
        break;

      STATUS = __d_kappa(d_kappa_k, eps_k, H, theta, K_0, K_inf, delta);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __d_kappa" RESET "\n");
        return EXIT_FAILURE;
      }

      d_PHI = __d_yield_function(d_kappa_k, G);

      d_gamma_k += -PHI / d_PHI;

      eps_k = eps_n + sqrt(2. / 3.) * d_gamma_k;

      STATUS = __kappa(kappa_k, sigma_y, eps_k, H, theta, K_0, K_inf, delta);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __kappa" RESET "\n");
        return EXIT_FAILURE;
      }

      PHI = __yield_function(kappa_k, kappa_n, J2, d_gamma_k, eps_k, G);
    }

    double d_K_kin = kappa_k[1] - kappa_n[1];

    STATUS = __update_internal_variables_plastic(
        Increment_E_plastic, IO_State.Stress, IO_State.Back_stress,
        IO_State.EPS, IO_State.D_phi, T_tr_vol, T_tr_dev, eigvec_b_e_tr, n,
        d_gamma_k, G, eps_k, d_K_kin);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in __update_internal_variables_plastic" RESET "\n");
      return EXIT_FAILURE;
    }
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

  printf("Number of iterations: %i \n", Iter);

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __compute_trial_b_e(double *eigval_b_e_tr, double *eigvec_b_e_tr,
                               const double *b_e, const double *d_phi) {

  unsigned Ndim = NumberDimensions;

  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      for (unsigned k = 0; k < Ndim; k++) {
        for (unsigned l = 0; l < Ndim; l++) {
          eigvec_b_e_tr[i * Ndim + j] +=
              d_phi[i * Ndim + k] * b_e[k * Ndim + l] * d_phi[j * Ndim + l];
        }
      }
    }
  }

#if NumberDimensions == 2
  /* Locals */
  int n = 2;
  int lda = 2;
  int ldvl = 2;
  int ldvr = 2;
  int info;
  int lwork;
  double wkopt;
  double *work;

  /* Local arrays */
  int IPIV[2] = {0, 0};
  double wi[2];
  double vl[4];

#else

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

#endif

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

#if NumberDimensions == 2
  eigval_b_e_tr[2] = b_e[4];
#endif

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __corrector_b_e(double *b_e, const double *eigvec_b_e_tr,
                           const double *E_hencky_trial) {

  int Ndim = NumberDimensions;

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

#else

  b_e[0] = 0.0;
  b_e[1] = 0.0;
  b_e[2] = 0.0;
  b_e[3] = 0.0;
  b_e[4] = 0.0;
  b_e[5] = 0.0;
  b_e[6] = 0.0;
  b_e[7] = 0.0;
  b_e[8] = 0.0;

#endif

  for (unsigned A = 0; A < Ndim; A++) {
    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        b_e[i * Ndim + j] += eigval_b_e[A] * eigvec_b_e_tr[A + i * Ndim] *
                             eigvec_b_e_tr[A + j * Ndim];
      }
    }
  }

#if NumberDimensions == 2
  b_e[4] = eigval_b_e[2];
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __trial_elastic(double *T_tr_vol, double *T_tr_dev,
                           const double *T_back, double *J2,
                           const double *E_hencky_trial, double K, double G) {

  double E_hencky_trial_vol[3];
  double tr_E_hencky_trial =
      E_hencky_trial[0] + E_hencky_trial[1] + E_hencky_trial[2];

  E_hencky_trial_vol[0] = (1.0 / 3.0) * tr_E_hencky_trial;
  E_hencky_trial_vol[1] = (1.0 / 3.0) * tr_E_hencky_trial;
  E_hencky_trial_vol[2] = (1.0 / 3.0) * tr_E_hencky_trial;

  T_tr_vol[0] = K * E_hencky_trial_vol[0];
  T_tr_vol[1] = K * E_hencky_trial_vol[1];
  T_tr_vol[2] = K * E_hencky_trial_vol[2];

  T_tr_dev[0] = 2 * G * (E_hencky_trial[0] - E_hencky_trial_vol[0]) - T_back[0];
  T_tr_dev[1] = 2 * G * (E_hencky_trial[1] - E_hencky_trial_vol[1]) - T_back[1];
  T_tr_dev[2] = 2 * G * (E_hencky_trial[2] - E_hencky_trial_vol[2]) - T_back[2];

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

  int Ndim = NumberDimensions;

  // Compute the transpose of D_phi

#if NumberDimensions == 2

  double D_phi_mT[4] = {0.0, 0.0, 0.0, 0.0};
  D_phi_mT[0] = D_phi[0];
  D_phi_mT[1] = D_phi[2];
  D_phi_mT[2] = D_phi[1];
  D_phi_mT[3] = D_phi[3];

  // Parameters for dgetrf_ and dgetri_
  int INFO;
  int N = 2;
  int LDA = 2;
  int LWORK = 2;
  int IPIV[2] = {0, 0};
  double WORK[2] = {0, 0};

#else
  double D_phi_mT[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  D_phi_mT[0] = D_phi[0];
  D_phi_mT[1] = D_phi[3];
  D_phi_mT[2] = D_phi[6];
  D_phi_mT[3] = D_phi[1];
  D_phi_mT[4] = D_phi[4];
  D_phi_mT[5] = D_phi[7];
  D_phi_mT[6] = D_phi[2];
  D_phi_mT[7] = D_phi[5];
  D_phi_mT[8] = D_phi[8];

  // Parameters for dgetrf_ and dgetri_
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};

#endif

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
#if NumberDimensions == 2
  printf("%f %f %f \n", D_phi_mT[0], D_phi_mT[1], 0.0);
  printf("%f %f %f \n", D_phi_mT[2], D_phi_mT[3], 0.0);
  printf("%f %f %f \n", 0.0, 0.0, 1.0);
#else
  printf("%f %f %f \n", D_phi_mT[0], D_phi_mT[1], D_phi_mT[2]);
  printf("%f %f %f \n", D_phi_mT[3], D_phi_mT[4], D_phi_mT[5]);
  printf("%f %f %f \n", D_phi_mT[6], D_phi_mT[7], D_phi_mT[8]);
#endif

#endif
#endif

#if NumberDimensions == 2
  double T_aux[4] = {0.0, 0.0, 0.0, 0.0};
#else
  double T_aux[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif
  double T_A = 0.0;

  for (unsigned A = 0; A < Ndim; A++) {

    T_A = T_tr_vol[A] + T_tr_dev[A];

    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        T_aux[i * Ndim + j] +=
            T_A * eigvec_b_e_tr[A + i * Ndim] * eigvec_b_e_tr[A + j * Ndim];
      }
    }
  }

  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      Stress[i * Ndim + j] = 0.0;

      for (unsigned k = 0; k < Ndim; k++) {
        Stress[i * Ndim + j] += T_aux[i * Ndim + k] * D_phi_mT[k * Ndim + j];
      }
    }
  }

#if NumberDimensions == 2
  Stress[4] = T_tr_vol[2] + T_tr_dev[2];
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

static int __compute_plastic_flow_direction(double *n, const double *T_tr_dev,
                                            double J2) {
  n[0] = T_tr_dev[0] / J2;
  n[1] = T_tr_dev[1] / J2;
  n[2] = T_tr_dev[2] / J2;

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __kappa(double *kappa_k, double sigma_y, double eps_k, double H,
                   double theta, double K_0, double K_inf, double delta) {
  if (eps_k < 0.0)
    return EXIT_FAILURE;

  // Isotropic hardening
  kappa_k[0] =
      sigma_y + theta * H * eps_k + (K_inf - K_0) * (1 - exp(-delta * eps_k));

  // Kinematic hardening
  kappa_k[1] = (1 - theta) * H * eps_k;

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __d_kappa(double *d_kappa, double eps_k, double H, double theta,
                     double K_0, double K_inf, double delta) {
  if (eps_k < 0.0)
    return EXIT_FAILURE;

  // Isotropic hardening
  d_kappa[0] = (theta * H + delta * (K_inf - K_0) * exp(-delta * eps_k));

  // Kinematic hardening
  d_kappa[1] = (1 - theta) * H;

  return EXIT_SUCCESS;
}

/**************************************************************/

static double __yield_function(const double *kappa_k, const double *kappa_n,
                               double J2, double d_gamma_k, double eps_k,
                               double G) {
  double K_iso_k = kappa_k[0];
  double K_kin_k = kappa_k[1];
  double K_kin_n = kappa_n[1];

  return J2 - sqrt(2. / 3.) * (K_iso_k + K_kin_k - K_kin_n) -
         2.0 * G * d_gamma_k;
}

/**************************************************************/

static double __d_yield_function(const double *d_kappa_k, double G) {
  double d_K_iso_k = d_kappa_k[0];
  double d_K_kin_k = d_kappa_k[1];

  return -2.0 * G * (1.0 + (d_K_iso_k + d_K_kin_k) / (3 * G));
}

/**************************************************************/

static int __update_internal_variables_plastic(
    double *Increment_E_plastic, double *Stress, double *T_back, double *eps_n1,
    const double *D_phi, const double *T_tr_vol, const double *T_tr_dev,
    const double *eigvec_b_e_tr, const double *n, double d_gamma_k, double G,
    double eps_k, double d_K_kin) {

  int Ndim = NumberDimensions;

  // Update equivalent plastic strain
  *eps_n1 = eps_k;

  // Update the increment of the plastic strain tensor
  Increment_E_plastic[0] = d_gamma_k * n[0];
  Increment_E_plastic[1] = d_gamma_k * n[1];
  Increment_E_plastic[2] = d_gamma_k * n[2];

  // Compute the transpose of D_phi
#if NumberDimensions == 2

  double D_phi_mT[4] = {0.0, 0.0, 0.0, 0.0};
  D_phi_mT[0] = D_phi[0];
  D_phi_mT[1] = D_phi[2];
  D_phi_mT[2] = D_phi[1];
  D_phi_mT[3] = D_phi[3];

  // Parameters for dgetrf_ and dgetri_
  int INFO;
  int N = 2;
  int LDA = 2;
  int LWORK = 2;
  int IPIV[2] = {0, 0};
  double WORK[2] = {0, 0};

#else
  double D_phi_mT[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  D_phi_mT[0] = D_phi[0];
  D_phi_mT[1] = D_phi[3];
  D_phi_mT[2] = D_phi[6];
  D_phi_mT[3] = D_phi[1];
  D_phi_mT[4] = D_phi[4];
  D_phi_mT[5] = D_phi[7];
  D_phi_mT[6] = D_phi[2];
  D_phi_mT[7] = D_phi[5];
  D_phi_mT[8] = D_phi[8];

  // Parameters for dgetrf_ and dgetri_
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};

#endif

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
#if NumberDimensions == 2
  printf("%f %f %f \n", D_phi_mT[0], D_phi_mT[1], 0.0);
  printf("%f %f %f \n", D_phi_mT[2], D_phi_mT[3], 0.0);
  printf("%f %f %f \n", 0.0, 0.0, 1.0);
#else
  printf("%f %f %f \n", D_phi_mT[0], D_phi_mT[1], D_phi_mT[2]);
  printf("%f %f %f \n", D_phi_mT[3], D_phi_mT[4], D_phi_mT[5]);
  printf("%f %f %f \n", D_phi_mT[6], D_phi_mT[7], D_phi_mT[8]);
#endif

#endif
#endif

#if NumberDimensions == 2
  double T_aux[4] = {0.0, 0.0, 0.0, 0.0};
#else
  double T_aux[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif
  double T_A = 0.0;

  for (unsigned A = 0; A < Ndim; A++) {

    T_A = -T_tr_vol[A] + T_tr_dev[A] + T_back[A] - d_gamma_k * 2 * G * n[A];

    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        T_aux[i * Ndim + j] +=
            T_A * eigvec_b_e_tr[A * Ndim + i] * eigvec_b_e_tr[A * Ndim + j];
      }
    }
  }

  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      Stress[i * Ndim + j] = 0.0;

      for (unsigned k = 0; k < Ndim; k++) {
        Stress[i * Ndim + j] += T_aux[i * Ndim + k] * D_phi_mT[k * Ndim + j];
      }
    }
  }

#if NumberDimensions == 2
  Stress[4] = -T_tr_vol[2] + T_tr_dev[2] + T_back[2] - d_gamma_k * 2 * G * n[2];
#endif

  // Update back stress
  T_back[0] += sqrt(2. / 3.) * d_K_kin * n[0];
  T_back[1] += sqrt(2. / 3.) * d_K_kin * n[1];
  T_back[2] += sqrt(2. / 3.) * d_K_kin * n[2];

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
