#include <math.h>
#include "nl-partsol.h"

static int __compute_trial_b_e(
    double *eigval_b_e_tr /**< [out] Eigenvalues of b elastic trial. */,
    double *eigvec_b_e_tr /**< [out] Eigenvector of b elastic trial. */,
    const double *b_e /**< [in] (n) Elastic left Cauchy-Green.*/,
    const double *d_phi /**< [in] Incremental deformation gradient. */);

static int __corrector_b_e(
    double *b_e /**< [out] (n+1) Elastic deformation gradient. */,
    const double *eigval_b_e_tr /**< [in] Eigenvalues of b elastic trial. */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */,
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
  double kirchhoff_tr_vol[3] = {0.0, 0.0, 0.0};
  double kirchhoff_tr_dev[3] = {0.0, 0.0, 0.0};

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
  if (STATUS) {
    ERROR();
    fprintf(stderr, "__compute_trial_b_e\n");
    RESET_ERROR();
    return STATUS;
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

  STATUS = __corrector_b_e(IO_State.b_e, eigval_b_e_tr, eigvec_b_e_tr,
                           IO_State.Increment_E_plastic);
  if (STATUS) {
    ERROR();
    fprintf(stderr, "__corrector_b_e\n");
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
  puts("t = n + 1 elastic left Cauchy-Green tensor");
  printf("%e, %e, %e \n", IO_State.b_e[0], IO_State.b_e[1], 0.0);
  printf("%e, %e, %e \n", IO_State.b_e[2], IO_State.b_e[3], 0.0);
  printf("%e, %e, %e \n", 0.0, 0.0, IO_State.b_e[4]);
#endif
#endif

  return STATUS;
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
    printf("Error in Eigen_analysis__TensorLib__() : The algorithm failed to "
           "compute eigenvalues.\n");
    return EXIT_FAILURE;
  }

  /* Solve eigenproblem */
  dgeev_("N", "V", &n, b_e_tr, &lda, eigval_b_e_tr, wi, vl, &ldvl,
         eigvec_b_e_tr, &ldvr, work, &lwork, &info);

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

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __corrector_b_e(double *b_e, const double *eigval_b_e_tr,
                           const double *eigvec_b_e_tr,
                           const double *Increment_E_plastic) {

  double eigval_b_e[3] = {0.0, 0.0, 0.0};

  eigval_b_e[0] = eigval_b_e_tr[0] * exp(-2 * Increment_E_plastic[0]);
  eigval_b_e[1] = eigval_b_e_tr[1] * exp(-2 * Increment_E_plastic[1]);
  eigval_b_e[2] = eigval_b_e_tr[2] * exp(-2 * Increment_E_plastic[2]);

#if NumberDimensions == 2

  b_e[0] = 0.0;
  b_e[1] = 0.0;
  b_e[2] = 0.0;
  b_e[3] = 0.0;
  b_e[4] = 0.0;

  for (unsigned i = 0; i < 3; i++) {
    b_e[0] +=
        eigval_b_e_tr[i] * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 0];
    b_e[1] +=
        eigval_b_e_tr[i] * eigvec_b_e_tr[i * 3 + 0] * eigvec_b_e_tr[i * 3 + 1];
    b_e[2] +=
        eigval_b_e_tr[i] * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 0];
    b_e[3] +=
        eigval_b_e_tr[i] * eigvec_b_e_tr[i * 3 + 1] * eigvec_b_e_tr[i * 3 + 1];
    b_e[4] +=
        eigval_b_e_tr[i] * eigvec_b_e_tr[i * 3 + 2] * eigvec_b_e_tr[i * 3 + 2];
  }

#else
  No esta implementado
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