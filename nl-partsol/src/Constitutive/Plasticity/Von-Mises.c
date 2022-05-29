
/**
 * @file Von-Mises.c
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2022-05-28
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Constitutive/Plasticity/Von-Mises.h"
#include "Globals.h"

/**************************************************************/

/**
 * @brief
 *
 * @param eigval_b_e_tr Eigenvalues of b elastic trial
 * @param eigvec_b_e_tr igenvector of b elastic trial
 * @param b_e Elastic left Cauchy-Green (n)
 * @param d_phi Incremental deformation gradient
 * @return int
 */
static int __compute_trial_b_e(double *eigval_b_e_tr, double *eigvec_b_e_tr,
                               const double *b_e, const double *d_phi);
/*******************************************************/

/**
 * @brief
 *
 * @param b_e (n+1) Elastic deformation gradient
 * @param eigvec_b_e_tr Eigenvector of b elastic trial
 * @param E_hencky_trial Corrected Henky strain
 * @return int
 */
static int __corrector_b_e(double *b_e, const double *eigvec_b_e_tr,
                           const double *E_hencky_trial);
/*******************************************************/

/**
 * @brief
 *
 * @param T_tr_vol Volumetric elastic stress tensor
 * @param T_tr_dev Deviatoric elastic stress tensor
 * @param T_back Back-stress (kinematic hardening)
 * @param J2 Second invariant of the deviatoric stress tensor
 * @param E_hencky_trial Trial elastic strain tensor
 * @param K First Lamé invariant
 * @param G Second Lamé invariant
 * @return int
 */
static int __trial_elastic(double *T_tr_vol, double *T_tr_dev,
                           const double *T_back, double *J2,
                           const double *E_hencky_trial, double K, double G);
/*******************************************************/

/**
 * @brief
 *
 * @param T_ppal Stress tensor in the ppal coordinates
 * @param T_tr_vol Volumetric elastic stress tensor
 * @param T_tr_dev Deviatoric elastic stress tensor
 * @return int
 */
static int __elastic_stress_ppal(double *T_ppal, const double *T_tr_vol,
                                 const double *T_tr_dev);
/**************************************************************/

/**
 * @brief
 *
 * @param T_xyz Nominal stress tensor
 * @param T_ppal Stress tensor in the ppal coordinates
 * @param eigvec_b_e_tr Eigenvector of b elastic trial
 * @return int
 */
static int __update_internal_variables_elastic(double *T_xyz,
                                               const double *T_ppal,
                                               const double *eigvec_b_e_tr);
/*******************************************************/

/**
 * @brief
 *
 * @param n Plastic flow direction
 * @param T_tr_dev Deviatoric elastic stress tensor
 * @param J2 Second invariant of the deviatoric stress tensor
 * @return int
 */
static int __compute_plastic_flow_direction(double *n, const double *T_tr_dev,
                                            double J2);
/*******************************************************/

/**
 * @brief
 *
 * @param kappa_k Hardening function (isotropic,kinematic)
 * @param sigma_y Reference hardening
 * @param eps_k Equivalent plastic strain
 * @param H Hardening modulus
 * @param theta Ratio isotropic/kinematic hardening
 * @param K_0 Reference non-linear hardening (saturation parameter)
 * @param K_inf Saturation non-linear hardening (saturation parameter)
 * @param delta Saturation eps (saturation parameter)
 * @return int
 */
static int __kappa(double *kappa_k, double sigma_y, double eps_k, double H,
                   double theta, double K_0, double K_inf, double delta);
/*******************************************************/

/**
 * @brief
 *
 * @param d_kappa Derivative of the hardening function (isotropic,kinematic)
 * @param eps_k Equivalent plastic strain
 * @param H Hardening modulus
 * @param theta Ratio isotropic/kinematic hardening
 * @param K_0 Reference non-linear hardening (saturation parameter)
 * @param K_inf Saturation non-linear hardening (saturation parameter)
 * @param delta Saturation eps (saturation parameter)
 * @return int
 */
static int __d_kappa(double *d_kappa, double eps_k, double H, double theta,
                     double K_0, double K_inf, double delta);
/*******************************************************/

/**
 * @brief
 *
 * @param kappa_k Hardening function (isotropic,kinematic), t = k
 * @param kappa_n Hardening function (isotropic,kinematic), t = n
 * @param J2 Second invariant of the deviatoric stress tensor
 * @param d_gamma_k Increment of the discrete plastic multiplier
 * @param eps_k Equivalent plastic strain
 * @param G Second Lamé invariant
 * @return double
 */
static double __yield_function(const double *kappa_k, const double *kappa_n,
                               double J2, double d_gamma_k, double eps_k,
                               double G);
/*******************************************************/

/**
 * @brief
 *
 * @param d_kappa_k Hardening function derivative (isotropic,kinematic), t = k
 * @param G Second Lamé invariant
 * @return double
 */
static double __d_yield_function(const double *d_kappa_k, double G);
/*******************************************************/

/**
 * @brief
 *
 * @param T_ppal Stress tensor in the ppal coordinates
 * @param T_tr_vol Volumetric elastic stress tensor
 * @param T_tr_dev Deviatoric elastic stress tensor
 * @param T_back Back-stress (kinematic hardening)
 * @param n Plastic flow direction
 * @param d_gamma_k Discrete plastic multiplier
 * @param G Second Lamé invariant
 * @return int
 */
static int __elastoplastic_stress_ppal(double *T_ppal, const double *T_tr_vol,
                                       const double *T_tr_dev,
                                       const double *T_back, const double *n,
                                       double d_gamma_k, double G);
/*******************************************************/

/**
 * @brief
 *
 * @param Increment_E_plastic Increment plastic strain
 * @param T_xyz Nominal stress tensor
 * @param T_ppal Stress tensor in the ppal coordinates
 * @param T_back Back-stress (kinematic hardening)
 * @param eps_n1 Equivalent plastic strain
 * @param eigvec_b_e_tr Eigenvector of b elastic trial
 * @param n Plastic flow direction
 * @param d_gamma_k Discrete plastic multiplier
 * @param eps_k quivalent plastic strain
 * @param d_K_kin Increment of the kinematic hardening
 * @return int
 */
static int __update_internal_variables_plastic(
    double *Increment_E_plastic, double *T_xyz, const double *T_ppal,
    double *T_back, double *eps_n1, const double *eigvec_b_e_tr,
    const double *n, double d_gamma_k, double eps_k, double d_K_kin);
/*******************************************************/

/**
 * @brief
 *
 * @param C_ep Elastoplastic tangent moduli
 * @param n Plastic flow direction
 * @param kappa_k Hardening function (isotropic,kinematic)
 * @param d_gamma_k Discrete plastic multiplier
 * @param J2 Second invariant of the deviatoric stress tensor
 * @param K First Lamé invariant
 * @param G Second Lamé invariant
 * @return int
 */
static int __tangent_moduli(double *C_ep, const double *n,
                            const double *kappa_k, double d_gamma_k, double J2,
                            double K, double G);
/*******************************************************/

int compute_Kirchhoff_Stress_Von_Mises__Constitutive__(
    State_Parameters IO_State, Material MatProp) {

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
  double T_ppal[3] = {0.0, 0.0, 0.0};

  STATUS = __compute_trial_b_e(eigval_b_e_tr, eigvec_b_e_tr, IO_State.b_e,
                               IO_State.d_phi);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __compute_trial_b_e" RESET "\n");
    return EXIT_FAILURE;
  }

  E_hencky_trial[0] = 0.5 * log(eigval_b_e_tr[0]);
  E_hencky_trial[1] = 0.5 * log(eigval_b_e_tr[1]);
  E_hencky_trial[2] = 0.5 * log(eigval_b_e_tr[2]);

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

    //! Compute elastic stress tensor in the ppal directions
    __elastic_stress_ppal(T_ppal, T_tr_vol, T_tr_dev);

    STATUS = __update_internal_variables_elastic(IO_State.Stress, T_ppal,
                                                 eigvec_b_e_tr);
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

    __elastoplastic_stress_ppal(T_ppal, T_tr_vol, T_tr_dev, T_back, n,
                                d_gamma_k, G);

    STATUS = __update_internal_variables_plastic(
        Increment_E_plastic, IO_State.Stress, T_ppal, IO_State.Back_stress,
        IO_State.EPS, eigvec_b_e_tr, n, d_gamma_k, eps_k, d_K_kin);
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

  //! Update elastic left Cauchy-Green tensor
  STATUS = __corrector_b_e(IO_State.b_e, eigvec_b_e_tr, E_hencky_trial);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __corrector_b_e" RESET "\n");
    return EXIT_FAILURE;
  }

  //! Compute Elastoplastic tangent moduli
  if (IO_State.compute_C_ep) {
    STATUS = __tangent_moduli(IO_State.C_ep, n, kappa_k, d_gamma_k, J2, K, G);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __tangent_moduli" RESET "\n");
      return EXIT_FAILURE;
    }
  }

  //! Compute deformation energy
  *(IO_State.W) =
      0.5 * (T_ppal[0] * E_hencky_trial[0] + T_ppal[1] * E_hencky_trial[1] +
             T_ppal[2] * E_hencky_trial[2]);

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __compute_trial_b_e(double *eigval_b_e_tr, double *eigvec_b_e_tr,
                               const double *b_e, const double *d_phi) {

  unsigned Ndim = NumberDimensions;
  lapack_int n = NumberDimensions;
  lapack_int lda = NumberDimensions;

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

  lapack_int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec_b_e_tr,
                                  lda, eigval_b_e_tr);

  /* Check for convergence */
  if (info > 0) {
    fprintf(stderr,
            "" RED "Error in LAPACKE_dsyev(): %s\n %s; \n %i+1:N \n %s " RESET
            "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    fprintf(stderr,
            "" RED "Error in LAPACKE_dsyev(): the %i-th argument had an "
            "illegal value." RESET "\n",
            abs(info));
    return EXIT_FAILURE;
  }

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

/**************************************************************/

static int __elastic_stress_ppal(double *T_ppal, const double *T_tr_vol,
                                 const double *T_tr_dev) {

  T_ppal[0] = T_tr_vol[0] + T_tr_dev[0];
  T_ppal[1] = T_tr_vol[1] + T_tr_dev[1];
  T_ppal[2] = T_tr_vol[2] + T_tr_dev[2];

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __update_internal_variables_elastic(double *T_xyz,
                                               const double *T_ppal,
                                               const double *eigvec_b_e_tr) {

  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  double T_aux[4] = {0.0, 0.0, 0.0, 0.0};
#else
  double T_aux[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif

  for (unsigned A = 0; A < Ndim; A++) {

    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        T_aux[i * Ndim + j] += T_ppal[A] * eigvec_b_e_tr[A + i * Ndim] *
                               eigvec_b_e_tr[A + j * Ndim];
      }
    }
  }

#if NumberDimensions == 2
  T_xyz[0] = T_aux[0];
  T_xyz[1] = T_aux[1];
  T_xyz[2] = T_aux[2];
  T_xyz[3] = T_aux[3];
  T_xyz[4] = T_ppal[2];
#else
  T_xyz[0] = T_aux[0];
  T_xyz[1] = T_aux[1];
  T_xyz[2] = T_aux[2];
  T_xyz[3] = T_aux[3];
  T_xyz[4] = T_aux[4];
  T_xyz[5] = T_aux[5];
  T_xyz[6] = T_aux[6];
  T_xyz[7] = T_aux[7];
  T_xyz[8] = T_aux[8];
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

static int __elastoplastic_stress_ppal(double *T_ppal, const double *T_tr_vol,
                                       const double *T_tr_dev,
                                       const double *T_back, const double *n,
                                       double d_gamma_k, double G) {

  T_ppal[0] = T_tr_vol[0] + T_tr_dev[0] + T_back[0] - d_gamma_k * 2 * G * n[0];
  T_ppal[1] = T_tr_vol[1] + T_tr_dev[1] + T_back[1] - d_gamma_k * 2 * G * n[1];
  T_ppal[2] = T_tr_vol[2] + T_tr_dev[2] + T_back[2] - d_gamma_k * 2 * G * n[2];
}

/**************************************************************/

static int __update_internal_variables_plastic(
    double *Increment_E_plastic, double *T_xyz, const double *T_ppal,
    double *T_back, double *eps_n1, const double *eigvec_b_e_tr,
    const double *n, double d_gamma_k, double eps_k, double d_K_kin) {

  int Ndim = NumberDimensions;

  // Update equivalent plastic strain
  *eps_n1 = eps_k;

  // Update the increment of the plastic strain tensor
  Increment_E_plastic[0] = d_gamma_k * n[0];
  Increment_E_plastic[1] = d_gamma_k * n[1];
  Increment_E_plastic[2] = d_gamma_k * n[2];

#if NumberDimensions == 2
  double T_aux[4] = {0.0, 0.0, 0.0, 0.0};
#else
  double T_aux[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif

  for (unsigned A = 0; A < Ndim; A++) {

    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        T_aux[i * Ndim + j] += T_ppal[A] * eigvec_b_e_tr[A * Ndim + i] *
                               eigvec_b_e_tr[A * Ndim + j];
      }
    }
  }

#if NumberDimensions == 2
  T_xyz[0] = T_aux[0];
  T_xyz[1] = T_aux[1];
  T_xyz[2] = T_aux[2];
  T_xyz[3] = T_aux[3];
  T_xyz[4] = T_ppal[2];
#else
  T_xyz[0] = T_aux[0];
  T_xyz[1] = T_aux[1];
  T_xyz[2] = T_aux[2];
  T_xyz[3] = T_aux[3];
  T_xyz[4] = T_aux[4];
  T_xyz[5] = T_aux[5];
  T_xyz[6] = T_aux[6];
  T_xyz[7] = T_aux[7];
  T_xyz[8] = T_aux[8];
#endif

  // Update back stress
  T_back[0] += sqrt(2. / 3.) * d_K_kin * n[0];
  T_back[1] += sqrt(2. / 3.) * d_K_kin * n[1];
  T_back[2] += sqrt(2. / 3.) * d_K_kin * n[2];

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __tangent_moduli(double *C_ep, const double *n,
                            const double *kappa_k, double d_gamma_k, double J2,
                            double K, double G) {

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;

#if NumberDimensions == 2
  double R2_Identity[2] = {1.0, 1.0};

  double R4_Identity[2][2] = {{1.0, 0.0}, {0.0, 1.0}};
#else
  double R2_Identity[3] = {1.0, 1.0, 1.0};

  double R4_Identity[3][3] = {
      {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
#endif

  double K_iso_k = kappa_k[0];
  double K_kin_k = kappa_k[1];
  double theta = 0.0;
  if (J2 > TOL_NR) {
    theta = 1.0 - 2.0 * G * d_gamma_k / J2;
  }

  double theta_bar =
      1.0 / (1.0 + (K_iso_k + K_kin_k) / (3.0 * G)) - (1.0 - theta);

  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      C_ep[i * Ndim + j] = K * R2_Identity[i] * R2_Identity[j] +
                           2.0 * G * theta *
                               (R4_Identity[i][j] -
                                (1.0 / 3.0) * R2_Identity[i] * R2_Identity[j]) -
                           2.0 * G * theta_bar * n[i] * n[j];
    }
  }

  return STATUS;
}

/**************************************************************/