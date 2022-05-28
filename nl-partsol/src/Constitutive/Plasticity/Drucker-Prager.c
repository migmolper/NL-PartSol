/**
 * @file Drucker-Prager.c
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2022-05-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Constitutive/Plasticity/Drucker-Prager.h"

/**************************************************************/
/*!
    \param[out] eigval_b_e_tr Eigenvalues of b elastic trial.
    \param[out] eigvec_b_e_tr Eigenvector of b elastic trial.
    \param[in] b_e (n) Elastic left Cauchy-Green.
    \param[in] d_phi Incremental deformation gradient.
*/
static int __compute_trial_b_e(double *eigval_b_e_tr, double *eigvec_b_e_tr,
                               const double *b_e, const double *d_phi);
/**************************************************************/

/*!
    \param[out] b_e (n+1) Elastic deformation gradient.
    \param[in] eigvec_b_e_tr Eigenvector of b elastic trial.
    \param[in] E_hencky_trial Corrected Henky strain
*/
static int __corrector_b_e(double *b_e, const double *eigvec_b_e_tr,
                           const double *E_hencky_trial);
/**************************************************************/

/*!
    \param[out] T_tr_vol Volumetric elastic stress tensor.
    \param[out] T_tr_dev Deviatoric elastic stress tensor.
    \param[out] pressure First invariant of the stress tensor
    \param[out] J2 Second invariant of the deviatoric stress tensor
    \param[in] E_hencky_trial, Trial elastic strain tensor.
    \param[in] K First Lamé invariant.
    \param[in] G Second Lamé invariant.
    \param[in] p_ref Reference pressure.
*/
static int __trial_elastic(double *T_tr_vol, double *T_tr_dev, double *pressure,
                           double *J2, const double *E_hencky_trial, double K,
                           double G, double p_ref);
/**************************************************************/

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

/*!
    \param[out] T_xyz Nominal stress tensor
    \param[in] T_ppal Stress tensor in the ppal coordinates
    \param[in] eigvec_b_e_tr Eigenvector of b elastic trial.
*/
static int __update_internal_variables_elastic(double *T_xyz,
                                               const double *T_ppal,
                                               const double *eigvec_b_e_tr);
/**************************************************************/

/*!
    \param[out] C_ep Elastoplastic tanget moduli
    \param[in] K First Lamé invariant.
    \param[in] G Second Lamé invariant.
*/
static int __tangent_moduli_elastic(double *C_ep, double K, double G);
/**************************************************************/

/*!
    \param[out] n Plastic flow direction
    \param[in] T_tr_dev Deviatoric elastic stress tensor
    \param[in] J2 Second invariant of the deviatoric stress tensor
*/
static int __compute_plastic_flow_direction(double *n, const double *T_tr_dev,
                                            double J2);
/**************************************************************/

/*!
    \param[out] eps_k Equivalent plastic strain
    \param[in] d_gamma_k Discrete plastic multiplier
    \param[in] eps_n Equivalent plastic strain in the last step
    \param[in] alpha_Q Plastic potential parameter
*/
static int __eps(double *eps_k, double d_gamma_k, double eps_n, double alpha_Q);
/**************************************************************/

/*!
    \param[out] kappa_k Hardening function.
    \param[in] kappa_0 Reference hardening
    \param[in] exp_param Hardening exponential
    \param[in] eps_k Equivalent plastic strain
    \param[in] eps_0 Reference plastic strain
*/
static int __kappa(double *kappa_k, double kappa_0, double exp_param,
                   double eps_k, double eps_0);
/**************************************************************/

/*!
    \param[out] d_kappa Derivative of the hardening function.
    \param[in] kappa_0 Reference hardening
    \param[in] eps_k Equivalent plastic strain
    \param[in] eps_0 Reference plastic strain
    \param[in] exp_param Hardening exponential
*/
static int __d_kappa(double *d_kappa, double kappa_0, double eps_k,
                     double eps_0, double exp_param);
/**************************************************************/

/*!
    \param[out] pressure_limit Limit for the apex region.
    \param[in] J2 Second invariant of the deviatoric stress tensor
    \param[in] kappa_n Hardening function.
    \param[in] d_kappa Derivative of the hardening function.
    \param[in] K First Lamé invariant.
    \param[in] G Second Lamé invariant.
    \param[in] alpha_F Yield surface parameter I.
    \param[in] alpha_Q Plastic potential parameter.
    \param[in] beta Yield surface parameter II.
*/
static int __compute_pressure_limit(double *pressure_limit, double J2,
                                    double kappa_n, double d_kappa, double K,
                                    double G, double alpha_F, double alpha_Q,
                                    double beta);
/**************************************************************/

/*!
    \param[in] pressure First invariant of the stress tensor
    \param[in] J2 Second invariant of the deviatoric stress tensor
    \param[in] d_gamma_k Discrete plastic multiplier
    \param[in] kappa_k Hardening function.
    \param[in] alpha_F Yield surface parameter I.
    \param[in] alpha_Q Plastic potential parameter.
    \param[in] beta Yield surface parameter II.
    \param[in] K First Lamé invariant.
    \param[in] G Second Lamé invariant.
*/
static double __yield_function_classical(double pressure, double J2,
                                         double d_gamma_k, double kappa_k,
                                         double alpha_F, double alpha_Q,
                                         double beta, double K, double G);
/**************************************************************/

/*!
    \param[in] d_kappa_k Derivative of the hardening function.
    \param[in] K First Lamé invariant.
    \param[in] G Second Lamé invariant.
    \param[in] alpha_F Yield surface parameter I.
    \param[in] alpha_Q Plastic potential parameter.
    \param[in] beta Yield surface parameter II.
*/
static double __d_yield_function_classical(double d_kappa_k, double K, double G,
                                           double alpha_F, double alpha_Q,
                                           double beta);
/**************************************************************/

/**
 * @brief
 *
 * @param T_ppal Stress tensor in the ppal coordinates
 * @param T_tr_vol Volumetric elastic stress tensor
 * @param T_tr_dev Deviatoric elastic stress tensor
 * @param n Plastic flow direction
 * @param d_gamma_k Discrete plastic multiplier
 * @param K First Lamé invariant
 * @param G Second Lamé invariant
 * @param alpha_Q Plastic potential parameter
 * @return int
 */
static int __elastoplastic_stress_ppal_classical(
    double *T_ppal, const double *T_tr_vol, const double *T_tr_dev,
    const double *n, double d_gamma_k, double K, double G, double alpha_Q);
/**************************************************************/

/**
 * @brief
 *
 * @param Increment_E_plastic Increment plastic strain
 * @param T_xyz Nominal stress tensor
 * @param eps_n1 Equivalent plastic strain
 * @param kappa_n1 Hardening function
 * @param T_ppal Stress tensor in the ppal coordinates
 * @param eigvec_b_e_tr Eigenvector of b elastic trial
 * @param n Plastic flow direction
 * @param d_gamma_k Discrete plastic multiplier
 * @param alpha_Q Plastic potential parameter
 * @param eps_k Equivalent plastic strain
 * @param kappa_k Hardening function
 * @return int
 */
static int __update_internal_variables_classical(
    double *Increment_E_plastic, double *T_xyz, double *eps_n1,
    double *kappa_n1, const double *T_ppal, const double *eigvec_b_e_tr,
    const double *n, double d_gamma_k, double alpha_Q, double eps_k,
    double kappa_k);

/**************************************************************/

/*!
  \param[out] C_ep  Elastoplastic tanget moduli
  \param[in] n  Plastic flow direction.
  \param[in] d_gamma_k  Derivative of the hardening function.
  \param[in] J2  Second invariant of the deviatoric stress tensor
  \param[in] d_kappa_k  Discrete plastic multiplier
  \param[in] K  First Lamé invariant.
  \param[in] G  Second Lamé invariant.
  \param[in] beta Yield surface parameter II.
  \param[in] alpha_F  Yield surface parameter I.
  \param[in] alpha_Q  Plastic potential parameter.
*/
static int __tangent_moduli_classical(double *C_ep, const double *n,
                                      double d_gamma_k, double J2,
                                      double d_kappa_k, double K, double G,
                                      double beta, double alpha_F,
                                      double alpha_Q);
/**************************************************************/

/*!
    \param[in] pressure First invariant of the stress tensor
    \param[in] d_gamma_k Discrete plastic multiplier
    \param[in] d_gamma_1 Discrete plastic multiplier I
    \param[in] kappa_k Hardening function.
    \param[in] d_kappa_k Discrete plastic multiplier
    \param[in] K First Lamé invariant.
    \param[in] alpha_F Yield surface parameter I.
    \param[in] alpha_Q Plastic potential parameter.
    \param[in] beta Yield surface parameter II.
*/
static double __yield_function_apex(double pressure, double d_gamma_k,
                                    double d_gamma_1, double kappa_k,
                                    double d_kappa_k, double K, double alpha_F,
                                    double alpha_Q, double beta);
/**************************************************************/

/*!
    \param[in] d_gamma_k Discrete plastic multiplier
    \param[in] d_gamma_1 Discrete plastic multiplier I
    \param[in] d_kappa_k Derivative of the hardening function
    \param[in] K First Lamé invariant.
    \param[in] alpha_F Yield surface parameter I.
    \param[in] alpha_Q Plastic potential parameter.
    \param[in] beta Yield surface parameter II.
*/
static double __d_yield_function_apex(double d_gamma_k, double d_gamma_1,
                                      double d_kappa_k, double K,
                                      double alpha_F, double alpha_Q,
                                      double beta);
/**************************************************************/

/**
 * @brief
 *
 * @param T_ppal
 * @param T_tr_vol Volumetric elastic stress tensor.
 * @param d_gamma_k Discrete plastic multiplier
 * @param K First Lamé invariant.
 * @param alpha_Q Plastic potential parameter.
 * @return int
 */
static int __elastoplastic_stress_ppal_apex(double *T_ppal,
                                            const double *T_tr_vol,
                                            double d_gamma_k, double K,
                                            double alpha_Q);
/**************************************************************/

/**
 * @brief
 *
 * @param Increment_E_plastic Increment plastic strain
 * @param T_xyz Nominal stress tensor
 * @param eps_n1 Equivalent plastic strain
 * @param kappa_n1 Hardening function
 * @param T_ppal Stress tensor in the ppal coordinates
 * @param eigvec_b_e_tr Eigenvector of b elastic trial
 * @param n Plastic flow direction
 * @param d_gamma_k Discrete plastic multiplier
 * @param d_gamma_1 Discrete plastic multiplier I
 * @param alpha_Q Plastic potential parameter
 * @param eps_k Equivalent plastic strain
 * @param kappa_k Hardening function.
 * @return int
 */
static int __update_internal_variables_apex(
    double *Increment_E_plastic, double *T_xyz, double *eps_n1,
    double *kappa_n1, const double *T_ppal, const double *eigvec_b_e_tr,
    const double *n, double d_gamma_k, double d_gamma_1, double alpha_Q,
    double eps_k, double kappa_k);
/**************************************************************/

/*!
  \param[out] C_ep Elastoplastic tanget moduli
  \param[in] n Plastic flow direction
  \param[in] d_gamma_k Discrete plastic multiplier
  \param[in] d_gamma_1 Discrete plastic multiplier I
  \param[in] d_kappa_k Derivative of the hardening function
  \param[in] K First Lamé invariant
  \param[in] G Second Lamé invariant
  \param[in] beta Yield surface parameter II
  \param[in] alpha_F Yield surface parameter I
  \param[in] alpha_Q Plastic potential parameter
*/
static int __tangent_moduli_apex(double *C_ep, const double *n,
                                 double d_gamma_k, double d_gamma_1,
                                 double d_kappa_k, double K, double G,
                                 double beta, double alpha_F, double alpha_Q);

/**************************************************************/

int compute_Kirchhoff_Stress_Drucker_Prager__Constitutive__(
    State_Parameters IO_State, Material MatProp)
/*
  Backward Euler algorithm for the Drucker-Prager (Lorenzo Sanavia)
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
  double p_ref = MatProp.ReferencePressure;
  double rad_friction_angle =
      (PI__MatrixLib__ / 180.0) * MatProp.phi_Frictional;
  double rad_dilatancy_angle =
      (PI__MatrixLib__ / 180.0) * MatProp.psi_Frictional;
  double exp_param = MatProp.Exponent_Hardening_Ortiz;
  double kappa_0 = MatProp.kappa_0;
  double eps_0 = MatProp.Plastic_Strain_0;

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
  double eps_n = *IO_State.EPS;
  double eps_k = eps_n;
  double kappa_k = *IO_State.Kappa;
  double d_kappa_k = 0.0;

  // Initialise solver parameters
  double TOL = TOL_Radial_Returning;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;

  STATUS = __trial_elastic(T_tr_vol, T_tr_dev, &pressure, &J2, E_hencky_trial,
                           K, G, p_ref);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __trial_elastic" RESET "\n");
    return EXIT_FAILURE;
  }

  PHI = PHI_0 = __yield_function_classical(pressure, J2, d_gamma_k, kappa_k,
                                           alpha_F, alpha_Q, beta, K, G);

  // Elastic
  if (PHI_0 <= TOL_NR) {

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

    if (IO_State.compute_C_ep) {
      STATUS = __tangent_moduli_elastic(IO_State.C_ep, K, G);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __tangent_moduli_elastic()" RESET "\n");
        return EXIT_FAILURE;
      }
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

    STATUS = __d_kappa(&d_kappa_k, kappa_0, eps_n, eps_0, exp_param);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __d_kappa" RESET "\n");
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
        if (fabs(d_PHI) < TOL) {
          fprintf(stderr,
                  "" RED "|d_PHI| = %f < TOL (classical loop)" RESET "\n",
                  fabs(d_PHI));
          return EXIT_FAILURE;
        }

        d_gamma_k += -PHI / d_PHI;
        if (d_gamma_k < 0.0) {
          fprintf(stderr,
                  "" RED "d_gamma_k = %f < 0 (classical loop)" RESET "\n",
                  d_gamma_k);
          return EXIT_FAILURE;
        }

        STATUS = __eps(&eps_k, d_gamma_k, eps_n, alpha_Q);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr, "" RED "Error in __eps (classical loop)" RESET "\n");
          return EXIT_FAILURE;
        }

        STATUS = __kappa(&kappa_k, kappa_0, exp_param, eps_k, eps_0);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr, "" RED "Error in __kappa" RESET "\n");
          return EXIT_FAILURE;
        }

        STATUS = __d_kappa(&d_kappa_k, kappa_0, eps_k, eps_0, exp_param);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr, "" RED "Error in __d_kappa" RESET "\n");
          return EXIT_FAILURE;
        }

        PHI = __yield_function_classical(pressure, J2, d_gamma_k, kappa_k,
                                         alpha_F, alpha_Q, beta, K, G);
      }

      __elastoplastic_stress_ppal_classical(T_ppal, T_tr_vol, T_tr_dev, n,
                                            d_gamma_k, K, G, alpha_Q);

      STATUS = __update_internal_variables_classical(
          Increment_E_plastic, IO_State.Stress, IO_State.EPS, IO_State.Kappa,
          T_ppal, eigvec_b_e_tr, n, d_gamma_k, alpha_Q, eps_k, kappa_k);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __update_internal_variables_classical" RESET
                "\n");
        return EXIT_FAILURE;
      }

      if (IO_State.compute_C_ep) {
        STATUS =
            __tangent_moduli_classical(IO_State.C_ep, n, d_gamma_k, J2,
                                       d_kappa_k, K, G, beta, alpha_F, alpha_Q);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr,
                  "" RED "Error in __tangent_moduli_classical" RESET "\n");
          return EXIT_FAILURE;
        }
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
        if (fabs(d_PHI) < TOL) {
          break;
        }

        d_gamma_2_k += -PHI / d_PHI;

        if (d_gamma_2_k < 0.0) {
          d_gamma_k = 0.0;
          d_gamma_2_k = 0.0;
          break;
        } else {
          d_gamma_k = d_gamma_1 + d_gamma_2_k;
        }

        PHI = __yield_function_apex(pressure, d_gamma_k, d_gamma_1, kappa_k,
                                    d_kappa_k, K, alpha_F, alpha_Q, beta);
      }

      STATUS = __eps(&eps_k, d_gamma_k, eps_n, alpha_Q);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __eps (apex loop)" RESET "\n");
        return EXIT_FAILURE;
      }

      __elastoplastic_stress_ppal_apex(T_ppal, T_tr_vol, d_gamma_k, K, alpha_Q);

      STATUS = __update_internal_variables_apex(
          Increment_E_plastic, IO_State.Stress, IO_State.EPS, IO_State.Kappa,
          T_ppal, eigvec_b_e_tr, n, d_gamma_k, d_gamma_1, alpha_Q, eps_k,
          kappa_k);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in__update_internal_variables_apex" RESET "\n");
        return EXIT_FAILURE;
      }

      if (IO_State.compute_C_ep) {
        STATUS = __tangent_moduli_apex(IO_State.C_ep, n, d_gamma_k, d_gamma_1,
                                       d_kappa_k, K, G, beta, alpha_F, alpha_Q);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr, "" RED "Error in __tangent_moduli_apex" RESET "\n");
          return EXIT_FAILURE;
        }
      }
    }
  }

  // Plastic corrector step for the left Cauchy-Green tensor
  E_hencky_trial[0] -= Increment_E_plastic[0];
  E_hencky_trial[1] -= Increment_E_plastic[1];
  E_hencky_trial[2] -= Increment_E_plastic[2];

  //! Compute deformation energy
  *(IO_State.W) =
      0.5 * (T_ppal[0] * E_hencky_trial[0] + T_ppal[1] * E_hencky_trial[1] +
             T_ppal[2] * E_hencky_trial[2]);

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

  T_tr_vol[0] = -p_ref - K * E_hencky_trial_vol[0];
  T_tr_vol[1] = -p_ref - K * E_hencky_trial_vol[1];
  T_tr_vol[2] = -p_ref - K * E_hencky_trial_vol[2];

  T_tr_dev[0] = 2 * G * (E_hencky_trial[0] - E_hencky_trial_vol[0]);
  T_tr_dev[1] = 2 * G * (E_hencky_trial[1] - E_hencky_trial_vol[1]);
  T_tr_dev[2] = 2 * G * (E_hencky_trial[2] - E_hencky_trial_vol[2]);

  *pressure = (T_tr_vol[0] + T_tr_vol[1] + T_tr_vol[2]) / 3.0;
  *J2 = sqrt(T_tr_dev[0] * T_tr_dev[0] + T_tr_dev[1] * T_tr_dev[1] +
             T_tr_dev[2] * T_tr_dev[2]);

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __elastic_stress_ppal(double *T_ppal, const double *T_tr_vol,
                                 const double *T_tr_dev) {

  T_ppal[0] = -T_tr_vol[0] + T_tr_dev[0];
  T_ppal[1] = -T_tr_vol[1] + T_tr_dev[1];
  T_ppal[2] = -T_tr_vol[2] + T_tr_dev[2];

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __update_internal_variables_elastic(double *T_xyz,
                                               const double *T_ppal,
                                               const double *eigvec_b_e_tr) {

  int Ndim = NumberDimensions;

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

/***************************************************************************/

static int __compute_plastic_flow_direction(double *n, const double *T_tr_dev,
                                            double J2) {

  if (J2 > TOL_NR) {
    n[0] = T_tr_dev[0] / J2;
    n[1] = T_tr_dev[1] / J2;
    n[2] = T_tr_dev[2] / J2;
  } else {
    n[0] = 0.0;
    n[1] = 0.0;
    n[2] = 0.0;
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __eps(double *eps_k, double d_gamma_k, double eps_n,
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

static int __kappa(double *kappa_k, double kappa_0, double exp_param,
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

static int __d_kappa(double *d_kappa, double kappa_0, double eps_k,
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

static int __elastoplastic_stress_ppal_classical(
    double *T_ppal, const double *T_tr_vol, const double *T_tr_dev,
    const double *n, double d_gamma_k, double K, double G, double alpha_Q) {

  T_ppal[0] =
      -T_tr_vol[0] + T_tr_dev[0] + d_gamma_k * (3 * K * alpha_Q - 2 * G * n[0]);
  T_ppal[1] =
      -T_tr_vol[1] + T_tr_dev[1] + d_gamma_k * (3 * K * alpha_Q - 2 * G * n[1]);
  T_ppal[2] =
      -T_tr_vol[2] + T_tr_dev[2] + d_gamma_k * (3 * K * alpha_Q - 2 * G * n[2]);

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __update_internal_variables_classical(
    double *Increment_E_plastic, double *T_xyz, double *eps_n1,
    double *kappa_n1, const double *T_ppal, const double *eigvec_b_e_tr,
    const double *n, double d_gamma_k, double alpha_Q, double eps_k,
    double kappa_k) {

  int Ndim = NumberDimensions;

  // Update hardening parameters
  *eps_n1 = eps_k;
  *kappa_n1 = kappa_k;

  // Update the increment of the plastic strain tensor
  Increment_E_plastic[0] = d_gamma_k * (alpha_Q + n[0]);
  Increment_E_plastic[1] = d_gamma_k * (alpha_Q + n[1]);
  Increment_E_plastic[2] = d_gamma_k * (alpha_Q + n[2]);

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

static int __elastoplastic_stress_ppal_apex(double *T_ppal,
                                            const double *T_tr_vol,
                                            double d_gamma_k, double K,
                                            double alpha_Q) {
  T_ppal[0] = -T_tr_vol[0] + d_gamma_k * 3 * K * alpha_Q;
  T_ppal[1] = -T_tr_vol[1] + d_gamma_k * 3 * K * alpha_Q;
  T_ppal[2] = -T_tr_vol[2] + d_gamma_k * 3 * K * alpha_Q;

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __update_internal_variables_apex(
    double *Increment_E_plastic, double *T_xyz, double *eps_n1,
    double *kappa_n1, const double *T_ppal, const double *eigvec_b_e_tr,
    const double *n, double d_gamma_k, double d_gamma_1, double alpha_Q,
    double eps_k, double kappa_k) {

  int Ndim = NumberDimensions;

  // Update hardening parameters
  *eps_n1 = eps_k;
  *kappa_n1 = kappa_k;

  // Update the increment of the plastic strain tensor
  Increment_E_plastic[0] = d_gamma_k * alpha_Q + d_gamma_1 * n[0];
  Increment_E_plastic[1] = d_gamma_k * alpha_Q + d_gamma_1 * n[1];
  Increment_E_plastic[2] = d_gamma_k * alpha_Q + d_gamma_1 * n[2];

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

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __tangent_moduli_elastic(double *C_ep, double K, double G) {

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

  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      C_ep[i * Ndim + j] = K * R2_Identity[i] * R2_Identity[j] +
                           2.0 * G *
                               (R4_Identity[i][j] -
                                (1.0 / 3.0) * R2_Identity[i] * R2_Identity[j]);
    }
  }

  return STATUS;
}

/**************************************************************/

static int __tangent_moduli_classical(double *C_ep, const double *n,
                                      double d_gamma_k, double J2,
                                      double d_kappa_k, double K, double G,
                                      double beta, double alpha_F,
                                      double alpha_Q) {

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

  double c0 = 9 * alpha_F * alpha_Q * K + 2 * G +
              beta * d_kappa_k * sqrt(2. / 3. * (1 + 3 * alpha_Q * alpha_Q));
  double c1 = 1.0 - 9.0 * alpha_F * alpha_Q * K / c0;
  double c2 = 0.0;
  if (J2 > TOL_NR) {
    c2 = d_gamma_k / J2;
  }

  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      C_ep[i * Ndim + j] =
          c1 * K * R2_Identity[i] * R2_Identity[j] +
          2 * G *
              (R4_Identity[i][j] - (1. / 3.) * (1.0 - 2.0 * G * c2) *
                                       R2_Identity[i] * R2_Identity[j]) -
          (6.0 * alpha_Q * K * G / c0) * R2_Identity[i] * n[j] -
          (6.0 * alpha_Q * K * G / c0) * n[i] * R2_Identity[j] -
          4 * G * G * (1.0 / c0 - c2) * n[i] * n[j];
    }
  }

  return STATUS;
}

/***************************************************************************/

static int __tangent_moduli_apex(double *C_ep, const double *n,
                                 double d_gamma_k, double d_gamma_1,
                                 double d_kappa_k, double K, double G,
                                 double beta, double alpha_F, double alpha_Q) {

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;

#if NumberDimensions == 2
  double R2_Identity[2] = {1.0, 1.0};
#else
  double R2_Identity[3] = {1.0, 1.0, 1.0};
#endif

  double c0 = 0.0;
  double c1 = 0.0;

  if (d_gamma_k > 0.0) {
    c0 = (alpha_Q * beta * sqrt(2. / 3.) * d_kappa_k * d_gamma_k) /
         (3.0 * alpha_F * K *
              sqrt(d_gamma_1 * d_gamma_1 +
                   3.0 * alpha_Q * alpha_Q * d_gamma_k * d_gamma_k) +
          alpha_Q * beta * sqrt(2. / 3.) * d_kappa_k * d_gamma_k);

    c1 = c0 * K / (2.0 * alpha_Q * G * d_gamma_k);
  }

  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      C_ep[i * Ndim + j] =
          c0 * K * R2_Identity[i] * R2_Identity[j] + c1 * R2_Identity[i] * n[j];
    }
  }

  return STATUS;
}

/**************************************************************/