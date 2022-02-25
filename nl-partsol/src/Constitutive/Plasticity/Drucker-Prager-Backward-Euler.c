#include "nl-partsol.h"

static int __trial_elastic(
    double *T_trial_vol /**< [in/out] Volumetric elastic stress tensor. */,
    double *T_trial_dev /**< [in/out] Deviatoric elastic stress tensor. */,
    double *pressure /**< [out] First invariant of the stress tensor */,
    double *J2 /**< [out] Second invariant of the deviatoric stress tensor */,
    const double *E_trial, /**< [in] Trial elastic strain tensor. */
    double K /**< [in] First Lamé invariant. */,
    double G /**< [in] Second Lamé invariant. */,
    double p_ref /**< [in] Reference pressure. */);

static int __update_internal_variables_elastic(
    State_Parameters IO_State /**< [in/out] List of internal variables */,
    const double *T_trial_vol /**< [in] Volumetric elastic stress tensor */,
    const double *T_trial_dev /**< [in] Deviatoric elastic stress tensor */);

static int __compute_plastic_flow_direction(
    double *n /**< [out] Plastic flow direction */,
    const double *T_trial_dev /**< [in] Deviatoric elastic stress tensor */,
    double J2 /**< [in] Second invariant of the deviatoric stress tensor */);

static double __compute_eps(
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
    const double *T_trial_vol /**< [in] Volumetric elastic stress tensor. */,
    const double *T_trial_dev /**< [in] Deviatoric elastic stress tensor. */,
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
    const double *T_trial_vol /**< [in] Volumetric elastic stress tensor. */,
    const double *n /**< [out] Plastic flow direction. */,
    double d_gamma_k /**< [in] Discrete plastic multiplier */,
    double d_gamma_1 /**< [in] Discrete plastic multiplier I */,
    double alpha_Q /**< [in] Plastic potential parameter. */,
    double K /**< [in] First Lamé invariant. */,
    double eps_k /**< [in] Equivalent plastic strain*/,
    double kappa_k /**< [in] Hardening function. */);

/**************************************************************/

int Drucker_Prager_backward_euler(State_Parameters IO_State, Material MatProp)
/*
  Backward Euler algorithm for the Drucker-Prager (Lorenzo Sanavia)
*/
{

  int Ndim = NumberDimensions;

  // Read input/output parameters
  double *E_trial = IO_State.Strain;
  double T_trial_vol[3] = {0.0, 0.0, 0.0};
  double T_trial_dev[3] = {0.0, 0.0, 0.0};

  // Material parameters
  double K = MatProp.E / (3 * (1 - 2 * MatProp.nu));
  double G = MatProp.E / (2 * (1 + MatProp.nu));
  double p_ref = MatProp.ReferencePressure;
  double rad_friction_angle = (PI__MatrixLib__ / 180) * MatProp.phi_Frictional;
  double rad_dilatancy_angle = (PI__MatrixLib__ / 180) * MatProp.psi_Frictional;
  double exp_param = MatProp.Exponent_Hardening_Ortiz;
  double kappa_0 = MatProp.yield_stress_0;
  double eps_0 = MatProp.Reference_Plastic_Strain_Ortiz;

#if NumberDimensions == 2
  double alpha_F = sqrt(2 / 3.) * tan(rad_friction_angle) /
                   sqrt(3 + 4 * DSQR(tan(rad_friction_angle)));
  double alpha_Q = sqrt(2 / 3.) * tan(rad_dilatancy_angle) /
                   sqrt(3 + 4 * DSQR(tan(rad_dilatancy_angle)));
  double beta = sqrt(2 / 3.) * 3 / sqrt(3 + 4 * DSQR(tan(rad_friction_angle)));
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

  // Initialise solver parameters
  double TOL = TOL_Radial_Returning;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;
  bool Convergence = false;

  __trial_elastic(T_trial_vol, T_trial_dev, &pressure, &J2, E_trial, K, G,
                  p_ref);

  PHI_0 = __yield_function_classical(pressure, J2, d_gamma_k, kappa_k, alpha_F,
                                     alpha_Q, beta, K, G);
  // Elastic
  if (PHI <= 0.0) {

    __update_internal_variables_elastic(IO_State, T_trial_vol, T_trial_dev);
  }
  // Plastic (check yield condition)
  else {

    PHI = PHI_0;

    __compute_plastic_flow_direction(n, T_trial_dev, J2);

    __compute_d_kappa(&d_kappa_k, kappa_0, eps_n, eps_0, exp_param);

    __compute_pressure_limit(&pressure_limit, J2, kappa_k, d_kappa_k, K, G,
                             alpha_F, alpha_Q, beta);

    // Classic radial returning algorithm
    if (-pressure < pressure_limit) {

      while (fabs(PHI / PHI_0) >= TOL) {

        Iter++;

        if (Iter == MaxIter)
          break;

        d_PHI = __d_yield_function_classical(d_kappa_k, K, G, alpha_F, alpha_Q,
                                             beta);

        d_gamma_k += -PHI / d_PHI;

        eps_k = __compute_eps(d_gamma_k, eps_n, alpha_Q);

        __compute_kappa(&kappa_k, kappa_0, exp_param, eps_k, eps_0);

        if (kappa_k > 0) {
          __compute_d_kappa(&d_kappa_k, kappa_0, eps_k, eps_0, exp_param);
        } else {
          kappa_k = 0.0;
          d_kappa_k = 0.0;
        }

        PHI = __yield_function_classical(pressure, J2, d_gamma_k, kappa_k,
                                         alpha_F, alpha_Q, beta, K, G);
      }

      __update_internal_variables_classical(IO_State, T_trial_vol, T_trial_dev,
                                            n, d_gamma_k, alpha_Q, K, G, eps_k,
                                            kappa_k);

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

      eps_k = __compute_eps(d_gamma_k, eps_n, alpha_Q);

      __update_internal_variables_apex(IO_State, T_trial_vol, n, d_gamma_k,
                                       d_gamma_1, alpha_Q, K, eps_k, kappa_k);
    }
  }

  return EXIT_SUCCESS;
}
/***************************************************************************/

static int __trial_elastic(double *T_trial_vol, double *T_trial_dev,
                           double *pressure, double *J2, const double *E_trial,
                           double K, double G, double p_ref) {

  double E_trial_vol[3];
  double tr_E_trial = E_trial[0] + E_trial[1] + E_trial[2];

  E_trial_vol[0] = (1.0 / 3.0) * tr_E_trial;
  E_trial_vol[0] = (1.0 / 3.0) * tr_E_trial;
  E_trial_vol[0] = (1.0 / 3.0) * tr_E_trial;

  T_trial_vol[0] = -K * E_trial_vol[0] - p_ref;
  T_trial_vol[1] = -K * E_trial_vol[1] - p_ref;
  T_trial_vol[2] = -K * E_trial_vol[2] - p_ref;

  T_trial_dev[0] = 2 * G * (E_trial[0] - E_trial_vol[0]);
  T_trial_dev[1] = 2 * G * (E_trial[1] - E_trial_vol[1]);
  T_trial_dev[2] = 2 * G * (E_trial[2] - E_trial_vol[2]);

  *pressure = (T_trial_vol[0] + T_trial_vol[1] + T_trial_vol[2]) / 3.0;
  *J2 = sqrt(T_trial_dev[0] * T_trial_dev[0] + T_trial_dev[1] * T_trial_dev[1] +
             T_trial_dev[2] * T_trial_dev[2]);

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __update_internal_variables_elastic(State_Parameters IO_State,
                                               const double *T_trial_vol,
                                               const double *T_trial_dev) {

  IO_State.Increment_E_plastic[0] = 0.0;
  IO_State.Increment_E_plastic[1] = 0.0;
  IO_State.Increment_E_plastic[2] = 0.0;

  IO_State.Stress[0] = -T_trial_vol[0] + T_trial_dev[0];
  IO_State.Stress[1] = -T_trial_vol[1] + T_trial_dev[1];
  IO_State.Stress[2] = -T_trial_vol[2] + T_trial_dev[2];

  return EXIT_SUCCESS;
}

/***************************************************************************/

static int __compute_plastic_flow_direction(double *n,
                                            const double *T_trial_dev,
                                            double J2) {

  n[0] = T_trial_dev[0] / J2;
  n[1] = T_trial_dev[1] / J2;
  n[2] = T_trial_dev[2] / J2;

  return EXIT_SUCCESS;
}

/**************************************************************/

static double __compute_eps(double d_gamma_k, double eps_n, double alpha_Q) {

  double eps_k = eps_n + d_gamma_k * sqrt(3 * alpha_Q * alpha_Q + 1);

  return eps_k;
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
    State_Parameters IO_State, const double *T_trial_vol,
    const double *T_trial_dev, const double *n, double d_gamma_k,
    double alpha_Q, double K, double G, double eps_k, double kappa_k) {

  IO_State.Increment_E_plastic[0] = d_gamma_k * (alpha_Q + n[0]);
  IO_State.Increment_E_plastic[1] = d_gamma_k * (alpha_Q + n[1]);
  IO_State.Increment_E_plastic[2] = d_gamma_k * (alpha_Q + n[2]);

  IO_State.Stress[0] = -T_trial_vol[0] + T_trial_dev[0] -
                       d_gamma_k * (3 * K * alpha_Q + 2 * G * n[0]);
  IO_State.Stress[1] = -T_trial_vol[1] + T_trial_dev[1] -
                       d_gamma_k * (3 * K * alpha_Q + 2 * G * n[1]);
  IO_State.Stress[2] = -T_trial_vol[2] + T_trial_dev[2] -
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
                                            const double *T_trial_vol,
                                            const double *n, double d_gamma_k,
                                            double d_gamma_1, double alpha_Q,
                                            double K, double eps_k,
                                            double kappa_k) {

  IO_State.Increment_E_plastic[0] = d_gamma_k * alpha_Q + d_gamma_1 * n[0];
  IO_State.Increment_E_plastic[1] = d_gamma_k * alpha_Q + d_gamma_1 * n[1];
  IO_State.Increment_E_plastic[2] = d_gamma_k * alpha_Q + d_gamma_1 * n[2];

  IO_State.Stress[0] = -T_trial_vol[0] - d_gamma_k * 3 * K * alpha_Q;
  IO_State.Stress[1] = -T_trial_vol[1] - d_gamma_k * 3 * K * alpha_Q;
  IO_State.Stress[2] = -T_trial_vol[2] - d_gamma_k * 3 * K * alpha_Q;

  *IO_State.Equiv_Plast_Str = eps_k;
  *IO_State.Kappa = kappa_k;

  return EXIT_SUCCESS;
}

/***************************************************************************/
