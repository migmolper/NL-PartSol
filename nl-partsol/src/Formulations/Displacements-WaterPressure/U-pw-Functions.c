/**
 * @file U-pw-Functions.c
 * @author Miguel Molinos (@migmolper)
 * @brief Functions to compute the residual and the jacobian
 * for the equations of balance of momentum of the mixture
 * and the darcy equation
 * @version 0.1
 * @date 2022-05-15
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Formulations/Displacements-WaterPressure/U-pw-Functions.h"

/**************************************************************/

void compute_L0__upw__(double *L0, const double *dN_alpha_u_n1,
                       const double *kichhoff_stress, double V0) {

#if NumberDimensions == 2
  L0[0] = (kichhoff_stress[0] * dN_alpha_u_n1[0] +
           kichhoff_stress[1] * dN_alpha_u_n1[1]) *
          V0;
  L0[1] = (kichhoff_stress[2] * dN_alpha_u_n1[0] +
           kichhoff_stress[3] * dN_alpha_u_n1[1]) *
          V0;
#else
  L0[0] = (kichhoff_stress[0] * dN_alpha_u_n1[0] +
           kichhoff_stress[1] * dN_alpha_u_n1[1] +
           kichhoff_stress[2] * dN_alpha_u_n1[2]) *
          V0;
  L0[1] = (kichhoff_stress[3] * dN_alpha_u_n1[0] +
           kichhoff_stress[4] * dN_alpha_u_n1[1] +
           kichhoff_stress[5] * dN_alpha_u_n1[2]) *
          V0;
  L0[2] = (kichhoff_stress[6] * dN_alpha_u_n1[0] +
           kichhoff_stress[7] * dN_alpha_u_n1[1] +
           kichhoff_stress[8] * dN_alpha_u_n1[2]) *
          V0;
#endif
}

/**************************************************************/

void compute_L1__upw__(double *L1, const double *dN_alpha_u_n1,
                       double kichhoff_pressure, double V0) {

#if NumberDimensions == 2
  L1[0] = -kichhoff_pressure * dN_alpha_u_n1[0] * V0;
  L1[1] = -kichhoff_pressure * dN_alpha_u_n1[1] * V0;
#else
  L1[0] = -kichhoff_pressure * dN_alpha_u_n1[0] * V0;
  L1[1] = -kichhoff_pressure * dN_alpha_u_n1[1] * V0;
  L1[2] = -kichhoff_pressure * dN_alpha_u_n1[2] * V0;
#endif
}

/**************************************************************/

void compute_DL1_Du__upw__(double *DL1_Du, const double *dN_alpha_u_n1,
                           const double *dN_beta_u_n1, double kichhoff_pressure,
                           double V0) {

#if NumberDimensions == 2
  DL1_Du[0] = kichhoff_pressure * dN_beta_u_n1[0] * dN_alpha_u_n1[0] * V0;
  DL1_Du[1] = kichhoff_pressure * dN_beta_u_n1[0] * dN_alpha_u_n1[1] * V0;

  DL1_Du[2] = kichhoff_pressure * dN_beta_u_n1[1] * dN_alpha_u_n1[0] * V0;
  DL1_Du[3] = kichhoff_pressure * dN_beta_u_n1[1] * dN_alpha_u_n1[1] * V0;
#else
  DL1_Du[0] = kichhoff_pressure * dN_beta_u_n1[0] * dN_alpha_u_n1[0] * V0;
  DL1_Du[1] = kichhoff_pressure * dN_beta_u_n1[0] * dN_alpha_u_n1[1] * V0;
  DL1_Du[2] = kichhoff_pressure * dN_beta_u_n1[0] * dN_alpha_u_n1[2] * V0;

  DL1_Du[3] = kichhoff_pressure * dN_beta_u_n1[1] * dN_alpha_u_n1[0] * V0;
  DL1_Du[4] = kichhoff_pressure * dN_beta_u_n1[1] * dN_alpha_u_n1[1] * V0;
  DL1_Du[5] = kichhoff_pressure * dN_beta_u_n1[1] * dN_alpha_u_n1[2] * V0;

  DL1_Du[6] = kichhoff_pressure * dN_beta_u_n1[2] * dN_alpha_u_n1[0] * V0;
  DL1_Du[7] = kichhoff_pressure * dN_beta_u_n1[2] * dN_alpha_u_n1[1] * V0;
  DL1_Du[8] = kichhoff_pressure * dN_beta_u_n1[2] * dN_alpha_u_n1[2] * V0;
#endif
}

/**************************************************************/

void compute_DL1_Dpw__upw__(double *DL1_Dpw, const double *dN_alpha_u_n1,
                            double N_beta_pw_n1, double V0) {

#if NumberDimensions == 2
  DL1_Dpw[0] = -dN_alpha_u_n1[0] * N_beta_pw_n1 * V0;
  DL1_Dpw[1] = -dN_alpha_u_n1[1] * N_beta_pw_n1 * V0;
#else
  DL1_Dpw[0] = -dN_alpha_u_n1[0] * N_beta_pw_n1 * V0;
  DL1_Dpw[1] = -dN_alpha_u_n1[1] * N_beta_pw_n1 * V0;
  DL1_Dpw[2] = -dN_alpha_u_n1[2] * N_beta_pw_n1 * V0;
#endif
}

/**************************************************************/

void compute_L2__upw__(double *L2, double N_alpha_u_n1, const double *b,
                       double m) {
#if NumberDimensions == 2
  L2[0] = -N_alpha_u_n1 * m * b[0];
  L2[1] = -N_alpha_u_n1 * m * b[1];
#else
  L2[0] = -N_alpha_u_n1 * m * b[0];
  L2[1] = -N_alpha_u_n1 * m * b[1];
  L2[2] = -N_alpha_u_n1 * m * b[2];
#endif
}

/**************************************************************/

void compute_DL2_Du__upw__(double *DL2_Du, double N_alpha_u_n1,
                           double N_beta_u_n1, const double *dN_beta_u_n1,
                           const double *b, double kichhoff_pressure,
                           double phi_f, double intrinsic_rho_f, double kappa_f,
                           double Jacobian, double m, double alpha_1,
                           double V0) {

  double c0 =
      intrinsic_rho_f * (Jacobian - kichhoff_pressure * phi_f / kappa_f);

#if NumberDimensions == 2
  DL2_Du[0] = -c0 * N_alpha_u_n1 * b[0] * dN_beta_u_n1[0] * V0 -
              alpha_1 * N_alpha_u_n1 * N_beta_u_n1 * m;
  DL2_Du[1] = -c0 * N_alpha_u_n1 * b[0] * dN_beta_u_n1[1] * V0;

  DL2_Du[2] = -c0 * N_alpha_u_n1 * b[1] * dN_beta_u_n1[0] * V0;
  DL2_Du[3] = -c0 * N_alpha_u_n1 * b[1] * dN_beta_u_n1[1] * V0 -
              alpha_1 * N_alpha_u_n1 * N_beta_u_n1 * m;
#else
  DL2_Du[0] = -c0 * N_alpha_u_n1 * b[0] * dN_beta_u_n1[0] * V0 -
              alpha_1 * N_alpha_u_n1 * N_beta_u_n1 * m;
  DL2_Du[1] = -c0 * N_alpha_u_n1 * b[0] * dN_beta_u_n1[1] * V0;
  DL2_Du[2] = -c0 * N_alpha_u_n1 * b[0] * dN_beta_u_n1[2] * V0;

  DL2_Du[3] = -c0 * N_alpha_u_n1 * b[1] * dN_beta_u_n1[0] * V0;
  DL2_Du[4] = -c0 * N_alpha_u_n1 * b[1] * dN_beta_u_n1[1] * V0 -
              alpha_1 * N_alpha_u_n1 * N_beta_u_n1 * m;
  DL2_Du[5] = -c0 * N_alpha_u_n1 * b[1] * dN_beta_u_n1[2] * V0;

  DL2_Du[6] = -c0 * N_alpha_u_n1 * b[2] * dN_beta_u_n1[0] * V0;
  DL2_Du[7] = -c0 * N_alpha_u_n1 * b[2] * dN_beta_u_n1[1] * V0;
  DL2_Du[8] = -c0 * N_alpha_u_n1 * b[2] * dN_beta_u_n1[2] * V0 -
              alpha_1 * N_alpha_u_n1 * N_beta_u_n1 * m;
#endif
}

/**************************************************************/

void compute_DL2_Dpw__upw__(double *DL2_Dpw, double N_alpha_u_n1,
                            double N_beta_pw_n1, const double *b, double phi_f,
                            double intrinsic_rho_f, double kappa_f, double V0) {

  double c0 = phi_f * intrinsic_rho_f / kappa_f;

#if NumberDimensions == 2
  DL2_Dpw[0] = -c0 * b[0] * N_alpha_u_n1 * N_beta_pw_n1 * V0;
  DL2_Dpw[1] = -c0 * b[1] * N_alpha_u_n1 * N_beta_pw_n1 * V0;
#else
  DL2_Dpw[0] = -c0 * b[0] * N_alpha_u_n1 * N_beta_pw_n1 * V0;
  DL2_Dpw[1] = -c0 * b[1] * N_alpha_u_n1 * N_beta_pw_n1 * V0;
  DL2_Dpw[2] = -c0 * b[2] * N_alpha_u_n1 * N_beta_pw_n1 * V0;
#endif
}

/**************************************************************/

void compute_G0__upw__(double *G0, double N_alpha_pw_n1, double intrinsic_rho_f,
                       double relative_rho_f_p, double rate_kichhoff_pressure,
                       double rate_Jacobian, double kappa_f, double V0) {
  *G0 = N_alpha_pw_n1 *
        (relative_rho_f_p * rate_kichhoff_pressure / kappa_f +
         intrinsic_rho_f * rate_Jacobian) *
        V0;
}

/**************************************************************/

void compute_DG0_Du__upw__(double *DG0_Du, double N_alpha_pw_n1,
                           const double *dN_beta_u_n1, const double *grad_v,
                           double kichhoff_pressure,
                           double rate_kichhoff_pressure, double Jacobian,
                           double rate_Jacobian, double phi_s, double phi_f,
                           double intrinsic_rho_f, double kappa_f,
                           double alpha_4, double V0) {

  double c0 = (rate_kichhoff_pressure * intrinsic_rho_f / kappa_f) *
              (phi_s - kichhoff_pressure * phi_f / (Jacobian * kappa_f));
  double c1 = rate_Jacobian * intrinsic_rho_f * kichhoff_pressure /
              (Jacobian * kappa_f);
  double c2 = intrinsic_rho_f * (rate_Jacobian + alpha_4 * Jacobian);
  double c3 = intrinsic_rho_f * Jacobian;

#if NumberDimensions == 2
  DG0_Du[0] =
      N_alpha_pw_n1 *
      ((c0 - c1 + c2) * dN_beta_u_n1[0] -
       c3 * (grad_v[0] * dN_beta_u_n1[0] + grad_v[1] * dN_beta_u_n1[1])) *
      V0;

  DG0_Du[1] =
      N_alpha_pw_n1 *
      ((c0 - c1 + c2) * dN_beta_u_n1[1] -
       c3 * (grad_v[2] * dN_beta_u_n1[0] + grad_v[3] * dN_beta_u_n1[1])) *
      V0;

#else
  DG0_Du[0] = N_alpha_pw_n1 *
              ((c0 - c1 + c2) * dN_beta_u_n1[0] -
               c3 * (grad_v[0] * dN_beta_u_n1[0] + grad_v[1] * dN_beta_u_n1[1] +
                     grad_v[2] * dN_beta_u_n1[2])) *
              V0;
  DG0_Du[1] = N_alpha_pw_n1 *
              ((c0 - c1 + c2) * dN_beta_u_n1[1] -
               c3 * (grad_v[3] * dN_beta_u_n1[0] + grad_v[4] * dN_beta_u_n1[1] +
                     grad_v[5] * dN_beta_u_n1[2])) *
              V0;
  DG0_Du[2] = N_alpha_pw_n1 *
              ((c0 - c1 + c2) * dN_beta_u_n1[2] -
               c3 * (grad_v[6] * dN_beta_u_n1[0] + grad_v[7] * dN_beta_u_n1[1] +
                     grad_v[8] * dN_beta_u_n1[2])) *
              V0;
#endif
}

/**************************************************************/

void compute_DG0_Dpw__upw__(double *DG0_Dpw, double N_alpha_pw_n1,
                            double N_beta_pw_n1, double rate_kichhoff_pressure,
                            double Jacobian, double rate_Jacobian, double phi_f,
                            double intrinsic_rho_f, double relative_rho_f,
                            double kappa_f, double alpha_4, double V0) {

  double c0 = (rate_kichhoff_pressure * intrinsic_rho_f * phi_f) /
              (kappa_f * kappa_f * Jacobian);

  double c1 = relative_rho_f * alpha_4 / kappa_f;

  double c2 = rate_Jacobian * intrinsic_rho_f / (Jacobian * kappa_f);

  *DG0_Dpw = N_alpha_pw_n1 * (c0 + c1 + c2) * N_beta_pw_n1 * V0;
}

/**************************************************************/

void compute_G1__upw__(double *G1, const double *dN_alpha_pw_n1,
                       const double *K, const double *grad_kichhoff_pressure,
                       double g, double V0) {
#if NumberDimensions == 2
  *G1 = -(1.0 / g) *
        (dN_alpha_pw_n1[0] * (K[0] * grad_kichhoff_pressure[0] +
                              K[1] * grad_kichhoff_pressure[1]) +
         dN_alpha_pw_n1[1] * (K[2] * grad_kichhoff_pressure[0] +
                              K[3] * grad_kichhoff_pressure[1])) *
        V0;
#else
  *G1 = -(1.0 / g) *
        (dN_alpha_pw_n1[0] * (K[0] * grad_kichhoff_pressure[0] +
                              K[1] * grad_kichhoff_pressure[1] +
                              K[2] * grad_kichhoff_pressure[2]) +
         dN_alpha_pw_n1[1] * (K[3] * grad_kichhoff_pressure[0] +
                              K[4] * grad_kichhoff_pressure[1] +
                              K[5] * grad_kichhoff_pressure[2]) +
         dN_alpha_pw_n1[2] * (K[6] * grad_kichhoff_pressure[0] +
                              K[7] * grad_kichhoff_pressure[1] +
                              K[8] * grad_kichhoff_pressure[2])) *
        V0;
#endif
}

/**************************************************************/

void compute_DG1_Du__upw__(double *DG1_Du, const double *dN_alpha_pw_n1,
                           const double *dN_beta_u_n1,
                           const double *grad_kichhoff_pressure,
                           const double *K, double g, double V0) {
#if NumberDimensions == 2
  DG1_Du[0] = (1.0 / g) * ();

  dN_alpha_pw_n1[0] * dN_beta_u_n1[0] * +dN_alpha_pw_n1[0] * dN_beta_u_n1[1] *

      DG1_Du[1] = (1.0 / g) * ();
#else
  DG1_Du[0] = (1.0 / g) *;
  DG1_Du[1] = (1.0 / g) *;
  DG1_Du[2] = (1.0 / g) *;
#endif
}

/**************************************************************/

void compute_G2__upw__(double *G2, const double *dN_alpha_pw_n1,
                       const double *K, const double *b, double Jacobian,
                       double intrinsic_rho_f, double g, double V0) {
#if NumberDimensions == 2
  *G2 = -(Jacobian * intrinsic_rho_f / g) *
        (dN_alpha_pw_n1[0] * (K[0] * b[0] + K[1] * b[1]) +
         dN_alpha_pw_n1[1] * (K[2] * b[0] + K[3] * b[1])) *
        V0;
#else
  *G2 = -(Jacobian * intrinsic_rho_f / g) *
        (dN_alpha_pw_n1[0] * (K[0] * b[0] + K[1] * b[1] + K[2] * b[2]) +
         dN_alpha_pw_n1[1] * (K[3] * b[0] + K[4] * b[1] + K[5] * b[2]) +
         dN_alpha_pw_n1[2] * (K[6] * b[0] + K[7] * b[1] + K[8] * b[2])) *
        V0;
#endif
}

/**************************************************************/