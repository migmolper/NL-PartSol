#include <math.h>
#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

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

static int __eps(
    double *eps_k /**< [out] Equivalent plastic strain*/,
    double d_gamma_k /**< [in] Discrete plastic multiplier */,
    double eps_n /**< [in] Equivalent plastic strain in the last step */,
    double alpha_Q /**< [in] Plastic potential parameter */);

static int __kappa(double *kappa_k /**< [out] Hardening function. */,
                           double kappa_0 /**< [in] Reference hardening */,
                           double exp_param /**< [in] Hardening exponential*/,
                           double eps_k /**< [in] Equivalent plastic strain*/,
                           double eps_0 /**< [in] Reference plastic strain */);

static int __d_kappa(
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

static int __tangent_moduli_classical(
  double *  C_ep /**< [out] Elastoplastic tanget moduli */, 
  const double * n /**< [in] Plastic flow direction. */,
  double d_gamma_k /**< [in] Derivative of the hardening function. */, 
  double J2 /**< [in] Second invariant of the deviatoric stress tensor */, 
  double d_kappa_k /**< [in] Discrete plastic multiplier */,
  double K /**< [in] First Lamé invariant. */, 
  double G /**< [in] Second Lamé invariant. */, 
  double beta/**< [in] Yield surface parameter II. */, 
  double alpha_F /**< [in] Yield surface parameter I. */, 
  double alpha_Q /**< [in] Plastic potential parameter. */);

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
    double d_kappa_k /**< [in] Derivative of the hardening function */,
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

static int __tangent_moduli_apex(
  double *C_ep /**< [out] Elastoplastic tanget moduli */, 
  const double *n /**< [in] Plastic flow direction. */, 
  double d_gamma_k /**< [in] Discrete plastic multiplier */, 
  double d_gamma_1 /**< [in] Discrete plastic multiplier I */, 
  double d_kappa_k /**< [in] Derivative of the hardening function */,
  double K /**< [in] First Lamé invariant. */, 
  double G /**< [in] Second Lamé invariant. */, 
  double beta /**< [in] Yield surface parameter II. */, 
  double alpha_F /**< [in] Yield surface parameter I. */, 
  double alpha_Q /**< [in] Plastic potential parameter. */);    

/**************************************************************/ 

int compute_1PK_Drucker_Prager(State_Parameters IO_State, Material MatProp)
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
        IO_State.Stress, IO_State.D_phi_n1, T_tr_vol, T_tr_dev, eigvec_b_e_tr);
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
          fprintf(stderr, "" RED "|d_PHI| = %f < TOL (classical loop)" RESET "\n",fabs(d_PHI));
          return EXIT_FAILURE;
        }

        d_gamma_k += -PHI / d_PHI;
        if (d_gamma_k < 0.0) {
          fprintf(stderr, "" RED "d_gamma_k = %f < 0 (classical loop)" RESET "\n",d_gamma_k);
          return EXIT_FAILURE;
        }

        STATUS = __eps(&eps_k, d_gamma_k, eps_n, alpha_Q);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr,
                  "" RED "Error in __eps (classical loop)" RESET "\n");
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

      STATUS = __update_internal_variables_classical(
          Increment_E_plastic, IO_State.Stress, IO_State.EPS,
          IO_State.Kappa, IO_State.D_phi_n1, T_tr_vol, T_tr_dev, eigvec_b_e_tr, n,
          d_gamma_k, alpha_Q, K, G, eps_k, kappa_k);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __update_internal_variables_classical" RESET
                "\n");
        return EXIT_FAILURE;
      }

      if(IO_State.compute_C_ep)
      {
        STATUS = __tangent_moduli_classical(
          IO_State.C_ep,n,d_gamma_k,J2,d_kappa_k,
          K,G,beta,alpha_F,alpha_Q);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr, "" RED "Error in __tangent_moduli_classical" RESET "\n");
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
          fprintf(stderr, "" RED "|d_PHI| = %e < TOL (apex loop)" RESET "\n",fabs(d_PHI));
          break;
        }

        d_gamma_2_k += -PHI / d_PHI;
        
        if (d_gamma_2_k < 0.0) {
          fprintf(stderr, "" RED "Breackage (apex loop)" RESET "\n");
          d_gamma_k = 0.0;
          d_gamma_1 = 0.0;
        }
        else
        {
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

      STATUS = __update_internal_variables_apex(
          Increment_E_plastic, IO_State.Stress, IO_State.EPS,
          IO_State.Kappa, IO_State.D_phi_n1, T_tr_vol, T_tr_dev, eigvec_b_e_tr, n,
          d_gamma_k, d_gamma_1, alpha_Q, K, G, eps_k, kappa_k);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in__update_internal_variables_apex" RESET "\n");
        return EXIT_FAILURE;
      }

      if(IO_State.compute_C_ep)
      {
        STATUS = __tangent_moduli_apex(
          IO_State.C_ep, n, d_gamma_k, d_gamma_1, d_kappa_k,
          K, G, beta, alpha_F, alpha_Q);
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

    T_A = -T_tr_vol[A] + T_tr_dev[A];

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
  Stress[4] = -T_tr_vol[2] + T_tr_dev[2];
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

static int __update_internal_variables_classical(
    double *Increment_E_plastic, double *Stress, double *eps_n1,
    double *kappa_n1, const double *D_phi, const double *T_tr_vol,
    const double *T_tr_dev, const double *eigvec_b_e_tr, const double *n,
    double d_gamma_k, double alpha_Q, double K, double G, double eps_k,
    double kappa_k) {

  int Ndim = NumberDimensions;

  // Update hardening parameters
  *eps_n1 = eps_k;
  *kappa_n1 = kappa_k;

  // Update the increment of the plastic strain tensor
  Increment_E_plastic[0] = d_gamma_k * (alpha_Q + n[0]);
  Increment_E_plastic[1] = d_gamma_k * (alpha_Q + n[1]);
  Increment_E_plastic[2] = d_gamma_k * (alpha_Q + n[2]);

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

    T_A = -T_tr_vol[A] + T_tr_dev[A] +
          d_gamma_k * (3 * K * alpha_Q - 2 * G * n[A]);

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
  Stress[4] =
      -T_tr_vol[2] + T_tr_dev[2] + d_gamma_k * (3 * K * alpha_Q - 2 * G * n[2]);
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

  int Ndim = NumberDimensions;

  // Update hardening parameters
  *eps_n1 = eps_k;
  *kappa_n1 = kappa_k;

  // Update the increment of the plastic strain tensor
  Increment_E_plastic[0] = d_gamma_k * alpha_Q + d_gamma_1 * n[0];
  Increment_E_plastic[1] = d_gamma_k * alpha_Q + d_gamma_1 * n[1];
  Increment_E_plastic[2] = d_gamma_k * alpha_Q + d_gamma_1 * n[2];

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

    T_A = -T_tr_vol[A] + d_gamma_k * 3 * K * alpha_Q;

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
  Stress[4] = -T_tr_vol[2] + d_gamma_k * 3 * K * alpha_Q;
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

static int __tangent_moduli_classical(
          double *  C_ep, 
          const double * n, 
          double d_gamma_k, 
          double J2, 
          double d_kappa_k,
          double K, 
          double G, 
          double beta, 
          double alpha_F, 
          double alpha_Q)
{

  int STATUS = EXIT_SUCCESS;  
  int Ndim = NumberDimensions;

#if NumberDimensions == 2
  double R2_Identity[2] = {1.0,1.0};

  double R4_Identity[2][2] = {
      {1.0,0.0},
      {0.0,1.0}
  };
#else
  double R2_Identity[3] = {1.0,1.0,1.0};

  double R4_Identity[3][3] = {
    {1.0,0.0,0.0},
    {0.0,1.0,0.0},
    {0.0,0.0,1.0}
  };
#endif


  double c0 = 9*alpha_F*alpha_Q*K + 2*G + beta*d_kappa_k*sqrt(2./3. * (1 + 3*alpha_Q*alpha_Q));
  double c1 = 1.0 - 9.0*alpha_F*alpha_Q*K/c0; 

  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
      C_ep[i*Ndim + j] = 
      c1*K*R2_Identity[i]*R2_Identity[j] 
      + 2*G*(R4_Identity[i][j] - (1./3.)*(1.0 - 2.0*G*d_gamma_k/J2)*R2_Identity[i]*R2_Identity[j]) 
      - (6.0*alpha_Q*K*G/c0)*R2_Identity[i]*n[j] 
      - (6.0*alpha_Q*K*G/c0)*n[i]*R2_Identity[j] 
      - 4*G*G*(1.0/c0 - d_gamma_k/J2)*n[i]*n[j];
    }
  }
  
  return STATUS;
}

/***************************************************************************/

static int __tangent_moduli_apex(
          double *  C_ep, 
          const double * n, 
          double d_gamma_k, 
          double d_gamma_1, 
          double d_kappa_k,
          double K, 
          double G, 
          double beta, 
          double alpha_F, 
          double alpha_Q)
{

  int STATUS = EXIT_SUCCESS;  
  int Ndim = NumberDimensions;

#if NumberDimensions == 2
    double R2_Identity[2] = {1.0,1.0};
#else
    double R2_Identity[3] = {1.0,1.0,1.0};
#endif

  double c0 = (alpha_Q*beta*sqrt(2./3.)*d_kappa_k*d_gamma_k)/
  (3.0*alpha_F*K*sqrt(d_gamma_1*d_gamma_1 + 3.0*alpha_Q*alpha_Q*d_gamma_k*d_gamma_k) + 
  alpha_Q*beta*sqrt(2./3.)*d_kappa_k*d_gamma_k);

  double c1 = c0*K/(2.0*alpha_Q*G*d_gamma_k);

  for (unsigned i = 0; i < Ndim; i++)
  {
    for (unsigned j = 0; j < Ndim; j++)
    {
      C_ep[i*Ndim + j] =c0*K*R2_Identity[i]*R2_Identity[j] + c1*R2_Identity[i]*n[j];
    }
  }
  
  return STATUS;
}

/**************************************************************/