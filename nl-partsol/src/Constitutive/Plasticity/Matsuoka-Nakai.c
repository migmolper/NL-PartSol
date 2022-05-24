
#include "Constitutive/Plasticity/Matsuoka-Nakai.h"

/**************************************************************/
static int __compute_trial_b_e(
    double *eigval_b_e_tr /**< [out] Eigenvalues of b elastic trial. */,
    double *eigvec_b_e_tr /**< [out] Eigenvector of b elastic trial. */,
    const double *b_e /**< [in] (n) Elastic left Cauchy-Green.*/,
    const double *d_phi /**< [in] Incremental deformation gradient. */);
/**************************************************************/

static void __corrector_b_e(
    double *b_e /**< [out] (n+1) Elastic deformation gradient. */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */,
    const double *E_hencky_trial /**< [in] Corrected Henky strain */);
/**************************************************************/

static void __elastic_tangent(
    double * CC /**< [out] Elastic compliance */, 
    double * AA /**< [out] Elastic matrix */,
    double E /**< [in] Young modulus */, 
    double nu /**< [in] Poisson ratio */,
    double Lame /**< [in] Lamé parameter */, 
    double G /**< [in] Shear modulus */);
/**************************************************************/

static void __trial_elastic(
    double *T_tr /**< [out] Trial elastic stress tensor*/, 
    const double * E_hencky_trial /**< [in] Henky strain (trial) */, 
    const double * AA /**< [in] Elastic matrix */,
    double c_cotphi /**< [in] Cohesion parameter */); 
/**************************************************************/

static void __E_hencky(
    double * E_hencky_k /**< [out] Henky strain (iter k) */, 
    const double * T_k  /**< [in] Local stress tensor (iter k) */, 
    const double * CC /**< [in] Elastic compliance */,
    double c_cotphi /**< [in] Cohesion parameter */);
/**************************************************************/

static int __update_internal_variables_elastic(
    double *Stress /**< [in/out] Nominal stress tensor */,
    const double *D_phi /**< [in] Total deformation gradient. */,
    const double *T_tr /**< [in] Elastic stress tensor */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */,
    double c_cotphi /**< [in] Cohesion parameter */);
/**************************************************************/

static void __elastic_tangent_moduli(
  double * C_ep /**< [in/out] Elastoplastic tangent matrix */, 
  const double * AA /**< [in] Elastic matrix */);
/**************************************************************/

static void __kappa(
    double *kappa /**< [out] Hardening vector */, 
    const double * a /**< [in] Vector with fit parameters (hardening) */, 
    double Lambda /**< [in] Total plastic multiplier */, 
    double I1 /**< [in] First invariant of the stress tensor */, 
    double alpha /**< [in] Dilatance parameter*/);
/**************************************************************/

static void __d_kappa_phi_d_stress(
    double *d_kappa_phi_d_stress /**< [out] Stress derivative of kappa[0] */, 
    const double * a /**< [in] Vector with fit parameters (hardening) */,
    double Lambda /**< [in] Total plastic multiplier */, 
    double I1 /**< [in] First invariant of the stress tensor */);
/**************************************************************/

static void __d_kappa_phi_d_lambda(
    double * d_kappa_phi_d_lambda /**< [out] Lambda derivative of kappa[0] */, 
    const double * a /**< [in] Vector with fit parameters (hardening) */,
    double Lambda /**< [in] Total plastic multiplier */, 
    double I1 /**< [in] First invariant of the stress tensor */);
/**************************************************************/

static void __F(
    double * F /**< [out] Yield function evaluation */,
    double kappa_phi /**< [in] Friction angle hardening */, 
    double I1 /**< [in] First invariant of the stress tensor */, 
    double I2 /**< [in] Second invariant of the stress tensor */, 
    double I3 /**< [in] Third invariant of the stress tensor */);
/**************************************************************/

static void __d_F_d_stress(
    double *d_F_d_stress /**< [out] Yield function derivative (stress) */, 
    const double *T_k /**< [in] Local stress tensor (iter k) */,
    double I1 /**< [in] First invariant of the stress tensor */, 
    double I2 /**< [in] Second invariant of the stress tensor */, 
    double I3 /**< [in] Third invariant of the stress tensor */,
    double kappa_phi /**< [in] Friction angle hardening */);
/**************************************************************/

static void __d_F_d_kappa_phi(
    double * d_F_d_kappa_phi /**< [out] Yield function derivative (kappa[0]) */,
    double I1 /**< [in] First invariant of the stress tensor */,
    double I3 /**< [in] Third invariant of the stress tensor */, 
    double kappa_phi /**< [in] Friction angle hardening */);
/**************************************************************/

static void __d_G_d_stress(
    double *d_G_d_stress /**< [out] Plastic potential function derivative (stress) */,
    const double *T_k /**< [in] Local stress tensor (iter k) */, 
    double I1 /**< [in] First invariant of the stress tensor */,
    double I2 /**< [in] Second invariant of the stress tensor */, 
    double I3 /**< [in] Third invariant of the stress tensor */, 
    double kappa_psi /**< [in] Dilatance angle hardening */);
/**************************************************************/

static void __dd_G_dd_stress(
    double *dd_G_dd_stress /**< [out] Plastic potential hessian (stress) */, 
    const double *T_k /**< [in] Local stress tensor (iter k) */,
    double kappa_psi /**< [in] Dilatance angle hardening */, 
    double I1 /**< [in] First invariant of the stress tensor */,
    double I2 /**< [in] Second invariant of the stress tensor */, 
    double I3 /**< [in] Third invariant of the stress tensor */);
/**************************************************************/

static void __dd_G_d_stress_d_kappa_psi(
    double *dd_G_d_stress_d_kappa_psi /**< [out] Plastic potential deriv */, 
    const double *T_k /**< [in] Local stress tensor (iter k) */, 
    double I1 /**< [in] First invariant of the stress tensor */, 
    double I3 /**< [in] Third invariant of the stress tensor */,
    double kappa_psi /**< [in] Dilatance angle hardening */); 
/**************************************************************/

static int __residual(
    double *Residual /**< [out] Residual of the problem */, 
    double * Error_k /**< [out] Norm of the residual */,
    const double *E_hencky_trial /**< [in] Henky strain (trial) */,
    const double * E_hencky_k /**< [in] Henky strain (iter k) */,
    const double *d_G_d_stress /**< [in] Plastic potential function derivative (stress) */,
    const double *kappa_k /**< [in] Hardening vector (iter k) */,
    const double *kappa_hat /**< [in] Hardening vector (eval) */, 
    double F_k /**< [in] Yield function evaluation (iter k) */,
    double delta_lambda_k /**< [in] Discrete plastic multiplier (iter k) */); 
/**************************************************************/

static int __update_internal_variables_plastic(
    double *Stress /**< [out] Nominal stress tensor */, 
    double *eps_n1 /**< [out] Equivalent plastic strain */,
    double *kappa_n1 /**< [out] Friction angle hardening */, 
    const double *D_phi /**< [in] Total deformation gradient. */, 
    const double *T_tr_k /**< [in] Stress tensor (iter k). */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */, 
    double Lambda_k /**< [in] Total plastic multiplier (iter k) */,
    double kappa_phi_k /**< [out] Friction angle hardening (iter k)*/,
    double c_cotphi /**< [in] Cohesion parameter */);
/**************************************************************/

static int __elastoplastic_tangent_moduli(
  double * C_ep /**< [in/out] Elastoplastic tangent matrix */, 
  const double * CC /**< [in] Elastic compliance */,
  const double * dd_G_dd_stress /**< [in] Plastic potential hessian (stress) */,
  double delta_lambda_k /**< [in] Discrete plastic multiplier (iter k) */);
/**************************************************************/

static int __solver(
  double *Tangent_Matrix /**< [in/out] Tangent matrix of the problem */,
  double *Residual /**< [in/out] Residual of the problem */) ;
/**************************************************************/

static int __reciprocal_condition_number(
    double *RCOND /**< [out] Condition number of the tangent matrix */,
    double *Tangent_Matrix /**< [in/out] Tangent matrix of the problem */);
/**************************************************************/


int compute_1PK_Matsuoka_Nakai(State_Parameters IO_State, Material MatProp)
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
  double E = MatProp.E;
  double nu = MatProp.nu;
  double Lame = E*nu / ((1.0 + nu)*(1.0 - 2.0 * nu));
  double G = E / (2.0 * (1.0 + nu));
  double friction_angle = MatProp.phi_Frictional;
  double rad_friction_angle = (PI__MatrixLib__ / 180.0) * friction_angle;  
  double c_cotphi = rad_friction_angle > 0.0? MatProp.Cohesion/tan(rad_friction_angle) : 0.0;
  double alpha = MatProp.alpha_Hardening_Borja;
  double a[3] = {0.0, 0.0, 0.0};
  a[0] = MatProp.a_Hardening_Borja[0];
  a[1] = MatProp.a_Hardening_Borja[1];
  a[2] = MatProp.a_Hardening_Borja[2];
  
  bool Activate_CutOff = false;
  double CutOff =  0.0;
  
  double CC[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double AA[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  __elastic_tangent(CC, AA ,E,nu,Lame,G);

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

  // Initialize Newton-Raphson solver
  double TOL = 1E-10;// TOL_Radial_Returning;
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

  __trial_elastic(T_tr, E_hencky_trial, AA, c_cotphi);

  I1 = T_tr[0] + T_tr[1] + T_tr[2];
  I2 = T_tr[0] * T_tr[1] + T_tr[1] * T_tr[2] + T_tr[0] * T_tr[2];
  I3 = T_tr[0] * T_tr[1] * T_tr[2];

  // Check yield
  __F(&F_0, kappa_n[0], I1, I2, I3);

  // Assign trial (elastic) stress-strain values
  T_k1[0] = T_tr[0];
  T_k1[1] = T_tr[1];
  T_k1[2] = T_tr[2];  

  // Elastic
  if (F_0 <= 0.0) {

    STATUS = __update_internal_variables_elastic(
        IO_State.Stress, IO_State.D_phi_n1, T_k1, eigvec_b_e_tr, c_cotphi);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in __update_internal_variables_elastic()" RESET
              "\n");
      return EXIT_FAILURE;
    }

    if(IO_State.compute_C_ep)
    {        
      __elastic_tangent_moduli(IO_State.C_ep, AA);
    }    

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

      // Introduce a preconditioner
      Tangent_Matrix[0] += Residual_k1[0];
      Tangent_Matrix[6] += Residual_k1[1];
      Tangent_Matrix[12] += Residual_k1[2];
      Tangent_Matrix[18] += Residual_k1[3];
      Tangent_Matrix[24] += Residual_k1[4];

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

        Lambda_k1 = Lambda_n;
        kappa_k1[0] = kappa_n[0];

        break;
      }
    }

    STATUS = __update_internal_variables_plastic(IO_State.Stress,IO_State.EPS,
     IO_State.Kappa,IO_State.D_phi_n1,T_k1,eigvec_b_e_tr,Lambda_k1, kappa_k1[0], c_cotphi);    
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __update_internal_variables_plastic" RESET "\n");
      return EXIT_FAILURE;
    }

    if(IO_State.compute_C_ep)
    {        
      STATUS = __elastoplastic_tangent_moduli(IO_State.C_ep, CC, dd_G_dd_stress, delta_lambda_k1);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __elastoplastic_tangent_moduli" RESET "\n");
        return EXIT_FAILURE;
      }
    }

  }
  
  __corrector_b_e(IO_State.b_e, eigvec_b_e_tr, E_hencky_k1);

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


  lapack_int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec_b_e_tr, lda, eigval_b_e_tr);

  if (info > 0) {
    fprintf(stderr,
            "" RED "Error in dsyev_(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    fprintf(stderr,
            "" RED "Error in dsyev_(): the %i-th argument had an "
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

static void __corrector_b_e(double *b_e, const double *eigvec_b_e_tr,
                           const double *E_hencky_trial) {
  unsigned Ndim = NumberDimensions;

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

  T_tr[0] = AA[0] * E_hencky_trial[0] 
          + AA[1] * E_hencky_trial[1] 
          + AA[2] * E_hencky_trial[2] 
          - c_cotphi;
  T_tr[1] = AA[3] * E_hencky_trial[0] 
          + AA[4] * E_hencky_trial[1] 
          + AA[5] * E_hencky_trial[2] 
          - c_cotphi;
  T_tr[2] = AA[6] * E_hencky_trial[0] 
          + AA[7] * E_hencky_trial[1] 
          + AA[8] * E_hencky_trial[2] 
          - c_cotphi;
}

/**************************************************************/

static void __E_hencky(double *E_hencky_k, 
                      const double *T_k, 
                      const double *CC, 
                      double c_cotphi) {
  E_hencky_k[0] = CC[0] * (T_k[0] + c_cotphi) 
                + CC[1] * (T_k[1] + c_cotphi) 
                + CC[2] * (T_k[2] + c_cotphi);
  E_hencky_k[1] = CC[3] * (T_k[0] + c_cotphi) 
                + CC[4] * (T_k[1] + c_cotphi) 
                + CC[5] * (T_k[2] + c_cotphi);
  E_hencky_k[2] = CC[6] * (T_k[0] + c_cotphi) 
                + CC[7] * (T_k[1] + c_cotphi) 
                + CC[8] * (T_k[2] + c_cotphi);
}

/**************************************************************/

static int __update_internal_variables_elastic(double *T,
                                               const double *D_phi,
                                               const double *T_tr,
                                               const double *eigvec_b_e_tr,
                                               double c_cotphi) {
  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;


#if NumberDimensions == 2
  double T_aux[4] = {0.0, 0.0, 0.0, 0.0};
#else
  double T_aux[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif
  double T_A = 0.0;

  for (unsigned A = 0; A < Ndim; A++) {

    T_A = T_tr[A] + c_cotphi;

    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        T_aux[i * Ndim + j] +=
            T_A * eigvec_b_e_tr[A + i * Ndim] * eigvec_b_e_tr[A + j * Ndim];
      }
    }
  }

#if NumberDimensions == 2
  T[0] = T_aux[0];
  T[1] = T_aux[1];
  T[2] = T_aux[2];
  T[3] = T_aux[3];
  T[4] = T_tr[2] + c_cotphi;
#else
  T[0] = T_aux[0];
  T[1] = T_aux[1];
  T[2] = T_aux[2];
  T[3] = T_aux[3];
  T[4] = T_aux[4];
  T[5] = T_aux[5];
  T[6] = T_aux[6];
  T[7] = T_aux[7];
  T[8] = T_aux[8];
#endif

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
#if NumberDimensions == 2
  puts("Nominal stress tensor");
  printf("%f %f %f \n", T[0], T[1], 0.0);
  printf("%f %f %f \n", T[2], T[3], 0.0);
  printf("%f %f %f \n", 0.0, 0.0, T[4]);
#endif
#endif
#endif

  return STATUS;
}

/**************************************************************/

static void __elastic_tangent_moduli(
  double * C_ep, 
  const double * AA)
{
#if NumberDimensions == 2
    C_ep[0] = AA[0];
    C_ep[1] = AA[1];
    C_ep[2] = AA[3];
    C_ep[3] = AA[4];
#else
    C_ep[0] = AA[0];
    C_ep[1] = AA[1];
    C_ep[2] = AA[2];
    C_ep[3] = AA[3];
    C_ep[4] = AA[4];
    C_ep[5] = AA[5];
    C_ep[6] = AA[6]; 
    C_ep[7] = AA[7];
    C_ep[8] = AA[8];           
#endif
}

/**************************************************************/

static void __kappa(double *kappa, const double *a, double Lambda, double I1,
                   double alpha) {
  kappa[0] = a[0] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  kappa[1] = alpha * kappa[0];
}

/**************************************************************/

static void __d_kappa_phi_d_stress(double *d_kappa_phi_d_stress, const double *a,
                                  double Lambda, double I1) {
  d_kappa_phi_d_stress[0] =
      a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  d_kappa_phi_d_stress[1] =
      a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  d_kappa_phi_d_stress[2] =
      a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
}

/**************************************************************/

static void __d_kappa_phi_d_lambda(double *d_kappa_phi_d_lambda, const double *a,
                                  double Lambda, double I1) {
  *d_kappa_phi_d_lambda =
      (1 - a[2] * Lambda) * a[0] * exp(a[1] * I1) * exp(-a[2] * Lambda);
}

/**************************************************************/

static void __F(double * F, double kappa_phi, double I1, double I2, double I3) {
  double K1 = 9.0 + kappa_phi;

  *F = cbrt(K1 * I3) - cbrt(I1 * I2);
}

/**************************************************************/

static void __d_F_d_stress(double *d_F_d_stress, const double *T_k, double I1,
                          double I2, double I3,double kappa_phi) {
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
                            double kappa_psi, double I1, double I2, double I3) {
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
      dd_G_dd_stress[A * 3 + B] =
          (1.0 / 3.0) * cbrt(K2 * I3) *
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

  dd_G_d_stress_d_kappa_psi[0] = (cbrt(I3) / (3.0 * T_k[0])) / (3.0 * pow(cbrt(K2), 2));
  dd_G_d_stress_d_kappa_psi[1] = (cbrt(I3) / (3.0 * T_k[1])) / (3.0 * pow(cbrt(K2), 2));
  dd_G_d_stress_d_kappa_psi[2] = (cbrt(I3) / (3.0 * T_k[2])) / (3.0 * pow(cbrt(K2), 2));
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
  lapack_int INFO;
  lapack_int N_rows = 5;
  lapack_int N_cols = 5;
  lapack_int LDA = 5;

  // Compute 1-norm
  ANORM = LAPACKE_dlange(LAPACK_ROW_MAJOR,'1',N_rows, N_cols, Tangent_Matrix,LDA); 	

  // Compute the Reciprocal condition number
  INFO = LAPACKE_dgecon(LAPACK_ROW_MAJOR,'1', N_rows, Tangent_Matrix, LDA, ANORM, RCOND);

  if (INFO < 0) {
    fprintf(stderr,""RED"Error in LAPACKE_dgecon() : the %i-th argument had an illegal value "RESET"\n",
           abs(INFO));
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __solver(double *Tangent_Matrix, double *Residual) {
  lapack_int Order = 5;
  lapack_int LDA = 5;
  lapack_int LDB = 5;
  char TRANS = 'T'; /* (Transpose) */
  lapack_int INFO = 3;
  lapack_int IPIV[5] = {0, 0, 0, 0, 0};
  lapack_int NRHS = 1;

  //  Compute the LU factorization
  INFO = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,Order,Order,Tangent_Matrix,LDA,IPIV);

  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(
          stderr,
          "" RED
          "Error in LAPACKE_dgetrf(): the %i-th argument had an illegal value" RESET
          "",
          abs(INFO));
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in LAPACKE_dgetrf(): D_phi_mT(%i,%i) %s \n %s \n %s \n %s" RESET
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
  INFO = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'T',Order,NRHS, Tangent_Matrix, LDA,IPIV,Residual,LDB);

  if (INFO) {
    fprintf(stderr, ""RED"Error in LAPACKE_dgetrs() "RESET"\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __update_internal_variables_plastic(
    double *T, double *eps_n1, double *kappa_n1, const double *D_phi,
    const double *T_tr_k, const double *eigvec_b_e_tr, double Lambda_k,
    double kappa_phi_k, double c_cotphi) {

  int Ndim = NumberDimensions;

  // Update hardening parameters
  *eps_n1 = Lambda_k;
  *kappa_n1 = kappa_phi_k;

#if NumberDimensions == 2
  double T_aux[4] = {0.0, 0.0, 0.0, 0.0};
#else
  double T_aux[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif
  double T_A = 0.0;

  for (unsigned A = 0; A < Ndim; A++) {

    T_A = T_tr_k[A] + c_cotphi;

    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        T_aux[i * Ndim + j] +=
            T_A * eigvec_b_e_tr[A * Ndim + i] * eigvec_b_e_tr[A * Ndim + j];
      }
    }
  }

#if NumberDimensions == 2
  T[0] = T_aux[0];
  T[1] = T_aux[1];
  T[2] = T_aux[2];
  T[3] = T_aux[3];
  T[4] = T_tr_k[2] + c_cotphi;
#else
  T[0] = T_aux[0];
  T[1] = T_aux[1];
  T[2] = T_aux[2];
  T[3] = T_aux[3];
  T[4] = T_aux[4];
  T[5] = T_aux[5];
  T[6] = T_aux[6];
  T[7] = T_aux[7];
  T[8] = T_aux[8];
#endif  

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0
#if NumberDimensions == 2
  puts("Nominal stress tensor");
  printf("%f %f %f \n", T[0], T[1], 0.0);
  printf("%f %f %f \n", T[2], T[3], 0.0);
  printf("%f %f %f \n", 0.0, 0.0, T[4]);
#endif
#endif
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __elastoplastic_tangent_moduli(
  double * C_ep, 
  const double * CC,
  const double * dd_G_dd_stress,
  double delta_lambda_k)
{
  double C_ep_aux[9];

  C_ep_aux[0] = CC[0] + delta_lambda_k * dd_G_dd_stress[0];
  C_ep_aux[1] = CC[1] + delta_lambda_k * dd_G_dd_stress[1];
  C_ep_aux[2] = CC[2] + delta_lambda_k * dd_G_dd_stress[2];
  C_ep_aux[3] = CC[3] + delta_lambda_k * dd_G_dd_stress[3];
  C_ep_aux[4] = CC[4] + delta_lambda_k * dd_G_dd_stress[4];
  C_ep_aux[5] = CC[5] + delta_lambda_k * dd_G_dd_stress[5];
  C_ep_aux[6] = CC[6] + delta_lambda_k * dd_G_dd_stress[6];
  C_ep_aux[7] = CC[7] + delta_lambda_k * dd_G_dd_stress[7];
  C_ep_aux[8] = CC[8] + delta_lambda_k * dd_G_dd_stress[8];
    
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N, &N, C_ep_aux, &LDA, IPIV, &INFO);
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
             "Error in dgetrf_(): C_ep_aux(%i,%i) %s \n %s \n %s \n %s " RESET
             "\n",
             INFO, INFO, "is exactly zero. The factorization",
             "has been completed, but the factor C_ep_aux is exactly",
             "singular, and division by zero will occur if it is used",
             "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  dgetri_(&N, C_ep_aux, &LDA, IPIV, WORK, &LWORK, &INFO);
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(stderr, "" RED "%s: the %i-th argument %s" RESET "\n",
              "Error in dgetri_()", abs(INFO), "had an illegal value");
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetri_(): C_ep_aux(%i,%i) %s \n %s \n %s \n %s " RESET
              "\n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor C_ep_aux is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

#if NumberDimensions == 2
    C_ep[0] = C_ep_aux[0];
    C_ep[1] = C_ep_aux[1];
    C_ep[2] = C_ep_aux[3];
    C_ep[3] = C_ep_aux[4];
#endif
  
  return EXIT_SUCCESS;
}


/**************************************************************/