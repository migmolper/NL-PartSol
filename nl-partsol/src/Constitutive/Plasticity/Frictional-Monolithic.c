
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

#include "Macros.h"
#include "Types.h"
#include "Globals.h"


/*
  Define local global variables
*/
int Particle_Idx;
double Error0;
bool Is_Matsuoka_Nakai;
bool Is_Lade_Duncan;
bool Is_Modified_Lade_Duncan;

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

static int __elastic_tangent(
    double * CC /**< [out] Elastic compliance */, 
    double * AA /**< [out] Elastic matrix */,
    double E /**< [in] Young modulus */, 
    double nu /**< [in] Poisson ratio */,
    double K /**< [in] Lamé parameter */, 
    double G /**< [in] Shear modulus */);

static int __trial_elastic(
    double *T_tr /**< [out] Trial elastic stress tensor*/, 
    double * I1 /**< [out] First invariant of T_tr */, 
    double * I2 /**< [out] Second invariant of T_tr */, 
    double * I3 /**< [out] Third invariant of T_tr */,
    const double * E_hencky_trial /**< [in] Henky strain (trial) */, 
    const double * AA /**< [in] Elastic matrix */); 

static int __E_hencky(
    double * E_hencky_k /**< [out] Henky strain (iter k) */, 
    const double * T_k  /**< [in] Local stress tensor (iter k) */, 
    const double * CC /**< [in] Elastic compliance */);

static int __update_internal_variables_elastic(
    double *Stress /**< [in/out] Nominal stress tensor */,
    const double *D_phi /**< [in] Total deformation gradient. */,
    const double *T_tr /**< [in] Elastic stress tensor */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */);

static int __kappa(
    double *kappa /**< [out] Hardening vector */, 
    const double * a /**< [in] Vector with fit parameters (hardening) */, 
    double Lambda /**< [in] Total plastic multiplier */, 
    double I1 /**< [in] First invariant of the stress tensor */, 
    double alpha /**< [in] Dilatance parameter*/);

static int __d_kappa_phi_d_stress(
    double *d_kappa_phi_d_stress /**< [out] Stress derivative of kappa[0] */, 
    const double * a /**< [in] Vector with fit parameters (hardening) */,
    double Lambda /**< [in] Total plastic multiplier */, 
    double I1 /**< [in] First invariant of the stress tensor */);

static int __d_kappa_phi_d_lambda(
    double * d_kappa_phi_d_lambda /**< [out] Lambda derivative of kappa[0] */, 
    const double * a /**< [in] Vector with fit parameters (hardening) */,
    double Lambda /**< [in] Total plastic multiplier */, 
    double I1 /**< [in] First invariant of the stress tensor */);

static double __F(
    double c0 /**< [in] Yield function fit parameter */, 
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
    double * d_F_d_kappa_phi /**< [out] Yield function derivative (kappa[0]) */,
    double I1 /**< [in] First invariant of the stress tensor */,
    double I3 /**< [in] Third invariant of the stress tensor */, 
    double c0 /**< [in] Yield function fit parameter */,
    double m /**< [in] Yield function fit parameter */, 
    double pa /**< [in] Atmospheric pressure */,
    double kappa_phi /**< [in] Friction angle hardening */);

static double __G(
    double c0 /**< [in] Yield function fit parameter */, 
    double kappa_psi /**< [in] Dilatance angle hardening */, 
    double pa /**< [in] Atmospheric pressure */,
    double I1 /**< [in] First invariant of the stress tensor */, 
    double I2 /**< [in] Second invariant of the stress tensor */, 
    double I3 /**< [in] Third invariant of the stress tensor */,
    double m /**< [in] Yield function fit parameter */);

static int __d_G_d_stress(
    double *d_G_d_stress /**< [out] Plastic potential function derivative (stress) */,
    const double *T_k /**< [in] Local stress tensor (iter k) */, 
    double I1 /**< [in] First invariant of the stress tensor */,
    double I2 /**< [in] Second invariant of the stress tensor */, 
    double I3 /**< [in] Third invariant of the stress tensor */, 
    double c0 /**< [in] Yield function fit parameter */,
    double kappa_psi /**< [in] Dilatance angle hardening */, 
    double pa /**< [in] Atmospheric pressure */, 
    double m /**< [in] Yield function fit parameter */);

static int __dd_G_dd_stress(
    double *dd_G_dd_stress /**< [out] Plastic potential hessian (stress) */, 
    const double *T_k /**< [in] Local stress tensor (iter k) */,
    double kappa_psi /**< [in] Dilatance angle hardening */, 
    double I1 /**< [in] First invariant of the stress tensor */,
    double I2 /**< [in] Second invariant of the stress tensor */, 
    double I3 /**< [in] Third invariant of the stress tensor */, 
    double m /**< [in] Yield function fit parameter */,
    double pa /**< [in] Atmospheric pressure */, 
    double c0 /**< [in] Yield function fit parameter */);

static int __dd_G_d_stress_d_kappa_psi(
    double *dd_G_d_stress_d_kappa_psi /**< [out] Plastic potential deriv */, 
    const double *T_k /**< [in] Local stress tensor (iter k) */, 
    double I1 /**< [in] First invariant of the stress tensor */, 
    double I3 /**< [in] Third invariant of the stress tensor */,
    double m /**< [in] Yield function fit parameter */, 
    double pa /**< [in] Atmospheric pressure */, 
    double c0 /**< [in] Yield function fit parameter */, 
    double kappa_psi /**< [in] Dilatance angle hardening */); 

static int __residual(
    double *Residual /**< [out] Residual of the problem */, 
    double * Error_k /**< [out] Norm of the residual */,
    const double *E_hencky_trial /**< [in] Henky strain (trial) */,
    const double * E_hencky_k /**< [in] Henky strain (iter k) */,
    const double *d_G_d_stress /**< [in] Plastic potential function derivative (stress) */,
    const double *kappa_k /**< [in] Hardening vector (iter k) */,
    const double *kappa_hat /**< [in] Hardening vector (eval) */, 
    double F_k /**< [in] Yield function evaluation (iter k) */,
    double delta_lambda_k); 

static void assemble_tangent_matrix(double *, double *, double *, double *,
                                    double, double, Model_Parameters);

static double compute_condition_number(double *);

static bool check_convergence(double, int, int);

static void solver(double *, double *, double *);
static void update_variables(double *, double *, double *, double *, double *,
                             double, double, double);
static State_Parameters fill_Outputs(double *, double *, double *, double,
                                     double, double);

/**************************************************************/

int Frictional_Monolithic(State_Parameters IO_State, Material MatProp)
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
  double rad_friction_angle = (PI__MatrixLib__ / 180.0) *MatProp.phi_Frictional;
  double m;
  double c0;  
 if (strcmp(MatProp.Yield_Function_Frictional, "Matsuoka-Nakai") == 0) {
    m = 0.0;
    c0 = 9;
    Is_Matsuoka_Nakai = true;
    Is_Lade_Duncan = false;
    Is_Modified_Lade_Duncan = false;
  } else if (strcmp(MatProp.Yield_Function_Frictional, "Lade-Duncan") == 0) {
    m = 0.0;
    c0 = 27;
    Is_Matsuoka_Nakai = false;
    Is_Lade_Duncan = true;
    Is_Modified_Lade_Duncan = false;
  } else if (strcmp(MatProp.Yield_Function_Frictional,
                    "Modified-Lade-Duncan") == 0) {
    m = MatProp.m_Frictional;
    c0 = 27;
    Is_Matsuoka_Nakai = false;
    Is_Lade_Duncan = false;
    Is_Modified_Lade_Duncan = true;
  }
  double a[3] = {0.0,0.0,0.0};
  a[0] = MatProp.a_Hardening_Borja[0];
  a[1] = MatProp.a_Hardening_Borja[1];
  a[2] = MatProp.a_Hardening_Borja[2];
  
  double CC[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double AA[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  __elastic_tangent(CC,AA ,E,nu,K,G);


  // Define scalar variables internal variables
  double F_k1, F_k2, F_0 = 0.0;
  double I1, I2, I3 = 0.0;
  double Lambda_n = *IO_State.Equiv_Plast_Str;
  double Lambda_k1, Lambda_k2 = Lambda_n;
  double delta_lambda_k1, delta_lambda_k2 =  0.0;

  // Define tensorial internal variables
  double T_tr[3] = {0.0, 0.0, 0.0}; // Trial elastic stress
  double T_k1[3] = {0.0, 0.0, 0.0}; // Stress iteration k
  double T_k2[3] = {0.0, 0.0, 0.0}; // Stress iteration k (line search)
  double kappa_k1[2] = {0.0,0.0}; // Hardening iteration k
  kappa_k1[0] = (*IO_State.Kappa);
  kappa_k1[1] = alpha * (*IO_State.Kappa);
  double kappa_k2[2]; // Hardening iteration k (line search)
  kappa_k2[0] = (*IO_State.Kappa);
  kappa_k2[1] = alpha * (*IO_State.Kappa);
  double kappa_hat[3] = {0.0,0.0,0.0};
  double d_G_d_stress[3] = {0.0, 0.0, 0.0}; // Plastic flow
  double d_F_d_stress[3] = {0.0,0.0,0.0};
  double dd_G_dd_stress[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double d_kappa_phi_d_stress[3] = {0.0,0.0,0.0};
  double d_kappa_phi_d_lambda = 0.0;
  double d_F_d_stress[3] = {0.0,0.0,0.0};
  double d_F_d_kappa1 = 0;

  /*
    Initialize Newton-Raphson solver
  */
  double Residual[5];
  double D_Residual[5];
  double dd_G_dd_stress[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double dd_G_d_stress_d_kappa_psi[3] = {0.0,0.0,0.0};
  double Tangent_Matrix[25];
  double Norm_Residual_1;
  double Norm_Residual_2;
  double delta = 1;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;
  bool Convergence = false;

  __trial_elastic(T_tr, &I1, &I2, &I3, E_hencky_trial, AA);

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
    
    if ((fabs(I1) < TOL_Radial_Returning) && (m != 0)) {
      fprintf(stderr, "" RED "%s: %s %i %s" RED " \n", "Error in Frictional_Monolithic()",
            "The stress state of particle", Particle_Idx, "has reach the apex");
      return EXIT_FAILURE;
    }

    STATUS = __E_hencky(E_hencky_k1,T_k1,CC);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __E_hencky" RESET "\n");
      return EXIT_FAILURE;
    }

    STATUS = __kappa(kappa_hat, a, Lambda_k1, I1, alpha);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __kappa" RESET "\n");
      return EXIT_FAILURE;
    }

    STATUS = __d_G_d_stress(d_G_d_stress, T_tr, I1, I2, I3, c0, kappa_k1[1], pa, m);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __d_G_d_stress" RESET "\n");
      return EXIT_FAILURE;
    }

    F_k1 = __F(c0,kappa_k1[0], pa, I1, I2, I3, m);


    STATUS = __residual(Residual, &Norm_Residual_1, E_hencky_trial,E_hencky_k1, d_G_d_stress, kappa_k1,kappa_hat, F_k1, delta_lambda_k1);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __residual" RESET "\n");
      return EXIT_FAILURE;
    }

    Convergence = check_convergence(Norm_Residual_1, Iter, MaxIter);

    while (Convergence == false) {

      if (Convergence == false) {

        delta = 1;

        Iter++;

        __d_kappa_phi_d_stress(d_kappa_phi_d_stress, a, Lambda_k1, I1);

        __d_kappa_phi_d_lambda(&d_kappa_phi_d_lambda, a, Lambda_k1, I1);

        __d_F_d_stress(d_F_d_stress, T_k1, I1, I2, I3, c0, kappa_k1[0], pa, m);

        __d_F_d_kappa_phi(&d_F_d_kappa1,I1,I3, c0,m, pa, kappa_k1[0]);

        __dd_G_dd_stress(dd_G_dd_stress, T_k1, kappa_k1[1], I1, I2, I3, m, pa, c0);

        __dd_G_d_stress_d_kappa_psi(dd_G_d_stress_d_kappa_psi, T_k1, I1, I3,m, pa, c0, kappa_k1[1]); 

        assemble_tangent_matrix(Tangent_Matrix, n, Stress_k,
                                kappa_k, delta_lambda_k, Lambda_k, Params);

        solver(Tangent_Matrix, Residual, D_Residual);

        T_k2[0] = T_k1[0];
        T_k2[1] = T_k1[1];
        T_k2[2] = T_k1[2];
        kappa_k2[0] = kappa_k1[0];
        kappa_k2[1] = kappa_k1[1];
        Lambda_k2 = Lambda_k1;
        delta_lambda_k2 = delta_lambda_k1;

        update_variables(D_Residual, Stress_k2, kappa_k2, &delta_lambda_k2,
                         &Lambda_k2, Lambda_n, alpha, delta);

        I1 = T_k1[0] + T_k1[1] + T_k1[2];
        I2 = T_k1[0] * T_k1[1] + T_k1[1] * T_k1[2] +
              T_k1[0] * T_k1[2];
        I3 = T_k1[0] * T_k1[1] * T_k1[2];

        Norm_Residual_2 =
            __residual(Residual, Stress_k2, E_hencky_trial, n,
                              kappa_k2, delta_lambda_k2, Lambda_k2, Params);

        if ((Norm_Residual_2 - Norm_Residual_1) > TOL_Radial_Returning) {
          while ((Norm_Residual_2 - Norm_Residual_1) > TOL_Radial_Returning) {

            delta =
                pow(delta, 2) * 0.5 * Norm_Residual_1 /
                (Norm_Residual_2 - delta * Norm_Residual_1 + Norm_Residual_1);

            if (delta < TOL_Radial_Returning) {
              Stress_k2[0] = Stress_k[0];
              Stress_k2[1] = Stress_k[1];
              Stress_k2[2] = Stress_k[2];
              kappa_k2[0] = kappa_k[0];
              kappa_k2[1] = kappa_k[1];
              Lambda_k2 = Lambda_k;
              delta_lambda_k2 = delta_lambda_k;
              break;
            } else {
              Stress_k2[0] = Stress_k[0];
              Stress_k2[1] = Stress_k[1];
              Stress_k2[2] = Stress_k[2];
              kappa_k2[0] = kappa_k[0];
              kappa_k2[1] = kappa_k[1];
              Lambda_k2 = Lambda_k;
              delta_lambda_k2 = delta_lambda_k;

              update_variables(D_Residual, Stress_k2, kappa_k2,
                               &delta_lambda_k2, &Lambda_k2, Lambda_n,
                               Params.alpha, delta);

              Norm_Residual_2 = __residual(
                  Residual, Stress_k2, Strain_e_tri, n, kappa_k2,
                  delta_lambda_k2, Lambda_k2, Params);
            }
          }
        }

        Stress_k[0] = Stress_k2[0];
        Stress_k[1] = Stress_k2[1];
        Stress_k[2] = Stress_k2[2];
        kappa_k[0] = kappa_k2[0];
        kappa_k[1] = kappa_k2[1];
        Lambda_k = Lambda_k2;
        delta_lambda_k = delta_lambda_k2;
        Norm_Residual_1 = Norm_Residual_2;

        Convergence = check_convergence(Norm_Residual_1, Iter, MaxIter);
      }
    }

    /*
      Update equivalent plastic strain and increment of plastic deformation
    */
    Outputs_VarCons =
        fill_Outputs(Increment_E_plastic, Stress_k, n, Lambda_k,
                     delta_lambda_k, kappa_k[0]);

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

  eigvec_b_e_tr[0] = d_phi[0] * b_e[0] * d_phi[0] + d_phi[0] * b_e[1] * d_phi[1] +
              d_phi[1] * b_e[2] * d_phi[0] + d_phi[1] * b_e[3] * d_phi[1];

  eigvec_b_e_tr[1] = d_phi[0] * b_e[0] * d_phi[2] + d_phi[0] * b_e[1] * d_phi[3] +
              d_phi[1] * b_e[2] * d_phi[2] + d_phi[1] * b_e[3] * d_phi[3];

  eigvec_b_e_tr[3] = d_phi[2] * b_e[0] * d_phi[0] + d_phi[2] * b_e[1] * d_phi[1] +
              d_phi[3] * b_e[2] * d_phi[0] + d_phi[3] * b_e[3] * d_phi[1];

  eigvec_b_e_tr[4] = d_phi[2] * b_e[0] * d_phi[2] + d_phi[2] * b_e[1] * d_phi[3] +
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
  dsyev_("V","L",&n, eigvec_b_e_tr, &lda, eigval_b_e_tr,&wkopt, &lwork,&info);
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

  dsyev_("V","L",&n, eigvec_b_e_tr, &lda, eigval_b_e_tr,work, &lwork,&info);
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

static int __elastic_tangent(double * CC, double * AA ,double E, double nu,double K, double G)
{
  CC[0] = 1.0 / E;
  CC[1] = -nu / E;
  CC[2] = -nu / E;

  CC[3] = -nu / E;
  CC[4] = 1.0 / E;
  CC[5] = -nu / E;

  CC[6] = -nu / E;
  CC[7] = -nu / E;
  CC[8] = 1.0 / E;

  AA[0] = K + 2*G;
  AA[1] = K;
  AA[2] = K;

  AA[3] = K;
  AA[4] = K + 2*G;
  AA[5] = K;

  AA[6] = K;
  AA[7] = K;
  AA[8] = K + 2*G;

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __trial_elastic(double *T_tr, double * I1, double * I2, double * I3,
                           const double * E_hencky_trial, const double * AA) {

#if NumberDimensions == 2
  T_tr[0] = AA[0]*E_hencky_trial[0] + AA[1]*E_hencky_trial[1] + AA[2]*E_hencky_trial[2];
  T_tr[1] = AA[3]*E_hencky_trial[0] + AA[4]*E_hencky_trial[1] + AA[5]*E_hencky_trial[2];
  T_tr[2] = AA[6]*E_hencky_trial[0] + AA[7]*E_hencky_trial[1] + AA[8]*E_hencky_trial[2];

#else
  No esta implementado
#endif

  *I1 = T_tr[0] + T_tr[1] + T_tr[2];
  *I2 = T_tr[0] * T_tr[1] + T_tr[1] * T_tr[2] +
              T_tr[0] * T_tr[2];
  *I3 = T_tr[0] * T_tr[1] * T_tr[2];

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __E_hencky(double * E_hencky_k, const double * T_k, const double * CC)
{
  #if NumberDimensions == 2
  
  E_hencky_k[0] = CC[0]*T_k[0] + CC[1]*T_k[1] + CC[2]*T_k[2];
  E_hencky_k[1] = CC[3]*T_k[0] + CC[4]*T_k[1] + CC[5]*T_k[2];
  E_hencky_k[2] = CC[6]*T_k[0] + CC[7]*T_k[1] + CC[8]*T_k[2];
  
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

static int __kappa(double * kappa, const double * a, double Lambda, 
                  double I1, double alpha) {

  kappa[0] = a[0] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  kappa[1] = alpha * kappa[0];

  if(kappa[0] < 0.0)
  {
    fprintf(stderr, "" RED "Negative value of kappa: %f " RESET "\n",
            kappa[0]);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __d_kappa_phi_d_stress
    (double *d_kappa_phi_d_stress, 
    const double * a,
    double Lambda, 
    double I1) {

#if NumberDimensions == 2

  d_kappa_phi_d_stress[0] = a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  d_kappa_phi_d_stress[1] = a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);
  d_kappa_phi_d_stress[2] = a[0] * a[1] * Lambda * exp(a[1] * I1) * exp(-a[2] * Lambda);

#else
  No esta implementado
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __d_kappa_phi_d_lambda(double * d_kappa_phi_d_lambda, 
                              const double * a,
                              double Lambda, double I1) {

  *d_kappa_phi_d_lambda = (1 - a[2] * Lambda) * a[0] * exp(a[1] * I1) * exp(-a[2] * Lambda);

  return EXIT_SUCCESS;
}

/**************************************************************/

static double __F(double c0, double kappa_phi, double pa,
                  double I1, double I2, double I3, double m) {

  double K1 = c0 + kappa_phi * pow(pa / I1, m);

  double PHI = 0.0;

  if (Is_Matsuoka_Nakai) {
    PHI = cbrt(K1 * I3) - cbrt(I1 * I2);
  } else if (Is_Lade_Duncan) {
    PHI = cbrt(K1 * I3) - I1;
  } else if (Is_Modified_Lade_Duncan) {
    PHI = cbrt(K1 * I3) - I1;
  } 

  return PHI;
}

/**************************************************************/

static int __d_F_d_stress(double *d_F_d_stress, const double *T_k,
                                      double I1, double I2, double I3,
                                      double c0, double kappa_phi, double pa,
                                      double m) {
  double Grad_f[3];

  double K1 = c0 + kappa_phi * pow(pa / I1, m);
  double b1 = m * kappa_phi * (pow(pa / I1, m)) * (cbrt(I3) / I1);

  if (Is_Matsuoka_Nakai) {
    Grad_f[0] =
        (I1 * (I1 - T_k[0]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2));
    Grad_f[1] =
        (I1 * (I1 - T_k[1]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2));
    Grad_f[2] =
        (I1 * (I1 - T_k[2]) + I2) / (3.0 * pow(cbrt(I1 * I2), 2));
  } else if (Is_Lade_Duncan) {
    Grad_f[0] = 1.0;
    Grad_f[1] = 1.0;
    Grad_f[2] = 1.0;
  } else if (Is_Modified_Lade_Duncan) {
    Grad_f[0] = 1.0;
    Grad_f[1] = 1.0;
    Grad_f[2] = 1.0;
  } 

  d_F_d_stress[0] =
      cbrt(K1 * I3) / (3.0 * T_k[0]) - b1 / (3.0 * pow(cbrt(K1), 2.0)) - Grad_f[0];
  d_F_d_stress[1] =
      cbrt(K1 * I3) / (3.0 * T_k[1]) - b1 / (3.0 * pow(cbrt(K1), 2.0)) - Grad_f[1];
  d_F_d_stress[2] =
      cbrt(K1 * I3) / (3.0 * T_k[2]) - b1 / (3.0 * pow(cbrt(K1), 2.0)) - Grad_f[2];

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __d_F_d_kappa_phi(double * d_F_d_kappa_phi,double I1, 
                             double I3, double c0,
                             double m, double pa,double kappa_phi) {

  double K1 = c0 + kappa_phi * pow(pa / I1, m);

  *d_F_d_kappa_phi = cbrt(I3) / (3.0 * pow(cbrt(K1), 2.0)) * pow(pa / I1, m);

  return EXIT_SUCCESS;
}

/**************************************************************/


static double __G(double c0, double kappa2, double pa,
                                  double I1, double I2, double I3,
                                  double m) {

  double K2 = c0 + kappa2 * pow(pa / I1, m);

  double G = 0.0;

  if (Is_Matsuoka_Nakai) {
    G = cbrt(K2 * I3) - cbrt(I1 * I2);
  } else if (Is_Lade_Duncan) {
    G = cbrt(K2 * I3) - I1;
  } else if (Is_Modified_Lade_Duncan) {
    G = cbrt(K2 * I3) - I1;
  } 

  return G;
}

/**************************************************************/

static int __d_G_d_stress(double *d_G_d_stress,
                                          const double * T_k, double I1,
                                          double I2, double I3, double c0,
                                          double kappa_psi, double pa, double m) {

  double Grad_g[3] = {0.0,0.0,0.0};
  double K2 = c0 + kappa_psi * pow(pa / I1, m);
  double b2 = m * kappa_psi * (pow(pa / I1, m)) * (cbrt(I3) / I1);

  if (Is_Matsuoka_Nakai) {
    Grad_g[0] = (I1 * (I1 -  T_k[0]) + I2) / (3 * pow(cbrt(I1 * I2), 2));
    Grad_g[1] = (I1 * (I1 -  T_k[1]) + I2) / (3 * pow(cbrt(I1 * I2), 2));
    Grad_g[2] = (I1 * (I1 -  T_k[2]) + I2) / (3 * pow(cbrt(I1 * I2), 2));
  } else if (Is_Lade_Duncan) {
    Grad_g[0] = 1.0;
    Grad_g[1] = 1.0;
    Grad_g[2] = 1.0;
  } else if (Is_Modified_Lade_Duncan) {
    Grad_g[0] = 1.0;
    Grad_g[1] = 1.0;
    Grad_g[2] = 1.0;
  }

  d_G_d_stress[0] =
      cbrt(K2 * I3) / (3 *  T_k[0]) - b2 / (3 * pow(cbrt(K2), 2)) - Grad_g[0];
  d_G_d_stress[1] =
      cbrt(K2 * I3) / (3 *  T_k[1]) - b2 / (3 * pow(cbrt(K2), 2)) - Grad_g[1];
  d_G_d_stress[2] =
      cbrt(K2 * I3) / (3 *  T_k[2]) - b2 / (3 * pow(cbrt(K2), 2)) - Grad_g[2];

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __dd_G_dd_stress(double *dd_G_dd_stress, const double *T_k,
                                            double kappa_psi, double I1,
                                            double I2, double I3, double m,
                                            double pa, double c0) {

  double K2 = c0 + kappa_psi * pow(pa / I1, m);
  
  double b2 = m * kappa_psi * (pow(pa / I1, m)) * (cbrt(I3) / I1);
  
  double Grad_K2[3] = {0.0,0.0,0.0};
  Grad_K2[0] = -(m * kappa_psi / I1) * pow(pa / I1, m);
  Grad_K2[1] = -(m * kappa_psi / I1) * pow(pa / I1, m);
  Grad_K2[2] = -(m * kappa_psi / I1) * pow(pa / I1, m);

  double Grad_b2[3] = {0.0,0.0,0.0};
  Grad_b2[0] = (b2 / I1) * (I1 / (3 * T_k[0]) - m - 1.0);
  Grad_b2[1] = (b2 / I1) * (I1 / (3 * T_k[1]) - m - 1.0);
  Grad_b2[2] = (b2 / I1) * (I1 / (3 * T_k[2]) - m - 1.0);

  double Hess_g[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  if (Is_Matsuoka_Nakai) {
    
    double d_g_d_stress[3];

    d_g_d_stress[0] =
        (I1 * (I1 - T_k[0]) + I2) / (3 * pow(cbrt(I1 * I2), 2));
    d_g_d_stress[1] =
        (I1 * (I1 - T_k[1]) + I2) / (3 * pow(cbrt(I1 * I2), 2));
    d_g_d_stress[2] =
        (I1 * (I1 - T_k[2]) + I2) / (3 * pow(cbrt(I1 * I2), 2));

    for (unsigned A = 0; A < 3; A++) {
      for (unsigned B = 0; B < 3; B++) {
        Hess_g[A * 3 + B] =
            (3.0 * I1 - T_k[A] - T_k[B] - I1 * (A == B)) /
                (3.0 * pow(cbrt(I1 * I2), 2)) -
            (2.0 / cbrt(I1 * I2)) * d_g_d_stress[A] * d_g_d_stress[B];
      }
    }
  } 


  for (unsigned A = 0; A < 3; A++) {
    for (unsigned B = 0; B < 3; B++) {
      dd_G_dd_stress[A * 3 + B] = (1.0 / 3.0) * cbrt(K2 * I3) *
                              (1.0 / (3 * T_k[A] * T_k[B]) -
                               1.0 * (A == B) / pow(T_k[A], 2)) +
                          (cbrt(I3) / T_k[A] + 2.0 * b2 / K2) * Grad_K2[B] /
                              (9.0 * pow(cbrt(K2), 2)) -
                          Grad_b2[B] / (3.0 * pow(cbrt(K2), 2)) -
                          Hess_g[A * 3 + B];
    }
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __dd_G_d_stress_d_kappa_psi(
    double *dd_G_d_stress_d_kappa_psi, const double *T_k, double I1, double I3,
    double m, double pa, double c0, double kappa_psi) {

  double K2 = c0 + kappa_psi * pow(pa / I1, m);
  double b2 = m * kappa_psi * (pow(pa / I1, m)) * (cbrt(I3) / I1);

  dd_G_d_stress_d_kappa_psi[0] = pow((pa / I1), m) *
                              (cbrt(I3) / (3.0 * T_k[0]) +
                               2.0 * b2 / (3.0 * K2) - m * cbrt(I3) / I1) /
                              (3.0 * pow(cbrt(K2), 2));
  dd_G_d_stress_d_kappa_psi[1] = pow((pa / I1), m) *
                              (cbrt(I3) / (3.0 * T_k[1]) +
                               2.0 * b2 / (3.0 * K2) - m * cbrt(I3) / I1) /
                              (3.0 * pow(cbrt(K2), 2));
  dd_G_d_stress_d_kappa_psi[2] = pow((pa / I1), m) *
                              (cbrt(I3) / (3.0 * T_k[2]) +
                               2.0 * b2 / (3.0 * K2) - m * cbrt(I3) / I1) /
                              (3.0 * pow(cbrt(K2), 2));
  return EXIT_SUCCESS;
}

/**************************************************************/

static int __residual(double *Residual, double * Error_k, const double *E_hencky_trial,
                      const double * E_hencky_k, const double *d_G_d_stress, const double *kappa_k,
                      const double *kappa_hat, double F_k, double delta_lambda_k) {
  
  Residual[0] = E_hencky_k[0] - E_hencky_trial[0] + delta_lambda_k * d_G_d_stress[0];
  Residual[1] = E_hencky_k[1] - E_hencky_trial[1] + delta_lambda_k * d_G_d_stress[1];
  Residual[2] = E_hencky_k[2] - E_hencky_trial[2] + delta_lambda_k * d_G_d_stress[2];
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

static void assemble_tangent_matrix(double *Tangent_Matrix,
                                    double *d_G_d_stress, double *T_k,
                                    double *kappa_k, double delta_lambda,
                                    double Lambda_k, double c0, const double * a,
                                    const double * dd_G_dd_stress,
                                    const double * CC) {

  /* First row */
  Tangent_Matrix[0] = CC[0] + delta_lambda * dd_G_dd_stress[0];
  Tangent_Matrix[1] = CC[1] + delta_lambda * dd_G_dd_stress[1];
  Tangent_Matrix[2] = CC[2] + delta_lambda * dd_G_dd_stress[2];
  Tangent_Matrix[3] = alpha * delta_lambda * dd_G_d_stress_d_kappa2[0];
  Tangent_Matrix[4] = d_G_d_stress[0];

  /* Second row */
  Tangent_Matrix[5] = CC[3] + delta_lambda * dd_G_dd_stress[3];
  Tangent_Matrix[6] = CC[4] + delta_lambda * dd_G_dd_stress[4];
  Tangent_Matrix[7] = CC[5] + delta_lambda * dd_G_dd_stress[5];
  Tangent_Matrix[8] = alpha * delta_lambda * dd_G_d_stress_d_kappa2[1];
  Tangent_Matrix[9] = d_G_d_stress[1];

  /* Third row */
  Tangent_Matrix[10] = CC[6] + delta_lambda * dd_G_dd_stress[6];
  Tangent_Matrix[11] = CC[7] + delta_lambda * dd_G_dd_stress[7];
  Tangent_Matrix[12] = CC[8] + delta_lambda * dd_G_dd_stress[8];
  Tangent_Matrix[13] = alpha * delta_lambda * dd_G_d_stress_d_kappa2[2];
  Tangent_Matrix[14] = d_G_d_stress[2];

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
  Tangent_Matrix[23] = d_F_d_kappa1;
  Tangent_Matrix[24] = 0.0;
}

/**************************************************************/

static bool check_convergence(double Error, int Iter, int MaxIter) {
  bool convergence;
  double Error_relative = 0.0;

  /*
    Compute relative error
  */
  if (Iter == 0) {
    Error0 = Error;
    Error_relative = Error / Error0;

    if (Error0 < TOL_Radial_Returning) {
      return true;
    }
  } else {
    Error_relative = Error / Error0;
  }

  if ((Error > TOL_Radial_Returning * 100) &&
      (Error_relative > TOL_Radial_Returning) && (Iter < MaxIter)) {
    return false;
  } else {

    if (Iter >= MaxIter) {
      fprintf(stderr, "%s: %s %i \n", "Error in Frictional_Monolithic()",
              "Maximm number of iteration reached in particle", Particle_Idx);
      fprintf(stderr, "Iter: %i | Error0: %e | Error: %e | Erro_rel: %e \n",
              Iter, Error0, Error, Error_relative);
    }

    return true;
  }
}

/**************************************************************/

static double compute_condition_number(double *Tangent_Matrix)
/*
  C = rcond(Tangent_Matrix) returns an estimate for the reciprocal condition of
  Tangent_Matrix in 1-norm.
*/
{

  double RCOND;
  double ANORM;
  int INFO;
  int N_rows = 5;
  int N_cols = 5;
  int LDA = IMAX(N_rows, N_cols);
  double *AUX_MEMORY = (double *)calloc(N_rows * N_cols, sizeof(double));
  double *WORK_ANORM = (double *)calloc(IMAX(1, N_rows), sizeof(double));
  double *WORK_RCOND = (double *)calloc(4 * N_rows, sizeof(double));
  int *IPIV = (int *)calloc(IMIN(N_rows, N_cols), sizeof(int));
  int *IWORK_RCOND = (int *)calloc(N_rows, sizeof(int));

  // Copy matrix because dgetrf_ is a destructive operation
  memcpy(AUX_MEMORY, Tangent_Matrix, N_rows * N_cols * sizeof(double));

  // Compute 1-norm
  ANORM = dlange_("1", &N_rows, &N_cols, AUX_MEMORY, &LDA, WORK_ANORM);

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N_rows, &N_cols, AUX_MEMORY, &LDA, IPIV, &INFO);

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
    exit(EXIT_FAILURE);
  }

  // Compute the Reciprocal condition number
  dgecon_("1", &N_rows, AUX_MEMORY, &LDA, &ANORM, &RCOND, WORK_RCOND,
          IWORK_RCOND, &INFO);

  if (INFO < 0) {
    printf("Error in Frictional_Monolithic() : the %i-th argument of dgecon_ "
           "had an illegal value\n",
           abs(INFO));
    exit(EXIT_FAILURE);
  }

  // Free auxiliar memory
  free(IPIV);
  free(WORK_ANORM);
  free(WORK_RCOND);
  free(IWORK_RCOND);
  free(AUX_MEMORY);

  return RCOND;
}

/**************************************************************/
static void solver(double *Tangent_Matrix, double *Residual,
                   double *D_Residual) {
  int Order = 5;
  int LDA = 5;
  int LDB = 5;
  char TRANS = 'T'; /* (Transpose) */
  int INFO = 3;
  int *IPIV = (int *)malloc(Order * sizeof(int));
  int NRHS = 1;

  double rcond = compute_condition_number(Tangent_Matrix);

  if (rcond < 1E-10) {
    fprintf(stderr, "%s: %s %i, %s: %e\n", "Error in Frictional_Monolithic",
            "Tangent_Matrix is near to singular matrix for particle",
            Particle_Idx, "rcond", rcond);
  }

  /*
    Generate auxiliar copy of the mass matrix to avoid destructive operations
  */
  memcpy(D_Residual, Residual, 5 * sizeof(double));

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
    exit(EXIT_FAILURE);
  }

  /*
    Solve
  */
  dgetrs_(&TRANS, &Order, &NRHS, Tangent_Matrix, &LDA, IPIV, D_Residual, &LDB,
          &INFO);
  free(IPIV);

  /*
    Check error messages in the LAPACK solver
  */
  if (INFO) {
    fprintf(stderr, "%s : %s %s %s \n", "Error in Frictional_Monolithic",
            "The function", "dgetrs_", "returned an error message !!!");
    exit(EXIT_FAILURE);
  }
}

/**************************************************************/

static void update_variables(double *D_Residual, double *Stress_k,
                             double *kappa_k, double *delta_lambda,
                             double *lambda_k, double lambda_n, double alpha,
                             double delta) {

  Stress_k[0] -= delta * D_Residual[0];
  Stress_k[1] -= delta * D_Residual[1];
  Stress_k[2] -= delta * D_Residual[2];
  kappa_k[0] -= delta * D_Residual[3];
  kappa_k[1] = alpha * kappa_k[0];
  *delta_lambda -= delta * D_Residual[4];
  *lambda_k = *delta_lambda + lambda_n;
}

/**************************************************************/

static State_Parameters fill_Outputs(double *Increment_E_plastic,
                                     double *Stress_k, double *Plastic_Flow,
                                     double Lambda_k, double delta_lambda,
                                     double kappa_k1) {
  State_Parameters Outputs_VarCons;

  *Outputs_VarCons.Equiv_Plast_Str = Lambda_k;
  *Outputs_VarCons.Kappa = kappa_k1;
  Outputs_VarCons.Stress = Stress_k;
  Outputs_VarCons.Increment_E_plastic = Increment_E_plastic;
  Outputs_VarCons.Increment_E_plastic[0] = delta_lambda * Plastic_Flow[0];
  Outputs_VarCons.Increment_E_plastic[1] = delta_lambda * Plastic_Flow[1];
  Outputs_VarCons.Increment_E_plastic[2] = delta_lambda * Plastic_Flow[2];

  return Outputs_VarCons;
}

/**************************************************************/
