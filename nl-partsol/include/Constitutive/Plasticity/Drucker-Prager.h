#ifndef _DRUCKER_PRAGER_CONSTITUTIVE_H_
#define _DRUCKER_PRAGER_CONSTITUTIVE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/*!
    double *eigval_b_e_tr /**< [out] Eigenvalues of b elastic trial.
    double *eigvec_b_e_tr /**< [out] Eigenvector of b elastic trial.
    const double *b_e /**< [in] (n) Elastic left Cauchy-Green.
    const double *d_phi /**< [in] Incremental deformation gradient.
*/
static int __compute_trial_b_e(
    double *eigval_b_e_tr /**< [out] Eigenvalues of b elastic trial. */,
    double *eigvec_b_e_tr /**< [out] Eigenvector of b elastic trial. */,
    const double *b_e /**< [in] (n) Elastic left Cauchy-Green.*/,
    const double *d_phi /**< [in] Incremental deformation gradient. */);
/**************************************************************/

static int __corrector_b_e(
    double *b_e /**< [out] (n+1) Elastic deformation gradient. */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */,
    const double *E_hencky_trial /**< [in] Corrected Henky strain */);
/**************************************************************/

static int __trial_elastic(
    double *T_tr_vol /**< [in/out] Volumetric elastic stress tensor. */,
    double *T_tr_dev /**< [in/out] Deviatoric elastic stress tensor. */,
    double *pressure /**< [out] First invariant of the stress tensor */,
    double *J2 /**< [out] Second invariant of the deviatoric stress tensor */,
    const double *E_hencky_trial, /**< [in] Trial elastic strain tensor. */
    double K /**< [in] First Lamé invariant. */,
    double G /**< [in] Second Lamé invariant. */,
    double p_ref /**< [in] Reference pressure. */);
/**************************************************************/

static int __update_internal_variables_elastic(
    double *Stress /**< [in/out] Nominal stress tensor */,
    const double *T_tr_vol /**< [in] Volumetric elastic stress tensor */,
    const double *T_tr_dev /**< [in] Deviatoric elastic stress tensor */,
    const double *eigvec_b_e_tr /**< [in] Eigenvector of b elastic trial. */);
/**************************************************************/

static int __tangent_moduli_elastic(
    double * C_ep /**< [out] Elastoplastic tanget moduli */,
    double K /**< [in] First Lamé invariant. */,
    double G /**< [in] Second Lamé invariant. */);
/**************************************************************/

static int __compute_plastic_flow_direction(
    double *n /**< [out] Plastic flow direction */,
    const double *T_tr_dev /**< [in] Deviatoric elastic stress tensor */,
    double J2 /**< [in] Second invariant of the deviatoric stress tensor */);
/**************************************************************/

static int __eps(
    double *eps_k /**< [out] Equivalent plastic strain*/,
    double d_gamma_k /**< [in] Discrete plastic multiplier */,
    double eps_n /**< [in] Equivalent plastic strain in the last step */,
    double alpha_Q /**< [in] Plastic potential parameter */);
/**************************************************************/    

static int __kappa(double *kappa_k /**< [out] Hardening function. */,
                           double kappa_0 /**< [in] Reference hardening */,
                           double exp_param /**< [in] Hardening exponential*/,
                           double eps_k /**< [in] Equivalent plastic strain*/,
                           double eps_0 /**< [in] Reference plastic strain */);
/**************************************************************/

/*!
    \param[out] d_kappa Derivative of the hardening function.
    \param[in] kappa_0 Reference hardening
    \param[in] eps_k Equivalent plastic strain
    \param[in] eps_0 Reference plastic strain
    \param[in] exp_param Hardening exponential
*/
static int __d_kappa(
    double *d_kappa /**< [out] Derivative of the hardening function. */,
    double kappa_0 /**< [in] Reference hardening */,
    double eps_k /**< [in] Equivalent plastic strain*/,
    double eps_0 /**< [in] Reference plastic strain */,
    double exp_param /**< [in] Hardening exponential*/);
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
static int __compute_pressure_limit(
    double *pressure_limit,
    double J2,
    double kappa_n,
    double d_kappa,
    double K,
    double G,
    double alpha_F,
    double alpha_Q,
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
static double __yield_function_classical(
    double pressure,
    double J2,
    double d_gamma_k,
    double kappa_k,
    double alpha_F,
    double alpha_Q,
    double beta,
    double K,
    double G);
/**************************************************************/

/*!
    \param[in] d_kappa_k Derivative of the hardening function.
    \param[in] K First Lamé invariant.
    \param[in] G Second Lamé invariant.
    \param[in] alpha_F Yield surface parameter I.
    \param[in] alpha_Q Plastic potential parameter.
    \param[in] beta Yield surface parameter II.
*/
static double __d_yield_function_classical(
    double d_kappa_k,
    double K,
    double G,
    double alpha_F,
    double alpha_Q,
    double beta);
/**************************************************************/

/*!
    \param[out] Increment_E_plastic  Increment plastic strain
    \param[out] Stress  Nominal stress tensor
    \param[out] eps_n1  Equivalent plastic strain
    \param[out] kappa_n1  Hardening function.
    \param[in] T_tr_vol Volumetric elastic stress tensor.
    \param[in] T_tr_dev Deviatoric elastic stress tensor.
    \param[in] eigvec_b_e_tr Eigenvector of b elastic trial.
    \param[in] n Plastic flow direction.
    \param[in] d_gamma_k Discrete plastic multiplier
    \param[in] alpha_Q Plastic potential parameter.
    \param[in] K First Lamé invariant.
    \param[in] G Second Lamé invariant.
    \param[in] eps_k Equivalent plastic strain
    \param[in] kappa_k Hardening function.
*/
static int __update_internal_variables_classical(
    double *Increment_E_plastic,
    double *Stress,
    double *eps_n1,
    double *kappa_n1,
    const double *T_tr_vol,
    const double *T_tr_dev,
    const double *eigvec_b_e_tr,
    const double *n,
    double d_gamma_k,
    double alpha_Q,
    double K,
    double G,
    double eps_k,
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
static double __yield_function_apex(
    double pressure,
    double d_gamma_k,
    double d_gamma_1,
    double kappa_k,
    double d_kappa_k,
    double K,
    double alpha_F,
    double alpha_Q,
    double beta);
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
static double __d_yield_function_apex(
    double d_gamma_k,
    double d_gamma_1,
    double d_kappa_k,
    double K,
    double alpha_F,
    double alpha_Q,
    double beta);
/**************************************************************/

/*!
    \param[out] Increment_E_plastic Increment plastic strain
    \param[out] Stress Nominal stress tensor
    \param[out] eps_n1 Equivalent plastic strain
    \param[out] kappa_n1 Hardening function.
    \param[in] T_tr_vol Volumetric elastic stress tensor.
    \param[in] T_tr_dev Deviatoric elastic stress tensor.
    \param[in] eigvec_b_e_tr Eigenvector of b elastic trial.
    \param[in] n Plastic flow direction.
    \param[in] d_gamma_k Discrete plastic multiplier
    \param[in] d_gamma_1 Discrete plastic multiplier I
    \param[in] alpha_Q Plastic potential parameter.
    \param[in] K First Lamé invariant.
    \param[in] G Second Lamé invariant.
    \param[in] eps_k Equivalent plastic strain*
    \param[in] kappa_k Hardening function.
*/
static int __update_internal_variables_apex(
    double *Increment_E_plastic,
    double *Stress,
    double *eps_n1,
    double *kappa_n1,
    const double *T_tr_vol,
    const double *T_tr_dev,
    const double *eigvec_b_e_tr,
    const double *n,
    double d_gamma_k,
    double d_gamma_1,
    double alpha_Q,
    double K,
    double G,
    double eps_k,
    double kappa_k);
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
static int __tangent_moduli_apex(
  double *C_ep, 
  const double *n, 
  double d_gamma_k, 
  double d_gamma_1, 
  double d_kappa_k,
  double K,
  double G, 
  double beta, 
  double alpha_F, 
  double alpha_Q);    
/**************************************************************/

/*!
  \fn int compute_1PK_Drucker_Prager(State_Parameters IO_State, Material MatProp);

  \brief Compute the behaviour of a Henky-Hyperelastic material \n
  using the Drucker-Prager yield criterium.

  \param[in] Input_SP : State parameters of the particle
  \param[in] MatProp : Material properties of the model
*/
int compute_1PK_Drucker_Prager(State_Parameters IO_State, Material MatProp);

#endif
