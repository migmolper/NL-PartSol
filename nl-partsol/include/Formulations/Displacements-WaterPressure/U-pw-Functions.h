/**
 * @file U-pw-Functions.h
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2022-05-15
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _UPW_FUNCTIONS_H_
#define _UPW_FUNCTIONS_H_

#include "Macros.h"
#include "Matlib.h"
#include "Types.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

// Shape functions and auxilar tools
#include "Nodes/Nodes-Tools.h"
#include "Nodes/Shape-Functions.h"

// HPC libs
#ifdef USE_OPENMP
#include <omp.h>
#endif

// Linear-Solver libs
#ifdef USE_PETSC
#include "Linear-Solvers/ksp-PETSC.h"
#else
#include "Linear-Solvers/dgetrs-LAPACK.h"
#endif

typedef struct {

#ifdef USE_PETSC
  Vec value;
  Vec d_value_dt;
  Vec d2_value_dt2;
#else
  double *value;
  double *d_value_dt;
  double *d2_value_dt2;
#endif

} Nodal_Field;

typedef struct {
  double alpha_1;
  double alpha_2;
  double alpha_3;
  double alpha_4;
  double alpha_5;
  double alpha_6;
  double epsilon;
  double DeltaTimeStep;
} Newmark_parameters;

// Material parameters
extern int Number_Soil_Water_Mixtures;
extern Mixture * Soil_Water_Mixtures;

/**************************************************************/

/**
 * @brief Update the local deformation of the particles \n
 * ensuring the local compatibility conditions
 *
 * @param D_U Increment of nodal kinetics
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @param STATUS Returns failure or success
 */
void local_compatibility_conditions__upw__(const Nodal_Field D_U,
                                           Mask ActiveNodes, Particle MPM_Mesh,
                                           Mesh FEM_Mesh, int *STATUS);
/**************************************************************/

/**
 * @brief 
 * 
 * @param MPM_Mesh 
 * @param FEM_Mesh 
 * @param STATUS 
 */
void constitutive_update__upw__(Particle MPM_Mesh, Mesh FEM_Mesh,
                                  int *STATUS);
/**************************************************************/

/**
 * @brief 
 * 
 * @param gradient_theta_n1_p 
 * @param dPw 
 * @param Pw_n 
 * @param d_shapefunction_n1_p 
 * @param Nnodes_p 
 */
void compute_gradient_kirchoff_pw__upw__(
  double * gradient_theta_n1_p,
  const double * dPw,
  const double * Pw_n,
  const double * d_shapefunction_n1_p,
  unsigned Nnodes_p);
/**************************************************************/  

/**
 * @brief 
 * 
 * @param b_n1_p 
 * @param a_n_p 
 * @param da_I 
 * @param shapefunction_n_p 
 * @param ListNodes_p 
 * @param ActiveNodes 
 */
void compute_total_acceleration__upw__(double *b_n1_p, const double *a_n_p,
                                         const double *da_I,
                                         const double *shapefunction_n_p,
                                         const ChainPtr ListNodes_p,
                                         Mask ActiveNodes);
/**************************************************************/

/**
 * @brief
 *
 * @param L0
 * @param dN_alpha_u_n1 Displacement shape function gradient evaluated at node
 * alpha
 * @param kichhoff_stress Kirchhoff pore water pressure
 * @param V0 Volume of the particle
 */
void compute_L0__upw__(double *L0, const double *dN_alpha_u_n1,
                       const double *kichhoff_stress, double V0);
/**************************************************************/

/**
 * @brief
 *
 * @param L1
 * @param dN_alpha_u_n1 Displacement shape function gradient evaluated at node
 * alpha
 * @param kichhoff_pressure Kirchhoff pore water pressure
 * @param V0 Volume of the particle
 */
void compute_L1__upw__(double *L1, const double *dN_alpha_u_n1,
                       double kichhoff_pressure, double V0);
/**************************************************************/

/**
 * @brief
 *
 * @param DL1_Du
 * @param dN_alpha_u_n1 Displacement shape function gradient evaluated at node
 * alpha
 * @param dN_beta_u_n1 Displacement shape function gradient evaluated at node
 * beta
 * @param kichhoff_pressure Kirchhoff pore water pressure
 * @param V0 Volume of the particle
 */
void compute_DL1_Du__upw__(double *DL1_Du, const double *dN_alpha_u_n1,
                           const double *dN_beta_u_n1, double kichhoff_pressure,
                           double V0);
/**************************************************************/

/**
 * @brief
 *
 * @param DL1_Dpw
 * @param dN_alpha_u_n1 Displacement shape function gradient evaluated at node
 * alpha
 * @param N_beta_pw_n1 Pore water pressure shape function evaluated at node beta
 * @param V0 Volume of the particle
 */
void compute_DL1_Dpw__upw__(double *DL1_Dpw, const double *dN_alpha_u_n1,
                            double N_beta_pw_n1, double V0);
/**************************************************************/

/**
 * @brief
 *
 * @param L2
 * @param N_alpha_u_n1 Displacement shape function evaluated at node alpha
 * @param b Particle acceleration
 * @param m Mass
 */
void compute_L2__upw__(double *L2, double N_alpha_u_n1, const double *b,
                       double m);
/**************************************************************/

/**
 * @brief Construct a new compute DL2 Du  upw object
 *
 * @param DL2_Du
 * @param N_alpha_u_n1 Displacement shape function evaluated at node alpha
 * @param N_beta_u_n1 Displacement shape function evaluated at node beta
 * @param dN_beta_u_n1 Displacement shape function gradient evaluated at node
 * beta
 * @param b Particle acceleration
 * @param kichhoff_pressure Kirchhoff pore water pressure
 * @param phi_f Volume fraction of the fluid
 * @param intrinsic_rho_f Intrinsic fluid density
 * @param kappa_f Fluid compressibility
 * @param Jacobian Deformatation gradiend determinant
 * @param m Mass
 * @param alpha_1 Newmark-beta time integration parameter
 * @param V0 Volume of the particle
 */
compute_DL2_Du__upw__(double *DL2_Du, double N_alpha_u_n1, double N_beta_u_n1,
                      const double *dN_beta_u_n1, const double *b,
                      double kichhoff_pressure, double phi_f,
                      double intrinsic_rho_f, double kappa_f, double Jacobian,
                      double m, double alpha_1, double V0);
/**************************************************************/

/**
 * @brief
 *
 * @param DL2_Dpw
 * @param N_alpha_u_n1 Displacement shape function evaluated at node alpha
 * @param N_beta_pw_n1 Pore water pressure shape function evaluated at node beta
 * @param b Particle acceleration
 * @param phi_f Volume fraction of the fluid
 * @param intrinsic_rho_f Intrinsic fluid density
 * @param kappa_f Fluid compressibility
 * @param V0 Volume of the particle
 */
void compute_DL2_Dpw__upw__(double *DL2_Dpw, double N_alpha_u_n1,
                            double N_beta_pw_n1, const double *b, double phi_f,
                            double intrinsic_rho_f, double kappa_f, double V0);
/**************************************************************/

/**
 * @brief
 *
 * @param G0
 * @param N_alpha_pw_n1 Pore water pressure shape function evaluated at node
 * alpha
 * @param intrinsic_rho_f Intrinsic fluid density
 * @param relative_rho_f_p Relative fluid density
 * @param rate_kichhoff_pressure Rate of the Kirchhoff pore water pressure
 * @param rate_Jacobian Rate of the deformatation gradiend determinant
 * @param kappa_f Fluid compressibility
 * @param V0 Volume of the particle
 */
void compute_G0__upw__(double *G0, double N_alpha_pw_n1, double intrinsic_rho_f,
                       double relative_rho_f_p, double rate_kichhoff_pressure,
                       double rate_Jacobian, double kappa_f, double V0);
/**************************************************************/

/**
 * @brief
 *
 * @param DG0_Du
 * @param N_alpha_pw_n1 Pore water pressure shape function evaluated at node
 * alpha
 * @param dN_beta_u_n1 Displacement shape function gradient evaluated at node
 * beta
 * @param grad_v Gradient of the velocity
 * @param kichhoff_pressure Kirchhoff pore water pressure
 * @param rate_kichhoff_pressure Rate of the Kirchhoff pore water pressure
 * @param Jacobian Deformatation gradiend determinant
 * @param rate_Jacobian Rate of the deformatation gradiend determinant
 * @param phi_s Volume fraction of the solid
 * @param phi_f Volume fraction of the fluid
 * @param intrinsic_rho_f Intrinsic fluid density
 * @param relative_rho_f_p Relative fluid density
 * @param kappa_f  Fluid compressibility
 * @param alpha_4 Newmark-beta time integration parameter
 * @param V0 Volume of the particle
 */
void compute_DG0_Du__upw__(double *DG0_Du, double N_alpha_pw_n1,
                           const double *dN_beta_u_n1, const double *grad_v,
                           double kichhoff_pressure,
                           double rate_kichhoff_pressure, double Jacobian,
                           double rate_Jacobian, double phi_s, double phi_f,
                           double intrinsic_rho_f, double relative_rho_f_p,
                           double kappa_f, double alpha_4, double V0);
/**************************************************************/

/**
 * @brief
 *
 * @param[out] DG0_Dpw
 * @param[in] N_alpha_pw_n1 Pore water pressure shape function evaluated at node
 * alpha
 * @param[in] N_beta_pw_n1 Pore water pressure shape function evaluated at node
 * beta
 * @param[in] rate_kichhoff_pressure Rate of the Kirchhoff pore water pressure
 * @param[in] Jacobian Deformatation gradiend determinant
 * @param[in] rate_Jacobian Rate of the deformatation gradiend determinant
 * @param[in] phi_f Volume fraction of the fluid
 * @param[in] intrinsic_rho_f Intrinsic fluid density
 * @param[in] kappa_f Fluid compressibility
 * @param[in] alpha_4 Newmark-beta time integration parameter
 * @param[in] V0 Volume of the particle
 */
void compute_DG0_Dpw__upw__(double *DG0_Dpw, double N_alpha_pw_n1,
                            double N_beta_pw_n1, double rate_kichhoff_pressure,
                            double Jacobian, double rate_Jacobian, double phi_f,
                            double intrinsic_rho_f, double kappa_f, double alpha_4, double V0);
/**************************************************************/

/**
 * @brief
 *
 * @param G1
 * @param dN_alpha_pw_n1 Pore water pressure shape function gradient evaluated
 * at node alpha
 * @param K Permeability tensor
 * @param grad_kichhoff_pressure Gradient of the Kirchhoff pore water pressure
 * @param g Gravity
 * @param V0 Volume of the particle
 */
void compute_G1__upw__(double *G1, const double *dN_alpha_pw_n1,
                       const double *K, const double *grad_kichhoff_pressure,
                       double g, double V0);
/**************************************************************/

/**
 * @brief
 *
 * @param G2
 * @param dN_alpha_pw_n1 Pore water pressure shape function gradient evaluated
 * at node alpha
 * @param K Permeability tensor
 * @param b Particle acceleration
 * @param Jacobian Deformatation gradiend determinant
 * @param intrinsic_rho_f Intrinsic fluid density
 * @param g Gravity
 * @param V0 Volume of the particle
 */
void compute_G2__upw__(double *G2, const double *dN_alpha_pw_n1,
                       const double *K, const double *b, double Jacobian,
                       double intrinsic_rho_f, double g, double V0);

#endif