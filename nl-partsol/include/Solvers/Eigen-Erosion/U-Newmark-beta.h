
#ifndef _U_NEWMARK_BETA_EIGEN_EROSION_H_
#define _U_NEWMARK_BETA_EIGEN_EROSION_H_

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Constitutive.h"
#include "Particles.h"
#include "Nodes.h"
#include "InOutFun.h"

#ifdef USE_PETSC

#else 
  #include "Linear-Solvers/dgetrs-LAPACK.h"
#endif

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/*
  Call global variables
*/
double Thickness_Plain_Stress;
Event *Out_nodal_path_csv;
Event *Out_particles_path_csv;
int Number_Out_nodal_path_csv;
int Number_Out_particles_path_csv;


// Global variuables
unsigned InitialStep;
unsigned TimeStep;
unsigned NumTimeStep;

typedef struct {

  double * value;
  double * d_value_dt;
  double * d2_value_dt2;

} Nodal_Field;

typedef struct {
  double alpha_1;
  double alpha_2;
  double alpha_3;
  double alpha_4;
  double alpha_5;
  double alpha_6;
  double epsilon;
} Newmark_parameters;


//  Auxiliar functions

static double __compute_deltat(
  Particle MPM_Mesh /**< */, 
  double h /**< */, 
  Time_Int_Params Parameters_Solver /**< */);

/*!
  \brief Finite strains Newmark-beta
 
  \param[in] beta: First Newmark-beta parameter
  \param[in] gamma: Second Newmark-beta parameter
  \param[in] DeltaTimeStep: Timestep increment
  \param[in] epsilon: Preconditioner parameter for the mass matrix
*/
static Newmark_parameters __compute_Newmark_parameters(
  double beta, 
  double gamma,
  double DeltaTimeStep,
  double epsilon);

/*
*/
static int __compute_nodal_effective_mass(
  double * Effective_MassMatrix /**< */,
  Particle MPM_Mesh /**< */, 
  Mesh FEM_Mesh /**< */,
  Mask ActiveNodes /**< */, 
  double epsilon /**< */);

static int __get_nodal_field_tn(
  Nodal_Field U_n /**< */,
  Particle MPM_Mesh /**< */,
  Mesh FEM_Mesh /**< */, 
  Mask ActiveNodes /**< */);

static void __initialise_nodal_increments(
  Nodal_Field DU,
  Nodal_Field U_n /**< */, 
  Mesh FEM_Mesh /**< */,
  Mask ActiveNodes /**< */,
  Newmark_parameters Params /**< */);

static int __local_deformation(
  const Nodal_Field D_U /**< */, 
  Mask ActiveNodes /**< */,
  Particle MPM_Mesh /**< */, 
  Mesh FEM_Mesh /**< */,
  double TimeStep /**< */);

/*!
  \brief Function to compute is a material point is or not eroded.
  Here the notation is the same as in \cite Pandolfi_2012

  \param[in] p Index of the particle
  \param[in] Phi Particles fields
  \param[in] Properties Define the material properties of the particle
  \param[in] B_eps Define the particles close to each particle
  \param[in] DeltaX Mesh size
*/
static void Eigenerosion(
  int p, 
  Fields Phi, 
  Material MatPro, 
  ChainPtr *Beps,
  double DeltaX);

static int __Nodal_Internal_Forces(
  double * Residual /**< */, 
  double * Reactions /**< */,
  Mask ActiveNodes /**< */, 
  Mask ActiveDOFs /**< */,
  Particle MPM_Mesh /**< */, 
  Mesh FEM_Mesh /**< */,
  double TimeStep /**< */);

static void __internal_force_density(
  double * InternalForcesDensity_Ap,
  const double * P_p,
  const double * F_n_p,
  const double * gradient_pA);

static void __Nodal_Traction_Forces(
  double * Residual /**< */, 
  double * Reactions /**< */,
  Mask ActiveNodes /**< */, 
  Mask ActiveDOFs /**< */,
  Particle MPM_Mesh /**< */,
  Mesh FEM_Mesh /**< */);

static void __Nodal_Body_Forces(
  double * Residual /**< */, 
  double * Reactions /**< */,
  Mask ActiveNodes /**< */, 
  Mask ActiveDOFs /**< */,
  Particle MPM_Mesh /**< */, 
  Mesh FEM_Mesh /**< */);

static void __Nodal_Inertial_Forces(
  double * Residual /**< */,
  double * Mass /**< */,  
  Nodal_Field U_n /**< */, 
  Nodal_Field D_U /**< */,
  Mask ActiveNodes /**< */,  
  Mask ActiveDOFs /**< */,  
  Newmark_parameters Params /**< */);

static double __error_residual(
  const double * Residual /**< */,
  int Total_dof /**< */);

static void compute_local_intertia(
  double * Inertia_density_p /**< */, 
  double Na_p /**< */,
  double Nb_p /**< */,
  double m_p /**< */, 
  double alpha_1 /**< */, 
  double epsilon /**< */,
  unsigned A /**< */, 
  unsigned B /**< */);  

#ifdef USE_PETSC

#else
static int __assemble_tangent_stiffness(
  double * Tangent_Stiffness /**< */,
  Mask ActiveNodes /**< */,
  Mask ActiveDOFs /**< */,
  Particle MPM_Mesh /**< */, 
  Mesh FEM_Mesh /**< */,
  Newmark_parameters Params /**< */);
#endif

static int __solve_equilibrium(
  double * Tangent_Stiffness /**< */,
  double * Residual /**< */,
  unsigned Nactivedofs /**< */);

static void __update_Nodal_Increments(
  const double * Residual /**< */,
  Nodal_Field D_U /**< */,
  Nodal_Field U_n /**< */,
  Mask ActiveDOFs /**< */,
  Newmark_parameters Params /**< */,
  unsigned Ntotaldofs /**< */);

static void __update_Particles(
  Nodal_Field D_U /**< */,
  Particle MPM_Mesh /**< */, 
  Mesh FEM_Mesh /**< */,
  Mask ActiveNodes /**< */);

static void __output(
  Particle MPM_Mesh /**< */, 
  Mesh FEM_Mesh /**< */, 
  Mask ActiveNodes /**< */,
  Nodal_Field D_U /**< */, 
  Nodal_Field U_n /**< */, 
  double * Reactions /**< */, 
  int TimeStep /**< */,
  int ResultsTimeStep /**< */);


/*
  \brief Finite strains Newmark-beta
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
*/
int U_Newmark_beta_Eigen_Erosion(
    Mesh FEM_Mesh, 
    Particle MPM_Mesh,
    Time_Int_Params Parameters_Solver);

#endif