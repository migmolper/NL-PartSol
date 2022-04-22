
#ifndef _U_NEWMARK_BETA_FINITE_STRAINS_H_
#define _U_NEWMARK_BETA_FINITE_STRAINS_H_

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

//#ifdef USE_PETSC

#include <petscksp.h>
//  #include "Linear-Solvers/"

//#else

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__ 
#include <Accelerate/Accelerate.h>
#endif

#include "Linear-Solvers/dgetrs-LAPACK.h"

//#endif


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

//#ifdef USE_PETSC
//  Vec value;
//  Vec d_value_dt;
//  Vec d2_value_dt2;
//#else 
  double * value;
  double * d_value_dt;
  double * d2_value_dt2;
//#endif

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

static void __preallocation_tangent_matrix(
  int * nnz /**< */,
  Mask ActiveNodes /**< */,
  Mask ActiveDOFs /**< */,
  Particle MPM_Mesh /**< */);

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
static int __assemble_tangent_stiffness(
  Mat Tangent_Stiffness /**< */,
  Mask ActiveNodes /**< */,
  Mask ActiveDOFs /**< */,
  Particle MPM_Mesh /**< */, 
  Mesh FEM_Mesh /**< */,
  Newmark_parameters Params /**< */);
#else
static int __assemble_tangent_stiffness(
  double * Tangent_Stiffness /**< */,
  Mask ActiveNodes /**< */,
  Mask ActiveDOFs /**< */,
  Particle MPM_Mesh /**< */, 
  Mesh FEM_Mesh /**< */,
  Newmark_parameters Params /**< */);
#endif


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
int U_Newmark_beta_Finite_Strains(
    Mesh FEM_Mesh, 
    Particle MPM_Mesh,
    Time_Int_Params Parameters_Solver);

#endif