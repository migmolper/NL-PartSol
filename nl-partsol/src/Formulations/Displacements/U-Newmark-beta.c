#include "Formulations/Displacements/U-Newmark-beta.h"
#include "Macros.h"
#include "Types.h"
#include "petscsnes.h"
#include <petscistypes.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <petscvec.h>
#include <stdio.h>
#include <stdlib.h>

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

typedef struct {

  Mask ActiveNodes;
  Mask ActiveDOFs;
  IS Dirichlet_dofs;
  Particle MPM_Mesh;
  Mesh FEM_Mesh;
  Vec Lumped_Mass;
  Vec U_n_dt;
  Vec U_n_dt2;
  Newmark_parameters Time_Integration_Params;

} Ctx;

/**************************************************************/
static double __compute_deltat(Particle MPM_Mesh, double h,
                               Time_Int_Params Parameters_Solver);

static Newmark_parameters __compute_Newmark_parameters(double beta,
                                                       double gamma,
                                                       double DeltaTimeStep,
                                                       double epsilon);

static int *__create_sparsity_pattern(Mask ActiveNodes, Particle MPM_Mesh);

static IS __get_dirichlet_list_dofs(Mask ActiveNodes, Mesh FEM_Mesh, int Step,
                                    int NumTimeStep);

static PetscErrorCode __compute_nodal_lumped_mass(Vec Lumped_MassMatrix,
                                                  Particle MPM_Mesh,
                                                  Mesh FEM_Mesh,
                                                  Mask ActiveNodes);

static PetscErrorCode __get_nodal_field_n(Vec U_n_dt, Vec U_n_dt2,
                                          Vec Lumped_Mass, Particle MPM_Mesh,
                                          Mesh FEM_Mesh, Mask ActiveNodes,
                                          Mask ActiveDOFs,
                                          Newmark_parameters Params);

static PetscErrorCode __form_initial_guess(Vec DU, Vec U_n_dt, Vec U_n_dt2,
                                           Mesh FEM_Mesh, Mask ActiveNodes,
                                           Newmark_parameters Params);

static PetscErrorCode __lagrangian_evaluation(SNES snes, Vec D_U, Vec Residual,
                                              void *ctx);

static PetscScalar *__compute_nodal_velocity_increments(
    const PetscScalar *dU, const PetscScalar *Un_dt, const PetscScalar *Un_dt2,
    Newmark_parameters Params, PetscInt Ntotaldofs);

static PetscErrorCode __local_compatibility_conditions(const PetscScalar *dU,
                                                       const PetscScalar *dU_dt,
                                                       Mask ActiveNodes,
                                                       Particle MPM_Mesh,
                                                       Mesh FEM_Mesh);

static PetscErrorCode __constitutive_update(Particle MPM_Mesh, Mesh FEM_Mesh);

static PetscErrorCode __nodal_internal_forces(PetscScalar *Lagrangian,
                                              Mask ActiveNodes, Mask ActiveDOFs,
                                              Particle MPM_Mesh, Mesh FEM_Mesh);

static void __nodal_traction_forces(PetscScalar *Lagrangian, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh);

static void __nodal_inertial_forces(PetscScalar *Lagrangian,
                                    const PetscScalar *M_II,
                                    const PetscScalar *dU,
                                    const PetscScalar *Un_dt,
                                    const PetscScalar *Un_dt2, Mask ActiveNodes,
                                    Mask ActiveDOFs, Newmark_parameters Params);

static PetscErrorCode __jacobian_evaluation(SNES snes, Vec dU, Mat Jacobian,
                                            Mat B, void *ctx);

static PetscErrorCode __monitor(PetscInt Time, PetscInt NumTimeStep,
                                PetscInt SNES_Iter, PetscInt KSP_Iter,
                                PetscInt SNES_MaxIter, PetscScalar KSP_Norm,
                                PetscScalar SNES_Norm,
                                SNESConvergedReason converged_reason);

static PetscErrorCode
__compute_nodal_kinetic_increments(Vec dU_dt, Vec dU_dt2, Vec dU, Vec Un_dt,
                                   Vec Un_dt2, Newmark_parameters Params,
                                   unsigned Ntotaldofs);

static PetscErrorCode __update_particles_internal_variables(Particle MPM_Mesh, Mesh FEM_Mesh);

static PetscErrorCode __update_particles_kinetics_FLIP_PIC(double alpha_blend, Vec dU,  Vec U_n_dt,  Vec dU_dt, Vec dU_dt2,
                                         Particle MPM_Mesh, Mesh FEM_Mesh,
                                         Mask ActiveNodes);

static PetscErrorCode __jacobian_variational_smoothing(const double * Lumped_Mass, Particle MPM_Mesh,
                                          Mesh FEM_Mesh, Mask ActiveNodes);

/**************************************************************/

// Global variables
unsigned InitialStep;
unsigned NumTimeStep;
unsigned TimeStep;
bool Use_explicit_trial;

/**************************************************************/

PetscErrorCode U_Newmark_Beta(Mesh FEM_Mesh, Particle MPM_Mesh,
                   Time_Int_Params Parameters_Solver) {

  PetscErrorCode STATUS = EXIT_SUCCESS;

  //  Auxiliar variables for the solver
  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes;
  unsigned Ntotaldofs;

  // Time integration variables
  InitialStep = Parameters_Solver.InitialTimeStep;
  NumTimeStep = Parameters_Solver.NumTimeStep;
  TimeStep = InitialStep;
  Use_explicit_trial = Parameters_Solver.Use_explicit_trial;
  double epsilon = Parameters_Solver.epsilon_Mass_Matrix;
  double beta = Parameters_Solver.beta_Newmark_beta;
  double gamma = Parameters_Solver.gamma_Newmark_beta;
  double alpha_blend = 0.95;
  double DeltaTimeStep;
  double DeltaX = FEM_Mesh.DeltaX;

  Mat Tangent_Stiffness;
  int *sparsity_pattern;
  Vec Lumped_Mass;
  Vec Residual;
  Vec DU;
  Vec U_n_dt, dU_dt;
  Vec U_n_dt2, dU_dt2;
  IS Dirichlet_dofs;
  Mask ActiveNodes;
  Mask ActiveDOFs;

  Newmark_parameters Time_Integration_Params;

  // Define variables for the non-linear solver
  SNES snes;
  KSP ksp;
  PC pc;
  Ctx AplicationCtx;
  unsigned SNES_Max_Iter = Parameters_Solver.MaxIter;
  double Relative_TOL = Parameters_Solver.TOL_Newmark_beta;
  double Absolute_TOL = 100 * Relative_TOL;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Time step is defined at the init of the simulation throught the
    CFL condition. Notice that for this kind of solver, CFL confition is
    not required to be satisfied. The only purpose of it is to use the existing
    software interfase.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DeltaTimeStep = __compute_deltat(MPM_Mesh, DeltaX, Parameters_Solver);

  if (Driver_EigenErosion) {
    compute_Beps__Constitutive__(MPM_Mesh, FEM_Mesh, true);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute time integration parameters
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Time_Integration_Params =
      __compute_Newmark_parameters(beta, gamma, DeltaTimeStep, epsilon);

  while (TimeStep < NumTimeStep) {

    DoProgress("Simulation:", TimeStep, NumTimeStep);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Local search and compute list of active nodes and dofs
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    STATUS = local_search__MeshTools__(MPM_Mesh, FEM_Mesh);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED " Error in " RESET "" BOLDRED
                      "local_search__MeshTools__() " RESET " \n");
      return EXIT_FAILURE;
    }

    ActiveNodes = get_active_nodes__MeshTools__(FEM_Mesh);
    Nactivenodes = ActiveNodes.Nactivenodes;
    Ntotaldofs = Ndim * Nactivenodes;
    ActiveDOFs = get_active_dofs__MeshTools__(ActiveNodes, FEM_Mesh, TimeStep,
                                              NumTimeStep);
    sparsity_pattern = __create_sparsity_pattern(ActiveNodes, MPM_Mesh);

    if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
      compute_Beps__Constitutive__(MPM_Mesh, FEM_Mesh, false);
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Define and allocate the effective mass matrix
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    VecCreate(PETSC_COMM_WORLD, &Lumped_Mass);
    VecSetSizes(Lumped_Mass, PETSC_DECIDE, Ntotaldofs);
    VecSetFromOptions(Lumped_Mass);

    PetscCall(__compute_nodal_lumped_mass(Lumped_Mass, MPM_Mesh, FEM_Mesh,
                                          ActiveNodes));

    VecAssemblyBegin(Lumped_Mass);
    VecAssemblyEnd(Lumped_Mass);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Get the previous converged nodal value
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(VecCreate(PETSC_COMM_WORLD, &U_n_dt));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &U_n_dt2));
    PetscCall(VecSetSizes(U_n_dt, PETSC_DECIDE, Ntotaldofs));
    PetscCall(VecSetSizes(U_n_dt2, PETSC_DECIDE, Ntotaldofs));
    PetscCall(VecSetFromOptions(U_n_dt));
    PetscCall(VecSetFromOptions(U_n_dt2));

    PetscCall(__get_nodal_field_n(U_n_dt, U_n_dt2, Lumped_Mass, MPM_Mesh,
                                  FEM_Mesh, ActiveNodes, ActiveDOFs,
                                  Time_Integration_Params));
    PetscCall(VecAssemblyBegin(U_n_dt));
    PetscCall(VecAssemblyEnd(U_n_dt));
    PetscCall(VecAssemblyBegin(U_n_dt2));
    PetscCall(VecAssemblyEnd(U_n_dt2));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create structure to store the dirchlet boudary conditions
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Dirichlet_dofs =
        __get_dirichlet_list_dofs(ActiveNodes, FEM_Mesh, TimeStep, NumTimeStep);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Define user parameters for the SNES context
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    AplicationCtx.ActiveNodes = ActiveNodes;
    AplicationCtx.ActiveDOFs = ActiveDOFs;
    AplicationCtx.Dirichlet_dofs = Dirichlet_dofs;
    AplicationCtx.MPM_Mesh = MPM_Mesh;
    AplicationCtx.FEM_Mesh = FEM_Mesh;
    AplicationCtx.Lumped_Mass = Lumped_Mass;
    AplicationCtx.U_n_dt = U_n_dt;
    AplicationCtx.U_n_dt2 = U_n_dt2;
    AplicationCtx.Time_Integration_Params = Time_Integration_Params;

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create nonlinear solver context
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
    PetscCall(SNESSetType(snes, SNESNEWTONLS));
    PetscCall(SNESSetOptionsPrefix(snes, "SolidLagragian_"));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create matrix and vector data structures; set corresponding routines
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(VecCreate(PETSC_COMM_WORLD, &DU));
    PetscCall(VecSetSizes(DU, PETSC_DECIDE, Ntotaldofs));
    PetscCall(VecSetFromOptions(DU));
    PetscCall(VecSetOption(DU, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE));
    PetscCall(VecDuplicate(DU, &Residual));
    PetscCall(VecSetFromOptions(Residual));
    PetscCall(VecSetOption(Residual, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create Jacobian matrix data structure
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(MatCreateSeqAIJ(PETSC_COMM_SELF, Ntotaldofs, Ntotaldofs, 0,
                              sparsity_pattern, &Tangent_Stiffness));
    PetscCall(
        MatSetOption(Tangent_Stiffness, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
    PetscCall(MatSetFromOptions(Tangent_Stiffness));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set function evaluation routine and vector.
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(SNESSetFunction(snes, Residual, __lagrangian_evaluation,
                              &AplicationCtx));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set Jacobian matrix data structure and Jacobian evaluation routine
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(SNESSetJacobian(snes, Tangent_Stiffness, Tangent_Stiffness,
                              __jacobian_evaluation, &AplicationCtx));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Customize nonlinear solver; set runtime options :
       Set linear solver defaults for this problem. By extracting the
       KSP and PC contexts from the SNES context, we can then
       directly call any KSP and PC routines to set various options.
       Optionally allow user-provided preconditioner
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Petsc_Direct_solver = false;
    Petsc_Iterative_solver = true;
    if(Petsc_Direct_solver)
    {
      PetscCall(SNESGetKSP(snes, &ksp));
      PetscCall(KSPGetPC(ksp, &pc));
      PetscCall(PCSetType(pc, PCCHOLESKY));
      PCFactorSetMatSolverType(pc,MATSOLVERCHOLMOD);
      PetscCall(SNESSetTolerances(snes, Absolute_TOL, Relative_TOL, PETSC_DEFAULT,
                                SNES_Max_Iter, PETSC_DEFAULT));
      PetscCall(KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT,
                               PETSC_DEFAULT));
    }
    else if(Petsc_Iterative_solver)
    {
      PetscCall(SNESGetKSP(snes, &ksp));
      PetscCall(KSPGetPC(ksp, &pc));
      PetscCall(PCSetType(pc, PCJACOBI));
      PetscCall(SNESSetTolerances(snes, Absolute_TOL, Relative_TOL, PETSC_DEFAULT,
                                SNES_Max_Iter, PETSC_DEFAULT));
      PetscCall(KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT,
                               PETSC_DEFAULT));
    }

   

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Set SNES/KSP/KSP/PC runtime options, e.g.,
          -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(SNESSetLagJacobian(snes, 1));
    PetscCall(SNESSetFromOptions(snes));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Evaluate initial guess; then solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    if (Use_explicit_trial == true) {
      PetscCall(__form_initial_guess(DU, U_n_dt, U_n_dt2, FEM_Mesh, ActiveNodes,
                                     Time_Integration_Params));
    } else {
      PetscCall(VecZeroEntries(DU));
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Run solver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(SNESSolve(snes, PETSC_NULL, DU));

    if (Flag_Print_Convergence) {
      SNESConvergedReason converged_reason;
      PetscInt SNES_Iter, KSP_Iter;
      PetscScalar KSP_Norm, SNES_Norm;
      VecNorm(Residual, NORM_2, &SNES_Norm);
      KSPGetResidualNorm(ksp, &KSP_Norm);
      SNESGetConvergedReason(snes, &converged_reason);
      SNESGetIterationNumber(snes, &SNES_Iter);
      SNESGetLinearSolveIterations(snes, &KSP_Iter);

      __monitor(TimeStep, NumTimeStep, SNES_Iter, KSP_Iter, SNES_Max_Iter,
                KSP_Norm, SNES_Norm, converged_reason);
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Free work space.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(MatDestroy(&Tangent_Stiffness));
    PetscCall(VecDestroy(&Residual));
    PetscCall(SNESDestroy(&snes));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Compute kinetic nodal increments
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(VecCreate(PETSC_COMM_WORLD, &dU_dt));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &dU_dt2));
    PetscCall(VecSetSizes(dU_dt, PETSC_DECIDE, Ntotaldofs));
    PetscCall(VecSetSizes(dU_dt2, PETSC_DECIDE, Ntotaldofs));
    PetscCall(VecSetFromOptions(dU_dt));
    PetscCall(VecSetFromOptions(dU_dt2));

    PetscCall(__compute_nodal_kinetic_increments(
        dU_dt, dU_dt2, DU, U_n_dt, U_n_dt2, Time_Integration_Params,
        Ntotaldofs));

    PetscCall(VecAssemblyBegin(dU_dt));
    PetscCall(VecAssemblyEnd(dU_dt));
    PetscCall(VecAssemblyBegin(dU_dt2));
    PetscCall(VecAssemblyEnd(dU_dt2));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Update particle information
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(__update_particles_internal_variables(MPM_Mesh, FEM_Mesh)); 

    PetscCall(
        __update_particles_kinetics_FLIP_PIC(alpha_blend, DU, U_n_dt, dU_dt, dU_dt2, MPM_Mesh, FEM_Mesh, ActiveNodes));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Outputs
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    if (TimeStep % ResultsTimeStep == 0) {
      particle_results_vtk__InOutFun__(MPM_Mesh, TimeStep, ResultsTimeStep);
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Free memory
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(VecDestroy(&Lumped_Mass));
    PetscCall(VecDestroy(&U_n_dt));
    PetscCall(VecDestroy(&U_n_dt2));
    PetscCall(VecDestroy(&DU));
    PetscCall(VecDestroy(&dU_dt));
    PetscCall(VecDestroy(&dU_dt2));
    PetscCall(ISDestroy(&Dirichlet_dofs));
    free(ActiveNodes.Nodes2Mask);
    free(ActiveDOFs.Nodes2Mask);
    free(sparsity_pattern);

    //! Update time step
    TimeStep++;
  }

  return EXIT_SUCCESS;
}

/*********************************************************************/

static double __compute_deltat(Particle MPM_Mesh, double h,
                               Time_Int_Params Parameters_Solver) {

  double DeltaT;
  double CEL_MAX = 0;
  double C[3] = {0, 0, 0};
  int Ndim = NumberDimensions;
  bool DynamicTimeStep = false;

  /*
    Read paramter solver
  */
  double CFL = Parameters_Solver.CFL;
  double CEL_MAT = Parameters_Solver.Cel;

  /*
    Consider the velocity of the particles for the courant. In
    some cases, for instance Fr>1 is important.
   */
  if (DynamicTimeStep) {
    /*
      Get the maximum wave speed in any direction
    */
    for (int i = 0; i < MPM_Mesh.NumGP; i++) {
      for (int j = 0; j < Ndim; j++) {
        C[j] = DMAX(C[j], CEL_MAT + fabs(MPM_Mesh.Phi.vel.nM[i][j]));
      }
    }

    /*
      Get the absolute maximum value of the celerity
    */
    for (int j = 0; j < Ndim; j++) {
      CEL_MAX = DMAX(CEL_MAX, C[j]);
    }

  } else {
    CEL_MAX = CEL_MAT;
  }

  /*
    Get the minimum value of the time step
  */
  DeltaT = CFL * h / CEL_MAX;

  /*
    Return new time step
  */
  return DeltaT;
}

/**************************************************************/

/*!
  \brief Finite strains Newmark-beta

  \param[in] beta: First Newmark-beta parameter
  \param[in] gamma: Second Newmark-beta parameter
  \param[in] DeltaTimeStep: Timestep increment
  \param[in] epsilon: Preconditioner parameter for the mass matrix
*/
static Newmark_parameters __compute_Newmark_parameters(double beta,
                                                       double gamma,
                                                       double DeltaTimeStep,
                                                       double epsilon) {

  Newmark_parameters Params;

  Params.alpha_1 = 1 / (beta * DSQR(DeltaTimeStep));
  Params.alpha_2 = 1 / (beta * DeltaTimeStep);
  Params.alpha_3 = (1 - 2 * beta) / (2 * beta);
  Params.alpha_4 = gamma / (beta * DeltaTimeStep);
  Params.alpha_5 = 1 - gamma / beta;
  Params.alpha_6 = (1 - gamma / (2 * beta)) * DeltaTimeStep;
  Params.epsilon = epsilon;
  Params.DeltaTimeStep = DeltaTimeStep;

  return Params;
}

/**************************************************************/

/**
 * @brief This function returns the lumped mass matrix in order to reduce
 * storage, this matrix is presented as a vector.
 *
 * @param Lumped_MassMatrix
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @param ActiveNodes List of nodes which takes place in the computation
 * @return PetscErrorCode, returns failure or success
 */
static PetscErrorCode __compute_nodal_lumped_mass(Vec Lumped_MassMatrix,
                                                  Particle MPM_Mesh,
                                                  Mesh FEM_Mesh,
                                                  Mask ActiveNodes) {

  PetscErrorCode STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumberNodes_p;
  unsigned p;

#if NumberDimensions == 2
  double Local_Mass_Matrix_p[2];
  int Mask_dofs_A[2];
#else
  double Local_Mass_Matrix_p[3];
  int Mask_dofs_A[3];
#endif

#pragma omp parallel private(NumberNodes_p, Local_Mass_Matrix_p, Mask_dofs_A)
  {

#pragma omp for private(p)
    for (p = 0; p < Np; p++) {

      //!  Define tributary nodes of the particle and shape function
      NumberNodes_p = MPM_Mesh.NumberNodes[p];
      Element Nodes_p =
          nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumberNodes_p);
      Matrix ShapeFunction_p =
          compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      //!  Get the mass of the particle
      double m_p = MPM_Mesh.Phi.mass.nV[p];

      for (unsigned A = 0; A < NumberNodes_p; A++) {

        /*
           Get the node in the mass matrix with the mask
        */
        int Ap = Nodes_p.Connectivity[A];
        int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];
        double Na_p = ShapeFunction_p.nV[A];
        double M_AB_p = Na_p * m_p;

        //! Compute the contribution of particle p to the lumped mass matrix
        //! correspoding to node A. Compute a mask with the position of the dofs
        //! with contributions to the lumped mass matrix for node A and particle
        //! p.
        for (unsigned i = 0; i < Ndim; i++) {
          Local_Mass_Matrix_p[i] = M_AB_p;
          Mask_dofs_A[i] = Mask_node_A * Ndim + i;
        }

#pragma omp critical
        {
          VecSetValues(Lumped_MassMatrix, Ndim, Mask_dofs_A,
                       Local_Mass_Matrix_p, ADD_VALUES);
        } // #pragma omp critical
      }   // for A

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);

    } // for p
  }   // #pragma omp parallel

  return STATUS;
}

/**************************************************************/

/**
 * @brief Project the kinematic information of the particles towards the nodes
 * \n of the mesh to get the nodal fields at t = n
 *
 * @param U_n_dt Nodal velocity field t = n
 * @param U_n_dt2 Nodal acceleration field t = n
 * @param Lumped_Mass Lumped mass matrix used for the projection
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param ActiveDOFs List of dofs which takes place in the computation
 * @param Params Time integration parameters
 * @return Returns failure or success
 */
static PetscErrorCode __get_nodal_field_n(Vec U_n_dt, Vec U_n_dt2,
                                          Vec Lumped_Mass, Particle MPM_Mesh,
                                          Mesh FEM_Mesh, Mask ActiveNodes,
                                          Mask ActiveDOFs,
                                          Newmark_parameters Params) {
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumberNodes_p;
  unsigned p;

#if NumberDimensions == 2
  int Mask_active_dofs_A[2];
  double V_N_m_IP[2];
  double A_N_m_IP[2];
#else
  int Mask_active_dofs_A[3];
  double V_N_m_IP[3];
  double A_N_m_IP[3];
#endif

  /*
    Use this option to impose dirichlet boundary conditions
  */
  VecSetOption(U_n_dt, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetOption(U_n_dt2, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

#pragma omp parallel private(NumberNodes_p)
  {

#pragma omp for private(p, V_N_m_IP, A_N_m_IP, Mask_active_dofs_A)
    for (p = 0; p < Np; p++) {

      /* Define element of the particle */
      NumberNodes_p = MPM_Mesh.NumberNodes[p];
      Element Nodes_p =
          nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumberNodes_p);
      Matrix ShapeFunction_p =
          compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      //! Get the mass of the GP
      double m_p = MPM_Mesh.Phi.mass.nV[p];
      const double *vel_p = &MPM_Mesh.Phi.vel.nV[p * Ndim];
      const double *acc_p = &MPM_Mesh.Phi.acc.nV[p * Ndim];

      for (unsigned A = 0; A < NumberNodes_p; A++) {

        //  Get the node in the nodal momentum with the mask
        unsigned Ap = Nodes_p.Connectivity[A];
        int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];
        double ShapeFunction_pA = ShapeFunction_p.nV[A];
        double m__x__ShapeFunction_pA = m_p * ShapeFunction_pA;

        /*
        Project particle velocity and acceleration to the nodes
        using a variational recovery approach
        and get the assembling locations nodal kinetics
        */
        for (unsigned i = 0; i < Ndim; i++) {
          V_N_m_IP[i] = m__x__ShapeFunction_pA * vel_p[i];
          A_N_m_IP[i] = m__x__ShapeFunction_pA * acc_p[i];
          Mask_active_dofs_A[i] =
              ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + i] == -1
                  ? -1
                  : Mask_node_A * Ndim + i;
        }

#pragma omp critical
        {
          VecSetValues(U_n_dt, Ndim, Mask_active_dofs_A, V_N_m_IP, ADD_VALUES);
          VecSetValues(U_n_dt2, Ndim, Mask_active_dofs_A, A_N_m_IP, ADD_VALUES);
        } // #pragma omp critical
      }   // for A

      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);

    } // for p

  } // #pragma omp parallel

  PetscCall(VecPointwiseDivide(U_n_dt, U_n_dt, Lumped_Mass));
  PetscCall(VecPointwiseDivide(U_n_dt2, U_n_dt2, Lumped_Mass));

  //! Apply boundary condition
  unsigned NumBounds = FEM_Mesh.Bounds.NumBounds;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

  if (TimeStep > 0) {

    unsigned TimeStep_n = TimeStep - 1;

    for (unsigned i = 0; i < NumBounds; i++) {

      /*
        Get the number of nodes of this boundarie
      */
      unsigned NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

      for (unsigned j = 0; j < NumNodesBound; j++) {
        /*
          Get the index of the node
        */
        int Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
        int Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

        /*
          The boundary condition is not affecting any active node,
          continue interating
        */
        if (Id_BCC_mask == -1) {
          continue;
        }

        /*
          Loop over the dimensions of the boundary condition
        */
        for (unsigned k = 0; k < Ndim; k++) {

          /*
            Apply only if the direction is active (1)
          */
          if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + TimeStep] == 1) {

            int Mask_restricted_dofs_A = Id_BCC_mask * Ndim + k;

            double D_U_value_It = 0.0;
            double U_value_In1 = 0.0;
            double V_value_In = 0.0;
            double V_value_In1 = 0.0;
            double A_value_In = 0.0;
            double A_value_In1 = 0.0;

            for (unsigned t = 0; t <= TimeStep_n; t++) {

              V_value_In = V_value_In1;
              A_value_In = A_value_In1;

              D_U_value_It = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[t];
              U_value_In1 += D_U_value_It;
              V_value_In1 = alpha_4 * D_U_value_It +
                            (alpha_5 - 1) * V_value_In + alpha_6 * A_value_In;
              A_value_In1 = alpha_1 * D_U_value_It - alpha_2 * V_value_In -
                            (alpha_3 + 1) * A_value_In;
            }

            VecSetValues(U_n_dt, 1, &Mask_restricted_dofs_A, &V_value_In1,
                         ADD_VALUES);
            VecSetValues(U_n_dt2, 1, &Mask_restricted_dofs_A, &A_value_In1,
                         ADD_VALUES);
          }
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static IS __get_dirichlet_list_dofs(Mask ActiveNodes, Mesh FEM_Mesh, int Step,
                                    int NumTimeStep) {

  /*
    Define auxilar variables
  */
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Nnodes_mask * Ndim;
  int Order_dirichlet = 0;
  int Number_of_BCC = FEM_Mesh.Bounds.NumBounds;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound;   /* Number of dimensions */
  int Id_BCC;        /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

  /*
    Generate mask for the static condensation.
  */
  PetscInt *List_active_dofs;
  PetscCalloc1(Order, &List_active_dofs);

  /*
    Loop over the the boundaries to find the constrained dofs
  */
  for (int i = 0; i < Number_of_BCC; i++) {

    /*
      Get the number of nodes of this boundary
    */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

    /*
      Get the number of dimensions where the BCC it is applied
    */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;

    for (int j = 0; j < NumNodesBound; j++) {
      /*
        Get the index of the node
      */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

      /*
        If the boundary condition is under an active node
      */
      if (Id_BCC_mask != -1) {
        /*
          Loop over the dimensions of the boundary condition
        */
        for (int k = 0; k < NumDimBound; k++) {

          Id_BCC_mask_k = Id_BCC_mask * Ndim + k;
          /*
            Apply only if the direction is active
          */
          if ((FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + Step] == 1) && 
          (List_active_dofs[Id_BCC_mask_k] == 0)) {
            List_active_dofs[Id_BCC_mask_k] = 1;
            Order_dirichlet++;
          }
        }
      }
    }
  }

  /*
    Generate mask using the location of the constrained dofs
  */
  PetscInt *List_active_dirchlet_dofs;
  PetscInt aux_idx = 0;

  PetscMalloc1(Order_dirichlet, &List_active_dirchlet_dofs);

  for (int A_i = 0; A_i < Order; A_i++) {
    if (List_active_dofs[A_i] == 1) {
      List_active_dirchlet_dofs[aux_idx] = A_i;
      aux_idx++;
    }
  }

  /*
    Output
  */
  IS Dirichlet_dofs;

  ISCreateGeneral(PETSC_COMM_WORLD, Order_dirichlet, List_active_dirchlet_dofs,
                  PETSC_USE_POINTER, &Dirichlet_dofs);

  /*
    Free auxiliar pointers
  */
  PetscFree(List_active_dofs);

  return Dirichlet_dofs;
}

/**************************************************************/

static PetscErrorCode __form_initial_guess(Vec DU, Vec U_n_dt, Vec U_n_dt2,
                                           Mesh FEM_Mesh, Mask ActiveNodes,
                                           Newmark_parameters Params) {

  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  PetscInt Ntotaldofs = Ndim * Nactivenodes;
  PetscInt Dof_Ai;
  PetscScalar DeltaTimeStep = Params.DeltaTimeStep;

  PetscScalar *DU_ptr;
  PetscCall(VecGetArray(DU, &DU_ptr));
  const PetscScalar *Un_dt_ptr;
  PetscCall(VecGetArrayRead(U_n_dt, &Un_dt_ptr));
  const PetscScalar *Un_dt2_ptr;
  PetscCall(VecGetArrayRead(U_n_dt2, &Un_dt2_ptr));

//! Compute trial
#pragma omp for private(Dof_Ai)
  for (Dof_Ai = 0; Dof_Ai < Ntotaldofs; Dof_Ai++) {

    DU_ptr[Dof_Ai] = DeltaTimeStep * Un_dt_ptr[Dof_Ai] +
                     0.5 * DSQR(DeltaTimeStep) * Un_dt2_ptr[Dof_Ai];

  } // #pragma omp for private (Dof_Ai)

  //! Apply boundary condition
  unsigned NumBounds = FEM_Mesh.Bounds.NumBounds;

  for (unsigned i = 0; i < NumBounds; i++) {

    /*
      Get the number of nodes of this boundarie
    */
    unsigned NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

    for (unsigned j = 0; j < NumNodesBound; j++) {
      /*
        Get the index of the node
      */
      int Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      int Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

      /*
        The boundary condition is not affecting any active node,
        continue interating
      */
      if (Id_BCC_mask == -1) {
        continue;
      }

      /*
        Loop over the dimensions of the boundary condition
      */
      for (unsigned k = 0; k < Ndim; k++) {

        /*
          Apply only if the direction is active (1)
        */
        if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + TimeStep] == 1) {

          Dof_Ai = Id_BCC_mask * Ndim + k;

          DU_ptr[Dof_Ai] = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep];
        }
      }
    }
  }

  PetscCall(VecRestoreArray(DU, &DU_ptr));
  PetscCall(VecRestoreArrayRead(U_n_dt, &Un_dt_ptr));
  PetscCall(VecRestoreArrayRead(U_n_dt2, &Un_dt2_ptr));

  return EXIT_SUCCESS;
}

/**************************************************************/

/**
 * @brief Evaluates nonlinear function: L(D_U)
 *
 * @param snes
 * @param D_U Increment of the nodal displacement
 * @param Residual Lagrangian in expressed in residual shape
 * @param ctx User-defined context for the Lagrangian evaluation
 * @return PetscErrorCode
 */
static PetscErrorCode __lagrangian_evaluation(SNES snes, Vec dU, Vec Lagrangian,
                                              void *ctx) {

  PetscErrorCode STATUS = EXIT_SUCCESS;

  /**
   * Read variables from user-defined structure
   * ctx
   */
  Mask ActiveNodes = ((Ctx *)ctx)->ActiveNodes;
  Mask ActiveDOFs = ((Ctx *)ctx)->ActiveDOFs;
  Particle MPM_Mesh = ((Ctx *)ctx)->MPM_Mesh;
  Mesh FEM_Mesh = ((Ctx *)ctx)->FEM_Mesh;
  Vec Lumped_Mass = ((Ctx *)ctx)->Lumped_Mass;
  Vec Un_dt = ((Ctx *)ctx)->U_n_dt;
  Vec Un_dt2 = ((Ctx *)ctx)->U_n_dt2;

  Newmark_parameters Time_Integration_Params =
      ((Ctx *)ctx)->Time_Integration_Params;

  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;

  PetscScalar *Lagrangian_ptr;
  const PetscScalar *Lumped_Mass_ptr;
  const PetscScalar *dU_ptr;
  const PetscScalar *Un_dt_ptr;
  const PetscScalar *Un_dt2_ptr;

  /*
    Initialize the lagrangian for a new evaluation
  */
  PetscCall(VecZeroEntries(Lagrangian));

  /*
   Get pointers to vector data.
      - For default PETSc vectors, VecGetArray() returns a pointer to
        the data array.  Otherwise, the routine is implementation dependent.
      - You MUST call VecRestoreArray() when you no longer need access to
        the array.
   */
  PetscCall(VecGetArray(Lagrangian, &Lagrangian_ptr));
  PetscCall(VecGetArrayRead(Lumped_Mass, &Lumped_Mass_ptr));
  PetscCall(VecGetArrayRead(dU, &dU_ptr));
  PetscCall(VecGetArrayRead(Un_dt, &Un_dt_ptr));
  PetscCall(VecGetArrayRead(Un_dt2, &Un_dt2_ptr));

  PetscScalar *dU_dt_ptr = __compute_nodal_velocity_increments(
      dU_ptr, Un_dt_ptr, Un_dt2_ptr, Time_Integration_Params, Ntotaldofs);

  PetscCall(__local_compatibility_conditions(dU_ptr, dU_dt_ptr, ActiveNodes,
                                             MPM_Mesh, FEM_Mesh));

//  __jacobian_variational_smoothing(Lumped_Mass_ptr, MPM_Mesh, FEM_Mesh, ActiveNodes);

  PetscCall(__constitutive_update(MPM_Mesh, FEM_Mesh));

  PetscCall(__nodal_internal_forces(Lagrangian_ptr, ActiveNodes, ActiveDOFs,
                                    MPM_Mesh, FEM_Mesh));

  __nodal_traction_forces(Lagrangian_ptr, ActiveNodes, ActiveDOFs, MPM_Mesh,
                          FEM_Mesh);

  __nodal_inertial_forces(Lagrangian_ptr, Lumped_Mass_ptr, dU_ptr, Un_dt_ptr,
                          Un_dt2_ptr, ActiveNodes, ActiveDOFs,
                          Time_Integration_Params);

  /**
   * Restore vectors
   *
   */
  PetscCall(VecRestoreArray(Lagrangian, &Lagrangian_ptr));
  PetscCall(VecRestoreArrayRead(Lumped_Mass, &Lumped_Mass_ptr));
  PetscCall(VecRestoreArrayRead(dU, &dU_ptr));
  PetscCall(VecRestoreArrayRead(Un_dt, &Un_dt_ptr));
  PetscCall(VecRestoreArrayRead(Un_dt2, &Un_dt2_ptr));
  PetscCall(PetscFree(dU_dt_ptr));

  return STATUS;
}

/**************************************************************/

/**
 * @brief
 *
 * @param dU Incremental nodal displacement field
 * @param dU_dt Incremental nodal velocity field
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @return PetscErrorCode
 */
static PetscErrorCode __local_compatibility_conditions(const PetscScalar *dU,
                                                       const PetscScalar *dU_dt,
                                                       Mask ActiveNodes,
                                                       Particle MPM_Mesh,
                                                       Mesh FEM_Mesh) {

  PetscErrorCode STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumberNodes_p;
  unsigned Order_p;
  unsigned p;
  int Idx_Element_p;
  int Idx_Patch_p;

/*
  Loop in the material point set
*/
#pragma omp parallel private(NumberNodes_p, Order_p)
  {

#pragma omp for private(p)
    for (p = 0; p < Np; p++) {

      //  Define tributary nodes of the particle
      NumberNodes_p = MPM_Mesh.NumberNodes[p];
      Element Nodes_p =
          nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumberNodes_p);
      Order_p = NumberNodes_p * Ndim;

      //  Get the nodal increment of displacement using the mask
      double *D_Displacement_Ap = (double *)calloc(Order_p, __SIZEOF_DOUBLE__);
      double *D_Velocity_Ap = (double *)calloc(Order_p, __SIZEOF_DOUBLE__);
      if ((D_Displacement_Ap == NULL) || (D_Velocity_Ap == NULL)) {
        fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
        STATUS = EXIT_FAILURE;
      }

      get_set_field__MeshTools__(D_Displacement_Ap, dU, Nodes_p, ActiveNodes);
      get_set_field__MeshTools__(D_Velocity_Ap, dU_dt, Nodes_p, ActiveNodes);

      /*
        Evaluate the shape function gradient in the coordinates of the
        particle
      */
      Matrix gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
        Take the values of the deformation gradient from the previous step
      */
      double *F_n_p = MPM_Mesh.Phi.F_n.nM[p];
      double *F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
      double *DF_p = MPM_Mesh.Phi.DF.nM[p];
      double *dFdt_n_p = MPM_Mesh.Phi.dt_F_n.nM[p];
      double *dFdt_n1_p = MPM_Mesh.Phi.dt_F_n1.nM[p];
      double *dt_DF_p = MPM_Mesh.Phi.dt_DF.nM[p];

      update_increment_Deformation_Gradient__Particles__(
          DF_p, D_Displacement_Ap, gradient_p.nV, NumberNodes_p);

      update_rate_increment_Deformation_Gradient__Particles__(
          dt_DF_p, D_Velocity_Ap, gradient_p.nV, NumberNodes_p);

      /*
        Update the deformation gradient in t = n + 1 with the information
        from t = n and the increment of deformation gradient.
      */
      update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, DF_p);
      update_rate_Deformation_Gradient_n1__Particles__(dFdt_n1_p, dt_DF_p,
                                                       F_n_p, DF_p, dFdt_n_p);

      //  Compute Jacobian of the deformation gradient
      MPM_Mesh.Phi.J_n1.nV[p] = I3__TensorLib__(F_n1_p);
      if (MPM_Mesh.Phi.J_n1.nV[p] <= 0.0) {
        fprintf(stderr,
                "" RED "Negative jacobian in particle %i: %e" RESET " \n", p,
                MPM_Mesh.Phi.J_n1.nV[p]);
        MPM_Mesh.Phi.J_n1.nV[p] = 0.0;
      }
      

      // F-bar  Update patch
      if (FEM_Mesh.Locking_Control_Fbar) {
        Idx_Element_p = MPM_Mesh.Element_p[p];
        Idx_Patch_p = FEM_Mesh.Idx_Patch[Idx_Element_p];
        FEM_Mesh.Vol_Patch_n[Idx_Patch_p] +=
            MPM_Mesh.Phi.J_n.nV[p] * MPM_Mesh.Phi.Vol_0.nV[p];
        FEM_Mesh.Vol_Patch_n1[Idx_Patch_p] +=
            MPM_Mesh.Phi.J_n1.nV[p] * MPM_Mesh.Phi.Vol_0.nV[p];
      }

      free(D_Displacement_Ap);
      free(D_Velocity_Ap);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);
    }

// F-bar
#pragma omp barrier
    if (FEM_Mesh.Locking_Control_Fbar) {

      double Vn_patch;
      double Vn1_patch;
      double J_patch;

#pragma omp for private(p, Vn_patch, Vn1_patch, J_patch, Idx_Element_p,        \
                        Idx_Patch_p)
      for (p = 0; p < Np; p++) {

        Idx_Element_p = MPM_Mesh.Element_p[p];
        Idx_Patch_p = FEM_Mesh.Idx_Patch[Idx_Element_p];

        Vn_patch = FEM_Mesh.Vol_Patch_n[Idx_Patch_p];
        Vn1_patch = FEM_Mesh.Vol_Patch_n1[Idx_Patch_p];
        J_patch = Vn1_patch / Vn_patch;

        STATUS = get_locking_free_Deformation_Gradient_n1__Particles__(
            p, J_patch, MPM_Mesh);
        if (STATUS == EXIT_FAILURE) {
          fprintf(
              stderr,
              "" RED "Error in "
              "get_locking_free_Deformation_Gradient_n1__Particles__()" RESET
              " \n");
          STATUS = EXIT_FAILURE;
        }

        MPM_Mesh.Phi.Jbar.nV[p] *= J_patch;
      }
    }
  }

  return STATUS;
}

/**************************************************************/

/**
 * @brief Update the stress tensor of the particle
 *
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @return PetscErrorCode
 */
static PetscErrorCode __constitutive_update(Particle MPM_Mesh, Mesh FEM_Mesh) {
  PetscErrorCode STATUS = EXIT_SUCCESS;
  int STATUS_p = EXIT_SUCCESS;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned MatIndx_p;
  unsigned p;

#pragma omp for private(p, MatIndx_p, STATUS_p)
  for (p = 0; p < Np; p++) {

    //! If the particle is failed, skip
    if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
      if (MPM_Mesh.Phi.Damage_n[p] == 1.0) {
        MPM_Mesh.Phi.W[p] = 0.0;
        continue;
      }
    }

    //  Update the Kirchhoff stress tensor with an apropiate
    //  integration rule.
    MatIndx_p = MPM_Mesh.MatIdx[p];
    Material MatProp_p = MPM_Mesh.Mat[MatIndx_p];

    STATUS_p = Stress_integration__Constitutive__(p, MPM_Mesh, MatProp_p);
    if (STATUS_p == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in Stress_integration__Constitutive__(%i,,)" RESET
              " \n",
              p);
      STATUS = STATUS_p;
    }
  }

  return STATUS;
}

/**************************************************************/

/**
 * @brief Function used to compute the contribution of the \n
 * internal forces to the Lagrangian
 *
 * @param Lagrangian Lagrangian vector
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param ActiveDOFs List of dofs which takes place in the computation
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @return PetscErrorCode
 */
static PetscErrorCode __nodal_internal_forces(PetscScalar *Lagrangian,
                                              Mask ActiveNodes, Mask ActiveDOFs,
                                              Particle MPM_Mesh,
                                              Mesh FEM_Mesh) {

  PetscErrorCode STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned p;

#if NumberDimensions == 2
  double InternalForcesDensity_Ap[2];
  int Mask_dofs_A[2];
#else
  double InternalForcesDensity_Ap[3];
  int Mask_dofs_A[3];
#endif

  double *Damage_field_n1 = MPM_Mesh.Phi.Damage_n1;
  const double *Vol_0 = MPM_Mesh.Phi.Vol_0.nV;

#pragma omp parallel private(NumNodes_p, InternalForcesDensity_Ap, Mask_dofs_A)
  {
#pragma omp for private(p)
    for (p = 0; p < Np; p++) {

      PetscErrorCode STATUS_p;

      //  Get the volume of the particle in the reference configuration
      double V0_p = Vol_0[p];

      // Get the incremental deformation gradient
      double *DF_p = MPM_Mesh.Phi.DF.nM[p];

      //  Define nodal connectivity for each particle
      //  and compute gradient of the shape function
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      Element Nodes_p =
          nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);
      Matrix d_shapefunction_n_p =
          compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      // Pushforward the shape functions
      double *d_shapefunction_n1_p = push_forward_dN__MeshTools__(
          d_shapefunction_n_p.nV, DF_p, NumNodes_p, &STATUS_p);
      if (STATUS_p == EXIT_FAILURE) {
        fprintf(stderr, "" RED " Error in " RESET "" BOLDRED
                        "push_forward_dN__MeshTools__() " RESET " \n");
        STATUS = EXIT_FAILURE;
      }

      // Get the Kirchhoff stress tensor pointer
      double *kirchhoff_p = MPM_Mesh.Phi.Stress.nM[p];

      // Compute damage parameter (eigenerosion/eigensoftening)
      if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {

        STATUS_p = compute_damage__Constitutive__(p, MPM_Mesh, FEM_Mesh.DeltaX);
        if (STATUS_p == EXIT_FAILURE) {
          fprintf(stderr, "" RED " Error in " RESET "" BOLDRED
                          "compute_damage__Constitutive__() " RESET " \n");
          STATUS = EXIT_FAILURE;
        }

#if NumberDimensions == 2
        for (unsigned i = 0; i < 5; i++) {
          kirchhoff_p[i] *= (1.0 - Damage_field_n1[p]);
        }
#else
        for (unsigned i = 0; i < 9; i++) {
          kirchhoff_p[i] *= (1.0 - Damage_field_n1[p]);
        }
#endif
      }

      for (unsigned A = 0; A < NumNodes_p; A++) {

        //! Get the gradient evaluation in node A \n
        //! and the masked index of the node A
        double *d_shapefunction_n1_pA = &d_shapefunction_n1_p[A * Ndim];
        int Ap = Nodes_p.Connectivity[A];
        int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

        //! Compute the contribution to the internal forces in the node A of the
        //! particle p (InternalForcesDensity_Ap). And, compute a mask
        //! (Mask_dofs_A) with the position of the computed values.
        for (unsigned i = 0; i < Ndim; i++) {
          InternalForcesDensity_Ap[i] = 0.0;
          for (unsigned j = 0; j < Ndim; j++) {
            InternalForcesDensity_Ap[i] +=
                kirchhoff_p[i * Ndim + j] * d_shapefunction_n1_pA[j];
          }
          Mask_dofs_A[i] = Mask_node_A * Ndim + i;
        }

        //  Asign the local internal forces (InternalForcesDensity_Ap)
        //  contribution to the node A using the volume of the particle as
        //  integration weight
#pragma omp critical
        {
          for (unsigned i = 0; i < Ndim; i++) {
            if (ActiveDOFs.Nodes2Mask[Mask_dofs_A[i]] != -1) {
              Lagrangian[Mask_dofs_A[i]] += InternalForcesDensity_Ap[i] * V0_p;
            }
          }
        } // #pragma omp critical
      }   // for unsigned A

      //   Free memory
      free__MatrixLib__(d_shapefunction_n_p);
      free(d_shapefunction_n1_p);
      free(Nodes_p.Connectivity);
    } // For unsigned p
  }   // #pragma omp parallel

  return STATUS;
}

/**************************************************************/

/**
 * @brief Function used to compute the contribution of the contact forces to the
 * residual
 *
 * @param Lagrangian Lagrangian vector
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param ActiveDOFs List of dofs which takes place in the computation
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 */
static void __nodal_traction_forces(PetscScalar *Lagrangian, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh) {

  unsigned Ndim = NumberDimensions;
  unsigned NumContactForces = MPM_Mesh.Neumann_Contours.NumBounds;
  unsigned NumNodesLoad;

  Load Load_i;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_p;      /* Nodal values of the sahpe function */
  double N_pa;
  double A0_p;

#if NumberDimensions == 2
  double T[2] = {0.0, 0.0};
#else
  double T[3] = {0.0, 0.0};
#endif

#if NumberDimensions == 2
  double LocalTractionForce_Ap[2];
  int Mask_dofs_A[2];
#else
  double LocalTractionForce_Ap[3];
  int Mask_dofs_A[3];
#endif

  unsigned p;
  int Ap, Mask_node_A;
  unsigned NumNodes_p; /* Number of nodes of each particle */

  for (unsigned cf_idx = 0; cf_idx < NumContactForces; cf_idx++) {

    /*
      Read load i
    */
    Load_i = MPM_Mesh.Neumann_Contours.BCC_i[cf_idx];

    NumNodesLoad = Load_i.NumNodes;

    for (unsigned nl_idx = 0; nl_idx < NumNodesLoad; nl_idx++) {

      /*
        Get the index of the particle
      */
      p = Load_i.Nodes[nl_idx];

      /*
        Get the area of each particle
      */
#if NumberDimensions == 2
      A0_p = MPM_Mesh.Phi.Vol_0.nV[p] / Thickness_Plain_Stress;
#else
      A0_p = MPM_Mesh.Phi.Area_0.nV[p];
#endif

      /*
        Define tributary nodes of the particle
      */
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

      /*
        Evaluate the shape function in the coordinates of the particle
      */
      N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
        Fill vector of contact forces
      */
      for (unsigned i = 0; i < Ndim; i++) {
        if (Load_i.Dir[i * NumTimeStep + TimeStep] == 1) {
          T[i] = Load_i.Value[i].Fx[TimeStep];
        }
      }

      /*
        Get the node of the mesh for the contribution
      */
      for (unsigned A = 0; A < NumNodes_p; A++) {

        //! Get the shape function evaluation in node A \n
        //! and the masked index of the node A
        N_pa = N_p.nV[A];
        Ap = Nodes_p.Connectivity[A];
        Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

        //! Compute the contribution to the local traction forces in the node A
        //! of the particle p (LocalTractionForce_Ap). And, compute a mask
        //! (Mask_dofs_A) with the position of the computed values.
        for (unsigned i = 0; i < Ndim; i++) {
          LocalTractionForce_Ap[i] = -N_pa * T[i];
          Mask_dofs_A[i] = Mask_node_A * Ndim + i;
        }

        //  Asign the nodal contact forces (LocalTractionForce_Ap) contribution
        //  to the node using the area of the particle as integration weight
#pragma omp critical
        {
          for (unsigned i = 0; i < Ndim; i++) {
            if (ActiveDOFs.Nodes2Mask[Mask_dofs_A[i]] != -1) {
              Lagrangian[Mask_dofs_A[i]] += LocalTractionForce_Ap[i] * A0_p;
            }
          }
        } // #pragma omp critical

      } // for unsigned A

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(N_p);
      free(Nodes_p.Connectivity);
    }
  }
}

/**************************************************************/

/**
 * @brief Function used to compute the contribution of the \n
 * inertial forces to the lagrangian (a - b)
 *
 * @param Lagrangian Lagrangian vector
 * @param M_II Lumped mass matrix
 * @param dU Nodal field of incremental displacements
 * @param Un_dt Nodal field of velocities at t = n
 * @param Un_dt2 Nodal field of accelerations at t = n
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param ActiveDOFs List of dofs which takes place in the computation
 * @param Params Time integration parameters
 */
static void
__nodal_inertial_forces(PetscScalar *Lagrangian, const PetscScalar *M_II,
                        const PetscScalar *dU, const PetscScalar *Un_dt,
                        const PetscScalar *Un_dt2, Mask ActiveNodes,
                        Mask ActiveDOFs, Newmark_parameters Params) {

  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;
  unsigned idx;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;

#if NumberDimensions == 2
  double b[2] = {0.0, 0.0};
#else
  double b[3] = {0.0, 0.0, 0.0};
#endif

  if (gravity_field.STATUS == true) {
    for (unsigned i = 0; i < Ndim; i++) {
      b[i] = gravity_field.Value[i].Fx[TimeStep];
    }
  }

#pragma omp for private(idx)
  for (idx = 0; idx < Ntotaldofs; idx++) {

#pragma omp critical
    {
      if (ActiveDOFs.Nodes2Mask[idx] != -1) {
        Lagrangian[idx] +=
            M_II[idx] * (alpha_1 * dU[idx] - alpha_2 * Un_dt[idx] -
                         alpha_3 * Un_dt2[idx] - b[idx % Ndim]);
      }
    }
  }
}

/**************************************************************/

/**
 * @brief Returns the sparsity pattern (non-zeros per row)
 *
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param MPM_Mesh Information of the particles
 * @return sparsity pattern
 */
static int *__create_sparsity_pattern(Mask ActiveNodes, Particle MPM_Mesh) {

  unsigned Ndim = NumberDimensions;
  unsigned Ntotaldofs = Ndim * ActiveNodes.Nactivenodes;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;

  int *Active_dof_Mat = (int *)calloc(Ntotaldofs * Ntotaldofs, __SIZEOF_INT__);

  // Spatial discretization variables
  Element Nodes_p;
  int Ap, Mask_node_A, Dof_Ai;
  int Bp, Mask_node_B, Dof_Bj;

  for (unsigned p = 0; p < Np; p++) {

    //  Define nodal connectivity for each particle
    //  and compute gradient of the shape function
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    for (unsigned A = 0; A < NumNodes_p; A++) {

      // Get the gradient evaluation in node A
      // and the masked index of the node A
      Ap = Nodes_p.Connectivity[A];
      Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

      for (unsigned B = 0; B < NumNodes_p; B++) {

        // Get the gradient evaluation in node B
        // and the masked index of the node B
        Bp = Nodes_p.Connectivity[B];
        Mask_node_B = ActiveNodes.Nodes2Mask[Bp];

        //  Assembling process
        for (unsigned i = 0; i < Ndim; i++) {

          Dof_Ai = Mask_node_A * Ndim + i;

          for (unsigned j = 0; j < Ndim; j++) {

            Dof_Bj = Mask_node_B * Ndim + j;

            Active_dof_Mat[Dof_Ai * Ntotaldofs + Dof_Bj] = 1;
          }
        }
      }
    }

    free(Nodes_p.Connectivity);
  }

  int *sparsity_pattern = (int *)calloc(Ntotaldofs, __SIZEOF_INT__);

  for (unsigned A = 0; A < Ntotaldofs; A++) {
    for (unsigned B = 0; B < Ntotaldofs; B++) {
      sparsity_pattern[A] += Active_dof_Mat[A * Ntotaldofs + B];
    }
  }

  free(Active_dof_Mat);

  return sparsity_pattern;
}

/**************************************************************/

/**
 * @brief Evaluates Jacobian matrix.
 *
 * @param snes the SNES context
 * @param DU Input vector
 * @param jac Jacobian matrix
 * @param Jacobian Optionally different preconditioning matrix
 * @param ctx User-defined context
 * @return PetscErrorCode
 */
static PetscErrorCode __jacobian_evaluation(SNES snes, Vec dU, Mat Jacobian,
                                            Mat Preconditioner, void *ctx) {

  /**
   * Read variables from user-defined structure
   * ctx
   */
  Mask ActiveNodes = ((Ctx *)ctx)->ActiveNodes;
  IS Dirichlet_dofs = ((Ctx *)ctx)->Dirichlet_dofs;
  Particle MPM_Mesh = ((Ctx *)ctx)->MPM_Mesh;
  Mesh FEM_Mesh = ((Ctx *)ctx)->FEM_Mesh;
  Vec Lumped_Mass = ((Ctx *)ctx)->Lumped_Mass;
  Newmark_parameters Time_Integration_Params =
      ((Ctx *)ctx)->Time_Integration_Params;

  PetscErrorCode STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned MatIndx_p;
  unsigned p;

#if NumberDimensions == 2
  double Jacobian_p[4];
  double Stiffness_density_p[4];
  int Mask_dofs_A[2];
  int Mask_dofs_B[2];
#else
  double Jacobian_p[9];
  double Stiffness_density_p[9];
  int Mask_dofs_A[3];
  int Mask_dofs_B[3];
#endif

  // Time integartion parameters
  double alpha_1 = Time_Integration_Params.alpha_1;
  double alpha_4 = Time_Integration_Params.alpha_4;
  double epsilon = Time_Integration_Params.epsilon;

  const PetscScalar *Lumped_Mass_ptr;
  PetscCall(VecGetArrayRead(Lumped_Mass, &Lumped_Mass_ptr));

  // Set to zero the Jacobian
  PetscCall(MatZeroEntries(Jacobian));

#pragma omp parallel private(NumNodes_p, MatIndx_p, Stiffness_density_p,       \
                             Jacobian_p, Mask_dofs_A, Mask_dofs_B)
  {
#pragma omp for private(p)
    for (p = 0; p < Np; p++) {

      //! Get mass and the volume of the particle in the reference
      //! configuration and the jacobian of the deformation gradient
      double m_p = MPM_Mesh.Phi.mass.nV[p];
      double V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

      //! Material properties of the particle
      MatIndx_p = MPM_Mesh.MatIdx[p];
      Material MatProp_p = MPM_Mesh.Mat[MatIndx_p];

      //! Pointer to the incremental deformation gradient
      double *DF_p = MPM_Mesh.Phi.DF.nM[p];

      //!  Define nodal connectivity for each particle
      //!  and compute gradient of the shape function
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      Element Nodes_p =
          nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);
      Matrix shapefunction_n_p =
          compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
      Matrix d_shapefunction_n_p =
          compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      //! Pushforward the shape function gradient
      double *d_shapefunction_n1_p = push_forward_dN__MeshTools__(
          d_shapefunction_n_p.nV, DF_p, NumNodes_p, &STATUS);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in push_forward_dN__MeshTools__()" RESET " \n");
        STATUS = EXIT_FAILURE;
      }

      for (unsigned A = 0; A < NumNodes_p; A++) {

        //! Get the gradient evaluation in node A \n
        //! and the masked index of the node A
        double shapefunction_n_pA = shapefunction_n_p.nV[A];
        double *d_shapefunction_n_pA = d_shapefunction_n_p.nM[A];
        int Ap = Nodes_p.Connectivity[A];
        int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

        for (unsigned B = 0; B < NumNodes_p; B++) {

          //! Get the gradient evaluation in node B \n
          //! and the masked index of the node B
          double shapefunction_n_pB = shapefunction_n_p.nV[B];
          double *d_shapefunction_n_pB = d_shapefunction_n_p.nM[B];
          int Bp = Nodes_p.Connectivity[B];
          int Mask_node_B = ActiveNodes.Nodes2Mask[Bp];

          //! Do the local and global assembly process for the tangent matrix
          STATUS = stiffness_density__Constitutive__(
              p, Stiffness_density_p, &d_shapefunction_n1_p[A * Ndim],
              &d_shapefunction_n1_p[B * Ndim], d_shapefunction_n_pA,
              d_shapefunction_n_pB, alpha_4, MPM_Mesh, MatProp_p);
          if (STATUS == EXIT_FAILURE) {
            fprintf(stderr,
                    "" RED "Error in stiffness_density__Constitutive__" RESET
                    "\n");
          }

          //! Damage contribution of the particle to the residual
          if ((Driver_EigenErosion == true) ||
              (Driver_EigenSoftening == true)) {
            double damage_p = MPM_Mesh.Phi.Damage_n1[p];
            for (unsigned i = 0; i < Ndim * Ndim; i++) {
              Stiffness_density_p[i] *= (1.0 - damage_p);
            }
          }

          /*
          Assemble the jacobian using the stiffness density derived from a
          component-free aproximation and simultaneasly assemble the mask with
          the dofs
          */
          for (unsigned i = 0; i < Ndim; i++) {
            for (unsigned j = 0; j < Ndim; j++) {

              /* Contribution of the internal forces to the tangent matrix */
              Jacobian_p[i * Ndim + j] =
                  Stiffness_density_p[i * Ndim + j] * V0_p;
            }

            Mask_dofs_A[i] = Mask_node_A * Ndim + i;
            Mask_dofs_B[i] = Mask_node_B * Ndim + i;
          }

#pragma omp critical
          {

            MatSetValues(Jacobian, Ndim, Mask_dofs_A, Ndim, Mask_dofs_B,
                         Jacobian_p, ADD_VALUES);

          } // #pragma omp critical
        }   // for B (node)
      }     // for A (node)

      // Free memory
      free__MatrixLib__(shapefunction_n_p);
      free__MatrixLib__(d_shapefunction_n_p);
      free(d_shapefunction_n1_p);
      free(Nodes_p.Connectivity);

    } // for p (particle)
  }   // #pragma omp parallel

  /*
    Add the contribution of the lumped mass matrix to the tangent matrix
  */
  PetscInt Dof_Ai;
  PetscInt TotalDofs = Ndim * ActiveNodes.Nactivenodes;
#pragma omp parallel
  {
#pragma omp for private(Dof_Ai)
    for (Dof_Ai = 0; Dof_Ai < TotalDofs; Dof_Ai++) {
      const PetscScalar Inertia_Ai = alpha_1 * Lumped_Mass_ptr[Dof_Ai];
      MatSetValues(Jacobian, 1, &Dof_Ai, 1, &Dof_Ai, &Inertia_Ai, ADD_VALUES);
    }
  }

  PetscCall(VecRestoreArrayRead(Lumped_Mass, &Lumped_Mass_ptr));

  /*
     Assemble matrix
  */
  PetscCall(MatAssemblyBegin(Jacobian, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Jacobian, MAT_FINAL_ASSEMBLY));

  /*
    Dirichlet boundary conditions
  */
  PetscCall(MatZeroRowsColumnsIS(Jacobian, Dirichlet_dofs, 1.0, NULL, NULL));


  return STATUS;
}

/**************************************************************/

static PetscScalar *__compute_nodal_velocity_increments(
    const PetscScalar *dU, const PetscScalar *Un_dt, const PetscScalar *Un_dt2,
    Newmark_parameters Params, PetscInt Ntotaldofs) {

  PetscInt Dof_Ai;
  PetscScalar alpha_4 = Params.alpha_4;
  PetscScalar alpha_5 = Params.alpha_5;
  PetscScalar alpha_6 = Params.alpha_6;

  PetscScalar *dU_dt;
  PetscMalloc(sizeof(PetscScalar) * Ntotaldofs, &dU_dt);

#pragma omp for private(Dof_Ai)
  for (Dof_Ai = 0; Dof_Ai < Ntotaldofs; Dof_Ai++) {

    dU_dt[Dof_Ai] = alpha_4 * dU[Dof_Ai] + (alpha_5 - 1) * Un_dt[Dof_Ai] +
                    alpha_6 * Un_dt2[Dof_Ai];

  } // #pragma omp for private (Dof_Ai)

  return dU_dt;
}

/**************************************************************/

static PetscErrorCode
__compute_nodal_kinetic_increments(Vec dU_dt, Vec dU_dt2, Vec dU, Vec Un_dt,
                                   Vec Un_dt2, Newmark_parameters Params,
                                   unsigned Ntotaldofs) {

  PetscErrorCode STATUS = EXIT_SUCCESS;

  unsigned Dof_Ai;
  PetscScalar alpha_1 = Params.alpha_1;
  PetscScalar alpha_2 = Params.alpha_2;
  PetscScalar alpha_3 = Params.alpha_3;
  PetscScalar alpha_4 = Params.alpha_4;
  PetscScalar alpha_5 = Params.alpha_5;
  PetscScalar alpha_6 = Params.alpha_6;

  PetscScalar *dU_dt_ptr;
  PetscCall(VecGetArray(dU_dt, &dU_dt_ptr));
  PetscScalar *dU_dt2_ptr;
  PetscCall(VecGetArray(dU_dt2, &dU_dt2_ptr));

  const PetscScalar *dU_ptr;
  PetscCall(VecGetArrayRead(dU, &dU_ptr));
  const PetscScalar *Un_dt_ptr;
  PetscCall(VecGetArrayRead(Un_dt, &Un_dt_ptr));
  const PetscScalar *Un_dt2_ptr;
  PetscCall(VecGetArrayRead(Un_dt2, &Un_dt2_ptr));

#pragma omp for private(Dof_Ai)
  for (Dof_Ai = 0; Dof_Ai < Ntotaldofs; Dof_Ai++) {

    dU_dt_ptr[Dof_Ai] = alpha_4 * dU_ptr[Dof_Ai] +
                        (alpha_5 - 1) * Un_dt_ptr[Dof_Ai] +
                        alpha_6 * Un_dt2_ptr[Dof_Ai];

    dU_dt2_ptr[Dof_Ai] = alpha_1 * dU_ptr[Dof_Ai] -
                         alpha_2 * Un_dt_ptr[Dof_Ai] -
                         (alpha_3 + 1) * Un_dt2_ptr[Dof_Ai];

  } // #pragma omp for private (Dof_Ai)

  PetscCall(VecRestoreArray(dU_dt, &dU_dt_ptr));
  PetscCall(VecRestoreArray(dU_dt2, &dU_dt2_ptr));
  PetscCall(VecRestoreArrayRead(dU, &dU_ptr));
  PetscCall(VecRestoreArrayRead(Un_dt, &Un_dt_ptr));
  PetscCall(VecRestoreArrayRead(Un_dt2, &Un_dt2_ptr));

  return STATUS;
}

/**************************************************************/

/**
 * @brief
 *
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @return PetscErrorCode
 */
static PetscErrorCode __update_particles_internal_variables(Particle MPM_Mesh, Mesh FEM_Mesh) {

  PetscErrorCode STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned p;

#pragma omp for private(p)
    for (p = 0; p < Np; p++) {
      
      // Update the determinant of the deformation gradient
      MPM_Mesh.Phi.J_n.nV[p] = MPM_Mesh.Phi.J_n1.nV[p];

      //! Update density
      MPM_Mesh.Phi.rho.nV[p] =
          MPM_Mesh.Phi.mass.nV[p] /
          (MPM_Mesh.Phi.Vol_0.nV[p] * MPM_Mesh.Phi.J_n.nV[p]);

      //! Update hardening
      MPM_Mesh.Phi.Kappa_n[p] = MPM_Mesh.Phi.Kappa_n1[p];

      //! Update equivalent plastic strains
      MPM_Mesh.Phi.EPS_n[p] = MPM_Mesh.Phi.EPS_n1[p];

      //! Update elastic left Cauchy-Green tensor
#if NumberDimensions == 2
      for (unsigned i = 0; i < 5; i++)
        MPM_Mesh.Phi.b_e_n.nM[p][i] = MPM_Mesh.Phi.b_e_n1.nM[p][i];
#else
      for (unsigned i = 0; i < 9; i++)
        MPM_Mesh.Phi.b_e_n.nM[p][i] = MPM_Mesh.Phi.b_e_n1.nM[p][i];
#endif

      //! Update damage variable
      if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
        MPM_Mesh.Phi.Damage_n[p] = MPM_Mesh.Phi.Damage_n1[p];
      }
      if (Driver_EigenSoftening == true) {
        MPM_Mesh.Phi.Strain_f_n[p] = MPM_Mesh.Phi.Strain_f_n1[p];
      }

      //! Update deformation gradient
#if NumberDimensions == 2
      for (unsigned i = 0; i < 5; i++)
        MPM_Mesh.Phi.F_n.nM[p][i] = MPM_Mesh.Phi.F_n1.nM[p][i];
#else
      for (unsigned i = 0; i < 9; i++)
        MPM_Mesh.Phi.F_n.nM[p][i] = MPM_Mesh.Phi.F_n1.nM[p][i];
#endif

        //! Update rate of deformation gradient
#if NumberDimensions == 2
      for (unsigned i = 0; i < 5; i++)
        MPM_Mesh.Phi.dt_F_n.nM[p][i] = MPM_Mesh.Phi.dt_F_n1.nM[p][i];
#else
      for (unsigned i = 0; i < 9; i++)
        MPM_Mesh.Phi.dt_F_n.nM[p][i] = MPM_Mesh.Phi.dt_F_n1.nM[p][i];
#endif
    } // #pragma omp for private(p)

  return STATUS;
}

/**************************************************************/

/**
 * @brief
 *
 * @param dU Incremental nodal displacement field
 * @param dU_dt Nodal velocity increment
 * @param dU_dt2 Nodal acceleration increment
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @param ActiveNodes List of nodes which takes place in the computation
 * @return PetscErrorCode
 */
static PetscErrorCode __update_particles_kinetics_FLIP_PIC(double alpha_blend,Vec dU, Vec U_n_dt, Vec dU_dt, Vec dU_dt2,
                                         Particle MPM_Mesh, Mesh FEM_Mesh,
                                         Mask ActiveNodes) {

  PetscErrorCode STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned p;

  double DU_pI;
  double V_n_pI;
  double D_V_pI;
  double D_A_pI;
  double beta_blend = 1 - alpha_blend;

  const PetscScalar *dU_ptr;
  PetscCall(VecGetArrayRead(dU, &dU_ptr));
  const PetscScalar *U_n_dt_ptr;
  PetscCall(VecGetArrayRead(U_n_dt, &U_n_dt_ptr));
  const PetscScalar *dU_dt_ptr;
  PetscCall(VecGetArrayRead(dU_dt, &dU_dt_ptr));
  const PetscScalar *dU_dt2_ptr;
  PetscCall(VecGetArrayRead(dU_dt2, &dU_dt2_ptr));

#pragma omp parallel private(NumNodes_p, DU_pI, V_n_pI, D_V_pI, D_A_pI)
  {
#pragma omp for private(p)
    for (p = 0; p < Np; p++) {
      
      //  Define nodal connectivity for each particle
      //  and compute the shape function
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      Element Nodes_p =
          nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);
      Matrix ShapeFunction_p =
          compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);


      //! First part of the contribution
      for (unsigned i = 0; i < Ndim; i++) {
        MPM_Mesh.Phi.vel.nM[p][i] = alpha_blend*MPM_Mesh.Phi.vel.nM[p][i];
      }


      //!  Update acceleration, velocity and position of the particles
      for (unsigned A = 0; A < NumNodes_p; A++) {

        // Get the shape function evaluation in node A
        // and the masked index of the node A
        double ShapeFunction_pI = ShapeFunction_p.nV[A];
        int Ap = Nodes_p.Connectivity[A];
        int A_mask = ActiveNodes.Nodes2Mask[Ap];

        for (unsigned i = 0; i < Ndim; i++) {
          DU_pI = ShapeFunction_pI * dU_ptr[A_mask * Ndim + i];
          D_V_pI = ShapeFunction_pI * dU_dt_ptr[A_mask * Ndim + i];
          V_n_pI = ShapeFunction_pI * U_n_dt_ptr[A_mask * Ndim + i];
          D_A_pI = ShapeFunction_pI * dU_dt2_ptr[A_mask * Ndim + i];

          MPM_Mesh.Phi.acc.nM[p][i] += D_A_pI;
          MPM_Mesh.Phi.vel.nM[p][i] += D_V_pI + beta_blend*V_n_pI;
          MPM_Mesh.Phi.dis.nM[p][i] += DU_pI;
          MPM_Mesh.Phi.x_GC.nM[p][i] += DU_pI;
        }
      } // for unsigned A

      free(Nodes_p.Connectivity);
      free__MatrixLib__(ShapeFunction_p);

    } // #pragma omp for private(p)
  }   // #pragma omp parallel

  PetscCall(VecRestoreArrayRead(dU, &dU_ptr));
  PetscCall(VecRestoreArrayRead(U_n_dt, &U_n_dt_ptr));
  PetscCall(VecRestoreArrayRead(dU_dt, &dU_dt_ptr));
  PetscCall(VecRestoreArrayRead(dU_dt2, &dU_dt2_ptr));

  return STATUS;
}

/**************************************************************/

static PetscErrorCode __monitor(PetscInt Time, PetscInt NumTimeStep,
                                PetscInt SNES_Iter, PetscInt KSP_Iter,
                                PetscInt SNES_MaxIter, PetscScalar KSP_Norm,
                                PetscScalar SNES_Norm,
                                SNESConvergedReason converged_reason) {

  if (NumTimeStep < 10) {
    PetscPrintf(PETSC_COMM_WORLD, "" GREEN "Step" RESET ": [%01d/%01d] | ",
                Time, NumTimeStep);
  } else if (NumTimeStep < 100) {
    PetscPrintf(PETSC_COMM_WORLD, "" GREEN "Step" RESET ": [%02d/%02d] | ",
                Time, NumTimeStep);
  } else if (NumTimeStep < 1000) {
    PetscPrintf(PETSC_COMM_WORLD, "" GREEN "Step" RESET ": [%03d/%03d] | ",
                Time, NumTimeStep);
  } else if (NumTimeStep < 10000) {
    PetscPrintf(PETSC_COMM_WORLD, "" GREEN "Step" RESET ": [%04d/%04d] | ",
                Time, NumTimeStep);
  } else if (NumTimeStep < 100000) {
    PetscPrintf(PETSC_COMM_WORLD, "" GREEN "Step" RESET ": [%05d/%05d] | ",
                Time, NumTimeStep);
  } else if (NumTimeStep < 1000000) {
    PetscPrintf(PETSC_COMM_WORLD, "" GREEN "Step" RESET ": [%i/%i] | ", Time,
                NumTimeStep);
  }

  PetscPrintf(PETSC_COMM_WORLD, "" GREEN "SNES L2-norm" RESET ": %1.4e | ",
              SNES_Norm);

  if (SNES_MaxIter < 10) {
    PetscPrintf(PETSC_COMM_WORLD,
                "" GREEN "SNES Iterations" RESET ": [%01d/%01d] | ", SNES_Iter,
                SNES_MaxIter);
  } else if (SNES_MaxIter < 100) {
    PetscPrintf(PETSC_COMM_WORLD,
                "" GREEN "SNES Iterations" RESET ": [%02d/%02d] | ", SNES_Iter,
                SNES_MaxIter);
  }

  PetscPrintf(PETSC_COMM_WORLD, "" GREEN "KSP L2-norm" RESET ": %1.4e | ",
              KSP_Norm);

  PetscPrintf(PETSC_COMM_WORLD, "" GREEN "KSP Iterations" RESET ": %02d | ",
              KSP_Iter);

  PetscPrintf(PETSC_COMM_WORLD, "" GREEN "Converged reason" RESET ": %s \n",
              SNESConvergedReasons[converged_reason]);

  FILE *Stats_Solver;
  char Name_file_t[10000];
  sprintf(Name_file_t, "%s/Stats_Solver.csv", OutputDir);
  Stats_Solver = fopen(Name_file_t, "a");

  if (Time == 0) {
    fprintf(Stats_Solver, "%s,%s,%s,%s,%s\n", "SNES Iterations",
            "KSP Iterations", "SNES L2-norm", "KSP L2-norm",
            "Converged reason");
  }
  fprintf(Stats_Solver, "%i,%i,%1.4e,%1.4e,%s\n", SNES_Iter, KSP_Iter,
          SNES_Norm, KSP_Norm, SNESConvergedReasons[converged_reason]);

  fclose(Stats_Solver);

  return EXIT_SUCCESS;
}

/*********************************************************************/

/**
 * @brief Project the kinematic information of the particles towards the nodes
 * \n of the mesh to get the nodal fields at t = n
 *
 * @param Lumped_Mass Lumped mass matrix used for the projection
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @param ActiveNodes List of nodes which takes place in the computation
 * @return Returns failure or success
 */
static PetscErrorCode __jacobian_variational_smoothing(const double * Lumped_Mass, Particle MPM_Mesh,
                                          Mesh FEM_Mesh, Mask ActiveNodes) {
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;
  unsigned NumberNodes_p;
  unsigned p;

#if NumberDimensions == 2
  int Mask_dofs_A[2];
  double local_jacobian_Ip[2];
#else
  int Mask_dofs_A[3];
  double local_jacobian_Ip[3];
#endif

  double * jacobian_I = (double *)calloc(Nactivenodes,sizeof(double));

#pragma omp parallel private(NumberNodes_p)
  {

#pragma omp for private(p)
    for (p = 0; p < Np; p++) {

      /* Define element of the particle */
      NumberNodes_p = MPM_Mesh.NumberNodes[p];
      Element Nodes_p =
          nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumberNodes_p);
      Matrix ShapeFunction_p =
          compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      //! Get the mass of the GP
      double m_p = MPM_Mesh.Phi.mass.nV[p];      
      double J_n1_p = MPM_Mesh.Phi.J_n1.nV[p];

      for (unsigned A = 0; A < NumberNodes_p; A++) {

        //  Get the node in the nodal momentum with the mask
        unsigned Ap = Nodes_p.Connectivity[A];
        int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];
        double ShapeFunction_pA = ShapeFunction_p.nV[A];
        double m__x__ShapeFunction_pA = m_p * ShapeFunction_pA;

#pragma omp critical
        {
          jacobian_I[Mask_node_A] += m__x__ShapeFunction_pA * J_n1_p;
        } // #pragma omp critical
      }   // for A

      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);

    } // for p

  } // #pragma omp parallel


#pragma omp for
  for (unsigned int idx = 0; idx < Nactivenodes; idx++) {

#pragma omp critical
    {
      jacobian_I[idx] = jacobian_I[idx]/Lumped_Mass[idx*Ndim];
    }
  }

#pragma omp for private(p)
    for (p = 0; p < Np; p++) {

      MPM_Mesh.Phi.J_n1.nV[p] = 0.0;

      //  Define nodal connectivity for each particle
      //  and compute the shape function
      unsigned int NumNodes_p = MPM_Mesh.NumberNodes[p];
      Element Nodes_p =
          nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);
      Matrix ShapeFunction_p =
          compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      //!  Update acceleration, velocity and position of the particles
      for (unsigned A = 0; A < NumNodes_p; A++) {

        // Get the shape function evaluation in node A
        // and the masked index of the node A
        double ShapeFunction_pI = ShapeFunction_p.nV[A];
        int Ap = Nodes_p.Connectivity[A];
        int A_mask = ActiveNodes.Nodes2Mask[Ap];

        MPM_Mesh.Phi.J_n1.nV[p] += ShapeFunction_pI * jacobian_I[A_mask];

      } // for unsigned A

    } // #pragma omp for private(p)


  free(jacobian_I);  

  return EXIT_SUCCESS;
}

/*********************************************************************/