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

  Mask ActiveNodes;
  Mask ActiveDOFs;
  IS Dirichlet_dofs;
  Particle MPM_Mesh;
  Mesh FEM_Mesh;
  Vec Lumped_Mass;

} Ctx;

/**************************************************************/
static double __compute_deltat(Particle MPM_Mesh, double h,
                               Time_Int_Params Parameters_Solver);

static int *__create_sparsity_pattern(Mask ActiveNodes, Particle MPM_Mesh);

static IS __get_dirichlet_list_dofs(Mask ActiveNodes, Mesh FEM_Mesh, int Step,
                                    int NumTimeStep);

static PetscErrorCode __compute_nodal_lumped_mass(Vec Lumped_MassMatrix,
                                                  Particle MPM_Mesh,
                                                  Mesh FEM_Mesh,
                                                  Mask ActiveNodes);

static PetscErrorCode __lagrangian_evaluation(SNES snes, Vec D_U, Vec Residual,
                                              void *ctx);

static PetscErrorCode __local_compatibility_conditions(const PetscScalar *dU,
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

static void
__nodal_inertial_forces(PetscScalar *Lagrangian, const PetscScalar *M_II,
                        Mask ActiveNodes,
                        Mask ActiveDOFs);

static PetscErrorCode __jacobian_evaluation(SNES snes, Vec dU, Mat Jacobian,
                                            Mat B, void *ctx);

static PetscErrorCode __monitor(PetscInt Time, PetscInt NumTimeStep,
                                PetscInt SNES_Iter, PetscInt KSP_Iter,
                                PetscInt SNES_MaxIter, PetscScalar KSP_Norm,
                                PetscScalar SNES_Norm,
                                SNESConvergedReason converged_reason);

static PetscErrorCode __update_Particles(Vec dU,
                                         Particle MPM_Mesh, Mesh FEM_Mesh,
                                         Mask ActiveNodes);

/**************************************************************/

// Global variables
static unsigned InitialStep;
static unsigned NumTimeStep;
static unsigned TimeStep;

/**************************************************************/

PetscErrorCode U_Static(Mesh FEM_Mesh, Particle MPM_Mesh,
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

  double DeltaTimeStep;
  double DeltaX = FEM_Mesh.DeltaX;

  Mat Tangent_Stiffness;
  int *sparsity_pattern;
  Vec Lumped_Mass;
  Vec Residual;
  Vec DU;
  IS Dirichlet_dofs;
  Mask ActiveNodes;
  Mask ActiveDOFs;

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
//    MatCreateAIJ(PETSC_COMM_SELF, Ntotaldofs, Ntotaldofs, 
//    PETSC_DETERMINE, PETSC_DETERMINE, 
//    PetscInt d_nz, const PetscInt d_nnz[], 
//    PetscInt o_nz, const PetscInt o_nnz[],&Tangent_Stiffness);

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
    Petsc_Direct_solver = true;
    Petsc_Iterative_solver = false;
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
    PetscCall(VecZeroEntries(DU));

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
       Update particle information
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(
        __update_Particles(DU, MPM_Mesh, FEM_Mesh, ActiveNodes));

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
    PetscCall(VecDestroy(&DU));
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

          /*
            Apply only if the direction is active
          */
          if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + Step] == 1) {
            Id_BCC_mask_k = Id_BCC_mask * Ndim + k;
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

  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;

  PetscScalar *Lagrangian_ptr;
  const PetscScalar *dU_ptr;
  const PetscScalar *Lumped_Mass_ptr;

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

  PetscCall(__local_compatibility_conditions(dU_ptr, ActiveNodes,
                                             MPM_Mesh, FEM_Mesh));

  PetscCall(__constitutive_update(MPM_Mesh, FEM_Mesh));

  PetscCall(__nodal_internal_forces(Lagrangian_ptr, ActiveNodes, ActiveDOFs,
                                    MPM_Mesh, FEM_Mesh));

  __nodal_traction_forces(Lagrangian_ptr, ActiveNodes, ActiveDOFs, MPM_Mesh,
                          FEM_Mesh);

  __nodal_inertial_forces(Lagrangian_ptr, Lumped_Mass_ptr, ActiveNodes, ActiveDOFs);

  /**
   * Restore vectors
   *
   */
  PetscCall(VecRestoreArray(Lagrangian, &Lagrangian_ptr));
  PetscCall(VecRestoreArrayRead(dU, &dU_ptr));

  return STATUS;
}

/**************************************************************/

/**
 * @brief
 *
 * @param dU Incremental nodal displacement field
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @return PetscErrorCode
 */
static PetscErrorCode __local_compatibility_conditions(const PetscScalar *dU,
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
      if (D_Displacement_Ap == NULL) {
        fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
        STATUS = EXIT_FAILURE;
      }

      get_set_field__MeshTools__(D_Displacement_Ap, dU, Nodes_p, ActiveNodes);

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
 
      update_increment_Deformation_Gradient__Particles__(
          DF_p, D_Displacement_Ap, gradient_p.nV, NumberNodes_p);

      /*
        Update the deformation gradient in t = n + 1 with the information
        from t = n and the increment of deformation gradient.
      */
      update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, DF_p);

      //  Compute Jacobian of the deformation gradient
      MPM_Mesh.Phi.J_n1.nV[p] = I3__TensorLib__(F_n1_p);
      if (MPM_Mesh.Phi.J_n1.nV[p] <= 0.0) {
        fprintf(stderr,
                "" RED "Negative jacobian in particle %i: %e" RESET " \n", p,
                MPM_Mesh.Phi.J_n1.nV[p]);
        STATUS = EXIT_FAILURE;
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
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param ActiveDOFs List of dofs which takes place in the computation
 */
static void
__nodal_inertial_forces(PetscScalar *Lagrangian, const PetscScalar *M_II,
                        Mask ActiveNodes,
                        Mask ActiveDOFs) {

  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;
  unsigned idx;

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
          - M_II[idx] * b[idx % Ndim];
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
              d_shapefunction_n_pB, 0.0, MPM_Mesh, MatProp_p);
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

/**
 * @brief
 *
 * @param dU Incremental nodal displacement field
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @param ActiveNodes List of nodes which takes place in the computation
 * @return PetscErrorCode
 */
static PetscErrorCode __update_Particles(Vec dU, 
                                         Particle MPM_Mesh, Mesh FEM_Mesh,
                                         Mask ActiveNodes) {

  PetscErrorCode STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned p;

  double DU_pI;
  const PetscScalar *dU_ptr;
  PetscCall(VecGetArrayRead(dU, &dU_ptr));

#pragma omp parallel private(NumNodes_p, DU_pI)
  {
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

      //  Define nodal connectivity for each particle
      //  and compute the shape function
      NumNodes_p = MPM_Mesh.NumberNodes[p];
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

        for (unsigned i = 0; i < Ndim; i++) {
          DU_pI = ShapeFunction_pI * dU_ptr[A_mask * Ndim + i];
          MPM_Mesh.Phi.dis.nM[p][i] += DU_pI;
          MPM_Mesh.Phi.x_GC.nM[p][i] += DU_pI;
        }
      } // for unsigned A

      free(Nodes_p.Connectivity);
      free__MatrixLib__(ShapeFunction_p);

    } // #pragma omp for private(p)
  }   // #pragma omp parallel

  PetscCall(VecRestoreArrayRead(dU, &dU_ptr));

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