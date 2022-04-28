
#include "Formulations/Displacements/U-Newmark-beta.h"

typedef struct {

  //#ifdef USE_PETSC
  //  Vec value;
  //  Vec d_value_dt;
  //  Vec d2_value_dt2;
  //#else
  double *value;
  double *d_value_dt;
  double *d2_value_dt2;
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
  double DeltaTimeStep;
} Newmark_parameters;

static double __compute_deltat(Particle MPM_Mesh /**< */, double h /**< */,
                               Time_Int_Params Parameters_Solver /**< */);
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
                                                       double epsilon);
/**************************************************************/

/*!
  \param[in,out] Effective_MassMatrix
*/
static int __compute_nodal_effective_mass(double *Effective_MassMatrix /**< */,
                                          Particle MPM_Mesh /**< */,
                                          Mesh FEM_Mesh /**< */,
                                          Mask ActiveNodes /**< */,
                                          double epsilon /**< */);
/**************************************************************/

static int __get_nodal_field_tn(Nodal_Field U_n /**< */,
                                Particle MPM_Mesh /**< */,
                                Mesh FEM_Mesh /**< */,
                                Mask ActiveNodes /**< */);
/**************************************************************/

/*!

*/
static void __initialise_nodal_increments(Nodal_Field DU,
                                          Nodal_Field U_n /**< */,
                                          Mesh FEM_Mesh /**< */,
                                          Mask ActiveNodes /**< */,
                                          Newmark_parameters Params /**< */);
/**************************************************************/

/*!
  \brief Update the local deformation of the particles \n
  ensuring the local compatibility conditions
  \param[in] D_U Increment of nodal kinetics
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[out] STATUS Returns failure or success
*/
static void __local_compatibility_conditions(const Nodal_Field D_U,
                                             Mask ActiveNodes,
                                             Particle MPM_Mesh, Mesh FEM_Mesh,
                                             int *STATUS);
/**************************************************************/

/*!
  \brief Function used to compute the equilibrium residual

  \param[in] U_n Nodal kinetics information from the step n
  \param[in] D_U Increment of nodal kinetics
  \param[in] Effective_Mass Effective mass matrix
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] Params Time integration parameters
  \param[in] Is_compute_Residual The function computes the residual
  \param[in] Is_compute_Reactions The function computes the reaction
  \param[out] STATUS Returns failure or success
  \return Residual vector
*/
#ifdef USE_PETSC
static Vec __assemble_residual(Nodal_Field U_n, Nodal_Field D_U,
                               double *Effective_Mass, Mask ActiveNodes,
                               Mask ActiveDOFs, Particle MPM_Mesh,
                               Mesh FEM_Mesh, Newmark_parameters Params,
                               bool Is_compute_Residual,
                               bool Is_compute_Reactions, int *STATUS);
#else
static double *__assemble_residual(Nodal_Field U_n, Nodal_Field D_U,
                                   double *Effective_Mass, Mask ActiveNodes,
                                   Mask ActiveDOFs, Particle MPM_Mesh,
                                   Mesh FEM_Mesh, Newmark_parameters Params,
                                   bool Is_compute_Residual,
                                   bool Is_compute_Reactions, int *STATUS);
#endif
/**************************************************************/

/*!
  \brief Function used to compute the contribution of the \n
  internal forces to the residual

  \param[in,out] Residual Residual/reaction vector
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] Is_compute_Residual The function computes the residual
  \param[in] Is_compute_Reactions The function computes the reaction
*/
#ifdef USE_PETSC
static int __Nodal_Internal_Forces(Vec Residual, Mask ActiveNodes,
                                   Mask ActiveDOFs, Particle MPM_Mesh,
                                   Mesh FEM_Mesh, bool Is_compute_Residual,
                                   bool Is_compute_Reactions);
#else
static int __Nodal_Internal_Forces(double *Residual, Mask ActiveNodes,
                                   Mask ActiveDOFs, Particle MPM_Mesh,
                                   Mesh FEM_Mesh, bool Is_compute_Residual,
                                   bool Is_compute_Reactions);
#endif
/**************************************************************/

static void __internal_force_density(double *InternalForcesDensity_Ap,
                                     const double *kirchhoff_p,
                                     const double *gradient_n1_pA);
/**************************************************************/

/*!
  \brief Function used to compute the contribution of the \n
  contact forces to the residual

  \param[in,out] Residual Residual vector
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] Is_compute_Residual The function computes the residual
  \param[in] Is_compute_Reactions The function computes the reaction
*/
#ifdef USE_PETSC
static void __Nodal_Traction_Forces(Vec Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions);
#else
static void __Nodal_Traction_Forces(double *Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions);
#endif
/**************************************************************/

/*!
  \brief Function used to compute the contribution of the \n
  body (distance) forces to the residual

  \param[in,out] Residual Residual vector
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] Is_compute_Residual The function computes the residual
  \param[in] Is_compute_Reactions The function computes the reaction
*/
#ifdef USE_PETSC
static void __Nodal_Body_Forces(Vec Residual, Mask ActiveNodes, Mask ActiveDOFs,
                                Particle MPM_Mesh, Mesh FEM_Mesh,
                                bool Is_compute_Residual,
                                bool Is_compute_Reactions);
#else
static void __Nodal_Body_Forces(double *Residual, Mask ActiveNodes,
                                Mask ActiveDOFs, Particle MPM_Mesh,
                                Mesh FEM_Mesh, bool Is_compute_Residual,
                                bool Is_compute_Reactions);
#endif
/**************************************************************/

/*!
  \brief Function used to compute the contribution of the \n
  inertial forces to the residual

  \param[in,out] Residual Residual vector
  \param[in] Mass Mass matrix
  \param[in] U_n Nodal kinetics information from the step n
  \param[in] D_U Increment of nodal kinetics
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] Params Time integration parameters
*/
#ifdef USE_PETSC
static void __Nodal_Inertial_Forces(Vec Residual, double *Mass, Nodal_Field U_n,
                                    Nodal_Field D_U, Mask ActiveNodes,
                                    Mask ActiveDOFs, Newmark_parameters Params);
#else
static void __Nodal_Inertial_Forces(double *Residual, double *Mass,
                                    Nodal_Field U_n, Nodal_Field D_U,
                                    Mask ActiveNodes, Mask ActiveDOFs,
                                    Newmark_parameters Params);
#endif
/**************************************************************/

/*!
  \brief Compute the 2-norm of the residual
  \param[in] Residual
  \param[in] Total_dof Size of the vector
  \return The norm of the residual
*/
#ifdef USE_PETSC
static double __error_residual(const Vec Residual, unsigned Total_dof);
#else
static double __error_residual(const double *Residual, unsigned Total_dof);
#endif
/**************************************************************/

static void __preallocation_tangent_matrix(int *nnz /**< */,
                                           Mask ActiveNodes /**< */,
                                           Mask ActiveDOFs /**< */,
                                           Particle MPM_Mesh /**< */);
/**************************************************************/

static void compute_local_intertia(double *Inertia_density_p /**< */,
                                   double Na_p /**< */, double Nb_p /**< */,
                                   double m_p /**< */, double alpha_1 /**< */,
                                   double epsilon /**< */, unsigned A /**< */,
                                   unsigned B /**< */);
/**************************************************************/

/*!
  \param[in] nnz Array for the Yale format storage
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] Params Time integration parameters
  \param[out] STATUS Returns failure or success
  \return The tangent matrix for the problem
*/
#ifdef USE_PETSC
static Mat __assemble_tangent_stiffness(int *nnz, Mask ActiveNodes,
                                        Mask ActiveDOFs, Particle MPM_Mesh,
                                        Mesh FEM_Mesh,
                                        Newmark_parameters Params, int *STATUS);
#else
static double *__assemble_tangent_stiffness(int *nnz, Mask ActiveNodes,
                                            Mask ActiveDOFs, Particle MPM_Mesh,
                                            Mesh FEM_Mesh,
                                            Newmark_parameters Params,
                                            int *STATUS);
#endif
/**************************************************************/

/*!
  \brief This function takes the nodal increments comming from \n
  the linear solver and computes the updated nodal increments \n
  for the nodal kinetics
  \param[in] Residual Increments during the k iteration
  \param[in,out] D_U Nodal incremented (displacement,velocity,acceleration)
  \param[in] U_n Previously converged nodal kinetics
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] Params Time integration parameters
  \param[in] Ntotaldofs Number of dofs
*/
#ifdef USE_PETSC
static void __update_Nodal_Increments(const Vec Residual, Nodal_Field D_U,
                                      Nodal_Field U_n, Mask ActiveDOFs,
                                      Newmark_parameters Params,
                                      unsigned Ntotaldofs);
#else
static void __update_Nodal_Increments(const double *Residual, Nodal_Field D_U,
                                      Nodal_Field U_n, Mask ActiveDOFs,
                                      Newmark_parameters Params,
                                      unsigned Ntotaldofs);
#endif
/**************************************************************/

static void __update_Particles(Nodal_Field D_U /**< */,
                               Particle MPM_Mesh /**< */, Mesh FEM_Mesh /**< */,
                               Mask ActiveNodes /**< */);
/**************************************************************/

int U_Newmark_Beta(Mesh FEM_Mesh, Particle MPM_Mesh,
                   Time_Int_Params Parameters_Solver) {

  int STATUS = EXIT_SUCCESS;

  /*
    Auxiliar variables for the solver
  */
  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes;
  unsigned Ntotaldofs;
  unsigned Nactivedofs;
  unsigned MaxIter = Parameters_Solver.MaxIter;
  unsigned Iter;

  // Time integration variables
  InitialStep = Parameters_Solver.InitialTimeStep;
  NumTimeStep = Parameters_Solver.NumTimeStep;
  TimeStep = InitialStep;

  double TOL = Parameters_Solver.TOL_Newmark_beta;
  double epsilon = Parameters_Solver.epsilon_Mass_Matrix;
  double beta = Parameters_Solver.beta_Newmark_beta;
  double gamma = Parameters_Solver.gamma_Newmark_beta;
  double DeltaTimeStep;
  double DeltaX = FEM_Mesh.DeltaX;

  double Error_0;
  double Error_i;
  double Error_relative;
  //  double Error_increment_i;

#ifdef USE_PETSC
  //  Mat Effective_Mass;
  Mat Tangent_Stiffness;
  Vec Residual;
  Vec Reactions;
  double *Effective_Mass;
  int *nnz;
#else
  double *Tangent_Stiffness;
  double *Residual;
  double *Reactions;
  double *Effective_Mass;
  int *nnz;
#endif

  Nodal_Field U_n;
  Nodal_Field D_U;

  Mask ActiveNodes;
  Mask ActiveDOFs;

  Newmark_parameters Time_Integration_Params;

// InoutParameters
#ifdef USE_PETSC
  PetscViewer viewer;
#else

#endif

  /*
    Time step is defined at the init of the simulation throught the
    CFL condition. Notice that for this kind of solver, CFL confition is
    not required to be satisfied. The only purpose of it is to use the existing
    software interfase.
  */
  DeltaTimeStep = __compute_deltat(MPM_Mesh, DeltaX, Parameters_Solver);

  //  Compute time integrator parameters
  Time_Integration_Params =
      __compute_Newmark_parameters(beta, gamma, DeltaTimeStep, epsilon);

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
  PetscLogEvent event_1;
  PetscLogDouble event_1_flops;
  PetscLogEvent event_2;
  PetscLogDouble event_2_flops;
  PetscLogEvent event_3;
  PetscLogDouble event_3_flops;
  PetscLogEvent event_4;
  PetscLogDouble event_4_flops;
  PetscLogEvent event_5;
  PetscLogDouble event_5_flops;
  PetscLogEvent event_6;
  PetscLogDouble event_6_flops;
#endif
#endif

  while (TimeStep < NumTimeStep) {
    print_Status("*************************************************", TimeStep);
    print_step(TimeStep, DeltaTimeStep);

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
    PetscLogEventRegister("Set shape functions", 0, &event_1);
    PetscLogEventBegin(event_1, 0, 0, 0, 0);
#endif
#endif

    local_search__MeshTools__(MPM_Mesh, FEM_Mesh);
    ActiveNodes = get_active_nodes__MeshTools__(FEM_Mesh);
    Nactivenodes = ActiveNodes.Nactivenodes;
    Ntotaldofs = Ndim * Nactivenodes;
    ActiveDOFs = get_active_dofs__MeshTools__(ActiveNodes, FEM_Mesh, TimeStep,
                                              NumTimeStep);
    Nactivedofs = ActiveDOFs.Nactivenodes;

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
    PetscLogFlops(event_1_flops);
    PetscLogEventEnd(event_1, 0, 0, 0, 0);
#endif
#endif

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
    PetscLogEventRegister("Set Initial nodal kinetics", 0, &event_2);
    PetscLogEventBegin(event_2, 0, 0, 0, 0);
#endif
#endif

    // Get the previous converged nodal value
    U_n.value = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    U_n.d_value_dt = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    U_n.d2_value_dt2 = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    if ((U_n.value == NULL) || (U_n.d_value_dt == NULL) ||
        (U_n.d2_value_dt2 == NULL)) {
      fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
      return EXIT_FAILURE;
    }
    STATUS = __get_nodal_field_tn(U_n, MPM_Mesh, FEM_Mesh, ActiveNodes);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __get_nodal_field_tn()" RESET " \n");
      return EXIT_FAILURE;
    }

    // Compute kinematic nodal values
    D_U.value = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    D_U.d_value_dt = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    D_U.d2_value_dt2 = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    if ((D_U.value == NULL) || (D_U.d_value_dt == NULL) ||
        (D_U.d2_value_dt2 == NULL)) {
      fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
      return EXIT_FAILURE;
    }
    __initialise_nodal_increments(D_U, U_n, FEM_Mesh, ActiveNodes,
                                  Time_Integration_Params);

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
    PetscLogFlops(event_2_flops);
    PetscLogEventEnd(event_2, 0, 0, 0, 0);
#endif
#endif

    // Define and allocate the effective mass matrix
    Effective_Mass =
        (double *)calloc(Ntotaldofs * Ntotaldofs, __SIZEOF_DOUBLE__);
    if (Effective_Mass == NULL) {
      fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
      return EXIT_FAILURE;
    }
    STATUS = __compute_nodal_effective_mass(Effective_Mass, MPM_Mesh, FEM_Mesh,
                                            ActiveNodes, epsilon);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in __compute_nodal_effective_mass()" RESET " \n");
      return EXIT_FAILURE;
    }

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
    PetscLogEventRegister("Set Initial Residual", 0, &event_3);
    PetscLogEventBegin(event_3, 0, 0, 0, 0);
#endif
#endif

    // Trial residual
    Residual = __assemble_residual(
        U_n, D_U, Effective_Mass, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
        Time_Integration_Params, true, false, &STATUS);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __assemble_residual()" RESET " \n");
      return EXIT_FAILURE;
    }

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
    PetscLogFlops(event_3_flops);
    PetscLogEventEnd(event_3, 0, 0, 0, 0);
#endif
#endif

    // Compute error
    Error_0 = Error_i = __error_residual(Residual, Nactivedofs);
    Error_relative = Error_i / Error_0;
    Iter = 0;

    nnz = (int *)calloc(Nactivedofs, __SIZEOF_INT__);
    __preallocation_tangent_matrix(nnz, ActiveNodes, ActiveDOFs, MPM_Mesh);

    while (Error_relative > TOL) {

      if ((Error_i < TOL * 100) || (Error_relative < TOL) || (Iter > MaxIter)) {
        break;
      }

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
      PetscLogEventRegister("Assemble tangent matrix", 0, &event_4);
      PetscLogEventBegin(event_4, 0, 0, 0, 0);
#endif
#endif

      Tangent_Stiffness = __assemble_tangent_stiffness(
          nnz, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
          Time_Integration_Params, &STATUS);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __assemble_tangent_stiffness()" RESET " \n");
        return EXIT_FAILURE;
      }

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
      PetscLogFlops(event_4_flops);
      PetscLogEventEnd(event_4, 0, 0, 0, 0);
#endif
#endif

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
      PetscLogEventRegister("Solve linear system", 0, &event_5);
      PetscLogEventBegin(event_5, 0, 0, 0, 0);
#endif
#endif

#ifdef USE_PETSC
      STATUS = krylov_PETSC(&Tangent_Stiffness, &Residual, Nactivedofs);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in krylov_PETSC()" RESET " \n");
        return EXIT_FAILURE;
      }
#else
      STATUS = dgetrs_LAPACK(Tangent_Stiffness, Residual, Nactivedofs);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in dgetrs_LAPACK()" RESET " \n");
        return EXIT_FAILURE;
      }
#endif

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
      PetscLogFlops(event_5_flops);
      PetscLogEventEnd(event_5, 0, 0, 0, 0);
#endif
#endif

      __update_Nodal_Increments(Residual, D_U, U_n, ActiveDOFs,
                                Time_Integration_Params, Ntotaldofs);

#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
      PetscLogEventRegister("Update local deformation", 0, &event_6);
      PetscLogEventBegin(event_6, 0, 0, 0, 0);
#endif
#endif


      __local_compatibility_conditions(D_U, ActiveNodes, MPM_Mesh, FEM_Mesh,
                                           &STATUS);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __local_compatibility_conditions()" RESET
                " \n");
        return EXIT_FAILURE;
      }
#ifdef USE_PETSC
#if defined(PETSC_USE_LOG)
      PetscLogFlops(event_6_flops);
      PetscLogEventEnd(event_6, 0, 0, 0, 0);
#endif
#endif

      // Free memory
#ifdef USE_PETSC
      VecDestroy(&Residual);
      MatDestroy(&Tangent_Stiffness);
#else
      free(Residual);
      free(Tangent_Stiffness);
#endif

      // Compute residual (NR-loop)
      Residual = __assemble_residual(
          U_n, D_U, Effective_Mass, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
          Time_Integration_Params, true, false, &STATUS);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __assemble_residual()" RESET " \n");
        return EXIT_FAILURE;
      }

      // Get stats for the convergence
      Error_i = __error_residual(Residual, Nactivedofs);
      Error_relative = Error_i / Error_0;
      Iter++;
      printf("Iter: [%i/%i]. Total Error: %e, Relative Error: %e \n", Iter,
             MaxIter, Error_i, Error_relative);
    }

    print_convergence_stats(TimeStep, Iter, Error_0, Error_i, Error_relative);

    if (Iter > MaxIter) {
      fprintf(
          stderr,
          "" RED
          "Convergence not reached in the maximum number of iterations" RESET
          " \n");
    }

    __update_Particles(D_U, MPM_Mesh, FEM_Mesh, ActiveNodes);

    if (TimeStep % ResultsTimeStep == 0) {
      Reactions = __assemble_residual(
          U_n, D_U, Effective_Mass, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
          Time_Integration_Params, false, true, &STATUS);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __assemble_residual()" RESET " \n");
        return EXIT_FAILURE;
      }

      //      Matrix Reactions_aux =
      //      memory_to_matrix__MatrixLib__(Nactivenodes,Ndim,Reactions);

      //      nodal_results_vtk__InOutFun__(FEM_Mesh, ActiveNodes,
      //      Reactions_aux, TimeStep, ResultsTimeStep);

      particle_results_vtk__InOutFun__(MPM_Mesh, TimeStep, ResultsTimeStep);

#ifdef USE_PETSC
      VecDestroy(&Reactions);
#else
      free(Reactions);
#endif
      //      free(Reactions_aux.nM);
    }

    TimeStep++;

#ifdef USE_PETSC
    VecDestroy(&Residual);
#else
    free(Residual);
#endif

    free(Effective_Mass);
    free(U_n.value);
    free(U_n.d_value_dt);
    free(U_n.d2_value_dt2);
    free(D_U.value);
    free(D_U.d_value_dt);
    free(D_U.d2_value_dt2);
    free(ActiveNodes.Nodes2Mask);
    free(ActiveDOFs.Nodes2Mask);
    free(nnz);
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

static int __compute_nodal_effective_mass(double *Effective_MassMatrix,
                                          Particle MPM_Mesh, Mesh FEM_Mesh,
                                          Mask ActiveNodes, double epsilon)
/*
  This function computes the effective mass matrix as a convex combination
  of the lumped mass matrix and the consistent mass matrix. Later assemble
  a total mass matrix with the contribution of each degree of freedom.
*/
{

  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned Order = Ndim * Nnodes_mask;
  unsigned NumNodes_p;
  int Ap;
  int Bp;
  int A_mask;
  int B_mask;
  int idx_A_mask_i;
  int idx_B_mask_i;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;

  /* Evaluation of the particle in the node */
  double ShapeFunction_pA, ShapeFunction_pB;
  /* Mass of the particle */
  double m_p;
  /* Nodal contribution A of the particle p */
  double m_A_p;
  /* Nodal contribution A-B of the particle p */
  double m_AB_p;
  /* Element for each particle */
  Element Nodes_p;

  // Define and allocate the lumped mass matrix
  double *Lumped_MassMatrix = (double *)calloc(Order, __SIZEOF_DOUBLE__);
  if (Lumped_MassMatrix == NULL) {
    fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
    return EXIT_FAILURE;
  }

  /*
    Iterate over the particles to get the nodal values
  */
  for (unsigned p = 0; p < Np; p++) {

    /*
      Define tributary nodes of the particle
    */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /*
       Evaluate the shape function in the coordinates of the particle
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Get the mass of the particle
    */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    for (unsigned A = 0; A < NumNodes_p; A++) {

      /*
         Get the node in the mass matrix with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /* Get the value of the shape function */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /*
        Compute the nodal A contribution of the particle p
      */
      m_A_p = m_p * ShapeFunction_pA;

      /*
         Fill the Lumped mass matrix considering the number of dofs
      */
      for (unsigned i = 0; i < Ndim; i++) {
        idx_A_mask_i = A_mask * Ndim + i;
        Lumped_MassMatrix[idx_A_mask_i] += m_A_p;
      }

      for (unsigned B = 0; B < NumNodes_p; B++) {
        /*
           Get the node in the mass matrix with the mask
        */
        Bp = Nodes_p.Connectivity[B];
        B_mask = ActiveNodes.Nodes2Mask[Bp];

        /* Get the value of the shape function */
        ShapeFunction_pB = ShapeFunction_p.nV[B];

        /*
          Compute the nodal AB contribution of the particle p
        */
        m_AB_p = m_p * ShapeFunction_pA * ShapeFunction_pB;

        /*
           Fill the effective mass matrix considering the number of dofs
        */
        for (unsigned i = 0; i < Ndim; i++) {
          /*
            Compute the vectorized index
          */
          idx_A_mask_i = A_mask * Ndim + i;
          idx_B_mask_i = B_mask * Ndim + i;
          Effective_MassMatrix[idx_A_mask_i * Order + idx_B_mask_i] += m_AB_p;
        }
      }
    }

    /* Free the value of the shape functions */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

  /*
    At this point the effective mass matrix coincides with the consistent mass
    matrix. We can tune it by a convecx combination with the lumped mass matrix
  */
  for (unsigned A = 0; A < Order; A++) {
    for (unsigned B = 0; B < Order; B++) {
      Effective_MassMatrix[A * Order + B] =
          (1 - epsilon) * Effective_MassMatrix[A * Order + B] +
          (A == B) * epsilon * Lumped_MassMatrix[A];
    }
  }

  free(Lumped_MassMatrix);

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __get_nodal_field_tn(Nodal_Field U_n, Particle MPM_Mesh,
                                Mesh FEM_Mesh, Mask ActiveNodes)
/*
  Call the LAPACK solver to compute the nodal velocity. The operation is
  linearized and all the dof split the velocity array in n components like : | M
  0 |   |V.x|   | p.x | | 0 M | * |V.y| = | p.y |

*/
{
  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned Ap;
  int A_mask;
  int idx_A_mask_i;
  double epsilon = 1;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;

  /* Evaluation of the particle in the node */
  double ShapeFunction_pA;
  /* Mass of the particle */
  double m_p;
  /* Element for each particle */
  Element Nodes_p;

  /* Iterate over the particles to get the nodal values */
  for (unsigned p = 0; p < Np; p++) {

    /* Define element of the particle */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /* Evaluate the shape function in the coordinates of the particle */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Get the mass of the GP */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* Get the nodal mommentum */
    for (unsigned A = 0; A < NumNodes_p; A++) {

      /*
        Get the node in the nodal momentum with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /* Evaluate the GP function in the node */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /* Nodal velocity and acceleration  */
      for (unsigned i = 0; i < Ndim; i++) {

        idx_A_mask_i = A_mask * Ndim + i;

        U_n.value[idx_A_mask_i] +=
            m_p * ShapeFunction_pA * MPM_Mesh.Phi.dis.nM[p][i];
        U_n.d_value_dt[idx_A_mask_i] +=
            m_p * ShapeFunction_pA * MPM_Mesh.Phi.vel.nM[p][i];
        U_n.d2_value_dt2[idx_A_mask_i] +=
            m_p * ShapeFunction_pA * MPM_Mesh.Phi.acc.nM[p][i];
      }
    }

    /* Free the value of the shape functions */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Call the LAPACK solver to compute the accelerations and velocities
  */
  int Order = Ndim * Nnodes_mask;
  int LDA = Order;
  int LDB = Order;
  char TRANS = 'N';
  int INFO = 3;
  int NRHS = 1;

  int *IPIV = (int *)calloc(Order, __SIZEOF_INT__);
  if (IPIV == NULL) {
    fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
    return EXIT_FAILURE;
  }

  double *Effective_Mass = (double *)calloc(Order * Order, __SIZEOF_DOUBLE__);
  if (Effective_Mass == NULL) {
    fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
    return EXIT_FAILURE;
  }
  STATUS = __compute_nodal_effective_mass(Effective_Mass, MPM_Mesh, FEM_Mesh,
                                          ActiveNodes, epsilon);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,
            "" RED "Error in __compute_nodal_effective_mass()" RESET " \n");
    return EXIT_FAILURE;
  }

  //  Compute the LU factorization for the mass matrix
  dgetrf_(&Order, &Order, Effective_Mass, &LDA, IPIV, &INFO);
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

  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass, &LDA, IPIV, U_n.value, &LDB,
          &INFO);
  if (INFO) {
    fprintf(stderr, "" RED "Error in dgetrs_() " RESET "\n");
    return EXIT_FAILURE;
  }

  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass, &LDA, IPIV, U_n.d_value_dt,
          &LDB, &INFO);
  if (INFO) {
    fprintf(stderr, "" RED "Error in dgetrs_() " RESET "\n");
    return EXIT_FAILURE;
  }

  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass, &LDA, IPIV, U_n.d2_value_dt2,
          &LDB, &INFO);
  if (INFO) {
    fprintf(stderr, "" RED "Error in dgetrs_() " RESET "\n");
    return EXIT_FAILURE;
  }

  //  Free auxiliar memory
  free(Effective_Mass);
  free(IPIV);

  return STATUS;
}

/**************************************************************/

static void __initialise_nodal_increments(Nodal_Field D_U, Nodal_Field U_n,
                                          Mesh FEM_Mesh, Mask ActiveNodes,
                                          Newmark_parameters Params)
/*
  Apply the boundary conditions over the nodes
*/
{
  unsigned NumNodesBound;
  unsigned Ndim = NumberDimensions;
  unsigned NumBounds = FEM_Mesh.Bounds.NumBounds;
  int Id_BCC;
  int Id_BCC_mask;

  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

  double D_U_value_It;

  /*
    Loop over the the boundaries to set boundary conditions
  */
  for (unsigned i = 0; i < NumBounds; i++) {

    /*
      Get the number of nodes of this boundarie
    */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

    for (unsigned j = 0; j < NumNodesBound; j++) {
      /*
        Get the index of the node
      */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

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
          /*
            Assign the boundary condition :
            First remove the value obtained during the projection and compute
            the value using the evolution of the boundary condition
          */
          U_n.value[Id_BCC_mask * Ndim + k] = 0.0;
          U_n.d_value_dt[Id_BCC_mask * Ndim + k] = 0.0;
          U_n.d2_value_dt2[Id_BCC_mask * Ndim + k] = 0.0;

          for (unsigned t = 0; t < TimeStep; t++) {
            D_U_value_It = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[t];
            U_n.value[Id_BCC_mask * Ndim + k] += D_U_value_It;
            U_n.d2_value_dt2[Id_BCC_mask * Ndim + k] +=
                alpha_1 * D_U_value_It -
                alpha_2 * U_n.d_value_dt[Id_BCC_mask * Ndim + k] -
                (alpha_3 + 1) * U_n.d2_value_dt2[Id_BCC_mask * Ndim + k];
            U_n.d_value_dt[Id_BCC_mask * Ndim + k] +=
                alpha_4 * D_U_value_It +
                (alpha_5 - 1) * U_n.d_value_dt[Id_BCC_mask * Ndim + k] +
                alpha_6 * U_n.d2_value_dt2[Id_BCC_mask * Ndim + k];
          }

          /*
            Initialise increments using newmark and the value of the boundary
            condition
          */
          D_U.value[Id_BCC_mask * Ndim + k] =
              FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep];
          D_U.d2_value_dt2[Id_BCC_mask * Ndim + k] =
              alpha_1 * D_U.value[Id_BCC_mask * Ndim + k] -
              alpha_2 * U_n.d_value_dt[Id_BCC_mask * Ndim + k] -
              (alpha_3 + 1) * U_n.d2_value_dt2[Id_BCC_mask * Ndim + k];
          D_U.d_value_dt[Id_BCC_mask * Ndim + k] =
              alpha_4 * D_U.value[Id_BCC_mask * Ndim + k] +
              (alpha_5 - 1) * U_n.d_value_dt[Id_BCC_mask * Ndim + k] +
              alpha_6 * U_n.d2_value_dt2[Id_BCC_mask * Ndim + k];
        }
      }
    }
  }
}

/**************************************************************/

static void __local_compatibility_conditions(Nodal_Field D_U, Mask ActiveNodes,
                                             Particle MPM_Mesh, Mesh FEM_Mesh,
                                             int *STATUS) {

  *STATUS = EXIT_SUCCESS;
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
  #pragma omp parallel private(NumberNodes_p,Order_p)
  {

    #pragma omp for private(p)
    for (p = 0; p < Np; p++) {

    //  Define tributary nodes of the particle
    NumberNodes_p = MPM_Mesh.NumberNodes[p];
    Element Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumberNodes_p);
    Order_p = NumberNodes_p * Ndim;

    //  Get the nodal increment of displacement using the mask
    double *D_Displacement_Ap = (double *)calloc(Order_p, __SIZEOF_DOUBLE__);
    double *D_Velocity_Ap = (double *)calloc(Order_p, __SIZEOF_DOUBLE__);
    if ((D_Displacement_Ap == NULL) || (D_Velocity_Ap == NULL)) {
      fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
      *STATUS = EXIT_FAILURE;
    }

    get_set_field__MeshTools__(D_Displacement_Ap, D_U.value, Nodes_p,
                               ActiveNodes);
    get_set_field__MeshTools__(D_Velocity_Ap, D_U.d_value_dt, Nodes_p,
                               ActiveNodes);

    /*
      Evaluate the shape function gradient in the coordinates of the particle
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
    update_rate_Deformation_Gradient_n1__Particles__(dFdt_n1_p, dt_DF_p, F_n_p,
                                                     DF_p, dFdt_n_p);

    //  Compute Jacobian of the deformation gradient
    MPM_Mesh.Phi.J_n1.nV[p] = I3__TensorLib__(F_n1_p);
    if (MPM_Mesh.Phi.J_n1.nV[p] <= 0.0) {
      fprintf(stderr, "" RED "Negative jacobian in particle %i" RESET " \n", p);
      *STATUS = EXIT_FAILURE;
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

    #pragma omp for private(p,Vn_patch,Vn1_patch,J_patch,Idx_Element_p,Idx_Patch_p)
    for (p = 0; p < Np; p++) {

      Idx_Element_p = MPM_Mesh.Element_p[p];
      Idx_Patch_p = FEM_Mesh.Idx_Patch[Idx_Element_p];

      Vn_patch = FEM_Mesh.Vol_Patch_n[Idx_Patch_p];
      Vn1_patch = FEM_Mesh.Vol_Patch_n1[Idx_Patch_p];
      J_patch = Vn1_patch / Vn_patch;

      *STATUS = get_locking_free_Deformation_Gradient_n1__Particles__(
          p, J_patch, MPM_Mesh);
      if (*STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in "
                "get_locking_free_Deformation_Gradient_n1__Particles__()" RESET
                " \n");
        *STATUS = EXIT_FAILURE;
      }

      MPM_Mesh.Phi.Jbar.nV[p] *= J_patch;
    }
  }

  }
}

/**************************************************************/

#ifdef USE_PETSC
static Vec __assemble_residual(Nodal_Field U_n, Nodal_Field D_U,
                               double *Effective_Mass, Mask ActiveNodes,
                               Mask ActiveDOFs, Particle MPM_Mesh,
                               Mesh FEM_Mesh, Newmark_parameters Params,
                               bool Is_compute_Residual,
                               bool Is_compute_Reactions, int *STATUS)
#else
static double *__assemble_residual(Nodal_Field U_n, Nodal_Field D_U,
                                   double *Effective_Mass, Mask ActiveNodes,
                                   Mask ActiveDOFs, Particle MPM_Mesh,
                                   Mesh FEM_Mesh, Newmark_parameters Params,
                                   bool Is_compute_Residual,
                                   bool Is_compute_Reactions, int *STATUS)
#endif
{

  unsigned Size;

  if ((Is_compute_Residual == true) && (Is_compute_Reactions == false)) {
    Size = ActiveDOFs.Nactivenodes;
  } else if ((Is_compute_Residual == false) && (Is_compute_Reactions == true)) {
    Size = ActiveNodes.Nactivenodes * NumberDimensions;
  }

#ifdef USE_PETSC
  Vec Residual;
  VecCreate(PETSC_COMM_WORLD, &Residual);
  VecSetSizes(Residual, PETSC_DECIDE, Size);
  VecSetFromOptions(Residual);
#else
  double *Residual = (double *)calloc(Size, __SIZEOF_DOUBLE__);
  if (Residual == NULL) {
    fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
    *STATUS = EXIT_FAILURE;
    return Residual;
  }
#endif

  *STATUS = __Nodal_Internal_Forces(Residual, ActiveNodes, ActiveDOFs, MPM_Mesh,
                                    FEM_Mesh, Is_compute_Residual,
                                    Is_compute_Reactions);
  if (*STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __Nodal_Internal_Forces()" RESET " \n");
    return Residual;
  }

  __Nodal_Traction_Forces(Residual, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
                          Is_compute_Residual, Is_compute_Reactions);

  __Nodal_Body_Forces(Residual, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
                      Is_compute_Residual, Is_compute_Reactions);

  if ((Is_compute_Residual == true) && (Is_compute_Reactions == false)) {
    __Nodal_Inertial_Forces(Residual, Effective_Mass, U_n, D_U, ActiveNodes,
                            ActiveDOFs, Params);
  }

#ifdef USE_PETSC
  VecAssemblyBegin(Residual);
  VecAssemblyEnd(Residual);
#endif
  return Residual;
}

/**************************************************************/

#ifdef USE_PETSC
static int __Nodal_Internal_Forces(Vec Residual, Mask ActiveNodes,
                                   Mask ActiveDOFs, Particle MPM_Mesh,
                                   Mesh FEM_Mesh, bool Is_compute_Residual,
                                   bool Is_compute_Reactions)
#else
static int __Nodal_Internal_Forces(double *Residual, Mask ActiveNodes,
                                   Mask ActiveDOFs, Particle MPM_Mesh,
                                   Mesh FEM_Mesh, bool Is_compute_Residual,
                                   bool Is_compute_Reactions)
#endif
{

  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned MatIndx_p;

  int Ap, Mask_node_A, Mask_total_dof_Ai, Mask_active_dof_Ai;

  double V0_p;         // Volume of the particle in the reference configuration
  double *kirchhoff_p; // Kirchhoff Stress tensor
  double *DF_p;

#if NumberDimensions == 2
  double InternalForcesDensity_Ap[2];
#else
  double InternalForcesDensity_Ap[3];
#endif
  double Residual_val;

  Element Nodes_p; /* List of nodes for particle */
  Material MatProp_p;
  Matrix d_shapefunction_n_p; /* Shape functions gradients */
  double *d_shapefunction_n1_p;
  double *d_shapefunction_n1_pA;

  for (unsigned p = 0; p < Np; p++) {

    //  Get the volume of the particle in the reference configuration
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    // Get the incremental deformation gradient
    DF_p = MPM_Mesh.Phi.DF.nM[p];

    //  Define nodal connectivity for each particle
    //  and compute gradient of the shape function
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);
    d_shapefunction_n_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    // Pushforward the shape functions
    double *d_shapefunction_n1_p = push_forward_dN__MeshTools__(
        d_shapefunction_n_p.nV, DF_p, NumNodes_p, &STATUS);
    if (d_shapefunction_n1_p == NULL) {
      fprintf(stderr, "" RED "Error in calloc()" RESET " \n");
      return EXIT_FAILURE;
    }
    STATUS = push_forward_dN__MeshTools__(
        d_shapefunction_n1_p, d_shapefunction_n_p.nV, DF_p, NumNodes_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in push_forward_dN__MeshTools__()" RESET " \n");
      return STATUS;
    }

    //  Update the Kirchhoff stress tensor with an apropiate
    //  integration rule.
    MatIndx_p = MPM_Mesh.MatIdx[p];
    MatProp_p = MPM_Mesh.Mat[MatIndx_p];
    STATUS = Stress_integration__Constitutive__(p, MPM_Mesh, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in Stress_integration__Constitutive__(,)" RESET
              " \n");
      return EXIT_FAILURE;
    }

    // Get the Kirchhoff stress tensor
    kirchhoff_p = MPM_Mesh.Phi.Stress.nM[p];

    for (unsigned A = 0; A < NumNodes_p; A++) {

      //  Compute the gradient in the reference configuration
      d_shapefunction_n1_pA = &d_shapefunction_n1_p[A * Ndim];

      //  Compute the nodal forces of the particle
      __internal_force_density(InternalForcesDensity_Ap, kirchhoff_p,
                               d_shapefunction_n1_pA);

      //  Get the node of the mesh for the contribution
      Ap = Nodes_p.Connectivity[A];
      Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

      //  Asign the nodal forces contribution to the node
      for (unsigned i = 0; i < Ndim; i++) {

        Mask_total_dof_Ai = Mask_node_A * Ndim + i;
        Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

        Residual_val = InternalForcesDensity_Ap[i] * V0_p;

        // Compute residual or reactions
        if ((Mask_active_dof_Ai != -1) && (Is_compute_Residual == true)) {
#ifdef USE_PETSC
          VecSetValues(Residual, 1, &Mask_active_dof_Ai, &Residual_val,
                       ADD_VALUES);
#else
          Residual[Mask_active_dof_Ai] += Residual_val;
#endif
        } else if ((Mask_active_dof_Ai == -1) &&
                   (Is_compute_Reactions == true)) {
#ifdef USE_PETSC
          VecSetValues(Residual, 1, &Mask_total_dof_Ai, &Residual_val,
                       ADD_VALUES);
#else
          Residual[Mask_total_dof_Ai] += Residual_val;
#endif
        }
      }
    }

    //   Free memory
    free__MatrixLib__(d_shapefunction_n_p);
    free(d_shapefunction_n1_p);
    free(Nodes_p.Connectivity);
  }

  return STATUS;
}

/**************************************************************/

static void __internal_force_density(double *InternalForcesDensity_Ap,
                                     const double *kirchhoff_p,
                                     const double *gradient_n1_pA) {
#if NumberDimensions == 2
  InternalForcesDensity_Ap[0] =
      kirchhoff_p[0] * gradient_n1_pA[0] + kirchhoff_p[1] * gradient_n1_pA[1];
  InternalForcesDensity_Ap[1] =
      kirchhoff_p[2] * gradient_n1_pA[0] + kirchhoff_p[3] * gradient_n1_pA[1];
#else
  InternalForcesDensity_Ap[0] = kirchhoff_p[0] * gradient_n1_pA[0] +
                                kirchhoff_p[1] * gradient_n1_pA[1] +
                                kirchhoff_p[2] * gradient_n1_pA[2];
  InternalForcesDensity_Ap[1] = kirchhoff_p[3] * gradient_n1_pA[0] +
                                kirchhoff_p[4] * gradient_n1_pA[1] +
                                kirchhoff_p[5] * gradient_n1_pA[2];
  InternalForcesDensity_Ap[2] = kirchhoff_p[6] * gradient_n1_pA[0] +
                                kirchhoff_p[7] * gradient_n1_pA[1] +
                                kirchhoff_p[8] * gradient_n1_pA[2];
#endif
}

/**************************************************************/

#ifdef USE_PETSC
static void __Nodal_Traction_Forces(Vec Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions)
#else
static void __Nodal_Traction_Forces(double *Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions)
#endif
{

  unsigned Ndim = NumberDimensions;
  unsigned NumContactForces = MPM_Mesh.Neumann_Contours.NumBounds;
  unsigned NumNodesLoad;

  Load Load_i;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_p;      /* Nodal values of the sahpe function */
  double N_pa;
  double A0_p;
  double Residual_val;

#if NumberDimensions == 2
  double T[2] = {0.0, 0.0};
#else
  double T[3] = {0.0, 0.0};
#endif

  unsigned p;
  int Ap, Mask_node_A, Mask_total_dof_Ai, Mask_active_dof_Ai;
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

        /*
          Pass the value of the nodal shape function to a scalar
        */
        N_pa = N_p.nV[A];

        //  Get the node of the mesh for the contribution of the contact force
        Ap = Nodes_p.Connectivity[A];
        Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

        //  Asign the nodal contact forces contribution to the node
        for (unsigned i = 0; i < Ndim; i++) {

          Mask_total_dof_Ai = Mask_node_A * Ndim + i;
          Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

          Residual_val = -N_pa * T[i] * A0_p;

          if ((Mask_active_dof_Ai != -1) && (Is_compute_Residual == true)) {
#ifdef USE_PETSC
            VecSetValues(Residual, 1, &Mask_active_dof_Ai, &Residual_val,
                         ADD_VALUES);
#else
            Residual[Mask_active_dof_Ai] += Residual_val;
#endif
          } else if ((Mask_active_dof_Ai == -1) &&
                     (Is_compute_Reactions == true)) {
#ifdef USE_PETSC
            VecSetValues(Residual, 1, &Mask_total_dof_Ai, &Residual_val,
                         ADD_VALUES);
#else
            Residual[Mask_total_dof_Ai] += Residual_val;
#endif
          }
        }
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(N_p);
      free(Nodes_p.Connectivity);
    }
  }
}

/**************************************************************/

#ifdef USE_PETSC
static void __Nodal_Body_Forces(Vec Residual, Mask ActiveNodes, Mask ActiveDOFs,
                                Particle MPM_Mesh, Mesh FEM_Mesh,
                                bool Is_compute_Residual,
                                bool Is_compute_Reactions)
#else
static void __Nodal_Body_Forces(double *Residual, Mask ActiveNodes,
                                Mask ActiveDOFs, Particle MPM_Mesh,
                                Mesh FEM_Mesh, bool Is_compute_Residual,
                                bool Is_compute_Reactions)
#endif
{

#if NumberDimensions == 2
  double b[2] = {0.0, 0.0};
#else
  double b[3] = {0.0, 0.0, 0.0};
#endif

  /* Define auxilar variables */
  unsigned Ndim = NumberDimensions;
  //  unsigned NumBodyForces = MPM_Mesh.NumberBodyForces;
  unsigned NumGP = MPM_Mesh.NumGP;
  unsigned NumNodes_p;

  int Ap, Mask_node_A, Mask_total_dof_Ai, Mask_active_dof_Ai;

  //  Load *B = MPM_Mesh.B;    /* List with the load cases */
  Element Nodes_p;         /* Element for each particle */
  Matrix ShapeFunction_p;  /* Nodal values of the sahpe function */
  double ShapeFunction_pA; /* Evaluation in the node I for the particle p */
  double m_p;              /* Mass of the particle */
  double Residual_val;

  //  b[1] = -9.81;

  for (unsigned p = 0; p < NumGP; p++) {

    /* Get the value of the mass */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* Define tributary nodes of the particle */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /* Compute shape functions */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Get the node of the mesh for the contribution */
    for (unsigned A = 0; A < NumNodes_p; A++) {

      /* Pass the value of the nodal shape function to a scalar */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /* Get the node of the mesh for the contribution */
      Ap = Nodes_p.Connectivity[A];
      Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

      /* Compute body forces */
      for (unsigned i = 0; i < Ndim; i++) {

        Mask_total_dof_Ai = Mask_node_A * Ndim + i;
        Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

        Residual_val = -ShapeFunction_pA * b[i] * m_p;

        if ((Mask_active_dof_Ai != -1) && (Is_compute_Residual == true)) {
#ifdef USE_PETSC
          VecSetValues(Residual, 1, &Mask_active_dof_Ai, &Residual_val,
                       ADD_VALUES);
#else
          Residual[Mask_active_dof_Ai] += Residual_val;
#endif
        } else if ((Mask_active_dof_Ai == -1) &&
                   (Is_compute_Reactions == true)) {
#ifdef USE_PETSC
          VecSetValues(Residual, 1, &Mask_total_dof_Ai, &Residual_val,
                       ADD_VALUES);
#else
          Residual[Mask_total_dof_Ai] += Residual_val;
#endif
        }
      }
    }

    /* Free the matrix with the nodal gradient of the element */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }
}

/**************************************************************/

#ifdef USE_PETSC
static void __Nodal_Inertial_Forces(Vec Residual, double *Mass, Nodal_Field U_n,
                                    Nodal_Field D_U, Mask ActiveNodes,
                                    Mask ActiveDOFs, Newmark_parameters Params)
#else
static void __Nodal_Inertial_Forces(double *Residual, double *Mass,
                                    Nodal_Field U_n, Nodal_Field D_U,
                                    Mask ActiveNodes, Mask ActiveDOFs,
                                    Newmark_parameters Params)
#endif
{

  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;
  int Mask_active_dof_Ai;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double Residual_val;

  for (unsigned Mask_total_dof_Ai = 0; Mask_total_dof_Ai < Ntotaldofs;
       Mask_total_dof_Ai++) {

    Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

    if (Mask_active_dof_Ai != -1) {

      for (unsigned Mask_total_dof_Bj = 0; Mask_total_dof_Bj < Ntotaldofs;
           Mask_total_dof_Bj++) {

        Residual_val =
            Mass[Mask_total_dof_Ai * Ntotaldofs + Mask_total_dof_Bj] *
            (alpha_1 * D_U.value[Mask_total_dof_Bj] -
             alpha_2 * U_n.d_value_dt[Mask_total_dof_Bj] -
             alpha_3 * U_n.d2_value_dt2[Mask_total_dof_Bj]);

#ifdef USE_PETSC
        VecSetValues(Residual, 1, &Mask_active_dof_Ai, &Residual_val,
                     ADD_VALUES);
#else
        Residual[Mask_active_dof_Ai] += Residual_val;
#endif
      }
    }
  }
}

/**************************************************************/
#ifdef USE_PETSC
static double __error_residual(const Vec Residual, unsigned Total_dof)
#else
static double __error_residual(const double *Residual, unsigned Total_dof)
#endif
{

  double Error = 0;

#ifdef USE_PETSC
  VecNorm(Residual, NORM_2, &Error);
#else
  for (unsigned A = 0; A < Total_dof; A++) {
    Error += DSQR(Residual[A]);
  }
  Error = pow(Error, 0.5);
#endif

  return Error;
}

/**************************************************************/

static void __preallocation_tangent_matrix(int *nnz, Mask ActiveNodes,
                                           Mask ActiveDOFs, Particle MPM_Mesh) {
  unsigned Ndim = NumberDimensions;
  unsigned Nactivedofs = ActiveDOFs.Nactivenodes;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;

  int *Active_dof_Mat =
      (int *)calloc(Nactivedofs * Nactivedofs, __SIZEOF_INT__);

  // Spatial discretization variables
  Element Nodes_p;
  int Ap, Mask_node_A, Mask_total_dof_Ai, Mask_active_dof_Ai;
  int Bp, Mask_node_B, Mask_total_dof_Bj, Mask_active_dof_Bj;

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

          Mask_total_dof_Ai = Mask_node_A * Ndim + i;
          Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

          for (unsigned j = 0; j < Ndim; j++) {

            Mask_total_dof_Bj = Mask_node_B * Ndim + j;
            Mask_active_dof_Bj = ActiveDOFs.Nodes2Mask[Mask_total_dof_Bj];

            if ((Mask_active_dof_Ai != -1) && (Mask_active_dof_Bj != -1)) {
              Active_dof_Mat[Mask_active_dof_Ai * Nactivedofs +
                             Mask_active_dof_Bj] = 1;
            }
          }
        }
      }
    }

    free(Nodes_p.Connectivity);
  }

  for (unsigned A = 0; A < Nactivedofs; A++) {
    for (unsigned B = 0; B < Nactivedofs; B++) {
      nnz[A] += Active_dof_Mat[A * Nactivedofs + B];
    }
  }

  free(Active_dof_Mat);
}

/**************************************************************/

static void compute_local_intertia(double *Inertia_density_p, double Na_p,
                                   double Nb_p, double m_p, double alpha_1,
                                   double epsilon, unsigned A, unsigned B) {

  double ID_AB_p =
      alpha_1 * ((1 - epsilon) * Na_p * Nb_p + (A == B) * epsilon * Na_p) * m_p;

#if NumberDimensions == 2
  Inertia_density_p[0] = ID_AB_p;
  Inertia_density_p[3] = ID_AB_p;
#else
  Inertia_density_p[0] = ID_AB_p;
  Inertia_density_p[4] = ID_AB_p;
  Inertia_density_p[8] = ID_AB_p;
#endif
}

/**************************************************************/

#ifdef USE_PETSC
static Mat __assemble_tangent_stiffness(int *nnz, Mask ActiveNodes,
                                        Mask ActiveDOFs, Particle MPM_Mesh,
                                        Mesh FEM_Mesh,
                                        Newmark_parameters Params, int *STATUS)
#else
static double *
__assemble_tangent_stiffness(int *nnz, Mask ActiveNodes, Mask ActiveDOFs,
                             Particle MPM_Mesh, Mesh FEM_Mesh,
                             Newmark_parameters Params, int *STATUS)
#endif
{

  *STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Nactivedofs = ActiveDOFs.Nactivenodes;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned MatIndx_p;

  int Ap, Mask_node_A, Mask_total_dof_Ai, Mask_active_dof_Ai;
  int Bp, Mask_node_B, Mask_total_dof_Bj, Mask_active_dof_Bj;

#ifdef USE_PETSC
  Mat Tangent_Stiffness;
  MatCreateSeqAIJ(PETSC_COMM_SELF, Nactivedofs, Nactivedofs, 0, nnz,
                  &Tangent_Stiffness);
  MatSetOption(Tangent_Stiffness, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
  PetscInt Istart, Iend;
  MatGetOwnershipRange(Tangent_Stiffness, &Istart, &Iend);
#else

  double *Tangent_Stiffness =
      (double *)calloc(Nactivedofs * Nactivedofs, __SIZEOF_DOUBLE__);
  if (Tangent_Stiffness == NULL) {
    fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
    *STATUS = EXIT_FAILURE;
    return Tangent_Stiffness;
  }
#endif

  // Spatial discretization variables
  Element Nodes_p;
  Matrix shapefunction_n_p;
  Matrix d_shapefunction_n_p;
  double *d_shapefunction_n1_p;

  double *d_shapefunction_n_pA;
  double *d_shapefunction_n_pB;

  double shapefunction_n_pA;
  double shapefunction_n_pB;

#if NumberDimensions == 2
  double Stiffness_density_p[4];
  double Inertia_density_p[4] = {0.0, 0.0, 0.0, 0.0};
#else
  double Stiffness_density_p[9];
  double Inertia_density_p[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif
  double Tangent_Stiffness_val;

  // Auxiliar pointers to tensors
  double *DF_p;

  Material MatProp_p;
  double m_p;  /* Mass of the particle */
  double V0_p; /* Volume of the particle in the reference configuration */
  double alpha_1 = Params.alpha_1;
  double alpha_4 = Params.alpha_4;
  double epsilon = Params.epsilon;

  for (unsigned p = 0; p < Np; p++) {

    // Get mass and the volume of the particle in the reference configuration
    // and the jacobian of the deformation gradient
    m_p = MPM_Mesh.Phi.mass.nV[p];
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    // Material properties of the particle
    MatIndx_p = MPM_Mesh.MatIdx[p];
    MatProp_p = MPM_Mesh.Mat[MatIndx_p];

    //
    DF_p = MPM_Mesh.Phi.DF.nM[p];

    //  Define nodal connectivity for each particle
    //  and compute gradient of the shape function
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);
    shapefunction_n_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    d_shapefunction_n_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    // Pushforward the shape function gradient
    d_shapefunction_n1_p =
        (double *)calloc(NumNodes_p * Ndim, __SIZEOF_DOUBLE__);
    if (d_shapefunction_n1_p == NULL) {
      fprintf(stderr, "" RED "Error in calloc()" RESET " \n");
      *STATUS = EXIT_FAILURE;
      return Tangent_Stiffness;
    }
    *STATUS = push_forward_dN__MeshTools__(
        d_shapefunction_n1_p, d_shapefunction_n_p.nV, DF_p, NumNodes_p);
    if (*STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in push_forward_dN__MeshTools__()" RESET " \n");
      *STATUS = EXIT_FAILURE;
      return Tangent_Stiffness;
    }

    for (unsigned A = 0; A < NumNodes_p; A++) {

      // Get the gradient evaluation in node A
      // and the masked index of the node A
      shapefunction_n_pA = shapefunction_n_p.nV[A];
      d_shapefunction_n_pA = d_shapefunction_n_p.nM[A];
      Ap = Nodes_p.Connectivity[A];
      Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

      for (unsigned B = 0; B < NumNodes_p; B++) {

        // Get the gradient evaluation in node B
        // and the masked index of the node B
        shapefunction_n_pB = shapefunction_n_p.nV[B];
        d_shapefunction_n_pB = d_shapefunction_n_p.nM[B];
        Bp = Nodes_p.Connectivity[B];
        Mask_node_B = ActiveNodes.Nodes2Mask[Bp];

        // Compute the stiffness density matrix
        *STATUS = stiffness_density__Constitutive__(
            p, Stiffness_density_p, &d_shapefunction_n1_p[A * Ndim],
            &d_shapefunction_n1_p[B * Ndim], d_shapefunction_n_pA,
            d_shapefunction_n_pB, alpha_4, MPM_Mesh, MatProp_p);
        if (*STATUS == EXIT_FAILURE) {
          fprintf(stderr,
                  "" RED "Error in stiffness_density__Constitutive__" RESET
                  "\n");
          return Tangent_Stiffness;
        }

        // Local mass matrix
        compute_local_intertia(Inertia_density_p, shapefunction_n_pA,
                               shapefunction_n_pB, m_p, alpha_1, epsilon,
                               Mask_node_A, Mask_node_B);

        //  Assembling process
        for (unsigned i = 0; i < Ndim; i++) {

          Mask_total_dof_Ai = Mask_node_A * Ndim + i;
          Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

          for (unsigned j = 0; j < Ndim; j++) {

            Mask_total_dof_Bj = Mask_node_B * Ndim + j;
            Mask_active_dof_Bj = ActiveDOFs.Nodes2Mask[Mask_total_dof_Bj];

            if ((Mask_active_dof_Ai != -1) && (Mask_active_dof_Bj != -1)) {
              Tangent_Stiffness_val = Stiffness_density_p[i * Ndim + j] * V0_p +
                                      Inertia_density_p[i * Ndim + j];
#ifdef USE_PETSC
              MatSetValues(Tangent_Stiffness, 1, &Mask_active_dof_Ai, 1,
                           &Mask_active_dof_Bj, &Tangent_Stiffness_val,
                           ADD_VALUES);
#else
              Tangent_Stiffness[Mask_active_dof_Ai * Nactivedofs +
                                Mask_active_dof_Bj] += Tangent_Stiffness_val;
#endif
            }
          }
        }
      }
    }

    free__MatrixLib__(shapefunction_n_p);
    free__MatrixLib__(d_shapefunction_n_p);
    free(d_shapefunction_n1_p);
    free(Nodes_p.Connectivity);
  }

#ifdef USE_PETSC
  MatAssemblyBegin(Tangent_Stiffness, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Tangent_Stiffness, MAT_FINAL_ASSEMBLY);
  MatSetOption(Tangent_Stiffness, MAT_SYMMETRIC, PETSC_TRUE);
#endif

  return Tangent_Stiffness;
}

/**************************************************************/

#ifdef USE_PETSC
static void __update_Nodal_Increments(const Vec Residual, Nodal_Field D_U,
                                      Nodal_Field U_n, Mask ActiveDOFs,
                                      Newmark_parameters Params,
                                      unsigned Ntotaldofs)
#else
static void __update_Nodal_Increments(const double *Residual, Nodal_Field D_U,
                                      Nodal_Field U_n, Mask ActiveDOFs,
                                      Newmark_parameters Params,
                                      unsigned Ntotaldofs)
#endif
{
  int Mask_active_dof_Ai;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

#ifdef USE_PETSC
  const PetscScalar *aux_Residual;
  VecGetArrayRead(Residual, &aux_Residual);
#endif

  /*
    Update nodal variables
  */
  for (unsigned Mask_total_dof_Ai = 0; Mask_total_dof_Ai < Ntotaldofs;
       Mask_total_dof_Ai++) {

    Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

    if (Mask_active_dof_Ai == -1) {

      D_U.d2_value_dt2[Mask_total_dof_Ai] = 0.0;

      D_U.d_value_dt[Mask_total_dof_Ai] =
          alpha_4 * D_U.value[Mask_total_dof_Ai] +
          (alpha_5 - 1) * U_n.d_value_dt[Mask_total_dof_Ai] +
          alpha_6 * U_n.d2_value_dt2[Mask_total_dof_Ai];

    } else {

#ifdef USE_PETSC
      D_U.value[Mask_total_dof_Ai] -= aux_Residual[Mask_active_dof_Ai];
#else
      D_U.value[Mask_total_dof_Ai] -= Residual[Mask_active_dof_Ai];
#endif

      D_U.d2_value_dt2[Mask_total_dof_Ai] =
          alpha_1 * D_U.value[Mask_total_dof_Ai] -
          alpha_2 * U_n.d_value_dt[Mask_total_dof_Ai] -
          (alpha_3 + 1) * U_n.d2_value_dt2[Mask_total_dof_Ai];

      D_U.d_value_dt[Mask_total_dof_Ai] =
          alpha_4 * D_U.value[Mask_total_dof_Ai] +
          (alpha_5 - 1) * U_n.d_value_dt[Mask_total_dof_Ai] +
          alpha_6 * U_n.d2_value_dt2[Mask_total_dof_Ai];
    }
  }

#ifdef USE_PETSC
  VecRestoreArrayRead(Residual, &aux_Residual);
#endif
}

/**************************************************************/

static void __update_Particles(Nodal_Field D_U, Particle MPM_Mesh,
                               Mesh FEM_Mesh, Mask ActiveNodes) {

  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;

  int Ap, A_mask;

  Matrix ShapeFunction_p;  /* Value of the shape-function in the particle */
  double ShapeFunction_pI; /* Nodal value for the particle */
  Element Nodes_p;         /* Element for each particle */
  double D_U_pI;
  double D_V_pI;
  double D_A_pI;

  /* iterate over the particles */
  for (unsigned p = 0; p < Np; p++) {

    // Update the determinant of the deformation gradient
    MPM_Mesh.Phi.J_n.nV[p] = MPM_Mesh.Phi.J_n1.nV[p];

    // Update density
    MPM_Mesh.Phi.rho.nV[p] =
        MPM_Mesh.Phi.mass.nV[p] /
        (MPM_Mesh.Phi.Vol_0.nV[p] * MPM_Mesh.Phi.J_n.nV[p]);

    // Update hardening
    MPM_Mesh.Phi.Kappa_n[p] = MPM_Mesh.Phi.Kappa_n1[p];

    // Update equivalent plastic strains
    MPM_Mesh.Phi.EPS_n[p] = MPM_Mesh.Phi.EPS_n1[p];

    // Update elastic left Cauchy-Green tensor
#if NumberDimensions == 2
    for (unsigned i = 0; i < 5; i++)
      MPM_Mesh.Phi.b_e_n.nM[p][i] = MPM_Mesh.Phi.b_e_n1.nM[p][i];
#else
    for (unsigned i = 0; i < 9; i++)
      MPM_Mesh.Phi.b_e_n.nM[p][i] = MPM_Mesh.Phi.b_e_n1.nM[p][i];
#endif

      // Update deformation gradient
#if NumberDimensions == 2
    for (unsigned i = 0; i < 5; i++)
      MPM_Mesh.Phi.F_n.nM[p][i] = MPM_Mesh.Phi.F_n1.nM[p][i];
#else
    for (unsigned i = 0; i < 9; i++)
      MPM_Mesh.Phi.F_n.nM[p][i] = MPM_Mesh.Phi.F_n1.nM[p][i];
#endif

      // Update rate of deformation gradient
#if NumberDimensions == 2
    for (unsigned i = 0; i < 5; i++)
      MPM_Mesh.Phi.dt_F_n.nM[p][i] = MPM_Mesh.Phi.dt_F_n1.nM[p][i];
#else
    for (unsigned i = 0; i < 9; i++)
      MPM_Mesh.Phi.dt_F_n.nM[p][i] = MPM_Mesh.Phi.dt_F_n1.nM[p][i];
#endif

    //  Define nodal connectivity for each particle
    //  and compute the shape function
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    for (unsigned A = 0; A < NumNodes_p; A++) {

      // Get the shape function evaluation in node A
      // and the masked index of the node A
      ShapeFunction_pI = ShapeFunction_p.nV[A];
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      //  Update acceleration, velocity and position of the particles
      for (unsigned i = 0; i < Ndim; i++) {
        D_U_pI = ShapeFunction_pI * D_U.value[A_mask * Ndim + i];
        D_V_pI = ShapeFunction_pI * D_U.d_value_dt[A_mask * Ndim + i];
        D_A_pI = ShapeFunction_pI * D_U.d2_value_dt2[A_mask * Ndim + i];

        MPM_Mesh.Phi.acc.nM[p][i] += D_A_pI;
        MPM_Mesh.Phi.vel.nM[p][i] += D_V_pI;
        MPM_Mesh.Phi.dis.nM[p][i] += D_U_pI;
        MPM_Mesh.Phi.x_GC.nM[p][i] += D_U_pI;
      }
    }

    free(Nodes_p.Connectivity);
    free__MatrixLib__(ShapeFunction_p);
  }
}

/**************************************************************/
