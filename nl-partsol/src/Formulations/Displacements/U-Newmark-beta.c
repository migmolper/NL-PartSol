#include "Formulations/Displacements/U-Newmark-beta.h"

typedef struct {

  Vec value;
  Vec d_value_dt;
  Vec d2_value_dt2;

} Nodal_Field;

typedef struct {

  Vec t;
  Vec t2;

} nodal_kinetics;

typedef struct {

  PetscScalar *t;
  PetscScalar *t2;

} ptr_nodal_kinetics;

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
  Particle MPM_Mesh;
  Mesh FEM_Mesh;
  Vec Lumped_Mass;
  Vec U_n;
  Vec U_n_dt;
  Vec U_n_dt2;
  Newmark_parameters Time_Integration_Params;

} Ctx;

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

/**
 * @brief Set local lagrangian values in to the global lagrangian additively
 *
 * @param lagrangian Lagrangian vector
 * @param Ndim Number of dimensions of the local lagrangian
 * @param Mask_dofs_A
 * @param local_lagrangian local lagrangian
 */
static void __set_local_values_lagragian(PetscScalar *lagrangian, unsigned Ndim,
                                         const int *Mask_dofs_A,
                                         const double *local_lagrangian);
/**************************************************************/

/*!
  \brief Compute the contribution of particle p to the lumped \n
  mass matrix correspoding to node A

  \param[out] Local_Mass_Matrix_p Local lumped Mass matrix
  \param[in] Na_p Shape function evaluation at node A
  \param[in] m_p Particle mass
*/
static void __compute_local_mass_matrix(double *Local_Mass_Matrix_p,
                                        double Na_p, double m_p);
/**************************************************************/

/*!
  \brief Returns a mask with the position of the dofs \n
  with contributions to the lumped mass matrix for node A \n
  and particle p.

  \param[out] Mask_active_dofs_A Global dofs positions for node A.
  \param[in] Mask_node_A Index of the node A with mask.
*/
static void __get_assembling_locations_lumped_mass(int *Mask_active_dofs_A,
                                                   int Mask_node_A);
/**************************************************************/

/*!
  \brief This function returns the lumped mass matrix \n
  in order to reduce storage, this matrix is presented as \n
  a vector.

  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[out] STATUS Returns failure or success
  \return The lumped mass matrix for the active system
*/
static Vec __compute_nodal_lumped_mass(Particle MPM_Mesh, Mesh FEM_Mesh,
                                       Mask ActiveNodes, int *STATUS);
/**************************************************************/

static void __local_projection_displacement(double *U_N_m_IP,
                                            const double *dis_p,
                                            double m__x__ShapeFunction_pA);
/**************************************************************/

static void __local_projection_velocity(double *V_N_m_IP, const double *vel_p,
                                        double m__x__ShapeFunction_pA);
/**************************************************************/

static void __local_projection_acceleration(double *A_N_m_IP,
                                            const double *acc_p,
                                            double m__x__ShapeFunction_pA);
/**************************************************************/

static void __local_projection_increment_displacement(
    double *DU_N_m_IP, const double *vel_p, const double *acc_p,
    double m__x__ShapeFunction_pA, double DeltaTimeStep);

/**************************************************************/

static void __get_assembling_locations_nodal_kinetics(int *Mask_active_dofs_A,
                                                      int Mask_node_A,
                                                      Mask ActiveDOFs);
/**************************************************************/

/*!
  \brief Project the kinematic information of the particles towards the nodes \n
  of the mesh to get the nodal fields at t = n
  \param[in] Lumped_Mass Lumped mass matrix used for the projection
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] Params Time integration parameters
  \param[out] STATUS Returns failure or success
  \return Nodal kinetics information from the step n
*/
static Nodal_Field __get_nodal_field_tn(Vec Lumped_Mass, Particle MPM_Mesh,
                                        Mesh FEM_Mesh, Mask ActiveNodes,
                                        Mask ActiveDOFs,
                                        Newmark_parameters Params, int *STATUS);
/**************************************************************/

/**
 * @brief
 *
 * @param Lumped_Mass
 * @param MPM_Mesh
 * @param FEM_Mesh
 * @param U_n
 * @param ActiveNodes
 * @param ActiveDOFs
 * @param Params
 * @param STATUS
 * @return Nodal_Field
 */
static Nodal_Field __trial_nodal_increments(Vec Lumped_Mass, Particle MPM_Mesh,
                                            Mesh FEM_Mesh, Nodal_Field U_n,
                                            Mask ActiveNodes, Mask ActiveDOFs,
                                            Newmark_parameters Params,
                                            int *STATUS);
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
static PetscErrorCode __Lagrangian_evaluation(SNES snes, Vec D_U, Vec Residual,
                                              void *ctx);
/**************************************************************/

/**
 * @brief
 *
 * @param dU Incremental nodal displacement field
 * @param dU_dt Incremental nodal velocity field
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @param STATUS Returns failure or success
 */
static void __local_compatibility_conditions(const PetscScalar *dU,
                                             const PetscScalar *dU_dt,
                                             Mask ActiveNodes,
                                             Particle MPM_Mesh, Mesh FEM_Mesh,
                                             int *STATUS);

/**************************************************************/

/*!
  \brief Update the stress tensor of the particle
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[out] STATUS Returns failure or success
*/
static int __constitutive_update(Particle MPM_Mesh, Mesh FEM_Mesh);
/**************************************************************/

/*!
  \brief Returns a mask with the position of the computed values.

  \param[out] Mask_active_dofs_A Global dofs positions for node A.
  \param[in] Mask_node_A Index of the node A with mask.
*/
static void __get_assembling_locations_lagrangian(int *Mask_active_dofs_A,
                                                  int Mask_node_A);

/**************************************************************/

/**
 * @brief Function used to compute the contribution of the \n
 * internal forces to the Lagrangian
 *
 * @param Lagrangian Lagrangian vector
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 * @param STATUS Returns failure or success
 */
static void __nodal_internal_forces(PetscScalar *Lagrangian, Mask ActiveNodes,
                                    Particle MPM_Mesh, Mesh FEM_Mesh,
                                    int *STATUS);

/**************************************************************/

/**
 * @brief Compute the contribution to the internal forces in the node A of
 * the particle p
 *
 * @param InternalForcesDensity_Ap contribution to the internal forces in
 * the node A of the particle p
 * @param kirchhoff_p Kirchhoff stress tensor
 * @param gradient_n1_pA Shape function gradient in the n+1 time step
 * @param V0_p Volume of the particle in the reference configuration
 */
static void __internal_force_density(double *InternalForcesDensity_Ap,
                                     const double *kirchhoff_p,
                                     const double *gradient_n1_pA, double V0_p);
/**************************************************************/

/**
 * @brief Function used to compute the contribution of the contact forces to the
 * residual
 *
 * @param Lagrangian Lagrangian vector
 * @param ActiveNodes List of nodes which takes place in the computation
 * @param MPM_Mesh Information of the particles
 * @param FEM_Mesh Information of the background nodes
 */
static void __nodal_traction_forces(PetscScalar *Lagrangian, Mask ActiveNodes,
                                    Particle MPM_Mesh, Mesh FEM_Mesh);
/**************************************************************/

/**
 * @brief Compute the contribution to the local traction forces in the node A of
 * the particle p
 *
 * @param LocalTractionForce_Ap contribution to the local traction forces in
 * the node A of the particle p
 * @param T Local traction force
 * @param N_n1_pA Shape function evaluation at node A
 * @param A0_p Reference area of the particle p
 */
static void __local_traction_force(double *LocalTractionForce_Ap,
                                   const double *T, double N_n1_pA,
                                   double A0_p);
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
 * @param Params Time integration parameters
 */
static void __nodal_inertial_forces(PetscScalar *Lagrangian,
                                    const PetscScalar *M_II,
                                    const PetscScalar *dU,
                                    const PetscScalar *Un_dt,
                                    const PetscScalar *Un_dt2, Mask ActiveNodes,
                                    Newmark_parameters Params);
/**************************************************************/

/*!
  \brief Compute the 2-norm of the residual
  \param[in] Residual
  \return The norm of the residual
*/
static double __error_residual(const Vec Residual);
/**************************************************************/

/*!

*/
static Mat __preallocation_tangent_matrix(Mask ActiveNodes, Mask ActiveDOFs,
                                          Particle MPM_Mesh, int *STATUS);
/**************************************************************/

static void compute_local_intertia(double *Inertia_density_p /**< */,
                                   double Na_p /**< */, double Nb_p /**< */,
                                   double m_p /**< */, double alpha_1 /**< */,
                                   double epsilon /**< */, unsigned A /**< */,
                                   unsigned B /**< */);
/**************************************************************/

/*!
  \brief Returns a mask with the position of the stifness in the global tangent
  \n matrix for an eeficient assembly process.

  \param[out] Mask_active_dofs_A Global dofs positions for node A.
  \param[in] Mask_node_A Index of the node A with mask.
  \param[out] Mask_active_dofs_B Global dofs positions for node B
  \param[in] Mask_node_B Index of the node B with mask
  \param[in] ActiveDOFs List of dofs which takes place in the computation
*/
static void __get_tangent_matrix_assembling_locations(int *Mask_active_dofs_A,
                                                      int Mask_node_A,
                                                      int *Mask_active_dofs_B,
                                                      int Mask_node_B,
                                                      Mask ActiveDOFs);
/**************************************************************/

/*!
  \brief Computes the local tangent stiffness as the addition of the \n
  local tangent stiffness density and the intertia density matrix

  \param[out] Tangent_Stiffness_p Local tangent stiffness
  \param[in] Stiffness_density_p Local stiffness density
  \param[in] Inertia_density_p Inertia density matrix, a.k.a mass matrix
  \param[in] V0_p Particle volume
*/
static void local_tangent_stiffness(double *Tangent_Stiffness_p,
                                    const double *Stiffness_density_p,
                                    const double *Inertia_density_p,
                                    double V0_p);
/**************************************************************/

/*!
  \param[in,out] Tangent_Stiffness The tangent matrix for the problem
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] Params Time integration parameters
  \param[in] Iter Current iteration of the solver
  \param[out] STATUS Returns failure or success
*/
static void __assemble_tangent_stiffness(Mat Tangent_Stiffness,
                                         Mask ActiveNodes, Mask ActiveDOFs,
                                         Particle MPM_Mesh, Mesh FEM_Mesh,
                                         Newmark_parameters Params,
                                         unsigned Iter, int *STATUS);
/**************************************************************/

/**
 * @brief Evaluates Jacobian matrix.
 *
 * @param snes the SNES context
 * @param[in] DU Input vector
 * @param[in,out] jac Jacobian matrix
 * @param[out] B Optionally different preconditioning matrix
 * @param dummy User-defined context
 * @return PetscErrorCode
 */
PetscErrorCode __Jacobian_evaluation(SNES snes, Vec DU, Mat jac, Mat B,
                                     void *dummy);
/**************************************************************/

/**
 * @brief
 *
 * @param D_U Nodal incremented (displacement,velocity,acceleration)
 * @param U_n Previously converged nodal kinetics
 * @param ActiveDOFs List of dofs which takes place in the computation
 * @param Params Time integration parameters
 * @param Ntotaldofs Number of dofs
 */
static void __trial_Nodal_Increments(Nodal_Field D_U, Nodal_Field U_n,
                                     Mask ActiveDOFs, Newmark_parameters Params,
                                     unsigned Ntotaldofs);
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

/**
 * @brief This function takes the nodal increments comming from \n
 * the linear solver and computes the updated nodal increments \n
 * for the nodal kinetics
 *
 * @param dU Incremental nodal displacement field
 * @param dU_dt Incremental nodal velocity field
 * @param dU_dt2 Incremental nodal acceleration field
 * @param Un_dt Nodal field of velocity at t = n
 * @param Un_dt2 Nodal field of acceleration at t = n
 * @param ActiveDOFs List of dofs which takes place in the computation
 * @param Params Time integration parameters
 * @param Ntotaldofs Number of dofs
 */
static void __compute_nodal_kinetics_increments(
    const PetscScalar *dU, PetscScalar *dU_dt, PetscScalar *dU_dt2,
    const PetscScalar *Un_dt, const PetscScalar *Un_dt2, Mask ActiveDOFs,
    Newmark_parameters Params, unsigned Ntotaldofs);
/**************************************************************/

/**
 * @brief
 *
 * @param D_U
 * @param MPM_Mesh
 * @param FEM_Mesh
 * @param ActiveNodes
 */
static void __update_Particles(Nodal_Field D_U, Particle MPM_Mesh,
                               Mesh FEM_Mesh, Mask ActiveNodes);
/**************************************************************/

// Global variables
unsigned InitialStep;
unsigned NumTimeStep;
unsigned TimeStep;
bool Use_explicit_trial;

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
  Use_explicit_trial = Parameters_Solver.Use_explicit_trial;

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

  Mat Tangent_Stiffness;
  Vec Lumped_Mass;
  Vec Residual;
  Vec Reactions;

  Vec U_n;
  Vec U_n_dt;
  Vec U_n_dt2;

  Vec DU;
  Vec DU_trial;
  Vec DU_dt;
  Vec DU_dt2;

  Mask ActiveNodes;
  Mask ActiveDOFs;

  Newmark_parameters Time_Integration_Params;

  // Define variables for the non-linear solver
  SNES snes;
  KSP ksp;
  PC pc;

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
    Nactivedofs = ActiveDOFs.Nactivenodes;

    if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
      compute_Beps__Constitutive__(MPM_Mesh, FEM_Mesh, false);
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Define and allocate the effective mass matrix
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Lumped_Mass =
        __compute_nodal_lumped_mass(MPM_Mesh, FEM_Mesh, ActiveNodes, &STATUS);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in __compute_nodal_lumped_mass()" RESET " \n");
      return EXIT_FAILURE;
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Get the previous converged nodal value
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    //    U_n = __get_nodal_field_tn(Lumped_Mass, MPM_Mesh, FEM_Mesh,
    //    ActiveNodes,
    //                               ActiveDOFs, Time_Integration_Params,
    //                               &STATUS);
    //    if (STATUS == EXIT_FAILURE) {
    //      fprintf(stderr, "" RED "Error in __get_nodal_field_tn()" RESET "
    //      \n"); return EXIT_FAILURE;
    //    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create nonlinear solver context
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    //    PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
    //    PetscCall(SNESSetType(snes, SNESNEWTONLS));
    //    PetscCall(SNESSetOptionsPrefix(snes, "mysolver_"));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create matrix and vector data structures; set corresponding routines
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
       Create vectors for solution and nonlinear function
    */
    //    PetscCall(VecCreate(PETSC_COMM_WORLD, &DU));
    //    PetscCall(VecSetSizes(DU, PETSC_DECIDE, Ntotaldofs));
    //    PetscCall(VecSetFromOptions(DU));
    //    PetscCall(VecDuplicate(DU, &Residual));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create Jacobian matrix data structure
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    //    PetscCall(MatCreate(PETSC_COMM_WORLD, &Tangent_Stiffness));
    //    PetscCall(MatSetSizes(Tangent_Stiffness, PETSC_DECIDE, PETSC_DECIDE,
    //                          Ntotaldofs, Ntotaldofs));
    //    PetscCall(MatSetFromOptions(Tangent_Stiffness));
    //    PetscCall(MatSetUp(Tangent_Stiffness));

    //    Tangent_Stiffness = __preallocation_tangent_matrix(ActiveNodes,
    //    ActiveDOFs,
    //                                                       MPM_Mesh, &STATUS);
    //    if (STATUS == EXIT_FAILURE) {
    //      fprintf(stderr,
    //              "" RED "Error in __preallocation_tangent_matrix()" RESET "
    //              \n");
    //      return EXIT_FAILURE;
    //    }

    /*
     Set function evaluation routine and vector.
    */
    //    PetscCall(SNESSetFunction(snes, Residual, __Lagrangian_evaluation,
    //    NULL));

    /*
     Set Jacobian matrix data structure and Jacobian evaluation routine
    */
    //    PetscCall(SNESSetJacobian(snes, Tangent_Stiffness, Tangent_Stiffness,
    //                              __Jacobian_evaluation, NULL));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Customize nonlinear solver; set runtime options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
       Set linear solver defaults for this problem. By extracting the
       KSP and PC contexts from the SNES context, we can then
       directly call any KSP and PC routines to set various options.
    */
    //    PetscCall(SNESGetKSP(snes, &ksp));
    //    PetscCall(KSPGetPC(ksp, &pc));
    //    PetscCall(PCSetType(pc, PCNONE));
    //    PetscCall(KSPSetTolerances(ksp, 1.e-4, PETSC_DEFAULT, PETSC_DEFAULT,
    //    20));

    /*
       Set SNES/KSP/KSP/PC runtime options, e.g.,
           -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
       These options will override those specified above as long as
       SNESSetFromOptions() is called _after_ any other customization
       routines.
    */
    //    PetscCall(SNESSetFromOptions(snes));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Evaluate initial guess; then solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    //    if (Use_explicit_trial == true) {
    //
    //      __trial_nodal_increments(DU, Lumped_Mass, MPM_Mesh, FEM_Mesh, U_n,
    //                               ActiveNodes, ActiveDOFs,
    //                               Time_Integration_Params, &STATUS);
    //      if (STATUS == EXIT_FAILURE) {
    //        fprintf(stderr,
    //                "" RED "Error in __trial_nodal_increments()" RESET " \n");
    //        return EXIT_FAILURE;
    //      }
    //    }
    //    else{
    //      PetscCall(VecZeroEntries(DU));
    //    }

    /*
       Note: The user should initialize the vector, x, with the initial guess
       for the nonlinear solver prior to calling SNESSolve().  In particular,
       to employ an initial guess of zero, the user should explicitly set
       this vector to zero by calling VecSet().
    */
    //    PetscCall(SNESSolve(snes, NULL, DU));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Free work space.  All PETSc objects should be destroyed when they
       are no longer needed.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    //    PetscCall(MatDestroy(&Tangent_Stiffness));
    //    PetscCall(SNESDestroy(&snes));

    //    print_convergence_stats(TimeStep, NumTimeStep, Iter, MaxIter, Error_0,
    //                            Error_i, Error_relative);

    //    if (Iter > MaxIter) {
    //      fprintf(
    //          stderr,
    //          "" RED
    //          "Convergence not reached in the maximum number of iterations"
    //          RESET " \n");
    //    }

    //    __update_Particles(D_U, MPM_Mesh, FEM_Mesh, ActiveNodes);

    //     particle_results_vtk__InOutFun__(MPM_Mesh, TimeStep,
    //     ResultsTimeStep);

    //      VecDestroy(&Reactions);
    //    }

    //! Update time step
    TimeStep++;

    PetscCall(VecDestroy(&Lumped_Mass));
    //    PetscCall(VecDestroy(&DU));
    //    PetscCall(VecDestroy(&Residual));

    //    VecDestroy(&U_n.value);
    //    VecDestroy(&U_n.d_value_dt);
    //    VecDestroy(&U_n.d2_value_dt2);
    //    VecDestroy(&D_U.value);
    //    VecDestroy(&D_U.d_value_dt);
    //    VecDestroy(&D_U.d2_value_dt2);

    free(ActiveNodes.Nodes2Mask);
    free(ActiveDOFs.Nodes2Mask);
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

static void __compute_local_mass_matrix(double *Local_Mass_Matrix_p,
                                        double Na_p, double m_p) {

  double M_AB_p = Na_p * m_p;

#if NumberDimensions == 2
  Local_Mass_Matrix_p[0] = M_AB_p;
  Local_Mass_Matrix_p[1] = M_AB_p;
#else
  Local_Mass_Matrix_p[0] = M_AB_p;
  Local_Mass_Matrix_p[1] = M_AB_p;
  Local_Mass_Matrix_p[2] = M_AB_p;
#endif
}

/**************************************************************/

static void __get_assembling_locations_lumped_mass(int *Mask_dofs_A,
                                                   int Mask_node_A) {
  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  Mask_dofs_A[0] = Mask_node_A * Ndim + 0;
  Mask_dofs_A[1] = Mask_node_A * Ndim + 1;
#else
  Mask_dofs_A[0] = Mask_node_A * Ndim + 0;
  Mask_dofs_A[1] = Mask_node_A * Ndim + 1;
  Mask_dofs_A[2] = Mask_node_A * Ndim + 2;
#endif
}

/**************************************************************/

static Vec __compute_nodal_lumped_mass(Particle MPM_Mesh, Mesh FEM_Mesh,
                                       Mask ActiveNodes, int *STATUS) {

  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned Ntotaldofs = Ndim * Nnodes_mask;
  unsigned NumberNodes_p;
  unsigned p;

#if NumberDimensions == 2
  double Local_Mass_Matrix_p[2];
  int Mask_dofs_A[2];
#else
  double Local_Mass_Matrix_p[3];
  int Mask_dofs_A[3];
#endif

  // Define and allocate the lumped mass matrix
  Vec Lumped_MassMatrix;
  VecCreate(PETSC_COMM_WORLD, &Lumped_MassMatrix);
  VecSetSizes(Lumped_MassMatrix, PETSC_DECIDE, Ntotaldofs);
  VecSetFromOptions(Lumped_MassMatrix);

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

        __compute_local_mass_matrix(Local_Mass_Matrix_p, Na_p, m_p);

        __get_assembling_locations_lumped_mass(Mask_dofs_A, Mask_node_A);

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

  VecAssemblyBegin(Lumped_MassMatrix);
  VecAssemblyEnd(Lumped_MassMatrix);

  return Lumped_MassMatrix;
}

/**************************************************************/

static void __local_projection_displacement(double *U_N_m_IP,
                                            const double *dis_p,
                                            double m__x__ShapeFunction_pA) {
#if NumberDimensions == 2
  U_N_m_IP[0] = m__x__ShapeFunction_pA * dis_p[0];
  U_N_m_IP[1] = m__x__ShapeFunction_pA * dis_p[1];
#else
  U_N_m_IP[0] = m__x__ShapeFunction_pA * dis_p[0];
  U_N_m_IP[1] = m__x__ShapeFunction_pA * dis_p[1];
  U_N_m_IP[2] = m__x__ShapeFunction_pA * dis_p[2];
#endif
}

/**************************************************************/

static void __local_projection_velocity(double *V_N_m_IP, const double *vel_p,
                                        double m__x__ShapeFunction_pA) {
#if NumberDimensions == 2
  V_N_m_IP[0] = m__x__ShapeFunction_pA * vel_p[0];
  V_N_m_IP[1] = m__x__ShapeFunction_pA * vel_p[1];
#else
  V_N_m_IP[0] = m__x__ShapeFunction_pA * vel_p[0];
  V_N_m_IP[1] = m__x__ShapeFunction_pA * vel_p[1];
  V_N_m_IP[2] = m__x__ShapeFunction_pA * vel_p[2];
#endif
}

/**************************************************************/

static void __local_projection_acceleration(double *A_N_m_IP,
                                            const double *acc_p,
                                            double m__x__ShapeFunction_pA) {
#if NumberDimensions == 2
  A_N_m_IP[0] = m__x__ShapeFunction_pA * acc_p[0];
  A_N_m_IP[1] = m__x__ShapeFunction_pA * acc_p[1];
#else
  A_N_m_IP[0] = m__x__ShapeFunction_pA * acc_p[0];
  A_N_m_IP[1] = m__x__ShapeFunction_pA * acc_p[1];
  A_N_m_IP[2] = m__x__ShapeFunction_pA * acc_p[2];
#endif
}

/**************************************************************/

static void __local_projection_increment_displacement(
    double *DU_N_m_IP, const double *vel_p, const double *acc_p,
    double m__x__ShapeFunction_pA, double DeltaTimeStep) {
#if NumberDimensions == 2
  DU_N_m_IP[0] =
      m__x__ShapeFunction_pA *
      (DeltaTimeStep * vel_p[0] + 0.5 * DSQR(DeltaTimeStep) * acc_p[0]);
  DU_N_m_IP[1] =
      m__x__ShapeFunction_pA *
      (DeltaTimeStep * vel_p[1] + 0.5 * DSQR(DeltaTimeStep) * acc_p[1]);
#else
  DU_N_m_IP[0] =
      m__x__ShapeFunction_pA *
      (DeltaTimeStep * vel_p[0] + 0.5 * DSQR(DeltaTimeStep) * acc_p[0]);
  DU_N_m_IP[1] =
      m__x__ShapeFunction_pA *
      (DeltaTimeStep * vel_p[1] + 0.5 * DSQR(DeltaTimeStep) * acc_p[1]);
  DU_N_m_IP[2] =
      m__x__ShapeFunction_pA *
      (DeltaTimeStep * vel_p[2] + 0.5 * DSQR(DeltaTimeStep) * acc_p[2]);
#endif
}

/**************************************************************/

static void __get_assembling_locations_nodal_kinetics(int *Mask_active_dofs_A,
                                                      int Mask_node_A,
                                                      Mask ActiveDOFs) {
  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  Mask_active_dofs_A[0] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 0] == -1
                              ? -1
                              : Mask_node_A * Ndim + 0;
  Mask_active_dofs_A[1] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 1] == -1
                              ? -1
                              : Mask_node_A * Ndim + 1;
#else
  Mask_active_dofs_A[0] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 0] == -1
                              ? -1
                              : Mask_node_A * Ndim + 0;
  Mask_active_dofs_A[1] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 1] == -1
                              ? -1
                              : Mask_node_A * Ndim + 1;
  Mask_active_dofs_A[2] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 2] == -1
                              ? -1
                              : Mask_node_A * Ndim + 2;
#endif
}

/**************************************************************/

static Nodal_Field __get_nodal_field_tn(Vec Lumped_Mass, Particle MPM_Mesh,
                                        Mesh FEM_Mesh, Mask ActiveNodes,
                                        Mask ActiveDOFs,
                                        Newmark_parameters Params,
                                        int *STATUS) {
  *STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumberNodes_p;
  unsigned p;

  Nodal_Field U_n;

  VecCreate(PETSC_COMM_WORLD, &U_n.value);
  VecCreate(PETSC_COMM_WORLD, &U_n.d_value_dt);
  VecCreate(PETSC_COMM_WORLD, &U_n.d2_value_dt2);
  VecSetSizes(U_n.value, PETSC_DECIDE, Ntotaldofs);
  VecSetSizes(U_n.d_value_dt, PETSC_DECIDE, Ntotaldofs);
  VecSetSizes(U_n.d2_value_dt2, PETSC_DECIDE, Ntotaldofs);
  VecSetFromOptions(U_n.value);
  VecSetFromOptions(U_n.d_value_dt);
  VecSetFromOptions(U_n.d2_value_dt2);
  VecSetOption(U_n.value, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetOption(U_n.d_value_dt, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetOption(U_n.d2_value_dt2, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

#if NumberDimensions == 2
  int Mask_active_dofs_A[2];
  double U_N_m_IP[2];
  double V_N_m_IP[2];
  double A_N_m_IP[2];
#else
  int Mask_active_dofs_A[3];
  double U_N_m_IP[3];
  double V_N_m_IP[3];
  double A_N_m_IP[3];
#endif

#pragma omp parallel private(NumberNodes_p)
  {

#pragma omp for private(p, U_N_m_IP, V_N_m_IP, A_N_m_IP, Mask_active_dofs_A)
    for (p = 0; p < Np; p++) {

      /* Define element of the particle */
      NumberNodes_p = MPM_Mesh.NumberNodes[p];
      Element Nodes_p =
          nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumberNodes_p);
      Matrix ShapeFunction_p =
          compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      //! Get the mass of the GP
      double m_p = MPM_Mesh.Phi.mass.nV[p];
      const double *dis_p = &MPM_Mesh.Phi.dis.nV[p * Ndim];
      const double *vel_p = &MPM_Mesh.Phi.vel.nV[p * Ndim];
      const double *acc_p = &MPM_Mesh.Phi.acc.nV[p * Ndim];

      for (unsigned A = 0; A < NumberNodes_p; A++) {

        //  Get the node in the nodal momentum with the mask
        unsigned Ap = Nodes_p.Connectivity[A];
        int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];
        double ShapeFunction_pA = ShapeFunction_p.nV[A];
        double m__x__ShapeFunction_pA = m_p * ShapeFunction_pA;

        __local_projection_displacement(U_N_m_IP, dis_p,
                                        m__x__ShapeFunction_pA);
        __local_projection_velocity(V_N_m_IP, vel_p, m__x__ShapeFunction_pA);
        __local_projection_acceleration(A_N_m_IP, acc_p,
                                        m__x__ShapeFunction_pA);

        __get_assembling_locations_nodal_kinetics(Mask_active_dofs_A,
                                                  Mask_node_A, ActiveDOFs);

#pragma omp critical
        {
          VecSetValues(U_n.value, Ndim, Mask_active_dofs_A, U_N_m_IP,
                       ADD_VALUES);
          VecSetValues(U_n.d_value_dt, Ndim, Mask_active_dofs_A, V_N_m_IP,
                       ADD_VALUES);
          VecSetValues(U_n.d2_value_dt2, Ndim, Mask_active_dofs_A, A_N_m_IP,
                       ADD_VALUES);
        } // #pragma omp critical
      }   // for A

      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);

    } // for p

  } // #pragma omp parallel

  VecPointwiseDivide(U_n.value, U_n.value, Lumped_Mass);
  VecPointwiseDivide(U_n.d_value_dt, U_n.d_value_dt, Lumped_Mass);
  VecPointwiseDivide(U_n.d2_value_dt2, U_n.d2_value_dt2, Lumped_Mass);

  //! Apply boundary condition
  unsigned NumBounds = FEM_Mesh.Bounds.NumBounds;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

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

          for (unsigned t = 0; t < TimeStep; t++) {

            V_value_In = V_value_In1;
            A_value_In = A_value_In1;

            D_U_value_It = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[t];
            U_value_In1 += D_U_value_It;
            V_value_In1 = alpha_4 * D_U_value_It + (alpha_5 - 1) * V_value_In +
                          alpha_6 * A_value_In;
            A_value_In1 = alpha_1 * D_U_value_It - alpha_2 * V_value_In -
                          (alpha_3 + 1) * A_value_In;
          }

          VecSetValues(U_n.value, 1, &Mask_restricted_dofs_A, &U_value_In1,
                       ADD_VALUES);
          VecSetValues(U_n.d_value_dt, 1, &Mask_restricted_dofs_A, &V_value_In1,
                       ADD_VALUES);
          VecSetValues(U_n.d2_value_dt2, 1, &Mask_restricted_dofs_A,
                       &A_value_In1, ADD_VALUES);
        }
      }
    }
  }

  VecAssemblyBegin(U_n.value);
  VecAssemblyEnd(U_n.value);
  VecAssemblyBegin(U_n.d_value_dt);
  VecAssemblyEnd(U_n.d_value_dt);
  VecAssemblyBegin(U_n.d2_value_dt2);
  VecAssemblyEnd(U_n.d2_value_dt2);

  return U_n;
}

/**************************************************************/

static Nodal_Field __trial_nodal_increments(Vec Lumped_Mass, Particle MPM_Mesh,
                                            Mesh FEM_Mesh, Nodal_Field U_n,
                                            Mask ActiveNodes, Mask ActiveDOFs,
                                            Newmark_parameters Params,
                                            int *STATUS) {
  unsigned NumNodesBound;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumberNodes_p;
  unsigned NumBounds = FEM_Mesh.Bounds.NumBounds;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;
  unsigned p;
  int Id_BCC;
  int Id_BCC_mask;

  Nodal_Field D_U;

  VecCreate(PETSC_COMM_WORLD, &D_U.value);
  VecCreate(PETSC_COMM_WORLD, &D_U.d_value_dt);
  VecCreate(PETSC_COMM_WORLD, &D_U.d2_value_dt2);
  VecSetSizes(D_U.value, PETSC_DECIDE, Ntotaldofs);
  VecSetSizes(D_U.d_value_dt, PETSC_DECIDE, Ntotaldofs);
  VecSetSizes(D_U.d2_value_dt2, PETSC_DECIDE, Ntotaldofs);
  VecSetFromOptions(D_U.value);
  VecSetFromOptions(D_U.d_value_dt);
  VecSetFromOptions(D_U.d2_value_dt2);
  VecSetOption(D_U.value, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  const PetscScalar *Un_dt;
  VecGetArrayRead(U_n.d_value_dt, &Un_dt);
  const PetscScalar *Un_dt2;
  VecGetArrayRead(U_n.d2_value_dt2, &Un_dt2);

  const double alpha_1 = Params.alpha_1;
  const double alpha_2 = Params.alpha_2;
  const double alpha_3 = Params.alpha_3;
  const double alpha_4 = Params.alpha_4;
  const double alpha_5 = Params.alpha_5;
  const double alpha_6 = Params.alpha_6;
  const double DeltaTimeStep = Params.DeltaTimeStep;

  /**
   * Initialize the values of the nodal increment using
   * an approximation based on the explicit-predictor, else
   * keep the vector as zero.
   */
  if (Use_explicit_trial == true) {

#if NumberDimensions == 2
    int Mask_active_dofs_A[2];
    double DU_N_m_IP[2];
#else
    int Mask_active_dofs_A[3];
    double DU_N_m_IP[3];
#endif

#pragma omp parallel private(NumberNodes_p)
    {

#pragma omp for private(p, DU_N_m_IP, Mask_active_dofs_A)
      for (p = 0; p < Np; p++) {

        /* Define element of the particle */
        NumberNodes_p = MPM_Mesh.NumberNodes[p];
        Element Nodes_p =
            nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumberNodes_p);
        Matrix ShapeFunction_p =
            compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

        //! Get the mass, velocity and acceleration of the GP
        double m_p = MPM_Mesh.Phi.mass.nV[p];
        const double *vel_p = &MPM_Mesh.Phi.vel.nV[p * Ndim];
        const double *acc_p = &MPM_Mesh.Phi.acc.nV[p * Ndim];

        for (unsigned A = 0; A < NumberNodes_p; A++) {

          //  Get the node in the nodal momentum with the mask
          unsigned Ap = Nodes_p.Connectivity[A];
          int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];
          double ShapeFunction_pA = ShapeFunction_p.nV[A];
          double m__x__ShapeFunction_pA = m_p * ShapeFunction_pA;

          __local_projection_increment_displacement(
              DU_N_m_IP, vel_p, acc_p, m__x__ShapeFunction_pA, DeltaTimeStep);

          __get_assembling_locations_nodal_kinetics(Mask_active_dofs_A,
                                                    Mask_node_A, ActiveDOFs);

#pragma omp critical
          {
            VecSetValues(D_U.value, Ndim, Mask_active_dofs_A, DU_N_m_IP,
                         ADD_VALUES);
          } // #pragma omp critical
        }   // for A

        free__MatrixLib__(ShapeFunction_p);
        free(Nodes_p.Connectivity);

      } // for p

    } // #pragma omp parallel

    VecPointwiseDivide(D_U.value, D_U.value, Lumped_Mass);

  } // Use_explicit_trial

  double D_U_value_It;

  //! Loop over the the boundaries to set boundary conditions
  for (unsigned i = 0; i < NumBounds; i++) {

    //! Get the number of nodes of this boundary
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

    for (unsigned j = 0; j < NumNodesBound; j++) {

      //!
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

      if (Id_BCC_mask == -1) {
        continue;
      }

      //!
      for (unsigned k = 0; k < Ndim; k++) {

        if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + TimeStep] == 1) {

          int idx = Id_BCC_mask * Ndim + k;
          double D_U_bcc = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep];
          double D_U_dt_bcc = alpha_1 * D_U_bcc - alpha_2 * Un_dt[idx] -
                              (alpha_3 + 1) * Un_dt2[idx];
          double D_U_dt2_bcc = alpha_4 * D_U_bcc + (alpha_5 - 1) * Un_dt[idx] +
                               alpha_6 * Un_dt2[idx];

          VecSetValues(D_U.value, 1, &idx, &D_U_bcc, ADD_VALUES);
          VecSetValues(D_U.d_value_dt, 1, &idx, &D_U_dt_bcc, ADD_VALUES);
          VecSetValues(D_U.d2_value_dt2, 1, &idx, &D_U_dt2_bcc, ADD_VALUES);
        }
      }
    }
  }

  VecAssemblyBegin(D_U.value);
  VecAssemblyEnd(D_U.value);
  VecAssemblyBegin(D_U.d_value_dt);
  VecAssemblyEnd(D_U.d_value_dt);
  VecAssemblyBegin(D_U.d2_value_dt2);
  VecAssemblyEnd(D_U.d2_value_dt2);
  VecRestoreArrayRead(U_n.d_value_dt, &Un_dt);
  VecRestoreArrayRead(U_n.d2_value_dt2, &Un_dt2);

  return D_U;
}

/**************************************************************/

static void __trial_Nodal_Increments(Nodal_Field D_U, Nodal_Field U_n,
                                     Mask ActiveDOFs, Newmark_parameters Params,
                                     unsigned Ntotaldofs) {
  unsigned Mask_total_dof_Ai;
  int Mask_active_dof_Ai;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

  const PetscScalar *Un_dt;
  VecGetArrayRead(U_n.d_value_dt, &Un_dt);
  const PetscScalar *Un_dt2;
  VecGetArrayRead(U_n.d2_value_dt2, &Un_dt2);
  PetscScalar *dU;
  VecGetArray(D_U.value, &dU);
  PetscScalar *dU_dt;
  VecGetArray(D_U.d_value_dt, &dU_dt);
  PetscScalar *dU_dt2;
  VecGetArray(D_U.d2_value_dt2, &dU_dt2);

#pragma omp for private(Mask_total_dof_Ai, Mask_active_dof_Ai)
  for (Mask_total_dof_Ai = 0; Mask_total_dof_Ai < Ntotaldofs;
       Mask_total_dof_Ai++) {

    Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

    if (Mask_active_dof_Ai == -1) {

      dU_dt2[Mask_total_dof_Ai] = 0.0;

      dU_dt[Mask_total_dof_Ai] = alpha_4 * dU[Mask_total_dof_Ai] +
                                 (alpha_5 - 1) * Un_dt[Mask_total_dof_Ai] +
                                 alpha_6 * Un_dt2[Mask_total_dof_Ai];

    } else {
      dU_dt2[Mask_total_dof_Ai] = alpha_1 * dU[Mask_total_dof_Ai] -
                                  alpha_2 * Un_dt[Mask_total_dof_Ai] -
                                  (alpha_3 + 1) * Un_dt2[Mask_total_dof_Ai];

      dU_dt[Mask_total_dof_Ai] = alpha_4 * dU[Mask_total_dof_Ai] +
                                 (alpha_5 - 1) * Un_dt[Mask_total_dof_Ai] +
                                 alpha_6 * Un_dt2[Mask_total_dof_Ai];
    } // if Mask_active_dof_Ai == -1)
  }   // #pragma omp for private (Mask_total_dof_Ai)

  VecRestoreArrayRead(U_n.d_value_dt, &Un_dt);
  VecRestoreArrayRead(U_n.d2_value_dt2, &Un_dt2);
  VecRestoreArray(D_U.value, &dU);
  VecRestoreArray(D_U.d_value_dt, &dU_dt);
  VecRestoreArray(D_U.d2_value_dt2, &dU_dt2);
}

/**************************************************************/

static PetscErrorCode __Lagrangian_evaluation(SNES snes, Vec dU, Vec Lagrangian,
                                              void *ctx) {

  int STATUS = EXIT_SUCCESS;

  /**
   * Read variables from user-defined structure
   * ctx
   */
  Mask ActiveDOFs = ((Ctx *)ctx)->ActiveDOFs;
  Mask ActiveNodes = ((Ctx *)ctx)->ActiveNodes;
  Particle MPM_Mesh = ((Ctx *)ctx)->MPM_Mesh;
  Mesh FEM_Mesh = ((Ctx *)ctx)->FEM_Mesh;
  Vec Lumped_Mass = ((Ctx *)ctx)->Lumped_Mass;
  Vec Un = ((Ctx *)ctx)->U_n;
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
  PetscScalar *dU_dt_ptr;
  PetscScalar *dU_dt2_ptr;
  const PetscScalar *Un_dt_ptr;
  const PetscScalar *Un_dt2_ptr;

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
  PetscCall(VecGetArray(dU_dt, &dU_dt_ptr));
  PetscCall(VecGetArray(dU_dt2, &dU_dt2_ptr));
  PetscCall(VecGetArrayRead(Un_dt, &Un_dt_ptr));
  PetscCall(VecGetArrayRead(Un_dt2, &Un_dt2_ptr));

  ptr_nodal_kinetics ptr_Un_d;
  ptr_Un_d.t = Un_dt_ptr;
  ptr_Un_d.t2 = Un_dt2_ptr;

  __compute_nodal_kinetics_increments(dU_ptr, dU_dt_ptr, dU_dt2_ptr, Un_dt_ptr,
                                      Un_dt2_ptr, ActiveDOFs,
                                      Time_Integration_Params, Ntotaldofs);

  __local_compatibility_conditions(dU_ptr, dU_dt_ptr, ActiveNodes, MPM_Mesh,
                                   FEM_Mesh, &STATUS);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,
            "" RED "Error in __local_compatibility_conditions()" RESET " \n");
    return EXIT_FAILURE;
  }

  STATUS = __constitutive_update(MPM_Mesh, FEM_Mesh);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __constitutive_update()" RESET " \n");
    return EXIT_FAILURE;
  }

  __nodal_internal_forces(Lagrangian_ptr, ActiveNodes, MPM_Mesh, FEM_Mesh,
                          STATUS);
  if (*STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __nodal_internal_forces()" RESET " \n");
    return EXIT_FAILURE;
  }

  __nodal_traction_forces(Lagrangian_ptr, ActiveNodes, MPM_Mesh, FEM_Mesh);

  __nodal_inertial_forces(Lagrangian_ptr, Lumped_Mass_ptr, dU_ptr, Un_dt_ptr,
                          Un_dt2_ptr, ActiveNodes, Params);

  /**
   * Restore vectors
   *
   */
  PetscCall(VecRestoreArray(Lagrangian, &Lagrangian_ptr));
  PetscCall(VecRestoreArrayRead(Lumped_Mass, &Lumped_Mass_ptr));
  PetscCall(VecRestoreArrayRead(dU, &dU_ptr));
  PetscCall(VecRestoreArray(dU_dt, &dU_dt_ptr));
  PetscCall(VecRestoreArray(dU_dt2, &dU_dt2_ptr));
  PetscCall(VecRestoreArrayRead(Un_dt, &Un_dt_ptr));
  PetscCall(VecRestoreArrayRead(Un_dt2, &Un_dt2_ptr));

  return EXIT_SUCCESS;
}

/**************************************************************/

static void __local_compatibility_conditions(const PetscScalar *dU,
                                             const PetscScalar *dU_dt,
                                             Mask ActiveNodes,
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
        *STATUS = EXIT_FAILURE;
      }

      get_set_field__MeshTools__(D_Displacement_Ap, dU, Nodes_p, ActiveNodes);
      get_set_field__MeshTools__(D_Velocity_Ap, dU_dt, Nodes_p, ActiveNodes);

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
      update_rate_Deformation_Gradient_n1__Particles__(dFdt_n1_p, dt_DF_p,
                                                       F_n_p, DF_p, dFdt_n_p);

      //  Compute Jacobian of the deformation gradient
      MPM_Mesh.Phi.J_n1.nV[p] = I3__TensorLib__(F_n1_p);
      if (MPM_Mesh.Phi.J_n1.nV[p] <= 0.0) {
        fprintf(stderr,
                "" RED "Negative jacobian in particle %i: %e" RESET " \n", p,
                MPM_Mesh.Phi.J_n1.nV[p]);
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

#pragma omp for private(p, Vn_patch, Vn1_patch, J_patch, Idx_Element_p,        \
                        Idx_Patch_p)
      for (p = 0; p < Np; p++) {

        Idx_Element_p = MPM_Mesh.Element_p[p];
        Idx_Patch_p = FEM_Mesh.Idx_Patch[Idx_Element_p];

        Vn_patch = FEM_Mesh.Vol_Patch_n[Idx_Patch_p];
        Vn1_patch = FEM_Mesh.Vol_Patch_n1[Idx_Patch_p];
        J_patch = Vn1_patch / Vn_patch;

        *STATUS = get_locking_free_Deformation_Gradient_n1__Particles__(
            p, J_patch, MPM_Mesh);
        if (*STATUS == EXIT_FAILURE) {
          fprintf(
              stderr,
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

static int __constitutive_update(Particle MPM_Mesh, Mesh FEM_Mesh) {
  int STATUS = EXIT_SUCCESS;
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

static void __get_assembling_locations_lagrangian(int *Mask_dofs_A,
                                                  int Mask_node_A) {
  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  Mask_dofs_A[0] = Mask_node_A * Ndim + 0;
  Mask_dofs_A[1] = Mask_node_A * Ndim + 1;
#else
  Mask_dofs_A[0] = Mask_node_A * Ndim + 0;
  Mask_dofs_A[1] = Mask_node_A * Ndim + 1;
  Mask_dofs_A[2] = Mask_node_A * Ndim + 2;
#endif
}

/**************************************************************/

static void __set_local_values_lagragian(PetscScalar *lagrangian, unsigned Ndim,
                                         const int *Mask_dofs_A,
                                         const double *local_lagrangian) {

  int Mask_dof_Ai;

  for (unsigned idx_A = 0; idx_A < Ndim; idx_A++) {

    Mask_dof_Ai = Mask_dofs_A[idx_A];

    lagrangian[Mask_dof_Ai] += local_lagrangian[idx_A];
  }
}

/**************************************************************/

static void __nodal_internal_forces(PetscScalar *Lagrangian, Mask ActiveNodes,
                                    Particle MPM_Mesh, Mesh FEM_Mesh,
                                    int *STATUS) {

  *STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned p;

  double DeltaX = FEM_Mesh.DeltaX;

#if NumberDimensions == 2
  double InternalForcesDensity_Ap[2];
#else
  double InternalForcesDensity_Ap[3];
#endif

  const double *Damage_field_n = MPM_Mesh.Phi.Damage_n;
  double *Damage_field_n1 = MPM_Mesh.Phi.Damage_n1;
  const double *Strain_Energy_field = MPM_Mesh.Phi.W;
  const double *J_n1 = MPM_Mesh.Phi.J_n1.nV;
  const double *Vol_0 = MPM_Mesh.Phi.Vol_0.nV;

#pragma omp parallel private(NumNodes_p, InternalForcesDensity_Ap)
  {
#pragma omp for private(p)
    for (p = 0; p < Np; p++) {

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
          d_shapefunction_n_p.nV, DF_p, NumNodes_p, STATUS);
      if (*STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED " Error in " RESET "" BOLDRED
                        "push_forward_dN__MeshTools__() " RESET " \n");
      }

      // Get the Kirchhoff stress tensor pointer
      double *kirchhoff_p = MPM_Mesh.Phi.Stress.nM[p];

      // Compute damage parameter (eigenerosion/eigensoftening)
      if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {

        *STATUS = compute_damage__Constitutive__(p, MPM_Mesh, FEM_Mesh.DeltaX);
        if (*STATUS == EXIT_FAILURE) {
          fprintf(stderr, "" RED " Error in " RESET "" BOLDRED
                          "compute_damage__Constitutive__() " RESET " \n");
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

        //! Compute the nodal forces of the particle
        __internal_force_density(InternalForcesDensity_Ap, kirchhoff_p,
                                 d_shapefunction_n1_pA, V0_p);

#if NumberDimensions == 2
        int Mask_dofs_A[2];
#else
        int Mask_dofs_A[3];
#endif

        __get_assembling_locations_lagrangian(Mask_dofs_A, Mask_node_A);

#pragma omp critical
        {
          __set_local_values_lagragian(Lagrangian, Ndim, Mask_dofs_A,
                                       InternalForcesDensity_Ap);
        } // #pragma omp critical
      }   // for unsigned A

      //   Free memory
      free__MatrixLib__(d_shapefunction_n_p);
      free(d_shapefunction_n1_p);
      free(Nodes_p.Connectivity);
    } // For unsigned p
  }   // #pragma omp parallel
}

/**************************************************************/

static void __internal_force_density(double *InternalForcesDensity_Ap,
                                     const double *kirchhoff_p,
                                     const double *gradient_n1_pA,
                                     double V0_p) {
#if NumberDimensions == 2
  InternalForcesDensity_Ap[0] = (kirchhoff_p[0] * gradient_n1_pA[0] +
                                 kirchhoff_p[1] * gradient_n1_pA[1]) *
                                V0_p;
  InternalForcesDensity_Ap[1] = (kirchhoff_p[2] * gradient_n1_pA[0] +
                                 kirchhoff_p[3] * gradient_n1_pA[1]) *
                                V0_p;
#else
  InternalForcesDensity_Ap[0] =
      (kirchhoff_p[0] * gradient_n1_pA[0] + kirchhoff_p[1] * gradient_n1_pA[1] +
       kirchhoff_p[2] * gradient_n1_pA[2]) *
      V0_p;
  InternalForcesDensity_Ap[1] =
      (kirchhoff_p[3] * gradient_n1_pA[0] + kirchhoff_p[4] * gradient_n1_pA[1] +
       kirchhoff_p[5] * gradient_n1_pA[2]) *
      V0_p;
  InternalForcesDensity_Ap[2] =
      (kirchhoff_p[6] * gradient_n1_pA[0] + kirchhoff_p[7] * gradient_n1_pA[1] +
       kirchhoff_p[8] * gradient_n1_pA[2]) *
      V0_p;
#endif
}

/**************************************************************/

static void __nodal_traction_forces(PetscScalar *Lagrangian, Mask ActiveNodes,
                                    Particle MPM_Mesh, Mesh FEM_Mesh) {

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
#else
  double LocalTractionForce_Ap[3];
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

        //! Get the shape function evaluation in node A \n
        //! and the masked index of the node A
        N_pa = N_p.nV[A];
        Ap = Nodes_p.Connectivity[A];
        Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

        //! Local contribution
        __local_traction_force(LocalTractionForce_Ap, T, N_pa, A0_p);

#if NumberDimensions == 2
        int Mask_dofs_A[2];
#else
        int Mask_dofs_A[3];
#endif

        __get_assembling_locations_lagrangian(Mask_dofs_A, Mask_node_A);

        //  Asign the nodal contact forces contribution to the node
#pragma omp critical
        {
          __set_local_values_lagragian(Lagrangian, Ndim, Mask_dofs_A,
                                       LocalTractionForce_Ap);
        } // #pragma omp critical

      } // for unsigned A

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(N_p);
      free(Nodes_p.Connectivity);
    }
  }
}

/**************************************************************/

static void __local_traction_force(double *LocalTractionForce_Ap,
                                   const double *T, double N_n1_pA,
                                   double A0_p) {
#if NumberDimensions == 2
  LocalTractionForce_Ap[0] = -N_n1_pA * T[0] * A0_p;
  LocalTractionForce_Ap[1] = -N_n1_pA * T[1] * A0_p;
#else
  LocalTractionForce_Ap[0] = -N_n1_pA * T[0] * A0_p;
  LocalTractionForce_Ap[1] = -N_n1_pA * T[1] * A0_p;
  LocalTractionForce_Ap[2] = -N_n1_pA * T[2] * A0_p;
#endif
}

/**************************************************************/

static void __nodal_inertial_forces(PetscScalar *Lagrangian,
                                    const PetscScalar *M_II,
                                    const PetscScalar *dU,
                                    const PetscScalar *Un_dt,
                                    const PetscScalar *Un_dt2, Mask ActiveNodes,
                                    Newmark_parameters Params) {

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
#if NumberDimensions == 2
    b[0] = gravity_field.Value[0].Fx[TimeStep];
    b[1] = gravity_field.Value[1].Fx[TimeStep];
#else
    b[0] = gravity_field.Value[0].Fx[TimeStep];
    b[1] = gravity_field.Value[1].Fx[TimeStep];
    b[2] = gravity_field.Value[2].Fx[TimeStep];
#endif
  }

#pragma omp for private(idx)
  for (idx = 0; idx < Ntotaldofs; idx++) {

#pragma omp critical
    {
      double R_Ai;

      R_Ai = M_II[idx] * (alpha_1 * dU[idx] - alpha_2 * Un_dt[idx] -
                          alpha_3 * Un_dt2[idx] - b[idx % Ndim]);

      __set_local_values_lagragian(Lagrangian, 1, &idx, &R_Ai);
    }
  }
}

/**************************************************************/
static double __error_residual(const Vec Residual) {

  double Error = 0;

  VecNorm(Residual, NORM_2, &Error);

  return Error;
}

/**************************************************************/
static Mat __preallocation_tangent_matrix(Mask ActiveNodes, Mask ActiveDOFs,
                                          Particle MPM_Mesh, int *STATUS) {
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

  int *nnz = (int *)calloc(Nactivedofs, __SIZEOF_INT__);

  for (unsigned A = 0; A < Nactivedofs; A++) {
    for (unsigned B = 0; B < Nactivedofs; B++) {
      nnz[A] += Active_dof_Mat[A * Nactivedofs + B];
    }
  }

  free(Active_dof_Mat);

  Mat Tangent_Stiffness;
  MatCreateSeqAIJ(PETSC_COMM_SELF, Nactivedofs, Nactivedofs, 0, nnz,
                  &Tangent_Stiffness);
  MatSetOption(Tangent_Stiffness, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
  PetscInt Istart, Iend;
  MatGetOwnershipRange(Tangent_Stiffness, &Istart, &Iend);

  free(nnz);

  return Tangent_Stiffness;
}

/**************************************************************/

static void compute_local_intertia(double *Inertia_density_p, double Na_p,
                                   double Nb_p, double m_p, double alpha_1,
                                   double epsilon, unsigned A, unsigned B) {

  double ID_AB_p =
      alpha_1 * ((1 - epsilon) * Na_p * Nb_p + (A == B) * epsilon * Na_p) * m_p;

#if NumberDimensions == 2
  Inertia_density_p[0] = ID_AB_p;
  Inertia_density_p[1] = ID_AB_p;
#else
  Inertia_density_p[0] = ID_AB_p;
  Inertia_density_p[1] = ID_AB_p;
  Inertia_density_p[2] = ID_AB_p;
#endif
}

/**************************************************************/

static void local_tangent_stiffness(double *Tangent_Stiffness_p,
                                    const double *Stiffness_density_p,
                                    const double *Inertia_density_p,
                                    double V0_p) {
#if NumberDimensions == 2
  Tangent_Stiffness_p[0] = Stiffness_density_p[0] * V0_p + Inertia_density_p[0];
  Tangent_Stiffness_p[1] = Stiffness_density_p[1] * V0_p;
  Tangent_Stiffness_p[2] = Stiffness_density_p[2] * V0_p;
  Tangent_Stiffness_p[3] = Stiffness_density_p[3] * V0_p + Inertia_density_p[1];
#else
  Tangent_Stiffness_p[0] = Stiffness_density_p[0] * V0_p + Inertia_density_p[0];
  Tangent_Stiffness_p[1] = Stiffness_density_p[1] * V0_p;
  Tangent_Stiffness_p[2] = Stiffness_density_p[2] * V0_p;
  Tangent_Stiffness_p[3] = Stiffness_density_p[3] * V0_p;
  Tangent_Stiffness_p[4] = Stiffness_density_p[4] * V0_p + Inertia_density_p[1];
  Tangent_Stiffness_p[5] = Stiffness_density_p[5] * V0_p;
  Tangent_Stiffness_p[6] = Stiffness_density_p[6] * V0_p;
  Tangent_Stiffness_p[7] = Stiffness_density_p[7] * V0_p;
  Tangent_Stiffness_p[8] = Stiffness_density_p[8] * V0_p + Inertia_density_p[2];
#endif
}

/**************************************************************/

static void __get_tangent_matrix_assembling_locations(int *Mask_active_dofs_A,
                                                      int Mask_node_A,
                                                      int *Mask_active_dofs_B,
                                                      int Mask_node_B,
                                                      Mask ActiveDOFs) {
  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  Mask_active_dofs_A[0] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 0];
  Mask_active_dofs_A[1] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 1];

  Mask_active_dofs_B[0] = ActiveDOFs.Nodes2Mask[Mask_node_B * Ndim + 0];
  Mask_active_dofs_B[1] = ActiveDOFs.Nodes2Mask[Mask_node_B * Ndim + 1];
#else
  Mask_active_dofs_A[0] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 0];
  Mask_active_dofs_A[1] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 1];
  Mask_active_dofs_A[2] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 2];

  Mask_active_dofs_B[0] = ActiveDOFs.Nodes2Mask[Mask_node_B * Ndim + 0];
  Mask_active_dofs_B[1] = ActiveDOFs.Nodes2Mask[Mask_node_B * Ndim + 1];
  Mask_active_dofs_B[2] = ActiveDOFs.Nodes2Mask[Mask_node_B * Ndim + 2];
#endif
}

/**************************************************************/

static void __assemble_tangent_stiffness(Mat Tangent_Stiffness,
                                         Mask ActiveNodes, Mask ActiveDOFs,
                                         Particle MPM_Mesh, Mesh FEM_Mesh,
                                         Newmark_parameters Params,
                                         unsigned Iter, int *STATUS) {

  *STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Nactivedofs = ActiveDOFs.Nactivenodes;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned MatIndx_p;

  unsigned p;

#if NumberDimensions == 2
  double Stiffness_density_p[4];
  double Inertia_density_p[2];
  double Tangent_Stiffness_p[4];
  int Mask_active_dofs_A[2];
  int Mask_active_dofs_B[2];
#else
  double Stiffness_density_p[9];
  double Inertia_density_p[3];
  double Tangent_Stiffness_p[9];
  int Mask_active_dofs_A[3];
  int Mask_active_dofs_B[3];
#endif

  // Set tangent matrix to zero if it is necessary
  if (Iter > 0) {
    MatZeroEntries(Tangent_Stiffness);
  }

  // Time integartion parameters
  double alpha_1 = Params.alpha_1;
  double alpha_4 = Params.alpha_4;
  double epsilon = Params.epsilon;

#pragma omp parallel private(NumNodes_p, MatIndx_p, Stiffness_density_p,       \
                             Inertia_density_p, Tangent_Stiffness_p,           \
                             Mask_active_dofs_A, Mask_active_dofs_B)
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
          d_shapefunction_n_p.nV, DF_p, NumNodes_p, STATUS);
      if (*STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in push_forward_dN__MeshTools__()" RESET " \n");
        *STATUS = EXIT_FAILURE;
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
          *STATUS = stiffness_density__Constitutive__(
              p, Stiffness_density_p, &d_shapefunction_n1_p[A * Ndim],
              &d_shapefunction_n1_p[B * Ndim], d_shapefunction_n_pA,
              d_shapefunction_n_pB, alpha_4, MPM_Mesh, MatProp_p);
          if (*STATUS == EXIT_FAILURE) {
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

          compute_local_intertia(Inertia_density_p, shapefunction_n_pA,
                                 shapefunction_n_pB, m_p, alpha_1, epsilon,
                                 Mask_node_A, Mask_node_B);

          local_tangent_stiffness(Tangent_Stiffness_p, Stiffness_density_p,
                                  Inertia_density_p, V0_p);

          __get_tangent_matrix_assembling_locations(
              Mask_active_dofs_A, Mask_node_A, Mask_active_dofs_B, Mask_node_B,
              ActiveDOFs);

#pragma omp critical
          {

            MatSetValues(Tangent_Stiffness, Ndim, Mask_active_dofs_A, Ndim,
                         Mask_active_dofs_B, Tangent_Stiffness_p, ADD_VALUES);

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

  MatAssemblyBegin(Tangent_Stiffness, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Tangent_Stiffness, MAT_FINAL_ASSEMBLY);
  MatSetOption(Tangent_Stiffness, MAT_SYMMETRIC, PETSC_TRUE);
}

/**************************************************************/

/*
   FormJacobian1 - Evaluates Jacobian matrix.

   Input Parameters:
.  snes - the SNES context
.  x - input vector
.  dummy - optional user-defined context (not used here)

   Output Parameters:
.  jac - Jacobian matrix
.  B - optionally different preconditioning matrix
.  flag - flag indicating matrix structure
*/
PetscErrorCode __Jacobian_evaluation(SNES snes, Vec x, Mat jac, Mat B,
                                     void *dummy) {
  const PetscScalar *xx;
  PetscScalar A[4];
  PetscInt idx[2] = {0, 1};

  /*
     Get pointer to vector data
  */
  PetscCall(VecGetArrayRead(x, &xx));

  /*
     Compute Jacobian entries and insert into matrix.
      - Since this is such a small problem, we set all entries for
        the matrix at once.
  */
  A[0] = 2.0 * xx[0] + xx[1];
  A[1] = xx[0];
  A[2] = xx[1];
  A[3] = xx[0] + 2.0 * xx[1];
  PetscCall(MatSetValues(B, 2, idx, 2, idx, A, INSERT_VALUES));

  /*
     Restore vector
  */
  PetscCall(VecRestoreArrayRead(x, &xx));

  /*
     Assemble matrix
  */
  PetscCall(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));

  __assemble_tangent_stiffness(Tangent_Stiffness, ActiveNodes, ActiveDOFs,
                               MPM_Mesh, FEM_Mesh, Time_Integration_Params,
                               Iter, &STATUS);

  if (jac != B) {
    PetscCall(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
  }
  return 0;
}

/**************************************************************/

static Vec __compute_nodal_velocity_increments(const PetscScalar *dU,
                                               const PetscScalar *Un_dt,
                                               const PetscScalar *Un_dt2,
                                               Newmark_parameters Params,
                                               PetscInt Ntotaldofs) {

  PetscInt Mask_total_dof_Ai;
  PetscScalar alpha_1 = Params.alpha_1;
  PetscScalar alpha_2 = Params.alpha_2;
  PetscScalar alpha_3 = Params.alpha_3;
  PetscScalar alpha_4 = Params.alpha_4;
  PetscScalar alpha_5 = Params.alpha_5;
  PetscScalar alpha_6 = Params.alpha_6;

  Vec dU_dt;
  VecCreate(PETSC_COMM_WORLD, &dU_dt);
  VecSetSizes(dU_dt, PETSC_DECIDE, Ntotaldofs);
  VecSetFromOptions(dU_dt);

#pragma omp for private(Mask_total_dof_Ai)
  for (Mask_total_dof_Ai = 0; Mask_total_dof_Ai < Ntotaldofs;
       Mask_total_dof_Ai++) {

    PetscScalar dU_dt_value;
    dU_dt_value = alpha_4 * dU[Mask_total_dof_Ai] +
                  (alpha_5 - 1) * Un_dt[Mask_total_dof_Ai] +
                  alpha_6 * Un_dt2[Mask_total_dof_Ai];

    VecSetValues(dU_dt, 1, &Mask_total_dof_Ai, &dU_dt_value, INSERT_VALUES);

  } // #pragma omp for private (Mask_total_dof_Ai)

  /**
   * Finalize assembling process
   * for the kinetic vector
   */
  VecAssemblyBegin(dU_dt);
  VecAssemblyEnd(dU_dt);

  return dU_d;
}

/**************************************************************/

static Vec __compute_nodal_acceleration_increments(
    const PetscScalar *dU, const ptr_nodal_kinetics Un_d, Mask ActiveDOFs,
    Newmark_parameters Params, unsigned Ntotaldofs) {

  unsigned Mask_total_dof_Ai;
  int Mask_active_dof_Ai;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

  nodal_kinetics dU_d;
  VecCreate(PETSC_COMM_WORLD, &dU_d.t);
  VecCreate(PETSC_COMM_WORLD, &dU_d.t2);
  VecSetSizes(dU_d.t, PETSC_DECIDE, Ntotaldofs);
  VecSetSizes(dU_d.t2, PETSC_DECIDE, Ntotaldofs);
  VecSetFromOptions(dU_d.t);
  VecSetFromOptions(dU_d.t2);

#pragma omp for private(Mask_total_dof_Ai, Mask_active_dof_Ai)
  for (Mask_total_dof_Ai = 0; Mask_total_dof_Ai < Ntotaldofs;
       Mask_total_dof_Ai++) {

    Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

    if (Mask_active_dof_Ai == -1) {

      dU_dt2[Mask_total_dof_Ai] = 0.0;

      dU_dt[Mask_total_dof_Ai] = alpha_4 * dU[Mask_total_dof_Ai] +
                                 (alpha_5 - 1) * Un_d.t[Mask_total_dof_Ai] +
                                 alpha_6 * Un_d.t2[Mask_total_dof_Ai];

    } else {

      dU_dt2[Mask_total_dof_Ai] = alpha_1 * dU[Mask_total_dof_Ai] -
                                  alpha_2 * Un_d.t[Mask_total_dof_Ai] -
                                  (alpha_3 + 1) * Un_d.t2[Mask_total_dof_Ai];

      dU_dt[Mask_total_dof_Ai] = alpha_4 * dU[Mask_total_dof_Ai] +
                                 (alpha_5 - 1) * Un_d.t[Mask_total_dof_Ai] +
                                 alpha_6 * Un_d.t2[Mask_total_dof_Ai];
    } // if Mask_active_dof_Ai == -1)
  }   // #pragma omp for private (Mask_total_dof_Ai)

  /**
   * Finalize assembling process
   * for the kinetic vector
   */
  VecAssemblyBegin(dU_d.t);
  VecAssemblyEnd(dU_d.t);
  VecAssemblyBegin(dU_d.t2);
  VecAssemblyEnd(dU_d.t2);

  return dU_d;
}

/**************************************************************/

static void __update_Particles(Nodal_Field D_U, Particle MPM_Mesh,
                               Mesh FEM_Mesh, Mask ActiveNodes) {

  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned p;

  double D_U_pI;
  double D_V_pI;
  double D_A_pI;

  const PetscScalar *dU;
  VecGetArrayRead(D_U.value, &dU);
  const PetscScalar *dU_dt;
  VecGetArrayRead(D_U.d_value_dt, &dU_dt);
  const PetscScalar *dU_dt2;
  VecGetArrayRead(D_U.d2_value_dt2, &dU_dt2);

#pragma omp parallel private(NumNodes_p, D_U_pI, D_V_pI, D_A_pI)
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

        //! Update rate of deformation gradient
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
          D_U_pI = ShapeFunction_pI * dU[A_mask * Ndim + i];
          D_V_pI = ShapeFunction_pI * dU_dt[A_mask * Ndim + i];
          D_A_pI = ShapeFunction_pI * dU_dt2[A_mask * Ndim + i];

          MPM_Mesh.Phi.acc.nM[p][i] += D_A_pI;
          MPM_Mesh.Phi.vel.nM[p][i] += D_V_pI;
          MPM_Mesh.Phi.dis.nM[p][i] += D_U_pI;
          MPM_Mesh.Phi.x_GC.nM[p][i] += D_U_pI;
        }
      } // for unsigned A

      free(Nodes_p.Connectivity);
      free__MatrixLib__(ShapeFunction_p);

    } // #pragma omp for private(p)
  }   // #pragma omp parallel

  VecRestoreArrayRead(D_U.value, &dU);
  VecRestoreArrayRead(D_U.d_value_dt, &dU_dt);
  VecRestoreArrayRead(D_U.d2_value_dt2, &dU_dt2);
}

/**************************************************************/
