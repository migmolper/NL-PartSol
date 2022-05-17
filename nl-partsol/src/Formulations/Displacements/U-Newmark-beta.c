
#include "Formulations/Displacements/U-Newmark-beta.h"
#include "Formulations/Displacements/U-Newmark-beta-aux.h"

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
  Mat Tangent_Stiffness;
  Vec Lumped_Mass;
  Vec Residual;
  Vec Reactions;
#else
  double *Tangent_Stiffness;
  double *Residual;
  double *Reactions;
  double *Lumped_Mass;
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

  while (TimeStep < NumTimeStep) {

    DoProgress("Simulation:", TimeStep, NumTimeStep);

    //! Local search and compute list of active nodes and dofs
    local_search__MeshTools__(MPM_Mesh, FEM_Mesh);
    ActiveNodes = get_active_nodes__MeshTools__(FEM_Mesh);
    Nactivenodes = ActiveNodes.Nactivenodes;
    Ntotaldofs = Ndim * Nactivenodes;
    ActiveDOFs = get_active_dofs__MeshTools__(ActiveNodes, FEM_Mesh, TimeStep,
                                              NumTimeStep);
    Nactivedofs = ActiveDOFs.Nactivenodes;

    //! Define and allocate the effective mass matrix
    Lumped_Mass =
        __compute_nodal_lumped_mass(MPM_Mesh, FEM_Mesh, ActiveNodes, &STATUS);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in __compute_nodal_lumped_mass()" RESET " \n");
      return EXIT_FAILURE;
    }

    //! Get the previous converged nodal value
    U_n = __get_nodal_field_tn(Lumped_Mass, MPM_Mesh, FEM_Mesh, ActiveNodes,
                               ActiveDOFs, Time_Integration_Params, &STATUS);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __get_nodal_field_tn()" RESET " \n");
      return EXIT_FAILURE;
    }

    //! Compute kinematic nodal values
    D_U = __initialise_nodal_increments(U_n, FEM_Mesh, ActiveNodes,
                                        Time_Integration_Params, &STATUS);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in __initialise_nodal_increments()" RESET " \n");
      return EXIT_FAILURE;
    }

    //! Trial residual
    Residual = __assemble_residual(
        U_n, D_U, Lumped_Mass, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
        Time_Integration_Params, TimeStep, true, false, &STATUS);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __assemble_residual()" RESET " \n");
      return EXIT_FAILURE;
    }

    //! Compute initial error
    Error_0 = Error_i = __error_residual(Residual, Nactivedofs);
    Iter = 0;

    if (Error_0 < TOL) {
#ifdef USE_PETSC
      VecDestroy(&Residual);
#else
      free(Residual);
#endif
      Error_relative = 0.0;
    } else {
      Error_relative = Error_i / Error_0;
      Tangent_Stiffness = __preallocation_tangent_matrix(
          ActiveNodes, ActiveDOFs, MPM_Mesh, &STATUS);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __preallocation_tangent_matrix()" RESET " \n");
        return EXIT_FAILURE;
      }
    }

    //! Start Newton-Raphson
    while (Error_relative > TOL) {

      if ((Error_i < TOL * 100) || (Error_relative < TOL) || (Iter > MaxIter)) {
        break;
      }

      __assemble_tangent_stiffness(Tangent_Stiffness, ActiveNodes, ActiveDOFs,
                                   MPM_Mesh, FEM_Mesh, Time_Integration_Params,
                                   Iter, &STATUS);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __assemble_tangent_stiffness()" RESET " \n");
        return EXIT_FAILURE;
      }

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

      __update_Nodal_Increments(Residual, D_U, U_n, ActiveDOFs,
                                Time_Integration_Params, Ntotaldofs);

      __local_compatibility_conditions(D_U, ActiveNodes, MPM_Mesh, FEM_Mesh,
                                       &STATUS);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in __local_compatibility_conditions()" RESET
                " \n");
        return EXIT_FAILURE;
      }

      __constitutive_update(MPM_Mesh, FEM_Mesh, &STATUS);

      // Free memory
#ifdef USE_PETSC
      VecDestroy(&Residual);
#else
      free(Residual);
#endif

      // Compute residual (NR-loop)
      Residual = __assemble_residual(
          U_n, D_U, Lumped_Mass, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
          Time_Integration_Params, TimeStep, true, false, &STATUS);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr, "" RED "Error in __assemble_residual()" RESET " \n");
        return EXIT_FAILURE;
      }

      // Get stats for the convergence
      Error_i = __error_residual(Residual, Nactivedofs);
      Error_relative = Error_i / Error_0;
      Iter++;
    }

    //! Free residual and tangent matrix
    if (Iter > 0) {
#ifdef USE_PETSC
      VecDestroy(&Residual);
      MatDestroy(&Tangent_Stiffness);
#else
      free(Residual);
      free(Tangent_Stiffness);
#endif
    }

    print_convergence_stats(TimeStep, NumTimeStep, Iter, Error_0, Error_i,
                            Error_relative);

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
          U_n, D_U, Lumped_Mass, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
          Time_Integration_Params, TimeStep, false, true, &STATUS);
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

    //! Update time step
    TimeStep++;

#ifdef USE_PETSC
    VecDestroy(&Lumped_Mass);
    VecDestroy(&U_n.value);
    VecDestroy(&U_n.d_value_dt);
    VecDestroy(&U_n.d2_value_dt2);
    VecDestroy(&D_U.value);
    VecDestroy(&D_U.d_value_dt);
    VecDestroy(&D_U.d2_value_dt2);
#else
    free(Lumped_Mass);
    free(U_n.value);
    free(U_n.d_value_dt);
    free(U_n.d2_value_dt2);
    free(D_U.value);
    free(D_U.d_value_dt);
    free(D_U.d2_value_dt2);
#endif

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

static void __get_assembling_locations_lumped_mass(int *Mask_active_dofs_A,
                                                   int Mask_node_A) {
  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  Mask_active_dofs_A[0] = Mask_node_A * Ndim + 0;
  Mask_active_dofs_A[1] = Mask_node_A * Ndim + 1;
#else
  Mask_active_dofs_A[0] = Mask_node_A * Ndim + 0;
  Mask_active_dofs_A[1] = Mask_node_A * Ndim + 1;
  Mask_active_dofs_A[2] = Mask_node_A * Ndim + 2;
#endif
}

/**************************************************************/

#ifdef USE_PETSC
static Vec __compute_nodal_lumped_mass(Particle MPM_Mesh, Mesh FEM_Mesh,
                                       Mask ActiveNodes, int *STATUS)
#else
static double *__compute_nodal_lumped_mass(Particle MPM_Mesh, Mesh FEM_Mesh,
                                           Mask ActiveNodes, int *STATUS)
#endif
{

  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned Ntotaldofs = Ndim * Nnodes_mask;
  unsigned NumberNodes_p;
  unsigned p;

#if NumberDimensions == 2
  double Local_Mass_Matrix_p[2];
  int Mask_active_dofs_A[2];
#else
  double Local_Mass_Matrix_p[3];
  int Mask_active_dofs_A[3];
#endif

  // Define and allocate the lumped mass matrix
#ifdef USE_PETSC
  Vec Lumped_MassMatrix;
  VecCreate(PETSC_COMM_WORLD, &Lumped_MassMatrix);
  VecSetSizes(Lumped_MassMatrix, PETSC_DECIDE, Ntotaldofs);
  VecSetFromOptions(Lumped_MassMatrix);
#else
  double *Lumped_MassMatrix = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
  if (Lumped_MassMatrix == NULL) {
    fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
    *STATUS = EXIT_FAILURE;
    return Lumped_MassMatrix;
  }
#endif

#pragma omp parallel private(NumberNodes_p, Local_Mass_Matrix_p,               \
                             Mask_active_dofs_A)
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

        __get_assembling_locations_lumped_mass(Mask_active_dofs_A, Mask_node_A);

#pragma omp critical
        {
#ifdef USE_PETSC
          VecSetValues(Lumped_MassMatrix, Ndim, Mask_active_dofs_A,
                       Local_Mass_Matrix_p, ADD_VALUES);
#else
          VecSetValues(Lumped_MassMatrix, Local_Mass_Matrix_p,
                       Mask_active_dofs_A);
#endif

        } // #pragma omp critical
      }   // for A

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);

    } // for p
  }   // #pragma omp parallel

#ifdef USE_PETSC
  VecAssemblyBegin(Lumped_MassMatrix);
  VecAssemblyEnd(Lumped_MassMatrix);
#endif

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

#ifdef USE_PETSC
static Nodal_Field __get_nodal_field_tn(Vec Lumped_Mass, Particle MPM_Mesh,
                                        Mesh FEM_Mesh, Mask ActiveNodes,
                                        Mask ActiveDOFs,
                                        Newmark_parameters Params, int *STATUS)
#else
static Nodal_Field __get_nodal_field_tn(double *Lumped_Mass, Particle MPM_Mesh,
                                        Mesh FEM_Mesh, Mask ActiveNodes,
                                        Mask ActiveDOFs,
                                        Newmark_parameters Params, int *STATUS)
#endif
{
  *STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumberNodes_p;
  unsigned p;

  Nodal_Field U_n;

#ifdef USE_PETSC
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
#else
  U_n.value = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
  U_n.d_value_dt = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
  U_n.d2_value_dt2 = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
  if ((U_n.value == NULL) || (U_n.d_value_dt == NULL) ||
      (U_n.d2_value_dt2 == NULL)) {
    fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
    *STATUS = EXIT_FAILURE;
  }
#endif

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
#ifdef USE_PETSC
          VecSetValues(U_n.value, Ndim, Mask_active_dofs_A, U_N_m_IP,
                       ADD_VALUES);
          VecSetValues(U_n.d_value_dt, Ndim, Mask_active_dofs_A, V_N_m_IP,
                       ADD_VALUES);
          VecSetValues(U_n.d2_value_dt2, Ndim, Mask_active_dofs_A, A_N_m_IP,
                       ADD_VALUES);
#else
          VecSetValues(U_n.value, U_N_m_IP, Mask_active_dofs_A);
          VecSetValues(U_n.d_value_dt, V_N_m_IP, Mask_active_dofs_A);
          VecSetValues(U_n.d2_value_dt2, A_N_m_IP, Mask_active_dofs_A);
#endif
        } // #pragma omp critical
      }   // for A

      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);

    } // for p

  } // #pragma omp parallel

#ifdef USE_PETSC
  VecPointwiseDivide(U_n.value, U_n.value, Lumped_Mass);
  VecPointwiseDivide(U_n.d_value_dt, U_n.d_value_dt, Lumped_Mass);
  VecPointwiseDivide(U_n.d2_value_dt2, U_n.d2_value_dt2, Lumped_Mass);
#else
  unsigned idx;
#pragma omp for private(idx)
  for (idx = 0; idx < Ntotaldofs; idx++) {
    if (ActiveDOFs.Nodes2Mask[idx] != -1) {
      U_n.value[idx] *= 1.0 / Lumped_Mass[idx];
      U_n.d_value_dt[idx] *= 1.0 / Lumped_Mass[idx];
      U_n.d2_value_dt2[idx] *= 1.0 / Lumped_Mass[idx];
    }
  } // for idx
#endif

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

#ifdef USE_PETSC
          VecSetValues(U_n.value, 1, &Mask_restricted_dofs_A, &U_value_In1,
                       ADD_VALUES);
          VecSetValues(U_n.d_value_dt, 1, &Mask_restricted_dofs_A, &V_value_In1,
                       ADD_VALUES);
          VecSetValues(U_n.d2_value_dt2, 1, &Mask_restricted_dofs_A,
                       &A_value_In1, ADD_VALUES);
#else
          U_n.value[Mask_restricted_dofs_A] += U_value_In1;
          U_n.d_value_dt[Mask_restricted_dofs_A] += V_value_In1;
          U_n.d2_value_dt2[Mask_restricted_dofs_A] += A_value_In1;
#endif
        }
      }
    }
  }

#ifdef USE_PETSC
  VecAssemblyBegin(U_n.value);
  VecAssemblyEnd(U_n.value);
  VecAssemblyBegin(U_n.d_value_dt);
  VecAssemblyEnd(U_n.d_value_dt);
  VecAssemblyBegin(U_n.d2_value_dt2);
  VecAssemblyEnd(U_n.d2_value_dt2);
#endif

  return U_n;
}

/**************************************************************/

static Nodal_Field
__initialise_nodal_increments(Nodal_Field U_n, Mesh FEM_Mesh, Mask ActiveNodes,
                              Newmark_parameters Params, int *STATUS)
/*
  Apply the boundary conditions over the nodes
*/
{
  unsigned NumNodesBound;
  unsigned Ndim = NumberDimensions;
  unsigned NumBounds = FEM_Mesh.Bounds.NumBounds;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;
  int Id_BCC;
  int Id_BCC_mask;

  Nodal_Field D_U;

#ifdef USE_PETSC
  VecCreate(PETSC_COMM_WORLD, &D_U.value);
  VecCreate(PETSC_COMM_WORLD, &D_U.d_value_dt);
  VecCreate(PETSC_COMM_WORLD, &D_U.d2_value_dt2);
  VecSetSizes(D_U.value, PETSC_DECIDE, Ntotaldofs);
  VecSetSizes(D_U.d_value_dt, PETSC_DECIDE, Ntotaldofs);
  VecSetSizes(D_U.d2_value_dt2, PETSC_DECIDE, Ntotaldofs);
  VecSetFromOptions(D_U.value);
  VecSetFromOptions(D_U.d_value_dt);
  VecSetFromOptions(D_U.d2_value_dt2);
  const PetscScalar *Un_dt;
  VecGetArrayRead(U_n.d_value_dt, &Un_dt);
  const PetscScalar *Un_dt2;
  VecGetArrayRead(U_n.d2_value_dt2, &Un_dt2);
#else
  D_U.value = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
  D_U.d_value_dt = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
  D_U.d2_value_dt2 = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
  if ((D_U.value == NULL) || (D_U.d_value_dt == NULL) ||
      (D_U.d2_value_dt2 == NULL)) {
    fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
    *STATUS = EXIT_FAILURE;
  }
  const double *Un_dt = U_n.d_value_dt;
  const double *Un_dt2 = U_n.d2_value_dt2;
#endif

  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

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

#ifdef USE_PETSC
          VecSetValues(D_U.value, 1, &idx, &D_U_bcc, ADD_VALUES);
          VecSetValues(D_U.d_value_dt, 1, &idx, &D_U_dt_bcc, ADD_VALUES);
          VecSetValues(D_U.d2_value_dt2, 1, &idx, &D_U_dt2_bcc, ADD_VALUES);
#else
          D_U.value[idx] += D_U_bcc;
          D_U.d2_value_dt2[idx] += D_U_dt_bcc;
          D_U.d_value_dt[idx] += D_U_dt2_bcc;
#endif
        }
      }
    }
  }

#ifdef USE_PETSC
  VecAssemblyBegin(D_U.value);
  VecAssemblyEnd(D_U.value);
  VecAssemblyBegin(D_U.d_value_dt);
  VecAssemblyEnd(D_U.d_value_dt);
  VecAssemblyBegin(D_U.d2_value_dt2);
  VecAssemblyEnd(D_U.d2_value_dt2);
//  VecRestoreArrayRead(U_n.d_value_dt, &Un_dt);
//  VecRestoreArrayRead(U_n.d2_value_dt2, &Un_dt2);
#endif

  return D_U;
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

#ifdef USE_PETSC
  const double *dU;
  VecGetArrayRead(D_U.value, &dU);
  const double *dU_dt;
  VecGetArrayRead(D_U.d_value_dt, &dU_dt);
#else
  const double *dU = D_U.value;
  const double *dU_dt = D_U.d_value_dt;
#endif

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
        fprintf(stderr, "" RED "Negative jacobian in particle %i" RESET " \n",
                p);
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

#ifdef USE_PETSC
  VecRestoreArrayRead(D_U.value, &dU);
  VecRestoreArrayRead(D_U.d_value_dt, &dU_dt);
#endif
}

/**************************************************************/

static void __constitutive_update(Particle MPM_Mesh, Mesh FEM_Mesh,
                                  int *STATUS) {
  *STATUS = EXIT_SUCCESS;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned MatIndx_p;
  unsigned p;

#pragma omp for private(p, MatIndx_p)
  for (p = 0; p < Np; p++) {

    //  Update the Kirchhoff stress tensor with an apropiate
    //  integration rule.
    MatIndx_p = MPM_Mesh.MatIdx[p];
    Material MatProp_p = MPM_Mesh.Mat[MatIndx_p];
    *STATUS = Stress_integration__Constitutive__(p, MPM_Mesh, MatProp_p);
    if (*STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in Stress_integration__Constitutive__(,)" RESET
              " \n");
    }
  }
}

/**************************************************************/

#ifdef USE_PETSC
static Vec __assemble_residual(Nodal_Field U_n, Nodal_Field D_U,
                               Vec Lumped_Mass, Mask ActiveNodes,
                               Mask ActiveDOFs, Particle MPM_Mesh,
                               Mesh FEM_Mesh, Newmark_parameters Params,
                               unsigned TimeStep, bool Is_compute_Residual,
                               bool Is_compute_Reactions, int *STATUS)
#else
static double *__assemble_residual(Nodal_Field U_n, Nodal_Field D_U,
                                   double *Lumped_Mass, Mask ActiveNodes,
                                   Mask ActiveDOFs, Particle MPM_Mesh,
                                   Mesh FEM_Mesh, Newmark_parameters Params,
                                   unsigned TimeStep, bool Is_compute_Residual,
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

  __Nodal_Internal_Forces(Residual, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
                          Is_compute_Residual, Is_compute_Reactions, STATUS);
  if (*STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED "Error in __Nodal_Internal_Forces()" RESET " \n");
    return Residual;
  }

  __Nodal_Traction_Forces(Residual, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
                          Is_compute_Residual, Is_compute_Reactions);

  __Nodal_Body_Forces(Residual, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh,
                      TimeStep, Is_compute_Residual, Is_compute_Reactions);

  if ((Is_compute_Residual == true) && (Is_compute_Reactions == false)) {
    __Nodal_Inertial_Forces(Residual, Lumped_Mass, U_n, D_U, ActiveNodes,
                            ActiveDOFs, Params);
  }

#ifdef USE_PETSC
  VecAssemblyBegin(Residual);
  VecAssemblyEnd(Residual);
#endif
  return Residual;
}

/**************************************************************/

static void __get_assembling_locations_residual(int *Mask_active_dofs_A,
                                                int Mask_node_A,
                                                Mask ActiveDOFs) {
  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  Mask_active_dofs_A[0] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 0];
  Mask_active_dofs_A[1] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 1];
#else
  Mask_active_dofs_A[0] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 0];
  Mask_active_dofs_A[1] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 1];
  Mask_active_dofs_A[2] = ActiveDOFs.Nodes2Mask[Mask_node_A * Ndim + 2];
#endif
}

/**************************************************************/

static void __get_assembling_locations_reactions(int *Mask_active_dofs_A,
                                                 int Mask_node_A) {
  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  Mask_active_dofs_A[0] = Mask_node_A * Ndim + 0;
  Mask_active_dofs_A[1] = Mask_node_A * Ndim + 1;
#else
  Mask_active_dofs_A[0] = Mask_node_A * Ndim + 0;
  Mask_active_dofs_A[1] = Mask_node_A * Ndim + 1;
  Mask_active_dofs_A[2] = Mask_node_A * Ndim + 2;
#endif
}

/**************************************************************/

#ifndef USE_PETSC
static void VecSetValues(double *Residual,
                         const double *InternalForcesDensity_Ap,
                         const int *Mask_active_dofs_A) {

  unsigned Ndim = NumberDimensions;
  int Mask_active_dof_Ai;
  int Mask_active_dof_Bj;

  for (unsigned idx_A = 0; idx_A < Ndim; idx_A++) {

    Mask_active_dof_Ai = Mask_active_dofs_A[idx_A];

    if (Mask_active_dof_Ai != -1) {
      Residual[Mask_active_dof_Ai] += InternalForcesDensity_Ap[idx_A];
    }
  }
}
#endif

/**************************************************************/

#ifdef USE_PETSC
static void __Nodal_Internal_Forces(Vec Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions, int *STATUS)
#else
static void __Nodal_Internal_Forces(double *Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions, int *STATUS)
#endif
{

  *STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;

  unsigned p;

#if NumberDimensions == 2
  double InternalForcesDensity_Ap[2];
#else
  double InternalForcesDensity_Ap[3];
#endif

#ifdef USE_PETSC
  VecSetOption(Residual, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
#endif

#pragma omp parallel private(NumNodes_p, InternalForcesDensity_Ap)
  {
#pragma omp for private(p)
    for (p = 0; p < Np; p++) {

      //  Get the volume of the particle in the reference configuration
      double V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

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
        fprintf(stderr,
                "" RED "Error in push_forward_dN__MeshTools__()" RESET " \n");
      }

      // Get the Kirchhoff stress tensor pointer
      double *kirchhoff_p = MPM_Mesh.Phi.Stress.nM[p];

      for (unsigned A = 0; A < NumNodes_p; A++) {

        //! Get the gradient evaluation in node A \n
        //! and the masked index of the node A
        double *d_shapefunction_n1_pA = &d_shapefunction_n1_p[A * Ndim];
        int Ap = Nodes_p.Connectivity[A];
        int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

        //!  Compute the nodal forces of the particle
        __internal_force_density(InternalForcesDensity_Ap, kirchhoff_p,
                                 d_shapefunction_n1_pA, V0_p);

#if NumberDimensions == 2
        int Mask_active_dofs_A[2];
#else
        int Mask_active_dofs_A[3];
#endif

        if (Is_compute_Residual == true) {
          __get_assembling_locations_residual(Mask_active_dofs_A, Mask_node_A,
                                              ActiveDOFs);
        } else if (Is_compute_Reactions == true) {
          __get_assembling_locations_reactions(Mask_active_dofs_A, Mask_node_A);
        }

#pragma omp critical
        {
#ifdef USE_PETSC
          VecSetValues(Residual, Ndim, Mask_active_dofs_A,
                       InternalForcesDensity_Ap, ADD_VALUES);
#else
          VecSetValues(Residual, InternalForcesDensity_Ap, Mask_active_dofs_A);
#endif
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
        int Mask_active_dofs_A[2];
#else
        int Mask_active_dofs_A[3];
#endif

        if (Is_compute_Residual == true) {
          __get_assembling_locations_residual(Mask_active_dofs_A, Mask_node_A,
                                              ActiveDOFs);
        } else if (Is_compute_Reactions == true) {
          __get_assembling_locations_reactions(Mask_active_dofs_A, Mask_node_A);
        }

        //  Asign the nodal contact forces contribution to the node
#pragma omp critical
        {
#ifdef USE_PETSC
          VecSetValues(Residual, Ndim, Mask_active_dofs_A,
                       LocalTractionForce_Ap, ADD_VALUES);
#else
          VecSetValues(Residual, LocalTractionForce_Ap, Mask_active_dofs_A);
#endif
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

#ifdef USE_PETSC
static int __Nodal_Body_Forces(Vec Residual, Mask ActiveNodes, Mask ActiveDOFs,
                               Particle MPM_Mesh, Mesh FEM_Mesh,
                               unsigned TimeStep, bool Is_compute_Residual,
                               bool Is_compute_Reactions)
#else
static int
__Nodal_Body_Forces(double *Residual, Mask ActiveNodes, Mask ActiveDOFs,
                    Particle MPM_Mesh, Mesh FEM_Mesh, unsigned TimeStep,
                    bool Is_compute_Residual, bool Is_compute_Reactions)
#endif
{

  /* Define auxilar variables */
  unsigned Ndim = NumberDimensions;
  //  unsigned NumBodyForces = MPM_Mesh.NumberBodyForces;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned p;

  if (gravity_field.STATUS == false)
    return EXIT_SUCCESS;

#if NumberDimensions == 2
  const double b[2] = {gravity_field.Value[0].Fx[TimeStep],
                       gravity_field.Value[1].Fx[TimeStep]};
  double LocalBodyForce_Ap[2];
#else
  const double b[3] = {gravity_field.Value[0].Fx[TimeStep],
                       gravity_field.Value[1].Fx[TimeStep],
                       gravity_field.Value[2].Fx[TimeStep]};
  double LocalBodyForce_Ap[3];
#endif

#ifdef USE_PETSC
  VecSetOption(Residual, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
#endif

#pragma omp parallel private(NumNodes_p, LocalBodyForce_Ap)
  {
#pragma omp for private(p)
    for (p = 0; p < Np; p++) {

      /* Get the value of the mass */
      double m_p = MPM_Mesh.Phi.mass.nV[p];

      /* Define tributary nodes of the particle */
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      Element Nodes_p =
          nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

      /* Compute shape functions */
      Matrix ShapeFunction_p =
          compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Get the node of the mesh for the contribution */
      for (unsigned A = 0; A < NumNodes_p; A++) {

        /* Pass the value of the nodal shape function to a scalar */
        double ShapeFunction_pA = ShapeFunction_p.nV[A];

        /* Get the node of the mesh for the contribution */
        int Ap = Nodes_p.Connectivity[A];
        int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

        //! Local contribution
        __local_body_force(LocalBodyForce_Ap, b, ShapeFunction_pA, m_p);

#if NumberDimensions == 2
        int Mask_active_dofs_A[2];
#else
        int Mask_active_dofs_A[3];
#endif

        if (Is_compute_Residual == true) {
          __get_assembling_locations_residual(Mask_active_dofs_A, Mask_node_A,
                                              ActiveDOFs);
        } else if (Is_compute_Reactions == true) {
          __get_assembling_locations_reactions(Mask_active_dofs_A, Mask_node_A);
        }

        //  Asign the nodal body forces contribution to the node
#pragma omp critical
        {
#ifdef USE_PETSC
          VecSetValues(Residual, Ndim, Mask_active_dofs_A, LocalBodyForce_Ap,
                       ADD_VALUES);
#else
          VecSetValues(Residual, LocalBodyForce_Ap, Mask_active_dofs_A);
#endif
        } // #pragma omp critical
      }   // for loop unsigned A

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    } // For unsigned p
  }   // #pragma omp parallel

  return EXIT_SUCCESS;
}

/**************************************************************/

static void __local_body_force(double *LocalBodyForce_Ap, const double *b,
                               double N_n1_pA, double m_p) {
#if NumberDimensions == 2
  LocalBodyForce_Ap[0] = -N_n1_pA * b[0] * m_p;
  LocalBodyForce_Ap[1] = -N_n1_pA * b[1] * m_p;
#else
  LocalBodyForce_Ap[0] = -N_n1_pA * b[0] * m_p;
  LocalBodyForce_Ap[1] = -N_n1_pA * b[1] * m_p;
  LocalBodyForce_Ap[2] = -N_n1_pA * b[2] * m_p;
#endif
}

/**************************************************************/

#ifdef USE_PETSC
static void __Nodal_Inertial_Forces(Vec Residual, Vec Lumped_Mass,
                                    Nodal_Field U_n, Nodal_Field D_U,
                                    Mask ActiveNodes, Mask ActiveDOFs,
                                    Newmark_parameters Params)
#else
static void __Nodal_Inertial_Forces(double *Residual, double *Lumped_Mass,
                                    Nodal_Field U_n, Nodal_Field D_U,
                                    Mask ActiveNodes, Mask ActiveDOFs,
                                    Newmark_parameters Params)
#endif
{

  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;
  unsigned idx;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;

#ifdef USE_PETSC
  const PetscScalar *M_II;
  VecGetArrayRead(Lumped_Mass, &M_II);
  const PetscScalar *dU;
  VecGetArrayRead(D_U.value, &dU);
  const PetscScalar *Un_dt;
  VecGetArrayRead(U_n.d_value_dt, &Un_dt);
  const PetscScalar *Un_dt2;
  VecGetArrayRead(U_n.d2_value_dt2, &Un_dt2);
#else
  const double *M_II = Lumped_Mass;
  const double *dU = D_U.value;
  const double *Un_dt = U_n.d_value_dt;
  const double *Un_dt2 = U_n.d2_value_dt2;
#endif

#pragma omp for private(idx)
  for (idx = 0; idx < Ntotaldofs; idx++) {

#pragma omp critical
    {
      int Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[idx];

      if (Mask_active_dof_Ai != -1) {

        double R_Ai = M_II[idx] * (alpha_1 * dU[idx] - alpha_2 * Un_dt[idx] -
                                   alpha_3 * Un_dt2[idx]);
#ifdef USE_PETSC
        VecSetValues(Residual, 1, &Mask_active_dof_Ai, &R_Ai, ADD_VALUES);
#else
        Residual[Mask_active_dof_Ai] += R_Ai;
#endif
      }
    }
  }

#ifdef USE_PETSC
  VecRestoreArrayRead(Lumped_Mass, &M_II);
  VecRestoreArrayRead(D_U.value, &dU);
  VecRestoreArrayRead(U_n.d_value_dt, &Un_dt);
  VecRestoreArrayRead(U_n.d2_value_dt2, &Un_dt2);
#endif
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

#ifdef USE_PETSC
static Mat __preallocation_tangent_matrix(Mask ActiveNodes, Mask ActiveDOFs,
                                          Particle MPM_Mesh, int *STATUS)
#else
static double *__preallocation_tangent_matrix(Mask ActiveNodes, Mask ActiveDOFs,
                                              Particle MPM_Mesh, int *STATUS)
#endif

{
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
  }
#endif

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
#ifndef USE_PETSC
static void MatSetValues(double *Tangent_Stiffness,
                         const double *Tangent_Stiffness_p,
                         const int *Mask_active_dofs_A,
                         const int *Mask_active_dofs_B, int Nactivedofs) {

  unsigned Ndim = NumberDimensions;
  int Mask_active_dof_Ai;
  int Mask_active_dof_Bj;

  for (unsigned idx_A = 0; idx_A < Ndim; idx_A++) {

    Mask_active_dof_Ai = Mask_active_dofs_A[idx_A];

    for (unsigned idx_B = 0; idx_B < Ndim; idx_B++) {

      Mask_active_dof_Bj = Mask_active_dofs_B[idx_B];

      if ((Mask_active_dof_Ai != -1) && (Mask_active_dof_Bj != -1)) {
        Tangent_Stiffness[Mask_active_dof_Ai * Nactivedofs +
                          Mask_active_dof_Bj] +=
            Tangent_Stiffness_p[idx_A * Ndim + idx_B];
      }
    }
  }
}
#endif

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

#ifdef USE_PETSC
static void __assemble_tangent_stiffness(Mat Tangent_Stiffness,
                                         Mask ActiveNodes, Mask ActiveDOFs,
                                         Particle MPM_Mesh, Mesh FEM_Mesh,
                                         Newmark_parameters Params,
                                         unsigned Iter, int *STATUS)
#else
static void __assemble_tangent_stiffness(double *Tangent_Stiffness,
                                         Mask ActiveNodes, Mask ActiveDOFs,
                                         Particle MPM_Mesh, Mesh FEM_Mesh,
                                         Newmark_parameters Params,
                                         unsigned Iter, int *STATUS)
#endif
{

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
#ifdef USE_PETSC
    MatZeroEntries(Tangent_Stiffness);
#else
#pragma omp for
    for (unsigned idx = 0; idx < Nactivedofs * Nactivedofs; idx++) {
      Tangent_Stiffness[idx] = 0.0;
    }
#endif
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

      //! Get mass and the volume of the particle in the reference configuration
      //! and the jacobian of the deformation gradient
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

#ifdef USE_PETSC
            MatSetValues(Tangent_Stiffness, Ndim, Mask_active_dofs_A, Ndim,
                         Mask_active_dofs_B, Tangent_Stiffness_p, ADD_VALUES);
#else

            MatSetValues(Tangent_Stiffness, Tangent_Stiffness_p,
                         Mask_active_dofs_A, Mask_active_dofs_B, Nactivedofs);
#endif

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

#ifdef USE_PETSC
  MatAssemblyBegin(Tangent_Stiffness, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Tangent_Stiffness, MAT_FINAL_ASSEMBLY);
  MatSetOption(Tangent_Stiffness, MAT_SYMMETRIC, PETSC_TRUE);
#endif
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
  unsigned Mask_total_dof_Ai;
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
#else
  const double *aux_Residual = Residual;
  const double *Un_dt = U_n.d_value_dt;
  const double *Un_dt2 = U_n.d2_value_dt2;
  double *dU = D_U.value;
  double *dU_dt = D_U.d_value_dt;
  double *dU_dt2 = D_U.d2_value_dt2;
#endif

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
      dU[Mask_total_dof_Ai] -= aux_Residual[Mask_active_dof_Ai];

      dU_dt2[Mask_total_dof_Ai] = alpha_1 * dU[Mask_total_dof_Ai] -
                                  alpha_2 * Un_dt[Mask_total_dof_Ai] -
                                  (alpha_3 + 1) * Un_dt2[Mask_total_dof_Ai];

      dU_dt[Mask_total_dof_Ai] = alpha_4 * dU[Mask_total_dof_Ai] +
                                 (alpha_5 - 1) * Un_dt[Mask_total_dof_Ai] +
                                 alpha_6 * Un_dt2[Mask_total_dof_Ai];
    } // if Mask_active_dof_Ai == -1)
  }   // #pragma omp for private (Mask_total_dof_Ai)

#ifdef USE_PETSC
  VecRestoreArrayRead(Residual, &aux_Residual);
  VecRestoreArrayRead(U_n.d_value_dt, &Un_dt);
  VecRestoreArrayRead(U_n.d2_value_dt2, &Un_dt2);
  VecRestoreArray(D_U.value, &dU);
  VecRestoreArray(D_U.d_value_dt, &dU_dt);
  VecRestoreArray(D_U.d2_value_dt2, &dU_dt2);
#endif
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

#ifdef USE_PETSC
  const PetscScalar *dU;
  VecGetArrayRead(D_U.value, &dU);
  const PetscScalar *dU_dt;
  VecGetArrayRead(D_U.d_value_dt, &dU_dt);
  const PetscScalar *dU_dt2;
  VecGetArrayRead(D_U.d2_value_dt2, &dU_dt2);
#else
  const double *dU = D_U.value;
  const double *dU_dt = D_U.d_value_dt;
  const double *dU_dt2 = D_U.d2_value_dt2;
#endif

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

#ifdef USE_PETSC
  VecRestoreArrayRead(D_U.value, &dU);
  VecRestoreArrayRead(D_U.d_value_dt, &dU_dt);
  VecRestoreArrayRead(D_U.d2_value_dt2, &dU_dt2);
#endif
}

/**************************************************************/
