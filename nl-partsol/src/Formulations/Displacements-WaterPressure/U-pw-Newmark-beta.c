/**
 * @file U-pw-Newmark-beta.c
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2022-05-15
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Formulations/Displacements-WaterPressure/U-pw-Newmark-beta.h"

/**************************************************************/

/**
 * @brief Finite strains Newmark-beta parameters
 *
 * @param[in] beta First Newmark-beta parameter
 * @param[in] gamma Second Newmark-beta parameter
 * @param[in] DeltaTimeStep Timestep increment
 * @return Newmark_parameters
 */
static Newmark_parameters
__compute_Newmark_parameters(double beta, double gamma, double DeltaTimeStep);
/**************************************************************/

static Matrix __compute_nodal_lumped_mass(Particle MPM_Mesh, Mesh FEM_Mesh,
                                          Mask ActiveNodes, double epsilon);
/**************************************************************/

static Nodal_Field __get_nodal_field_tn(Particle, Mesh, Mask);
/**************************************************************/

static Nodal_Field __initialise_nodal_incrementss(Nodal_Field, Mask, Mesh,
                                                  Newmark_parameters, int, int);
/**************************************************************/

static Matrix __assemble_residual(Nodal_Field, Nodal_Field, Mask, Particle,
                                  Mesh, Newmark_parameters, int, int);
/**************************************************************/

static void compute_nominal_traction_and_fluid_flux(Matrix, Mask, Particle,
                                                    Mesh, int, int);
/**************************************************************/

static Matrix __assemble_tangent_stiffness(Nodal_Field, Nodal_Field, Mask,
                                           Particle, Mesh, double,
                                           Newmark_parameters);
/**************************************************************/

/**
 * @brief
 *
 * @param D_upw
 * @param upw_n
 * @param Params
 */
static void __update_Nodal_Increments(Nodal_Field D_upw, Nodal_Field upw_n,
                                      Newmark_parameters Params);
/**************************************************************/

/**
 * @brief
 *
 * @param D_upw
 * @param MPM_Mesh
 * @param FEM_Mesh
 * @param ActiveNodes
 */
static void __update_Particles(Nodal_Field D_upw, Particle MPM_Mesh,
                               Mesh FEM_Mesh, Mask ActiveNodes);

/**************************************************************/

int Newmark_beta__upw__(Mesh FEM_Mesh, Particle MPM_Mesh,
                        Time_Int_Params Parameters_Solver) {

  int STATUS = EXIT_SUCCESS;

  //  Auxiliar variables for the solver
  unsigned Ndim = NumberDimensions;
  unsigned Ndof = Ndim + 1;
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
  double CFL = Parameters_Solver.CFL;
  double DeltaTimeStep;
  double DeltaX = FEM_Mesh.DeltaX;

  double Error_0;
  double Error_i;
  double Error_relative;

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

  Nodal_Field D_upw;
  Nodal_Field upw_n;

  Mask ActiveNodes;
  Mask ActiveDOFs;

  /*
    Alpha parameters for the Newmark-beta
  */
  Newmark_parameters Params;

  /*
    Time step is defined at the init of the simulation throught the
    CFL condition. Notice that for this kind of solver, CFL confition is
    not required to be satisfied. The only purpose of it is to use the existing
    software interfase.
  */
  DeltaTimeStep =
      0.01; // DeltaT_Coussy__SolversLib__(MPM_Mesh, DeltaX, 1.0, CFL);
  Params = __compute_Newmark_parameters(beta, gamma, DeltaTimeStep);

  while (TimeStep < NumTimeStep) {

    DoProgress("Simulation:", TimeStep, NumTimeStep);

    //! Local search and compute list of active nodes and dofs
    local_search__MeshTools__(MPM_Mesh, FEM_Mesh);
    ActiveNodes = get_active_nodes__MeshTools__(FEM_Mesh);
    Nactivenodes = ActiveNodes.Nactivenodes;
    Ntotaldofs = Ndof * Nactivenodes;
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
    upw_n = __get_nodal_field_tn(MPM_Mesh, FEM_Mesh, ActiveNodes);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __get_nodal_field_tn()" RESET " \n");
      return EXIT_FAILURE;
    }

    //! Compute kinematic nodal values
    D_upw = __initialise_nodal_increments(upw_n, ActiveNodes, FEM_Mesh, Params,
                                          TimeStep, NumTimeStep);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in __initialise_nodal_increments()" RESET " \n");
      return EXIT_FAILURE;
    }

    //! Trial residual
    Residual = __assemble_residual(upw_n, D_upw, ActiveNodes, MPM_Mesh,
                                   FEM_Mesh, Params, TimeStep, NumTimeStep);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in __assemble_residual()" RESET " \n");
      return EXIT_FAILURE;
    }

    //! Check initial error and preallocate tangent matrix if it is necessary
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

      Tangent_Stiffness = __assemble_tangent_stiffness(
          upw_n, D_upw, ActiveNodes, MPM_Mesh, FEM_Mesh, epsilon, Params);
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

      __update_Nodal_Increments(D_upw, upw_n, Params);

      local_compatibility_conditions__upw__(D_upw, ActiveNodes, MPM_Mesh,
                                            FEM_Mesh, &STATUS);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in local_compatibility_conditions__upw__()" RESET
                " \n");
        return EXIT_FAILURE;
      }

      constitutive_update__upw__(MPM_Mesh, FEM_Mesh, &STATUS);

      // Free memory
#ifdef USE_PETSC
      VecDestroy(&Residual);
#else
      free(Residual);
#endif

      // Compute residual (NR-loop)
      Residual = __assemble_residual(upw_n, D_upw, ActiveNodes, MPM_Mesh,
                                     FEM_Mesh, Params, TimeStep, NumTimeStep);
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

    __update_Particles(D_upw, MPM_Mesh, FEM_Mesh, ActiveNodes);

    //!  Outputs
    if (TimeStep % ResultsTimeStep == 0) {
      Residual = __assemble_residual(upw_n, D_upw, ActiveNodes, MPM_Mesh,
                                     FEM_Mesh, Params, TimeStep, NumTimeStep);
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

    //!  Free memory
#ifdef USE_PETSC
    VecDestroy(&Lumped_Mass);
    VecDestroy(&upw_n.value);
    VecDestroy(&upw_n.d_value_dt);
    VecDestroy(&upw_n.d2_value_dt2);
    VecDestroy(&D_upw.value);
    VecDestroy(&D_upw.d_value_dt);
    VecDestroy(&D_upw.d2_value_dt2);
#else
    free(Lumped_Mass);
    free(upw_n.value);
    free(upw_n.d_value_dt);
    free(upw_n.d2_value_dt2);
    free(D_upw.value);
    free(D_upw.d_value_dt);
    free(D_upw.d2_value_dt2);
#endif

    free(ActiveNodes.Nodes2Mask);
    free(ActiveDOFs.Nodes2Mask);
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static Newmark_parameters
__compute_Newmark_parameters(double beta, double gamma, double DeltaTimeStep) {

  Newmark_parameters Params;

  Params.alpha_1 = 1 / (beta * DSQR(DeltaTimeStep));
  Params.alpha_2 = 1 / (beta * DeltaTimeStep);
  Params.alpha_3 = (1 - 2 * beta) / (2 * beta);
  Params.alpha_4 = gamma / (beta * DeltaTimeStep);
  Params.alpha_5 = 1 - gamma / beta;
  Params.alpha_6 = (1 - gamma / (2 * beta)) * DeltaTimeStep;

  return Params;
}

/**************************************************************/

static void __compute_local_mass_matrix(double *Local_Mass_Matrix_p,
                                        double Na_p, double m_p) {

  double M_AB_p = Na_p * m_p;

#if NumberDimensions == 2
  Local_Mass_Matrix_p[0] = M_AB_p; // u.x
  Local_Mass_Matrix_p[1] = M_AB_p; // u.y
  Local_Mass_Matrix_p[2] = M_AB_p; // pw
#else
  Local_Mass_Matrix_p[0] = M_AB_p; // u.x
  Local_Mass_Matrix_p[1] = M_AB_p; // u.y
  Local_Mass_Matrix_p[2] = M_AB_p; // u.z
  Local_Mass_Matrix_p[2] = M_AB_p; // pw
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

static Matrix __compute_nodal_lumped_mass(Particle MPM_Mesh, Mesh FEM_Mesh,
                                          Mask ActiveNodes, double epsilon)
/*
  This function computes the effective mass matrix as a convex combination
  of the lumped mass matrix and the consistent mass matrix.
*/
{

  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Np = MPM_Mesh.NumGP;
  int Order = Ndof * Nnodes_mask;
  int Ap;
  int Bp;
  int A_mask;
  int B_mask;
  int idx_AB_mask_i;
  int idx_A_mask_i;

  /* Value of the shape-function */
  Matrix N_p;

  /* Evaluation of the particle in the node */
  double Na_p;
  double Nb_p;
  /* Mass of the particle */
  double m_p;
  /* Element for each particle */
  Element Nodes_p;

  /* Define and allocate the effective mass matrix */
  Matrix Effective_MassMatrix = allocZ__MatrixLib__(Order, Order);

  /*
    Iterate over the particles to get the nodal values
  */
  for (int p = 0; p < Np; p++) {

    /*
      Define tributary nodes of the particle
    */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],
                                     MPM_Mesh.NumberNodes[p]);

    /*
      Evaluate the shape function in the coordinates of the particle
    */
    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Get the mass of the particle
    */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    for (int A = 0; A < Nodes_p.NumberNodes; A++) {

      /*
        Get the node in the mass matrix with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Get the value of the shape function
      */
      Na_p = N_p.nV[A];

      for (int B = 0; B < Nodes_p.NumberNodes; B++) {
        /*
          Get the node in the mass matrix with the mask
        */
        Bp = Nodes_p.Connectivity[B];
        B_mask = ActiveNodes.Nodes2Mask[Bp];

        /* Get the value of the shape function */
        Nb_p = N_p.nV[B];

        /*
          Fill the effective mass matrix considering the number of dofs
        */
        for (int i = 0; i < Ndof; i++) {
          Effective_MassMatrix.nM[A_mask * Ndof + i][B_mask * Ndof + i] +=
              (1 - epsilon) * m_p * Na_p * Nb_p +
              (A_mask == B_mask) * epsilon * m_p * Na_p;
        }
      }
    }

    /*
      Free the value of the shape functions
    */
    free__MatrixLib__(N_p);
    free(Nodes_p.Connectivity);
  }

  return Effective_MassMatrix;
}

/**************************************************************/

static Nodal_Field compute_Nodal_Field(Particle MPM_Mesh, Mesh FEM_Mesh,
                                       Mask ActiveNodes)
/*
  Call the LAPACK solver to compute simultanesly :

  -> The nodal Acceleration/Second time derivative of the pore water pressure.
  The operation is linearized and all the dof split the solution array in n
  components like : | M 0 0 |   |d2 (u.x) dt2|   |m * d2 (u.x) dt2| | 0 M 0 | *
  |d2 (u.y) dt2| = |m * d2 (u.y) dt2| | 0 0 M |   |d2 (p.w) dt2|   |m * d2 (p.w)
  dt2|

  -> The nodal Velocity/First time derivative of the pore water pressure.
  The operation is linearized and all the dof split the solution array in n
  components like : | M 0 0 |   |d (u.x) dt|   |m * d (u.x) dt| | 0 M 0 | * |d
  (u.y) dt| = |m * d (u.y) dt| | 0 0 M |   |d (p.w) dt|   |m * d (p.w) dt|

  -> The nodal displacement/pore water pressure.
  The operation is linearized and all the dof split the solution array in n
  components like : | M 0 0 |   |u.x|   |m * u.x| | 0 M 0 | * |u.y| = |m * u.y|
  | 0 0 M |   |p.w|   |m * p.w|
*/
{

  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int idx_A_mask_i;

  /* Value of the shape-function */
  Matrix N_p;

  /* Evaluation of the particle in the node */
  double Na_p;
  /* Mass of the particle */
  double m_p;
  /* Element for each particle */
  Element Nodes_p;

  double epsilon = 1;
  Matrix Effective_Mass =
      __compute_nodal_lumped_mass(MPM_Mesh, FEM_Mesh, ActiveNodes, epsilon);

  /* Define and allocate the output vector */
  Nodal_Field upw;
  upw.value = allocZ__MatrixLib__(Nnodes_mask, Ndof);
  upw.d_value_dt = allocZ__MatrixLib__(Nnodes_mask, Ndof);
  upw.d2_value_dt2 = allocZ__MatrixLib__(Nnodes_mask, Ndof);

  /* Iterate over the particles to get the nodal values */
  for (int p = 0; p < Np; p++) {
    /*
      Define element of the particle
    */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],
                                     MPM_Mesh.NumberNodes[p]);

    /*
      Evaluate the shape function in the coordinates of the particle
    */
    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Get the mass of the GP
    */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    for (int A = 0; A < Nodes_p.NumberNodes; A++) {

      /*
        Get the node with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Evaluate the GP function in the node
      */
      Na_p = N_p.nV[A];

      /*
        Nodal displacement and pressure (and its rates)
      */
      for (int i = 0; i < Ndof; i++) {
        if (i < Ndim) {
          upw.value.nM[A_mask][i] += m_p * Na_p * MPM_Mesh.Phi.dis.nM[p][i];
          upw.d_value_dt.nM[A_mask][i] +=
              m_p * Na_p * MPM_Mesh.Phi.vel.nM[p][i];
          upw.d2_value_dt2.nM[A_mask][i] +=
              m_p * Na_p * MPM_Mesh.Phi.acc.nM[p][i];
        } else {
          upw.value.nM[A_mask][i] += m_p * Na_p * MPM_Mesh.Phi.Pw.nV[p];
          upw.d_value_dt.nM[A_mask][i] +=
              m_p * Na_p * MPM_Mesh.Phi.d_Pw_dt_n.nV[p];
          upw.d2_value_dt2.nM[A_mask][i] +=
              m_p * Na_p * MPM_Mesh.Phi.d2_Pw_dt2.nV[p];
        }
      }
    }

    /*
      Free the value of the shape functions
    */
    free__MatrixLib__(N_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Call the LAPACK solver to compute the accelerations and second derivative of
    the pore water pressure
  */
  int Order = Nnodes_mask * Ndof;
  int LDA = Order;
  int LDB = Order;
  char TRANS = 'T'; /* (Transpose) */
  int INFO = 3;
  int *IPIV = (int *)Allocate_Array(Order, sizeof(int));
  int NRHS = 1;

  /*
    Compute the LU factorization for the mass matrix
  */
  dgetrf_(&Order, &Order, Effective_Mass.nV, &LDA, IPIV, &INFO);

  if (INFO != 0) {
    if (INFO < 0) {
      printf("%s : \n", "Error in __get_nodal_field_tn()");
      printf("the %i-th argument had an illegal value", abs(INFO));
    } else if (INFO > 0) {
      printf("%s :\n", "Error in __get_nodal_field_tn()");
      printf(" M(%i,%i) %s \n %s \n %s \n %s \n", INFO, INFO,
             "is exactly zero. The factorization",
             "has been completed, but the factor M is exactly",
             "singular, and division by zero will occur if it is used",
             "to solve a system of equations.");
    }
    exit(EXIT_FAILURE);
  }

  /*
    Solve for the acceleration and second derivative of the pore water pressure
  */
  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass.nV, &LDA, IPIV, upw.value.nV,
          &LDB, &INFO);
  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass.nV, &LDA, IPIV,
          upw.d_value_dt.nV, &LDB, &INFO);
  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass.nV, &LDA, IPIV,
          upw.d2_value_dt2.nV, &LDB, &INFO);

  /*
    Free auxiliar memory
  */
  free(Effective_Mass.nV);
  free(IPIV);
  free__MatrixLib__(Effective_Mass);

  return upw;
}

/**************************************************************/

static Nodal_Field
__initialise_nodal_increments(Nodal_Field upw_n, Mask ActiveNodes,
                              Mesh FEM_Mesh, Newmark_parameters Params,
                              int TimeStep, int NumTimeStep)
/*
  Allocate the increments nodal values for the time integration scheme,
  and apply the boundary conditions over the nodes.
*/
{

  /* 1º Define auxilar variables */
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int NumNodesBound;           /* Number of nodes of the bound */
  int Ndim = NumberDimensions; /* Number of dimensions */
  int Ndof = NumberDOF;        /* Number of degree of freedom */
  int Id_BCC;                  /* Index of the node where we apply the BCC */
  int Id_BCC_mask;

  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

  double D_upw_value_It;
  Nodal_Field D_upw;
  D_upw.value = allocZ__MatrixLib__(Nnodes_mask, Ndof);
  D_upw.d_value_dt = allocZ__MatrixLib__(Nnodes_mask, Ndof);
  D_upw.d2_value_dt2 = allocZ__MatrixLib__(Nnodes_mask, Ndof);

  /*
    Loop over the the boundaries to set boundary conditions
  */
  for (int i = 0; i < FEM_Mesh.Bounds.NumBounds; i++) {

    /*
      Get the number of nodes of this boundarie
    */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

    /*
      Get the number of dimensions where the BCC it is applied
    */
    Ndof = FEM_Mesh.Bounds.BCC_i[i].Dim;

    for (int j = 0; j < NumNodesBound; j++) {
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
      for (int k = 0; k < Ndof; k++) {

        /*
          Apply only if the direction is active (1)
        */
        if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + TimeStep] == 1) {
          /*
            Assign the boundary condition :
            First remove the value obtained during the projection and compute
            the value using the evolution of the boundary condition
          */
          upw_n.value.nM[Id_BCC_mask][k] = 0.0;
          upw_n.d_value_dt.nM[Id_BCC_mask][k] = 0.0;
          upw_n.d2_value_dt2.nM[Id_BCC_mask][k] = 0.0;

          for (int t = 0; t < TimeStep; t++) {
            D_upw_value_It = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[t];
            upw_n.value.nM[Id_BCC_mask][k] += D_upw_value_It;
            upw_n.d2_value_dt2.nM[Id_BCC_mask][k] +=
                alpha_1 * D_upw_value_It -
                alpha_2 * upw_n.d_value_dt.nM[Id_BCC_mask][k] -
                (alpha_3 + 1) * upw_n.d2_value_dt2.nM[Id_BCC_mask][k];
            upw_n.d_value_dt.nM[Id_BCC_mask][k] +=
                alpha_4 * D_upw_value_It +
                (alpha_5 - 1) * upw_n.d_value_dt.nM[Id_BCC_mask][k] +
                alpha_6 * upw_n.d2_value_dt2.nM[Id_BCC_mask][k];
          }

          /*
            Initialise increments using newmark and the value of the boundary
            condition
          */
          D_upw.value.nM[Id_BCC_mask][k] =
              FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep];
          D_upw.d2_value_dt2.nM[Id_BCC_mask][k] =
              alpha_1 * D_upw.value.nM[Id_BCC_mask][k] -
              alpha_2 * upw_n.d_value_dt.nM[Id_BCC_mask][k] -
              (alpha_3 + 1) * upw_n.d2_value_dt2.nM[Id_BCC_mask][k];
          D_upw.d_value_dt.nM[Id_BCC_mask][k] =
              alpha_4 * D_upw.value.nM[Id_BCC_mask][k] +
              (alpha_5 - 1) * upw_n.d_value_dt.nM[Id_BCC_mask][k] +
              alpha_6 * upw_n.d2_value_dt2.nM[Id_BCC_mask][k];
        }
      }
    }
  }

  return D_upw;
}

/**************************************************************/

#ifdef USE_PETSC
static Vec __assemble_residual(Nodal_Field U_n, Nodal_Field D_U,
                               Vec Lumped_Mass, Mask ActiveNodes,
                               Mask ActiveDOFs, Particle MPM_Mesh,
                               Mesh FEM_Mesh, Newmark_parameters Params,
                               bool Is_compute_Residual,
                               bool Is_compute_Reactions, int *STATUS)
#else
static double *__assemble_residual(Nodal_Field U_n, Nodal_Field D_U,
                                   double *Lumped_Mass, Mask ActiveNodes,
                                   Mask ActiveDOFs, Particle MPM_Mesh,
                                   Mesh FEM_Mesh, Newmark_parameters Params,
                                   bool Is_compute_Residual,
                                   bool Is_compute_Reactions, int *STATUS)
#endif
{

  unsigned Size;
  unsigned Ndim = NumberDimensions;
  unsigned Ndof = Ndim + 1;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned p;
  unsigned MixtIdx_p;
  unsigned MatIndx_Water_p;

  if ((Is_compute_Residual == true) && (Is_compute_Reactions == false)) {
    Size = ActiveDOFs.Nactivenodes;
  } else if ((Is_compute_Residual == false) && (Is_compute_Reactions == true)) {
    Size = ActiveNodes.Nactivenodes * Ndof;
  }

#ifdef USE_PETSC
  Vec Residual;
  VecCreate(PETSC_COMM_WORLD, &Residual);
  VecSetSizes(Residual, PETSC_DECIDE, Size);
  VecSetOption(Residual, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  VecSetFromOptions(Residual);
#else
  double *Residual = (double *)calloc(Size, __SIZEOF_DOUBLE__);
  if (Residual == NULL) {
    fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
    *STATUS = EXIT_FAILURE;
    return Residual;
  }
#endif

#ifdef USE_PETSC
  const PetscScalar *dU;
  VecGetArrayRead(D_U.value, &dU);
  const PetscScalar *dU_dt;
  VecGetArrayRead(D_U.d_value_dt, &dU_dt);
  const PetscScalar *da_I;
  VecGetArrayRead(D_U.d2_value_dt2, &da_I);
#else
  const double *dU = D_U.value;
  const double *dU_dt = D_U.d_value_dt;
  const double *da_I = D_U.d2_value_dt2;
#endif

#if NumberDimensions == 2
  double L0[2], L1[2], L2[2];
#else
  double L0[3], L1[3], L2[3];
#endif
  double G0, G1, G2;
  double g = -10.0;

#pragma omp parallel private(NumNodes_p, L0)
  {
#pragma omp for private(p)
    for (p = 0; p < Np; p++) {

      // Get the volume of the particle in the reference configuration
      //     Get the current intrinsic density, volume fraction
      //     and compressibility for each material point (fluid).
      //     Compute relative density.
      //
      double V0_p = MPM_Mesh.Phi.Vol_0.nV[p];
      double m_p = MPM_Mesh.Phi.mass.nV[p];
      double Jacobian_p = MPM_Mesh.Phi.J_n1.nV[p];
      double rate_Jacobian_p = MPM_Mesh.Phi.dJ_dt.nV[p];
      double phi_f_p = MPM_Mesh.Phi.phi_f.nV[p];
      double intrinsic_rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];
      const double *DF_p = MPM_Mesh.Phi.DF.nM[p];
      const double *a_n_p = MPM_Mesh.Phi.acc.nM[p];
#if NumberDimensions == 2
      double b_n1[2];
#else
      double b_n1[3];
#endif

      /*
        Load intrinsic properties for the fluid phase to
        get the fluid compressibility
      */
      MixtIdx_p = MPM_Mesh.MixtIdx[p];
      MatIndx_Water_p = Soil_Water_Mixtures[MixtIdx_p].Water_Idx;
      Material MatProp_Water_p = MPM_Mesh.Mat[MatIndx_Water_p];
      double kappa_f = MatProp_Water_p.Compressibility;
      const double *k_p = Soil_Water_Mixtures[MixtIdx_p].Permeability;

      //  Define nodal connectivity for each particle
      //  and compute gradient of the shape function
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      const ChainPtr ListNodes = MPM_Mesh.ListNodes[p];
      Element Nodes_p = nodal_set__Particles__(p, ListNodes, NumNodes_p);
      Matrix ShapeFunction_p =
          compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
      Matrix d_shapefunction_n_p =
          compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      // Pushforward the shape functions
      const double *d_shapefunction_n1_p = push_forward_dN__MeshTools__(
          d_shapefunction_n_p.nV, DF_p, NumNodes_p, STATUS);
      if (*STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in push_forward_dN__MeshTools__()" RESET " \n");
      }

      // Get the Kirchhoff stress tensor pointer
      const double *kirchhoff_p = MPM_Mesh.Phi.Stress.nM[p];
      double kichhoff_pressure = MPM_Mesh.Phi.Pw_n1.nV[p];
      double rate_kichhoff_pressure = MPM_Mesh.Phi.d_Pw_dt_n1.nV[p];

      //!
      const double *grad_kichhoff_pressure = ;

      //! Compute the acceleration b_n1 = a_n1 - g_n1
      compute_total_acceleration__upw__(b_n1, a_n_p, da_I, ShapeFunction_p.nV,
                                        ListNodes, ActiveNodes);

      for (unsigned A = 0; A < NumNodes_p; A++) {

        //! Get the gradient evaluation in node A \n
        //! and the masked index of the node A
        int Ap = Nodes_p.Connectivity[A];
        int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];
        double ShapeFunction_pA = ShapeFunction_p.nV[A];
        const double *d_shapefunction_n1_pA = &d_shapefunction_n1_p[A * Ndim];

        //! Balance of mometum of the mixture
        compute_L0__upw__(L0, d_shapefunction_n1_pA, kirchhoff_p, V0_p);

        compute_L1__upw__(L1, d_shapefunction_n1_pA, kichhoff_pressure, V0_p);

        compute_L2__upw__(L2, ShapeFunction_pA, b_n1, m_p);

        //! Balance of momentum of the intersticial flow + darcy law
        compute_G0__upw__(&G0, ShapeFunction_pA, intrinsic_rho_f_p,
                          relative_rho_f_p, rate_kichhoff_pressure,
                          rate_Jacobian_p, kappa_f, V0_p);

        compute_G1__upw__(&G1, d_shapefunction_n1_pA, k_p,
                          grad_kichhoff_pressure, g, V0_p);

        compute_G2__upw__(&G2, d_shapefunction_n1_pA, k_p, b_n1, Jacobian_p,
                          intrinsic_rho_f_p, g, V0_p);

#if NumberDimensions == 2
        int Mask_active_dofs_L_A[2];
#else
        int Mask_active_dofs_L_A[3];
#endif
        int Mask_active_dofs_G_A;

        if (Is_compute_Residual == true) {
          __get_assembling_locations_residual(Mask_active_dofs_L_A, Mask_node_A,
                                              ActiveDOFs);
        } else if (Is_compute_Reactions == true) {
          __get_assembling_locations_reactions(Mask_active_dofs_L_A,
                                               Mask_node_A);
        }

#pragma omp critical
        {
#ifdef USE_PETSC
          VecSetValues(Residual, Ndim, Mask_active_dofs_L_A, L0, ADD_VALUES);
          VecSetValues(Residual, Ndim, Mask_active_dofs_L_A, L1, ADD_VALUES);
          VecSetValues(Residual, Ndim, Mask_active_dofs_L_A, L2, ADD_VALUES);

          VecSetValues(Residual, 1, Mask_active_dofs_G_A, &G0, ADD_VALUES);
          VecSetValues(Residual, 1, Mask_active_dofs_G_A, &G1, ADD_VALUES);
          VecSetValues(Residual, 1, Mask_active_dofs_G_A, &G2, ADD_VALUES);
#else
          VecSetValues(Residual, L0, Mask_active_dofs_L_A);
          VecSetValues(Residual, L1 Mask_active_dofs_L_A);
          VecSetValues(Residual, L2, Mask_active_dofs_L_A);
#endif
        } // #pragma omp critical
      }   // for unsigned A

      //   Free memory
      free__MatrixLib__(ShapeFunction_p);
      free__MatrixLib__(d_shapefunction_n_p);
      free(d_shapefunction_n1_p);
      free(Nodes_p.Connectivity);
    } // For unsigned p
  }   // #pragma omp parallel

#ifdef USE_PETSC
  VecRestoreArrayRead(D_U.value, &dU);
  VecRestoreArrayRead(D_U.d_value_dt, &dU_dt);
  VecRestoreArrayRead(D_U.d2_value_dt2, &da_I);
#endif

  return Residual;
}

/**************************************************************/

static void compute_nominal_traction_and_fluid_flux(Matrix Residual,
                                                    Mask ActiveNodes,
                                                    Particle MPM_Mesh,
                                                    Mesh FEM_Mesh, int TimeStep,
                                                    int NumTimeStep) {
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  Load Load_i;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_p;      /* Nodal values of the sahpe function */
  double Na_p;
  Tensor T = alloc__TensorLib__(1); // Nominal traction
  double Q;                         // Nominal flux
  double A0_p; /* Area of the particle in the reference configuration */

  int NumContactForces = MPM_Mesh.Neumann_Contours.NumBounds;
  int NumNodesLoad;
  int p;
  int Ap;
  int A_mask;
  int NumNodes_p; /* Number of nodes of each particle */

  for (int i = 0; i < NumContactForces; i++) {

    /*
      Read load i
    */
    Load_i = MPM_Mesh.Neumann_Contours.BCC_i[i];

    NumNodesLoad = Load_i.NumNodes;

    for (int j = 0; j < NumNodesLoad; j++) {

      /*
        Get the index of the particle
      */
      p = Load_i.Nodes[j];

      /*
        Get the area of each particle
      */
      if (Ndim == 2) {
        A0_p = MPM_Mesh.Phi.Vol_0.nV[p] / Thickness_Plain_Stress;
      } else if (Ndim == 3) {
        A0_p = MPM_Mesh.Phi.Area_0.nV[p];
      }

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
      for (int k = 0; k < Ndof; k++) {
        if (Load_i.Dir[k * NumTimeStep + TimeStep] == 1) {
          if (k < Ndim) {
            T.n[k] = Load_i.Value[k].Fx[TimeStep];
          } else if (k == Ndim) {
            Q = Load_i.Value[k].Fx[TimeStep];
          }

        } else {
          if (k < Ndim) {
            T.n[k] = 0.0;
          } else if (k == Ndim) {
            Q = 0.0;
          }
        }
      }

      /*
        Get the node of the mesh for the contribution
      */
      for (int A = 0; A < NumNodes_p; A++) {

        /*
          Pass the value of the nodal shape function to a scalar
        */
        Na_p = N_p.nV[A];

        /*
          Node for the contribution
        */
        Ap = Nodes_p.Connectivity[A];
        A_mask = ActiveNodes.Nodes2Mask[Ap];

        /*
          Compute Contact forces
        */
        for (int k = 0; k < Ndim; k++) {
          Residual.nM[A_mask][k] -= Na_p * T.n[k] * A0_p;
        }

        Residual.nM[A_mask][Ndim] -= Na_p * Q * A0_p;
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(N_p);
      free(Nodes_p.Connectivity);
    }
  }

  free__TensorLib__(T);
}

/**************************************************************/

static Matrix __assemble_tangent_stiffness(Nodal_Field upw_n, Nodal_Field D_upw,
                                           Mask ActiveNodes, Particle MPM_Mesh,
                                           Mesh FEM_Mesh, double epsilon,
                                           Newmark_parameters Params)
/*
  This function computes the tangent stiffness matrix as a combination
  of the geometrical stiffness matrix and the material stiffness matrix.
*/
{

  *STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Ndof * Nnodes_mask;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int Bp;
  int B_mask;
  int NumNodes_p;
  int Mixture_idx;
  int Material_Soil_idx;
  int Material_Water_idx;

  /*
    Nodal values
  */
  Matrix Nodal_D_Acceleration_p; // Nodal values of the acceleration increments
  Matrix Nodal_Pw_n_p;           // Nodal values of the pore-water pressure
  Matrix Nodal_D_Pw_p; // Nodal values of the pore-water pressure increments

#if NumberDimensions == 2
  double DL0_Du[4];
  double DL1_Du[4], DL1_Dpw[2];
  double DL2_Du[4], DL2_Dpw[2];
  double DG0_Du[2], DG0_Dpw;

  double Fm1_n1_p[4];
  double grad_v_p[4];
  double b_p[2];
  double gradient_theta_n1_p[2];

  int Mask_active_dofs_A[2];
  int Mask_active_dofs_B[2];
#else
  double DL0_Du[9];
  double DL1_Du[9], DL1_Dpw[3];
  double DL2_Du[9], DL2_Dpw[3];
  double DG0_Du[3], DG0_Dpw;

  double Fm1_n1_p[9];
  double grad_v_p[9];
  double b_p[3];
  double gradient_theta_n1_p[3];

  int Mask_active_dofs_A[3];
  int Mask_active_dofs_B[3];
#endif

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

  // Time integartion parameters
  double alpha_1 = Params.alpha_1;
  double alpha_4 = Params.alpha_4;

  /*
    Loop in the particles for the assembling process
  */
  for (unsigned p = 0; p < Np; p++) {

    //! Get mass and the volume of the particle in the reference configuration
    //! and the jacobian of the deformation gradient
    double m_p = MPM_Mesh.Phi.mass.nV[p];
    double V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    //! Material properties of the particle
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
    Material_Water_idx = Soil_Water_Mixtures[Mixture_idx].Water_Idx;
    Material MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];
    Material MatProp_Water_p = MPM_Mesh.Mat[Material_Water_idx];

    // Water compressibility
    double kappa_f_p = MatProp_Water_p.Compressibility;
    // Volume fraction of fluid for particle p
    double phi_f_p = MPM_Mesh.Phi.phi_f.nV[p];
    // Volume fraction of solid for particle p
    double phi_s_p = MPM_Mesh.Phi.phi_s.nV[p];
    // Intrinsic density of fluid for particle p
    double intrinsic_rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];
    // Intrinsic density of solid for particle p
    double intrinsic_rho_s_p = MPM_Mesh.Phi.rho_s.nV[p];
    // Particle permeability
    double *k_p = Soil_Water_Mixtures[Mixture_idx].Permeability;
    //! Pointer to the incremental deformation gradient
    double *DF_p = MPM_Mesh.Phi.DF.nM[p];
    //! Pointer to the deformation gradient at t = n + 1
    double *F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
    //! Pointer to the deformation gradient rate at t = n + 1
    double *dFdt_n1_p = MPM_Mesh.Phi.dt_F_n1.nM[p];

    // Compute the gradient of the velocity field using the definition
    // grad_v = dFdt_n1_p x Fm1_n1_p
    STATUS = compute_inverse__TensorLib__(Fm1_n1_p, F_n1_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in compute_inverse__TensorLib__()" RESET " \n");
      STATUS = EXIT_FAILURE;
    }
    matrix_product__TensorLib__(grad_v_p, dFdt_n1_p, Fm1_n1_p);

    // Get the jacobian of the deformation gradient and its rate
    double Jacobian = MPM_Mesh.Phi.J_n1.nV[p];
    double rate_Jacobian = MPM_Mesh.Phi.dJ_dt.nV[p];

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

    /*
      Get some nodal values
    */
    Nodal_D_Acceleration_p = get_U_set_field_upw__MeshTools__(
        D_upw.d2_value_dt2, Nodes_p, ActiveNodes);

    //  Get some usefull intermediate results

    //  Get Kirchhoff pore water pressure and its rate
    double kichhoff_pressure = MPM_Mesh.Phi.Pw_n1.nV[p];
    double rate_kichhoff_pressure = MPM_Mesh.Phi.d_Pw_dt_n1.nV[p];

    //  Compute particle pore water pressure gradient
    get_set_scalar_field__MeshTools__(dPw_Ap, dPw, Nodes_p, ActiveNodes);
    get_set_scalar_field__MeshTools__(Pw_n_Ap, Pw_n, Nodes_p, ActiveNodes);

    compute_gradient_kirchoff_pw__upw__(gradient_theta_n1_p, dPw_Ap, Pw_n_Ap,
                                        d_shapefunction_n1_p, NumNodes_p);

    //  Compute mixture/fluid dynamics in the next time step
    for (unsigned i = 0; i < Ndim; i++) {
      b_p[i] = MPM_Mesh.Phi.acc.nM[p][i] - gravity_field.Value[i].Fx[TimeStep];
    }
    for (unsigned A = 0; A < NumNodes_p; A++) {

      // Get the shape function evaluation in node A
      // and the masked index of the node A
      double shapefunction_n_pA = shapefunction_n_p.nV[A];
      int Ap = Nodes_p.Connectivity[A];
      int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

      for (unsigned i = 0; i < Ndim; i++) {
        b_p[i] += shapefunction_n_pA * dU_dt2[A_mask * Ndim + i];
      }

    } // for unsigned A

    // Assembly process
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
            p, DL0_Du, &d_shapefunction_n1_p[A * Ndim],
            &d_shapefunction_n1_p[B * Ndim], d_shapefunction_n_pA,
            d_shapefunction_n_pB, alpha_4, MPM_Mesh, MatProp_p);
        if (*STATUS == EXIT_FAILURE) {
          fprintf(stderr,
                  "" RED "Error in stiffness_density__Constitutive__" RESET
                  "\n");
        }

        compute_DL1_Du__upw__(DL1_Du, d_shapefunction_n_pA,
                              d_shapefunction_n_pB, kichhoff_pressure, V0_p);
        compute_DL1_Dpw__upw__(DL1_Dpw, d_shapefunction_n_pA,
                               shapefunction_n_pB, V0_p);

        compute_DL2_Du__upw__(DL2_Du, shapefunction_n_pA, shapefunction_n_pB,
                              d_shapefunction_n_pB, b_p, kichhoff_pressure,
                              phi_f_p, intrinsic_rho_f_p, kappa_f_p, Jacobian,
                              m_p, alpha_1, V0_p);
        compute_DL2_Dpw__upw__(DL2_Dpw, shapefunction_n_pA, shapefunction_n_pB,
                               b_p, phi_f_p, intrinsic_rho_f_p, kappa_f_p,
                               V0_p);

        compute_DG0_Du__upw__(
            DG0_Du, shapefunction_n_pA, d_shapefunction_n_pB, grad_v_p,
            kichhoff_pressure, rate_kichhoff_pressure, Jacobian, rate_Jacobian,
            phi_s_p, phi_f_p, intrinsic_rho_f_p, kappa_f_p, alpha_4, V0_p);
        compute_DG0_Dpw__upw__(&DG0_Dpw, shapefunction_n_pA, shapefunction_n_pB,
                               rate_kichhoff_pressure, Jacobian, rate_Jacobian,
                               phi_f_p, intrinsic_rho_f_p, kappa_f_p, alpha_4,
                               V0_p);

      } // for B (node)
    }   // for A (node)

    // Free memory
    free__MatrixLib__(shapefunction_n_p);
    free__MatrixLib__(d_shapefunction_n_p);
    free(d_shapefunction_n1_p);
    free(Nodes_p.Connectivity);
  } // for p (particle)

#ifdef USE_PETSC
  MatAssemblyBegin(Tangent_Stiffness, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Tangent_Stiffness, MAT_FINAL_ASSEMBLY);
#endif
}

/**************************************************************/

static void system_reduction(Matrix Tangent_Stiffness, Matrix Residual,
                             Mask ActiveNodes, Mesh FEM_Mesh, int TimeStep,
                             int NumTimeStep) {

  /*
  Define auxilar variables
*/
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Nnodes_mask * Ndof;
  int Number_of_BCC = FEM_Mesh.Bounds.NumBounds;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound;   /* Number of dimensions */
  int Id_BCC;        /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

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
          if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + TimeStep] == 1) {

            for (int A_mask = 0; A_mask < Nnodes_mask; A_mask++) {
              for (int l = 0; l < Ndof; l++) {
                Tangent_Stiffness
                    .nM[A_mask * Ndof + l][Id_BCC_mask * Ndof + k] = 0.0;
                Tangent_Stiffness
                    .nM[Id_BCC_mask * Ndof + k][A_mask * Ndof + l] = 0.0;
              }
            }

            Tangent_Stiffness
                .nM[Id_BCC_mask * Ndof + k][Id_BCC_mask * Ndof + k] = 1.0;
          }
        }
      }
    }
  }
}

/**************************************************************/

static void __update_Nodal_Increments(Nodal_Field D_upw, Nodal_Field upw_n,
                                      Newmark_parameters Params) {

  int Nnodes = upw_n.value.N_rows;
  int Ndof = NumberDOF;
  int Total_dof = Nnodes * NumberDOF;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

  /*
    Update nodal variables
  */
  for (int A = 0; A < Nnodes; A++) {
    for (int i = 0; i < Ndof; i++) {
      D_upw.d2_value_dt2.nM[A][i] = alpha_1 * D_upw.value.nM[A][i] -
                                    alpha_2 * upw_n.d_value_dt.nM[A][i] -
                                    (alpha_3 + 1) * upw_n.d2_value_dt2.nM[A][i];
      D_upw.d_value_dt.nM[A][i] = alpha_4 * D_upw.value.nM[A][i] +
                                  (alpha_5 - 1) * upw_n.d_value_dt.nM[A][i] +
                                  alpha_6 * upw_n.d2_value_dt2.nM[A][i];
    }
  }
}

/**************************************************************/

static void __update_Particles(Nodal_Field D_upw, Particle MPM_Mesh,
                               Mesh FEM_Mesh, Mask ActiveNodes) {

  unsigned Ndim = NumberDimensions;
  unsigned Ndof = NumberDOF;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned p;

  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ap;
  int A_mask;
  int idx_A_mask_i;
  int idx_ij;
  int Mixture_idx;
  int Material_Soil_idx;
  Element Nodes_p; /* Element for each particle */
  Material MatProp_Soil_p;
  Matrix N_p; /* Value of the shape-function in the particle */
  Matrix gradient_p;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor dFdt_n_p;
  Tensor dFdt_n1_p;
  double N_pI; /* Nodal value for the particle */
  double D_upw_pI;
  double D_dt_upw_pI;
  double D_dt2_upw_pI;
  double Vol_0_p;

  /* iterate over the particles */
  for (p = 0; p < Np; p++) {

    // Update the determinant of the deformation gradient
    MPM_Mesh.Phi.J_n.nV[p] = MPM_Mesh.Phi.J_n1.nV[p];

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

    /*
      Compute the deformation energy (reference volume + material properties
      (solid phase))
    */
    Vol_0_p = MPM_Mesh.Phi.Vol_0.nV[p];
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
    MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];

    /*
      Iterate over the nodes of the particle
    */
    for (unsigned A = 0; A < Nodes_p.NumberNodes; A++) {

      // Get the shape function evaluation in node A
      // and the masked index of the node A
      double ShapeFunction_pI = ShapeFunction_p.nV[A];
      int Ap = Nodes_p.Connectivity[A];
      int A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Update the particle primitive fields using nodal values of the
        increments
      */
      for (unsigned i = 0; i < Ndof; i++) {
        D_upw_pI = N_pI * D_upw.value.nM[A_mask][i];
        D_dt_upw_pI = N_pI * D_upw.d_value_dt.nM[A_mask][i];
        D_dt2_upw_pI = N_pI * D_upw.d2_value_dt2.nM[A_mask][i];

        if (i < Ndim) {
          MPM_Mesh.Phi.acc.nM[p][i] += D_dt2_upw_pI;
          MPM_Mesh.Phi.vel.nM[p][i] += D_dt_upw_pI;
          MPM_Mesh.Phi.dis.nM[p][i] += D_upw_pI;
          MPM_Mesh.Phi.x_GC.nM[p][i] += D_upw_pI;
        } else {
          MPM_Mesh.Phi.d2_Pw_dt2.nV[p] += D_dt2_upw_pI;
          MPM_Mesh.Phi.d_Pw_dt_n.nV[p] += D_dt_upw_pI;
          MPM_Mesh.Phi.Pw.nV[p] += D_upw_pI;
        }
      }
    }

    /*
      Free memory
    */
    free(Nodes_p.Connectivity);
    free__MatrixLib__(N_p);
    free__MatrixLib__(gradient_p);
  }
}

/**************************************************************/
