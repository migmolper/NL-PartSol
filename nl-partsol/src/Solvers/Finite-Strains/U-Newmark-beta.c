
#include "Solvers/Finite-Strains/U-Newmark-beta.h"

/**************************************************************/

int U_Newmark_beta_Finite_Strains(Mesh FEM_Mesh, Particle MPM_Mesh,
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
  double CFL = Parameters_Solver.CFL;
  double DeltaTimeStep;
  double DeltaX = FEM_Mesh.DeltaX;

  double Error_0;
  double Error_i;
  double Error_relative;


//#ifdef USE_PETSC
//  Mat Effective_Mass;
//  Mat Tangent_Stiffness;
//  Vec Residual;
//  Vec Reactions;
//#else
  double * Effective_Mass;
  double * Tangent_Stiffness;
  double * Residual;
  double * Reactions;
//#endif

  Nodal_Field U_n;
  Nodal_Field D_U;

  Mask ActiveNodes;
  Mask ActiveDOFs;

  Newmark_parameters Params;

  /*
    Time step is defined at the init of the simulation throught the
    CFL condition. Notice that for this kind of solver, CFL confition is
    not required to be satisfied. The only purpose of it is to use the existing
    software interfase.
  */
  DeltaTimeStep = __compute_deltat(MPM_Mesh, DeltaX, Parameters_Solver);

  //  Compute time integrator parameters
  Params = __compute_Newmark_parameters(beta, gamma, DeltaTimeStep,epsilon);

  while (TimeStep < NumTimeStep) {
    print_Status("*************************************************", TimeStep);
    print_step(TimeStep, DeltaTimeStep);

    local_search__MeshTools__(MPM_Mesh, FEM_Mesh);
    ActiveNodes = get_active_nodes__MeshTools__(FEM_Mesh);
    Nactivenodes = ActiveNodes.Nactivenodes;
    Ntotaldofs = Ndim * Nactivenodes;
    ActiveDOFs = get_active_dofs__MeshTools__(ActiveNodes, FEM_Mesh, TimeStep, NumTimeStep);
    Nactivedofs = ActiveDOFs.Nactivenodes;


    // Get the previous converged nodal value
    U_n.value = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    U_n.d_value_dt = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    U_n.d2_value_dt2 = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    if((U_n.value == NULL) 
    || (U_n.d_value_dt == NULL)
    || (U_n.d2_value_dt2 == NULL)){
      fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
      return EXIT_FAILURE;
    } 
    STATUS = __get_nodal_field_tn(U_n, MPM_Mesh, FEM_Mesh, ActiveNodes);
    if(STATUS == EXIT_FAILURE){
        fprintf(stderr, ""RED"Error in __get_nodal_field_tn()"RESET" \n");
        return EXIT_FAILURE;
    } 

    // Compute kinematic nodal values
    D_U.value = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    D_U.d_value_dt = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    D_U.d2_value_dt2 = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    if((D_U.value == NULL) 
    || (D_U.d_value_dt == NULL)
    || (D_U.d2_value_dt2 == NULL)){
      fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
      return EXIT_FAILURE;
    }
    __initialise_nodal_increments(D_U, U_n, FEM_Mesh, ActiveNodes, Params);
    
    
    // Define and allocate the effective mass matrix
    Effective_Mass = (double *)calloc(Ntotaldofs*Ntotaldofs, __SIZEOF_DOUBLE__);
    if(Effective_Mass == NULL){
      fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
      return EXIT_FAILURE;
    }
    STATUS = __compute_nodal_effective_mass(Effective_Mass, MPM_Mesh, FEM_Mesh, ActiveNodes, epsilon);
    if(STATUS == EXIT_FAILURE){
      fprintf(stderr, ""RED"Error in __compute_nodal_effective_mass()"RESET" \n");
      return EXIT_FAILURE;
    }

    // Trial residual
    Residual = (double *)calloc(Nactivedofs, __SIZEOF_DOUBLE__);
    Reactions = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
    if((Residual == NULL) 
    || (Reactions == NULL)){
      fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
      return EXIT_FAILURE;
    } 

    STATUS = __Nodal_Internal_Forces(Residual, Reactions, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
    if(STATUS == EXIT_FAILURE){
      fprintf(stderr, ""RED"Error in __Nodal_Internal_Forces()"RESET" \n");
      return EXIT_FAILURE;
    }

    __Nodal_Traction_Forces(Residual, Reactions, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh);

    __Nodal_Body_Forces(Residual, Reactions, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh);

    __Nodal_Inertial_Forces(Residual, Effective_Mass, U_n, D_U, ActiveNodes, ActiveDOFs, Params);

    // Compute error
    Error_0 = Error_i = __error_residual(Residual,Nactivedofs);  
    Error_relative = Error_i/Error_0;
    Iter = 0;

    while (Error_relative > TOL) {

      if ((Error_i < TOL * 100) 
      || (Error_relative < TOL) 
      || (Iter > MaxIter)) {
        break;
      }

      Tangent_Stiffness = (double *)calloc(Nactivedofs * Nactivedofs, __SIZEOF_DOUBLE__);
      if(Tangent_Stiffness == NULL){
        fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
        return EXIT_FAILURE;
      } 

//#ifdef USE_PETSC

//#else
      STATUS = __assemble_tangent_stiffness(
        Tangent_Stiffness, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh, Params);
      if(STATUS == EXIT_FAILURE){
          fprintf(stderr, ""RED"Error in __assemble_tangent_stiffness()"RESET" \n");
          return EXIT_FAILURE;
      } 
//#endif

//#ifdef USE_PETSC
//      STATUS = krylov_PETSC(Tangent_Stiffness, Residual, Nactivedofs);
//      if(STATUS == EXIT_FAILURE){
//        fprintf(stderr, ""RED"Error in dgetrs_LAPACK()"RESET" \n");
//        return EXIT_FAILURE;
//      }
//#else 
      STATUS = dgetrs_LAPACK(Tangent_Stiffness, Residual, Nactivedofs);
      if(STATUS == EXIT_FAILURE){
        fprintf(stderr, ""RED"Error in dgetrs_LAPACK()"RESET" \n");
        return EXIT_FAILURE;
      }
//#endif

      __update_Nodal_Increments(Residual, D_U, U_n, ActiveDOFs, Params, Ntotaldofs);

      STATUS = __local_deformation(D_U, ActiveNodes, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
      if(STATUS == EXIT_FAILURE){
        fprintf(stderr, ""RED"Error in __local_deformation()"RESET" \n");
        return EXIT_FAILURE;
      } 

      // Free memory
      free(Residual);
      free(Reactions);
      free(Tangent_Stiffness);

      // Compute residual (NR-loop)
      Residual = (double *)calloc(Nactivedofs, __SIZEOF_DOUBLE__);
      Reactions = (double *)calloc(Ntotaldofs, __SIZEOF_DOUBLE__);
      if((Residual == NULL) 
      || (Reactions == NULL)){
        fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
        return EXIT_FAILURE;
      } 

      STATUS = __Nodal_Internal_Forces(Residual, Reactions, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
      if(STATUS == EXIT_FAILURE){
        fprintf(stderr, ""RED"Error in __Nodal_Internal_Forces()"RESET" \n");
        return EXIT_FAILURE;
      }

      __Nodal_Traction_Forces(Residual, Reactions, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh);

      __Nodal_Body_Forces(Residual, Reactions, ActiveNodes, ActiveDOFs, MPM_Mesh, FEM_Mesh);

      __Nodal_Inertial_Forces(Residual, Effective_Mass, U_n, D_U, ActiveNodes, ActiveDOFs, Params);

      // Get stats for the convergence
      Error_i = __error_residual(Residual,Nactivedofs);
      Error_relative = Error_i/Error_0;
      Iter++;
      printf("Iter: [%i/%i]. Total Error: %e, Relative Error: %e \n", 
      Iter, MaxIter, Error_i, Error_relative);
    }

    print_convergence_stats(TimeStep, Iter, Error_0, Error_i, Error_relative);

    if(Iter > MaxIter)
    {
      fprintf(stderr, ""RED"Convergence not reached in the maximum number of iterations"RESET" \n");
    }

    __update_Particles(D_U, MPM_Mesh, FEM_Mesh, ActiveNodes);

    __output(MPM_Mesh, FEM_Mesh, ActiveNodes, D_U, U_n, Reactions, TimeStep, ResultsTimeStep);

    TimeStep++;

    free(Effective_Mass);
    free(U_n.value);
    free(U_n.d_value_dt);
    free(U_n.d2_value_dt2);
    free(D_U.value);
    free(D_U.d_value_dt);
    free(D_U.d2_value_dt2);
    free(Residual);
    free(Reactions);
    free(ActiveNodes.Nodes2Mask);
    free(ActiveDOFs.Nodes2Mask);

  }

  return EXIT_SUCCESS;
}

/*********************************************************************/

static double __compute_deltat(
  Particle MPM_Mesh, 
  double h, 
  Time_Int_Params Parameters_Solver) {

  double DeltaT;
  double CEL_MAX = 0;
  double C[3] = {0, 0, 0};
  int Ndim = NumberDimensions;
  int Nmat = MPM_Mesh.NumberMaterials;
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

static Newmark_parameters __compute_Newmark_parameters(
  double beta, 
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

  return Params;
}

/**************************************************************/

static int __compute_nodal_effective_mass(
  double * Effective_MassMatrix,
  Particle MPM_Mesh, 
  Mesh FEM_Mesh,
  Mask ActiveNodes, 
  double epsilon)
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
  int Ap;
  int Bp;
  int A_mask;
  int B_mask;
  int idx_AB_mask_i;
  int idx_A_mask_i;

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
  double * Lumped_MassMatrix = (double *)calloc(Order, __SIZEOF_DOUBLE__);
  if(Lumped_MassMatrix == NULL){
      fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
      return EXIT_FAILURE;
  }

  /*
    Iterate over the particles to get the nodal values
  */
  for (unsigned p = 0; p < Np; p++) {

    /*
      Define tributary nodes of the particle
    */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],MPM_Mesh.NumberNodes[p]);

    /*
       Evaluate the shape function in the coordinates of the particle
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Get the mass of the particle
    */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    for (unsigned A = 0; A < Nodes_p.NumberNodes; A++) {

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
        Lumped_MassMatrix[A_mask * Ndim + i] += m_A_p;
      }

      for (unsigned B = 0; B < Nodes_p.NumberNodes; B++) {
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
          Effective_MassMatrix[(A_mask * Ndim + i)*Order + (B_mask * Ndim + i)] += m_AB_p;
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
  for (int A = 0; A < Order; A++) {
    for (int B = 0; B < Order; B++) {
      Effective_MassMatrix[A*Order + B] =
          (1 - epsilon) * Effective_MassMatrix[A*Order + B] +
          (A == B) * epsilon * Lumped_MassMatrix[A];
    }
  }


  free(Lumped_MassMatrix);

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __get_nodal_field_tn(
  Nodal_Field U_n,
  Particle MPM_Mesh,
  Mesh FEM_Mesh, 
  Mask ActiveNodes)
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
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    /* Evaluate the shape function in the coordinates of the particle */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Get the mass of the GP */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* Get the nodal mommentum */
    for (unsigned A = 0; A < Nodes_p.NumberNodes; A++) {

      /*
        Get the node in the nodal momentum with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /* Evaluate the GP function in the node */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /* Nodal velocity and acceleration  */
      for (unsigned i = 0; i < Ndim; i++) {

        idx_A_mask_i = A_mask*Ndim + i;

        U_n.value[idx_A_mask_i] += m_p * ShapeFunction_pA * MPM_Mesh.Phi.dis.nM[p][i];
        U_n.d_value_dt[idx_A_mask_i] += m_p * ShapeFunction_pA * MPM_Mesh.Phi.vel.nM[p][i];
        U_n.d2_value_dt2[idx_A_mask_i] += m_p * ShapeFunction_pA * MPM_Mesh.Phi.acc.nM[p][i];
      }
    }

    /* Free the value of the shape functions */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Call the LAPACK solver to compute the accelerations and velocities
  */
  unsigned Order = Ndim*Nnodes_mask;
  unsigned LDA = Order;
  unsigned LDB = Order;
  char TRANS = 'N';
  int INFO = 3;
  int NRHS = 1;

  int *IPIV = (int *)calloc(Order, __SIZEOF_INT__);
  if(IPIV == NULL){
        fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
        return EXIT_FAILURE;
  }   
  
  double * Effective_Mass = (double *)calloc(Order*Order, __SIZEOF_DOUBLE__);
  if(Effective_Mass == NULL){
    fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
    return EXIT_FAILURE;
  }
  STATUS = __compute_nodal_effective_mass(Effective_Mass, MPM_Mesh, FEM_Mesh, ActiveNodes, epsilon);
  if(STATUS == EXIT_FAILURE){
    fprintf(stderr, ""RED"Error in __compute_nodal_effective_mass()"RESET" \n");
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


  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass, &LDA, IPIV, U_n.value, &LDB, &INFO);
  if (INFO) {
    fprintf(stderr, ""RED"Error in dgetrs_() "RESET"\n");
    return EXIT_FAILURE;
  }

  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass, &LDA, IPIV, U_n.d_value_dt, &LDB, &INFO);
  if (INFO) {
    fprintf(stderr, ""RED"Error in dgetrs_() "RESET"\n");
    return EXIT_FAILURE;
  }

  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass, &LDA, IPIV, U_n.d2_value_dt2, &LDB, &INFO);
  if (INFO) {
    fprintf(stderr, ""RED"Error in dgetrs_() "RESET"\n");
    return EXIT_FAILURE;
  }

  //  Free auxiliar memory
  free(Effective_Mass);
  free(IPIV);

  return STATUS;
}

/**************************************************************/

static void __initialise_nodal_increments(
  Nodal_Field D_U,
  Nodal_Field U_n, 
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Newmark_parameters Params)
/*
  Apply the boundary conditions over the nodes
*/
{
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned NumNodesBound;          
  unsigned Ndim = NumberDimensions;
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
  for (unsigned i = 0; i < FEM_Mesh.Bounds.NumBounds; i++) {

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
          U_n.value[Id_BCC_mask*Ndim + k] = 0.0;
          U_n.d_value_dt[Id_BCC_mask*Ndim + k] = 0.0;
          U_n.d2_value_dt2[Id_BCC_mask*Ndim + k] = 0.0;

          for (unsigned t = 0; t < TimeStep; t++) {
            D_U_value_It = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[t];
            U_n.value[Id_BCC_mask*Ndim + k] += D_U_value_It;
            U_n.d2_value_dt2[Id_BCC_mask*Ndim + k] +=
                alpha_1 * D_U_value_It -
                alpha_2 * U_n.d_value_dt[Id_BCC_mask*Ndim + k] -
                (alpha_3 + 1) * U_n.d2_value_dt2[Id_BCC_mask*Ndim + k];
            U_n.d_value_dt[Id_BCC_mask*Ndim + k] +=
                alpha_4 * D_U_value_It +
                (alpha_5 - 1) * U_n.d_value_dt[Id_BCC_mask*Ndim + k] +
                alpha_6 * U_n.d2_value_dt2[Id_BCC_mask*Ndim + k];
          }

          /*
            Initialise increments using newmark and the value of the boundary
            condition
          */
          D_U.value[Id_BCC_mask*Ndim + k] =
              FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep];
          D_U.d2_value_dt2[Id_BCC_mask*Ndim + k] =
              alpha_1 * D_U.value[Id_BCC_mask*Ndim + k] -
              alpha_2 * U_n.d_value_dt[Id_BCC_mask*Ndim + k] -
              (alpha_3 + 1) * U_n.d2_value_dt2[Id_BCC_mask*Ndim + k];
          D_U.d_value_dt[Id_BCC_mask*Ndim + k] =
              alpha_4 * D_U.value[Id_BCC_mask*Ndim + k] +
              (alpha_5 - 1) * U_n.d_value_dt[Id_BCC_mask*Ndim + k] +
              alpha_6 * U_n.d2_value_dt2[Id_BCC_mask*Ndim + k];
        }
      }
    }
  }

}

/**************************************************************/

static int __local_deformation(
  Nodal_Field D_U, 
  Mask ActiveNodes,
  Particle MPM_Mesh, 
  Mesh FEM_Mesh,
  double TimeStep) {

  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned NumberNodes_p;
  unsigned Order_p;
  int Idx_Element_p;
  int Idx_Patch_p;

  Element Nodes_p;
  Matrix gradient_p;
  double * D_Displacement_Ap;
  double * D_Velocity_Ap;
  double * F_n_p;
  double * F_n1_p;
  double * DF_p;
  double * dFdt_n_p;
  double * dFdt_n1_p;
  double * dt_DF_p;

  /*
    Loop in the material point set
  */
  for (unsigned p = 0; p < Np; p++) {
    
    //  Define tributary nodes of the particle
    NumberNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],NumberNodes_p);
    Order_p = NumberNodes_p*Ndim;

    //  Get the nodal increment of displacement using the mask
    D_Displacement_Ap = (double *)calloc(Order_p, __SIZEOF_DOUBLE__);
    D_Velocity_Ap = (double *)calloc(Order_p, __SIZEOF_DOUBLE__);    
      if((D_Displacement_Ap == NULL) 
      || (D_Velocity_Ap == NULL)){
        fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
        return EXIT_FAILURE;
      } 

    get_set_field__MeshTools__(D_Displacement_Ap, D_U.value, Nodes_p, ActiveNodes);   
    get_set_field__MeshTools__(D_Velocity_Ap, D_U.d_value_dt, Nodes_p, ActiveNodes);

    /*
      Evaluate the shape function gradient in the coordinates of the particle
    */
    gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Take the values of the deformation gradient from the previous step
    */
    F_n_p = MPM_Mesh.Phi.F_n.nM[p];
    F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
    DF_p = MPM_Mesh.Phi.DF.nM[p];
    dFdt_n_p = MPM_Mesh.Phi.dt_F_n.nM[p];
    dFdt_n1_p = MPM_Mesh.Phi.dt_F_n1.nM[p];
    dt_DF_p = MPM_Mesh.Phi.dt_DF.nM[p];

    update_increment_Deformation_Gradient__Particles__(DF_p, D_Displacement_Ap, gradient_p.nV,NumberNodes_p);

    update_rate_increment_Deformation_Gradient__Particles__(dt_DF_p, D_Velocity_Ap, gradient_p.nV, NumberNodes_p);

    /*
      Update the deformation gradient in t = n + 1 with the information
      from t = n and the increment of deformation gradient.
    */
    update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, DF_p);
    update_rate_Deformation_Gradient_n1__Particles__(dFdt_n1_p, dt_DF_p, F_n_p, DF_p, dFdt_n_p);

    //  Compute Jacobian of the deformation gradient
    MPM_Mesh.Phi.J_n1.nV[p] = I3__TensorLib__(F_n1_p);
    if (MPM_Mesh.Phi.J_n1.nV[p] <= 0.0) {
      fprintf(stderr, ""RED"Negative jacobian in particle %i"RESET" \n",p);
      return EXIT_FAILURE;
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
  if (FEM_Mesh.Locking_Control_Fbar) {

    double Vn_patch;
    double Vn1_patch;
    double J_patch;

    for (int p = 0; p < Np; p++) {

      Idx_Element_p = MPM_Mesh.Element_p[p];
      Idx_Patch_p = FEM_Mesh.Idx_Patch[Idx_Element_p];

      Vn_patch = FEM_Mesh.Vol_Patch_n[Idx_Patch_p];
      Vn1_patch = FEM_Mesh.Vol_Patch_n1[Idx_Patch_p];
      J_patch = Vn1_patch / Vn_patch;

      STATUS = get_locking_free_Deformation_Gradient_n1__Particles__(p, J_patch, MPM_Mesh);
      if(STATUS == EXIT_FAILURE){
        fprintf(stderr, ""RED"Error in get_locking_free_Deformation_Gradient_n1__Particles__()"RESET" \n");
        return EXIT_FAILURE;
      }      

      MPM_Mesh.Phi.Jbar.nV[p] *= J_patch;
    }
  }


  return EXIT_SUCCESS;
}

/**************************************************************/

static int __Nodal_Internal_Forces(
  double * Residual, 
  double * Reactions,
  Mask ActiveNodes, 
  Mask ActiveDOFs,
  Particle MPM_Mesh, 
  Mesh FEM_Mesh,
  double dt) {

  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned MatIndx_p;

  int Ap, Mask_node_A, Mask_total_dof_Ai, Mask_active_dof_Ai;

  double V0_p; // Volume of the particle in the reference configuration
  double W_n1, W_n; // Internal energy
  double * P_p; // First Piola-Kirchhoff Stress tensor
  double * F_n_p; // Deformation gradient t = n
  double * dFdt_n1_p; // Deformation gradient rate t = n + 1

  #if NumberDimensions == 2
  double InternalForcesDensity_Ap[2];
  #else
  double InternalForcesDensity_Ap[3];
  #endif  

  Element Nodes_p;   /* List of nodes for particle */
  Material MatProp_p;  
  Matrix gradient_p; /* Shape functions gradients */
  double * gradient_pA;

  
  for (unsigned p = 0; p < Np; p++) {

    //  Get the volume of the particle in the reference configuration
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    //  Define nodal connectivity for each particle 
    //  and compute gradient of the shape function
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);
    gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    //  Update the first Piola-Kirchhoff stress tensor with an apropiate
    //  integration rule.
    MatIndx_p = MPM_Mesh.MatIdx[p];
    MatProp_p = MPM_Mesh.Mat[MatIndx_p];   
    STATUS = Stress_integration__Particles__(p, MPM_Mesh, FEM_Mesh, MatProp_p);
    if(STATUS == EXIT_FAILURE){
      fprintf(stderr, ""RED"Error in Stress_integration__Particles__(,)"RESET" \n");
      return EXIT_FAILURE;
    }

    // Get the first Piola-Kirchhoff stress tensor
    P_p = MPM_Mesh.Phi.Stress.nM[p];
    F_n_p = MPM_Mesh.Phi.F_n.nM[p];
    dFdt_n1_p = MPM_Mesh.Phi.dt_F_n1.nM[p];

    // Update energy
//    W_n = MPM_Mesh.Phi.W.nV[p];
//    W_n1 = MPM_Mesh.Phi.W.nV[p];
//    W_n1 = W_n; 
//    finite_strains_internal_energy__Particles__(&W_n1, P_p, dFdt_n1_p, dt);

    for (unsigned A = 0; A < NumNodes_p; A++) {

      //  Compute the gradient in the reference configuration
      gradient_pA = gradient_p.nM[A];

      //  Compute the nodal forces of the particle
      __internal_force_density(InternalForcesDensity_Ap, P_p, F_n_p, gradient_pA);

      //  Get the node of the mesh for the contribution
      Ap = Nodes_p.Connectivity[A];
      Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

      //  Asign the nodal forces contribution to the node
      for (unsigned i = 0; i < Ndim; i++) {

        Mask_total_dof_Ai = Mask_node_A * Ndim + i;
        Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

        if(Mask_active_dof_Ai != -1)
        {
          Residual[Mask_active_dof_Ai] += InternalForcesDensity_Ap[i] * V0_p;
        }
        else
        {
          Reactions[Mask_total_dof_Ai] += InternalForcesDensity_Ap[i] * V0_p;
        }
      }

    }

    //   Free memory
    free__MatrixLib__(gradient_p);
    free(Nodes_p.Connectivity);
  }

  return STATUS;
}

/**************************************************************/

static void __internal_force_density(
  double * InternalForcesDensity_Ap,
  const double * P_p,
  const double * F_n_p,
  const double * gradient_pA)
{
#if NumberDimensions == 2 
  double P__x__FnT[4];

  P__x__FnT[0] = P_p[0]*F_n_p[0] + P_p[1]*F_n_p[1];
  P__x__FnT[1] = P_p[0]*F_n_p[2] + P_p[1]*F_n_p[3];
  P__x__FnT[2] = P_p[2]*F_n_p[0] + P_p[3]*F_n_p[1];
  P__x__FnT[3] = P_p[2]*F_n_p[2] + P_p[3]*F_n_p[3];

  InternalForcesDensity_Ap[0] = P__x__FnT[0]*gradient_pA[0] + P__x__FnT[1]*gradient_pA[1];
  InternalForcesDensity_Ap[1] = P__x__FnT[2]*gradient_pA[0] + P__x__FnT[3]*gradient_pA[1];

#else 
  double P__x__FnT[9];
  P__x__FnT[0] = P_p[0]*F_n_p[0] + P_p[1]*F_n_p[1] + P_p[2]*F_n_p[2];
  P__x__FnT[1] = P_p[0]*F_n_p[3] + P_p[1]*F_n_p[4] + P_p[2]*F_n_p[5];
  P__x__FnT[2] = P_p[0]*F_n_p[6] + P_p[1]*F_n_p[7] + P_p[2]*F_n_p[8];
  P__x__FnT[3] = P_p[3]*F_n_p[0] + P_p[4]*F_n_p[1] + P_p[5]*F_n_p[2];
  P__x__FnT[4] = P_p[3]*F_n_p[3] + P_p[4]*F_n_p[4] + P_p[5]*F_n_p[5];
  P__x__FnT[5] = P_p[3]*F_n_p[6] + P_p[4]*F_n_p[7] + P_p[5]*F_n_p[8];
  P__x__FnT[6] = P_p[6]*F_n_p[0] + P_p[7]*F_n_p[1] + P_p[8]*F_n_p[2];
  P__x__FnT[7] = P_p[6]*F_n_p[3] + P_p[7]*F_n_p[4] + P_p[8]*F_n_p[5];
  P__x__FnT[8] = P_p[6]*F_n_p[6] + P_p[7]*F_n_p[7] + P_p[8]*F_n_p[8];

  InternalForcesDensity_Ap[0] = P__x__FnT[0]*gradient_pA[0] + P__x__FnT[1]*gradient_pA[1] + P__x__FnT[2]*gradient_pA[2];
  InternalForcesDensity_Ap[1] = P__x__FnT[3]*gradient_pA[0] + P__x__FnT[4]*gradient_pA[1] + P__x__FnT[5]*gradient_pA[2];
  InternalForcesDensity_Ap[2] = P__x__FnT[6]*gradient_pA[0] + P__x__FnT[7]*gradient_pA[1] + P__x__FnT[8]*gradient_pA[2];

#endif
}


/**************************************************************/

static void __Nodal_Traction_Forces(
  double * Residual, 
  double * Reactions,
  Mask ActiveNodes, 
  Mask ActiveDOFs,
  Particle MPM_Mesh,
  Mesh FEM_Mesh) {

  unsigned Ndim = NumberDimensions;
  unsigned NumContactForces = MPM_Mesh.Neumann_Contours.NumBounds;
  unsigned NumNodesLoad;

  Load Load_i;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_p;      /* Nodal values of the sahpe function */
  double N_pa;
  double A0_p; /* Area of the particle in the reference configuration */

#if NumberDimensions == 2
  double T[2] = {0.0, 0.0};
#else
  double T[3] = {0.0, 0.0};
#endif

  int p;
  int Ap, Mask_node_A, Mask_total_dof_Ai, Mask_active_dof_Ai;
  int NumNodes_p; /* Number of nodes of each particle */

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

          if(Mask_active_dof_Ai != -1)
          {
            Residual[Mask_active_dof_Ai] -= N_pa * T[i] * A0_p;
          }
          else
          {
            Reactions[Mask_total_dof_Ai] -= N_pa * T[i] * A0_p;
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

static void __Nodal_Body_Forces(
  double * Residual, 
  double * Reactions,
  Mask ActiveNodes, 
  Mask ActiveDOFs,
  Particle MPM_Mesh, 
  Mesh FEM_Mesh) {
  
#if NumberDimensions == 2
  double b[2] = {0.0, 0.0};
#else
  double b[3] = {0.0, 0.0, 0.0};
#endif

  /* Define auxilar variables */
  unsigned Ndim = NumberDimensions;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned NumBodyForces = MPM_Mesh.NumberBodyForces;
  unsigned NumNodes_p;

  int Ap, Mask_node_A, Mask_total_dof_Ai, Mask_active_dof_Ai;

  Load *B = MPM_Mesh.B;    /* List with the load cases */
  Element Nodes_p;         /* Element for each particle */
  Matrix ShapeFunction_p;  /* Nodal values of the sahpe function */
  double ShapeFunction_pA; /* Evaluation in the node I for the particle p */
  double m_p;              /* Mass of the particle */

//  b[1] = -9.81;

  for (unsigned p = 0; p < MPM_Mesh.NumGP; p++) {

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

        if(Mask_active_dof_Ai != -1)
        {
          Residual[Mask_active_dof_Ai] -= ShapeFunction_pA * b[i] * m_p;
        }
        else
        {
          Reactions[Mask_total_dof_Ai] -= ShapeFunction_pA * b[i] * m_p;
        }
      }
    }

    /* Free the matrix with the nodal gradient of the element */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

}

/**************************************************************/

static void __Nodal_Inertial_Forces(
  double * Residual,
  double * Mass,  
  Nodal_Field U_n, 
  Nodal_Field D_U,
  Mask ActiveNodes,  
  Mask ActiveDOFs,  
  Newmark_parameters Params) {
  
  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim * Nactivenodes;
  int Mask_active_dof_Ai,Mask_active_dof_Bj;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;

  for (unsigned Mask_total_dof_Ai = 0; Mask_total_dof_Ai < Ntotaldofs; Mask_total_dof_Ai++) {

    Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

    if(Mask_active_dof_Ai != -1) {

      for (unsigned Mask_total_dof_Bj = 0; Mask_total_dof_Bj < Ntotaldofs; Mask_total_dof_Bj++) {

        Mask_active_dof_Bj = ActiveDOFs.Nodes2Mask[Mask_total_dof_Bj];

          Residual[Mask_active_dof_Ai] += Mass[Mask_total_dof_Ai*Ntotaldofs + Mask_total_dof_Bj] * 
          (alpha_1 * D_U.value[Mask_total_dof_Bj] -
          alpha_2 * U_n.d_value_dt[Mask_total_dof_Bj] -
          alpha_3 * U_n.d2_value_dt2[Mask_total_dof_Bj]);
      }
    }
  }
}

/**************************************************************/
static double __error_residual(const double * Residual, int Total_dof) {
  
  double Error = 0;
  
  for (unsigned A = 0; A < Total_dof; A++) {
      Error += DSQR(Residual[A]);
  }
  Error = pow(Error, 0.5);

  return Error;
}

/**************************************************************/

static void compute_local_intertia(
  double * Inertia_density_p, 
  double Na_p,
  double Nb_p,
  double m_p, 
  double alpha_1, 
  double epsilon,
  unsigned A, 
  unsigned B)
{

double ID_AB_p = alpha_1 * ((1 - epsilon) * Na_p * Nb_p + (A == B) * epsilon * Na_p) * m_p;

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

static int __assemble_tangent_stiffness(
  double * Tangent_Stiffness,
  Mask ActiveNodes,
  Mask ActiveDOFs,
  Particle MPM_Mesh, 
  Mesh FEM_Mesh,
  Newmark_parameters Params)
{
  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  unsigned Ntotaldofs = Ndim*Nactivenodes;
  unsigned Nactivedofs = ActiveDOFs.Nactivenodes;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned NumNodes_p;
  unsigned MatIndx_p;

  int Ap, Mask_node_A, Mask_total_dof_Ai, Mask_active_dof_Ai;
  int Bp, Mask_node_B, Mask_total_dof_Bj, Mask_active_dof_Bj;

  // Spatial discretization variables
  Element Nodes_p;
  Matrix shapefunction_p;
  Matrix d_shapefunction_p;
  double * d_shapefunction_pA;
  double * d_shapefunction_pB;    
  double shapefunction_pA;
  double shapefunction_pB;

#if NumberDimensions == 2
  double Stiffness_density_p[4];
  double Inertia_density_p[4] = {
    0.0,0.0,
    0.0,0.0};
#else
  double Stiffness_density_p[9];
  double Inertia_density_p[9] = {
    0.0,0.0,0.0,
    0.0,0.0,0.0,
    0.0,0.0,0.0};
#endif

  Material MatProp_p;
  State_Parameters IO_State_p;
  double m_p; /* Mass of the particle */
  double V0_p; /* Volume of the particle in the reference configuration */
  double J_p;  /* Jacobian of the deformation gradient */
  double alpha_1 = Params.alpha_1;
  double alpha_4 = Params.alpha_4;
  double epsilon = Params.epsilon;

  for (unsigned p = 0; p < Np; p++) {

    // Get mass and the volume of the particle in the reference configuration and
    // the jacobian of the deformation gradient
    m_p = MPM_Mesh.Phi.mass.nV[p];
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];
    J_p = MPM_Mesh.Phi.J_n1.nV[p];

    // Material properties of the particle
    MatIndx_p = MPM_Mesh.MatIdx[p];
    MatProp_p = MPM_Mesh.Mat[MatIndx_p];

    //  Define nodal connectivity for each particle 
    //  and compute gradient of the shape function
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);
    shapefunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    d_shapefunction_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);


    for (unsigned A = 0; A < NumNodes_p; A++) {
      
      // Get the gradient evaluation in node A
      // and the masked index of the node A
      shapefunction_pA = shapefunction_p.nV[A];
      d_shapefunction_pA = d_shapefunction_p.nM[A];
      Ap = Nodes_p.Connectivity[A];
      Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

      for (unsigned B = 0; B < NumNodes_p; B++) {

        // Get the gradient evaluation in node B
        // and the masked index of the node B
        shapefunction_pB = shapefunction_p.nV[B];
        d_shapefunction_pB = d_shapefunction_p.nM[B];
        Bp = Nodes_p.Connectivity[B];
        Mask_node_B = ActiveNodes.Nodes2Mask[Bp];

        // Compute local stiffness density
        if (strcmp(MatProp_p.Type, "Neo-Hookean-Wriggers") == 0) {
          IO_State_p.D_phi_n = MPM_Mesh.Phi.F_n.nM[p]; 
          IO_State_p.d_phi = MPM_Mesh.Phi.DF.nM[p];          
          IO_State_p.J = J_p;          
          STATUS = compute_stiffness_density_Neo_Hookean_Wriggers(Stiffness_density_p, d_shapefunction_pA, d_shapefunction_pB, IO_State_p, MatProp_p);
          if (STATUS == EXIT_FAILURE) {
            fprintf(stderr, "" RED "Error in compute_stiffness_density_Neo_Hookean_Wriggers" RESET "\n");
            return EXIT_FAILURE;
          }
        } 
        else if (strcmp(MatProp_p.Type, "Newtonian-Fluid-Compressible") == 0) {
          IO_State_p.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p]; 
          IO_State_p.D_phi_n = MPM_Mesh.Phi.F_n.nM[p]; 
          IO_State_p.d_phi = MPM_Mesh.Phi.DF.nM[p];          
          IO_State_p.dFdt = MPM_Mesh.Phi.dt_F_n1.nM[p];
          IO_State_p.J = J_p;   
          IO_State_p.alpha_4 = alpha_4;
          STATUS = compute_stiffness_density_Newtonian_Fluid(Stiffness_density_p, d_shapefunction_pA, d_shapefunction_pB,IO_State_p,MatProp_p);
          if (STATUS == EXIT_FAILURE) {
            fprintf(stderr, "" RED "Error in compute_stiffness_density_Newtonian_Fluid" RESET "\n");
            return EXIT_FAILURE;
          }          
        }
        else if (strcmp(MatProp_p.Type, "Drucker-Prager") == 0) {
          IO_State_p.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p]; 
          IO_State_p.D_phi_n = MPM_Mesh.Phi.F_n.nM[p];
          IO_State_p.b_e = MPM_Mesh.Phi.b_e_n1.nM[p];
          IO_State_p.Stress = MPM_Mesh.Phi.Stress.nM[p];
          IO_State_p.C_ep = MPM_Mesh.Phi.C_ep.nM[p];
          STATUS = compute_1PK_elastoplastic_tangent_matrix(Stiffness_density_p, d_shapefunction_pA, d_shapefunction_pB,IO_State_p);
          if (STATUS == EXIT_FAILURE) {
            fprintf(stderr, "" RED "Error in compute_1PK_elastoplastic_tangent_matrix" RESET "\n");
            return EXIT_FAILURE;
          }
        }  
        else if (strcmp(MatProp_p.Type, "Matsuoka-Nakai") == 0) {
          IO_State_p.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p]; 
          IO_State_p.D_phi_n = MPM_Mesh.Phi.F_n.nM[p];
          IO_State_p.b_e = MPM_Mesh.Phi.b_e_n1.nM[p];
          IO_State_p.Stress = MPM_Mesh.Phi.Stress.nM[p];
          IO_State_p.C_ep = MPM_Mesh.Phi.C_ep.nM[p];
          STATUS = compute_1PK_elastoplastic_tangent_matrix(Stiffness_density_p, d_shapefunction_pA, d_shapefunction_pB,IO_State_p);
          if (STATUS == EXIT_FAILURE) {
            fprintf(stderr, "" RED "Error in compute_1PK_elastoplastic_tangent_matrix" RESET "\n");
            return EXIT_FAILURE;
          }          
        }      
        else {
          fprintf(stderr, "%s : %s %s %s \n",
          "Error in __assemble_tangent_stiffness()", "The material",
          MatProp_p.Type, "has not been yet implemnented");
          return EXIT_FAILURE;
        }

        // Local mass matrix
        compute_local_intertia(
          Inertia_density_p, shapefunction_pA, shapefunction_pB,
          m_p, alpha_1, epsilon, Mask_node_A, Mask_node_B);

        //  Assembling process
        for (unsigned i = 0; i < Ndim; i++) {

          Mask_total_dof_Ai = Mask_node_A * Ndim + i;
          Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

          for (unsigned j = 0; j < Ndim; j++) {

            Mask_total_dof_Bj = Mask_node_B * Ndim + j;
            Mask_active_dof_Bj = ActiveDOFs.Nodes2Mask[Mask_total_dof_Bj];

            if((Mask_active_dof_Ai != -1) 
            && (Mask_active_dof_Bj != -1))
            {
              Tangent_Stiffness[Mask_active_dof_Ai*Nactivedofs + Mask_active_dof_Bj] +=  
              Stiffness_density_p[i*Ndim + j] * V0_p + Inertia_density_p[i*Ndim + j];
            }
          }
        }
      }
    }

    free__MatrixLib__(shapefunction_p);
    free__MatrixLib__(d_shapefunction_p);
    free(Nodes_p.Connectivity);
  }

  return EXIT_SUCCESS;
}



/**************************************************************/

static void __update_Nodal_Increments(
  const double * Residual,
  Nodal_Field D_U,
  Nodal_Field U_n,
  Mask ActiveDOFs,
  Newmark_parameters Params,
  unsigned Ntotaldofs) {

  unsigned Ndim = NumberDimensions;
  int Mask_active_dof_Ai;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

  /*
    Update nodal variables
  */
  for (unsigned Mask_total_dof_Ai = 0; Mask_total_dof_Ai < Ntotaldofs; Mask_total_dof_Ai++) {

      Mask_active_dof_Ai = ActiveDOFs.Nodes2Mask[Mask_total_dof_Ai];

      if (Mask_active_dof_Ai == -1) {

        D_U.d2_value_dt2[Mask_total_dof_Ai] = 0.0;

        D_U.d_value_dt[Mask_total_dof_Ai] = 
        alpha_4 * D_U.value[Mask_total_dof_Ai] 
        + (alpha_5 - 1) * U_n.d_value_dt[Mask_total_dof_Ai] 
        + alpha_6 * U_n.d2_value_dt2[Mask_total_dof_Ai];

      } else {

        D_U.value[Mask_total_dof_Ai] -= Residual[Mask_active_dof_Ai];

        D_U.d2_value_dt2[Mask_total_dof_Ai] = 
        alpha_1 * D_U.value[Mask_total_dof_Ai] 
        - alpha_2 * U_n.d_value_dt[Mask_total_dof_Ai] 
        - (alpha_3 + 1) * U_n.d2_value_dt2[Mask_total_dof_Ai];
        
        D_U.d_value_dt[Mask_total_dof_Ai] = 
        alpha_4 * D_U.value[Mask_total_dof_Ai] 
        + (alpha_5 - 1) * U_n.d_value_dt[Mask_total_dof_Ai] 
        + alpha_6 * U_n.d2_value_dt2[Mask_total_dof_Ai];
      }
    }
}

/**************************************************************/

static void __update_Particles(
  Nodal_Field D_U,
  Particle MPM_Mesh, 
  Mesh FEM_Mesh,
  Mask ActiveNodes) {
    
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned NumNodes_p;
  
  int Ap;
  int A_mask;
  int idx_A_mask_i;
  int idx_ij;

  Matrix D_Displacement_Ap;
  Matrix ShapeFunction_p; /* Value of the shape-function in the particle */
  double ShapeFunction_pI; /* Nodal value for the particle */
  Element Nodes_p; /* Element for each particle */
  double D_U_pI;
  double D_V_pI;
  double D_A_pI;

  /* iterate over the particles */
  for (unsigned p = 0; p < Np; p++) {

    // Update the determinant of the deformation gradient
    MPM_Mesh.Phi.J_n.nV[p] = MPM_Mesh.Phi.J_n1.nV[p];

    // Update density
    MPM_Mesh.Phi.rho.nV[p] = MPM_Mesh.Phi.mass.nV[p]/(MPM_Mesh.Phi.Vol_0.nV[p]*MPM_Mesh.Phi.J_n.nV[p]);

    // Update hardening
    MPM_Mesh.Phi.Kappa_n[p] = MPM_Mesh.Phi.Kappa_n1[p];
    
    // Update equivalent plastic strains
    MPM_Mesh.Phi.EPS_n[p] = MPM_Mesh.Phi.EPS_n1[p];

    // Update elastic left Cauchy-Green tensor
#if NumberDimensions == 2
    for (unsigned i = 0 ; i<5 ; i++) MPM_Mesh.Phi.b_e_n.nM[p][i] = MPM_Mesh.Phi.b_e_n1.nM[p][i]; 
#else
    for (unsigned i = 0 ; i<9 ; i++) MPM_Mesh.Phi.b_e_n.nM[p][i] = MPM_Mesh.Phi.b_e_n1.nM[p][i];
#endif

    // Update deformation gradient
#if NumberDimensions == 2
    for (unsigned i = 0 ; i<5 ; i++) MPM_Mesh.Phi.F_n.nM[p][i] = MPM_Mesh.Phi.F_n1.nM[p][i]; 
#else
    for (unsigned i = 0 ; i<9 ; i++) MPM_Mesh.Phi.F_n.nM[p][i] = MPM_Mesh.Phi.F_n1.nM[p][i];
#endif 

    // Update rate of deformation gradient
#if NumberDimensions == 2
    for (unsigned i = 0 ; i<5 ; i++) MPM_Mesh.Phi.dt_F_n.nM[p][i] = MPM_Mesh.Phi.dt_F_n1.nM[p][i]; 
#else
    for (unsigned i = 0 ; i<9 ; i++) MPM_Mesh.Phi.dt_F_n.nM[p][i] = MPM_Mesh.Phi.dt_F_n1.nM[p][i];
#endif 

    //  Define nodal connectivity for each particle 
    //  and compute the shape function
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],NumNodes_p);
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    for (unsigned A = 0; A < NumNodes_p; A++) {

      // Get the shape function evaluation in node A
      // and the masked index of the node A
      ShapeFunction_pI = ShapeFunction_p.nV[A];
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      //  Update acceleration, velocity and position of the particles
      for (unsigned i = 0; i < Ndim; i++) {
        D_U_pI = ShapeFunction_pI * D_U.value[A_mask*Ndim + i];
        D_V_pI = ShapeFunction_pI * D_U.d_value_dt[A_mask*Ndim + i];
        D_A_pI = ShapeFunction_pI * D_U.d2_value_dt2[A_mask*Ndim + i];

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

static void __output(
  Particle MPM_Mesh, 
  Mesh FEM_Mesh, 
  Mask ActiveNodes,
  Nodal_Field D_U, 
  Nodal_Field U_n, 
  double * Reactions, 
  int TimeStep,
  int ResultsTimeStep) {

  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes = ActiveNodes.Nactivenodes;
  Matrix ShapeFunction;

  /*
    vtk results
  */
  if (TimeStep % ResultsTimeStep == 0) {

    int Ngp = MPM_Mesh.NumGP;
    int NumNodes_p;
    Element Nodes_p;
    Matrix ShapeFunction_p;

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

      if (Out_Partition_Unity) {
      MPM_Mesh.Phi.PU = allocZ__MatrixLib__(Ngp, 1);
      for (int p_idx = 0; p_idx < Ngp; p_idx++) {
        NumNodes_p = MPM_Mesh.NumberNodes[p_idx];
        Nodes_p = nodal_set__Particles__(p_idx, MPM_Mesh.ListNodes[p_idx],
                                         NumNodes_p);
        ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

        for (int A = 0; A < NumNodes_p; A++) {
          MPM_Mesh.Phi.PU.nV[p_idx] += ShapeFunction_p.nV[A];
        }

        free(Nodes_p.Connectivity);
        free__MatrixLib__(ShapeFunction_p);
      }
    }

  #endif
#endif

    Matrix Reactions_aux = memory_to_matrix__MatrixLib__(Nactivenodes,Ndim,Reactions);

    particle_results_vtk__InOutFun__(MPM_Mesh, TimeStep, ResultsTimeStep);

    nodal_results_vtk__InOutFun__(FEM_Mesh, ActiveNodes, Reactions_aux, TimeStep,
                                  ResultsTimeStep);

    free(Reactions_aux.nM);

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  if (Out_Partition_Unity) {
      free__MatrixLib__(MPM_Mesh.Phi.PU);
    }

  #endif
#endif

  }

  int Backup_TimeStep = 100;

  if (TimeStep % Backup_TimeStep == 0) {
    particle_backup_vtk__InOutFun__(MPM_Mesh, TimeStep, Backup_TimeStep);
  }

  // /*
  //   csv results
  // */
  // for(int i = 0 ; i<Number_Out_nodal_path_csv ; i++)
  // {

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_Velocity)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(Velocity,
  //     FEM_Mesh.Coordinates,"Nodal_path_velocity_csv", ActiveNodes,
  //     Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_Acceleration)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(Acceleration,
  //     FEM_Mesh.Coordinates,"Nodal_path_acceleration_csv", ActiveNodes,
  //     Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_D_Displacement)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(D_Displacement,
  //     FEM_Mesh.Coordinates,"Nodal_path_displacement_csv", ActiveNodes,
  //     Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_Forces)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(Forces,
  //     FEM_Mesh.Coordinates,"Nodal_path_forces_csv", ActiveNodes,
  //     Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_Reactions)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(Reactions,
  //     FEM_Mesh.Coordinates,"Nodal_path_reactions_csv", ActiveNodes,
  //     Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_Residual)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(Residual,
  //     FEM_Mesh.Coordinates,"Nodal_path_residual_csv", ActiveNodes,
  //     Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }
  // }

  // for(int i = 0 ; i<Number_Out_particles_path_csv ; i++)
  // {
  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Damage)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.chi,
  //     MPM_Mesh.Phi.x_GC, "Particles_path_damage_csv",
  //     Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Velocity)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.vel,
  //     MPM_Mesh.Phi.x_GC, "Particles_path_velocity_csv",
  //     Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Acceleration)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.acc,
  //     MPM_Mesh.Phi.x_GC, "Particles_path_acceleration_csv",
  //     Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Displacement)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.dis,
  //     MPM_Mesh.Phi.x_GC, "Particles_path_displacement_csv",
  //     Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Stress)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.Stress,
  //     MPM_Mesh.Phi.x_GC, "Particles_path_stress_csv",
  //     Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Strain)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.Strain,
  //     MPM_Mesh.Phi.x_GC, "Particles_path_strain_csv",
  //     Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Deformation_gradient)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.F_n,
  //     MPM_Mesh.Phi.x_GC, "Particles_path_deformation_gradient_csv",
  //     Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  // }
}

/**************************************************************/
