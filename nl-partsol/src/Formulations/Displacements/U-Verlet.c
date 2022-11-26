#include "Formulations/Displacements/U-Verlet.h"

double DeltaTimeStep;

/*
  Auxiliar functions
*/
static int __mass_NODES(
    double * Lumped_MassMatrix /**< [out] Lumped mass matrix */,
    Particle MPM_Mesh /**< [in] Information of the particles */,
    Mesh FEM_Mesh /**< [in] Information of the nodes */,   
    Mask ActiveNodes /**< [in] Information of the active nodes */); 

static int __predictor_PARTICLES(
  double * D_dis /**< [in/out] Vector with particle increment of displacement */,
  double * vel /**< [in/out] Vector with particle velocity */,
  const double * acc /**< [in] Vector with particle acceleration */, 
  double gamma /**< [in] Parameter for the time integration scheme */,
  unsigned Np /**< [in] Number of particles */);

static int __gravity_NODES(
  double * Gravity_field /**< [in/out] Vector with the gravity field */,
  const Mask ActiveNodes /**< [in] Information of the active nodes */, 
  const Load * B /**< [in] List of distance forces  */,
  unsigned NumBodyForces /**< [in] Number of distance forces */,
  unsigned TimeStep /**< [in] Current time step */, 
  unsigned NumTimeStep /**< [in] Number of time steps of the simulation */);

static int __d_displacement_NODES(
  double * D_dis_NODES /**< [out] Vector with nodal displacement */,
  const double * Mass_NODES /**< [in] Vector with nodal mass */,
  const double * D_dis_PARTICLES /**< [in] Vector with particle displacement */,
  const double * Mass_PARTICLES /**< [in] Vector with particle mass */,  
  Particle MPM_Mesh /**< [in] Information of the particles */, 
  Mesh FEM_Mesh /**< [in] Information of the nodes */,
  const Mask ActiveNodes /**< [in] Information of the active nodes */);

static Matrix compute_Nodal_Velocity(Particle, Mesh, Mask, Matrix);

static void impose_Dirichlet_Boundary_Conditions(
  double * D_Displacement,
  Mesh FEM_Mesh, 
  Mask ActiveNodes,
  int TimeStep, 
  int NumTimeStep);

/* Step 4 */    
static int __update_Local_State(const double *, Mask, Particle, Mesh);
/* Step 5 */
static Matrix compute_Nodal_Forces(Mask, Particle, Mesh, int, int);
static void compute_Nodal_Internal_Forces(Matrix, Mask, Particle, Mesh);
static void compute_Nodal_Nominal_traction_Forces(Matrix, Mask, Particle, Mesh,
                                                  int, int);
static Matrix solve_Nodal_Equilibrium(Matrix, Matrix, Matrix, Matrix, Particle,
                                      Mesh, Mask, Mask);
/* Step 6 */
static void compute_Explicit_Newmark_Corrector(Particle, double);
/* Step 7 */
static void output_selector(Particle, Mesh, Mask, Matrix, Matrix, Matrix,
                            Matrix, double, int, int);

/**************************************************************/

int U_Verlet(
    Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver) {

  /*
    Auxiliar variables for the solver
  */
  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;
  int Nactivenodes;
  int InitialStep = Parameters_Solver.InitialTimeStep;
  int NumTimeStep = Parameters_Solver.NumTimeStep;

  double gamma = 0.5;
  double CFL = Parameters_Solver.CFL;
  double DeltaX = FEM_Mesh.DeltaX;

  double * Lumped_Mass;
  double * Gravity_field;
  double * D_Displacement;
  double * Forces;
  double * Reactions;

  Mask ActiveNodes;
  Mask Free_and_Restricted_Dofs;

  for (int TimeStep = InitialStep; TimeStep < NumTimeStep; TimeStep++) {

    DeltaTimeStep = U_DeltaT__SolversLib__(MPM_Mesh, DeltaX, Parameters_Solver);
    print_step(TimeStep, NumTimeStep, DeltaTimeStep);
    local_search__MeshTools__(MPM_Mesh, FEM_Mesh);
    ActiveNodes = get_active_nodes__MeshTools__(FEM_Mesh);
    Nactivenodes = ActiveNodes.Nactivenodes;
    Free_and_Restricted_Dofs =
        get_active_dofs__MeshTools__(
            ActiveNodes, FEM_Mesh, TimeStep, NumTimeStep);

    STATUS = __mass_NODES(Lumped_Mass, MPM_Mesh, FEM_Mesh, ActiveNodes);
    if(STATUS == EXIT_FAILURE)
    {
      fprintf(stderr,""RED"Error in __mass_NODES"RESET" \n");
      return EXIT_FAILURE;
    }

    STATUS = __predictor_PARTICLES(MPM_Mesh.Phi.D_dis.nV,
                         MPM_Mesh.Phi.vel.nV,
                         MPM_Mesh.Phi.acc.nV,
                         gamma, MPM_Mesh.NumGP);
    if(STATUS == EXIT_FAILURE)
    {
      fprintf(stderr,""RED"Error in __predictor_PARTICLES"RESET" \n");
      return EXIT_FAILURE;
    }

    STATUS = __d_displacement_NODES(
      D_Displacement, Lumped_Mass, 
      MPM_Mesh.Phi.D_dis.nV, MPM_Mesh.Phi.mass.nV,
      MPM_Mesh, FEM_Mesh, ActiveNodes);
    if(STATUS == EXIT_FAILURE)
    {
      fprintf(stderr,""RED"Error in __d_displacement_NODES"RESET" \n");
      return EXIT_FAILURE;
    }

    impose_Dirichlet_Boundary_Conditions(D_Displacement, FEM_Mesh,
    ActiveNodes, TimeStep, NumTimeStep);

    STATUS = __update_Local_State(D_Displacement, ActiveNodes, MPM_Mesh, FEM_Mesh);
    if(STATUS == EXIT_FAILURE)
    {
      fprintf(stderr,""RED"Error in __update_Local_State"RESET" \n");
      return EXIT_FAILURE;
    }

//    Forces = compute_Nodal_Forces(ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep,
//                                  NumTimeStep);

//    Reactions = solve_Nodal_Equilibrium(Lumped_Mass, Gravity_field, Forces,
//                                        D_Displacement, MPM_Mesh, FEM_Mesh,
//                                        ActiveNodes, Free_and_Restricted_Dofs);

    compute_Explicit_Newmark_Corrector(MPM_Mesh, gamma);


//    output_selector(MPM_Mesh, FEM_Mesh, ActiveNodes, D_Displacement,
//                    Forces, Reactions, DeltaTimeStep, TimeStep,
//                    ResultsTimeStep);

    free(Lumped_Mass);
    free(Gravity_field);
    free(D_Displacement);
    free(Forces);
    free(Reactions);
    free(ActiveNodes.Nodes2Mask);
    free(Free_and_Restricted_Dofs.Nodes2Mask);
    print_Status("DONE !!!", TimeStep);
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __mass_NODES(
    double * Lumped_MassMatrix,
    Particle MPM_Mesh,
    Mesh FEM_Mesh,
    Mask ActiveNodes)
{

  int STATUS = EXIT_FAILURE;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned Ndim = NumberDimensions;
  unsigned NumGP = MPM_Mesh.NumGP;
  unsigned Size_connectivity_p;
  unsigned Order = Ndim * Nnodes_mask;
  unsigned Ap;
  unsigned A_mask;
  unsigned idx_A_mask_i;
  Matrix N_p;
  Element Nodes_p;
  double N_pA;
  double m_p;
  double m_A_p;  


  Lumped_MassMatrix = (double *)calloc(Order, __SIZEOF_DOUBLE__);
  if (Lumped_MassMatrix == NULL) {
    fprintf(stderr, "" RED "Out of memory" RESET "\n");
    return EXIT_FAILURE;
  }

  for (unsigned p = 0; p < NumGP; p++) {

    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    m_p = MPM_Mesh.Phi.mass.nV[p];

    Size_connectivity_p = Nodes_p.NumberNodes;

    for (unsigned A = 0; A < Size_connectivity_p; A++) {

      Ap = Nodes_p.Connectivity[A];

      A_mask = ActiveNodes.Nodes2Mask[Ap];

      N_pA = N_p.nV[A];

      m_A_p = m_p * N_pA;

      for (unsigned i = 0; i < Ndim; i++) {
        Lumped_MassMatrix[A_mask * Ndim + i] += m_A_p;
      }
    }

    free__MatrixLib__(N_p);
    free(Nodes_p.Connectivity);
  }

  return STATUS;
}

/**************************************************************/

static int __predictor_PARTICLES(
  double * D_dis,
  double * vel,
  const double * acc, 
  double gamma,
  unsigned NumGP)
{

  unsigned Ndim = NumberDimensions;
  unsigned idx;

  for (unsigned p = 0; p < NumGP; p++) {

    for (unsigned i = 0; i < Ndim; i++) {

      idx = p*Ndim + i;

      D_dis[idx] = DeltaTimeStep * vel[idx] + 0.5 * DSQR(DeltaTimeStep) * acc[idx];

      vel[idx] += (1 - gamma) * DeltaTimeStep * acc[idx];
    }
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __gravity_NODES(double * Gravity_field,
                    const Mask ActiveNodes, 
                    const Load *B,
                    unsigned NumBodyForces,
                    unsigned TimeStep, 
                    unsigned NumTimeStep)
{
  // Define auxilar variables */
  unsigned Ndim = NumberDimensions;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  double m_p; 
#if NumberDimensions == 2
  double b[2] = {0.0,0.0};
#else
  double b[3] = {0.0,0.0,0.0};
#endif 

  Gravity_field = (double *)calloc(Nnodes_mask*Ndim,__SIZEOF_DOUBLE__);
  if (Gravity_field == NULL) {
    fprintf(stderr, "" RED "Out of memory" RESET "\n");
    return EXIT_FAILURE;
  }  

  // Fill vector b of body acclerations
  for (unsigned i = 0; i < NumBodyForces; i++) {
    for (unsigned k = 0; k < Ndim; k++) {
      if (B[i].Dir[k * NumTimeStep + TimeStep]) {
        b[k] += B[i].Value[k].Fx[TimeStep];
      }
    }
  }

  // Add nodal contribution
  for (unsigned A = 0; A < Nnodes_mask; A++) {
    for (unsigned k = 0; k < Ndim; k++) {
      Gravity_field[A*Ndim + k] = b[k];
    }
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static int __d_displacement_NODES(
  double * D_dis_NODES,
  const double * Mass_NODES,
  const double * D_dis_PARTICLES,
  const double * Mass_PARTICLES,  
  Particle MPM_Mesh, 
  Mesh FEM_Mesh,
  const Mask ActiveNodes)
{

  // Define auxilar variables */
  unsigned Ndim = NumberDimensions;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned NumGP = MPM_Mesh.NumGP;
  unsigned Size_connectivity_p;
  unsigned Order = Nnodes_mask * Ndim;
  unsigned Ap;
  unsigned A_mask;
  Element Nodes_p;         
  Matrix N_p;  
  double N_pA; 
  double m_p;              

  D_dis_NODES = (double *)calloc(Nnodes_mask*Ndim,__SIZEOF_DOUBLE__);
  if (D_dis_NODES == NULL) {
    fprintf(stderr, "" RED "Out of memory" RESET "\n");
    return EXIT_FAILURE;
  } 

  // Compute N_Ip * m_p * d_dis_p
  for (unsigned p = 0; p < NumGP; p++) {

    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    m_p = Mass_PARTICLES[p];

    Size_connectivity_p = Nodes_p.NumberNodes;

    for (unsigned A = 0; A < Size_connectivity_p; A++) {

      Ap = Nodes_p.Connectivity[A];

      A_mask = ActiveNodes.Nodes2Mask[Ap];

      N_pA = N_p.nV[A];

      for (unsigned i = 0; i < Ndim; i++) {
        D_dis_NODES[A_mask*Ndim + i] += m_p * N_pA * D_dis_PARTICLES[p*Ndim + i];
      }
    }

    free__MatrixLib__(N_p);
    free(Nodes_p.Connectivity);
  }

  // Solve M_IJ * d_dis_J = N_Ip * m_p * d_dis_p
  for (int A = 0; A < Nnodes_mask; A++) {
    for (int i = 0; i < Ndim; i++) {
      D_dis_NODES[A*Ndim + i] = D_dis_NODES[A*Ndim + i] / Mass_NODES[A*Ndim + i];
    }
  }


  return EXIT_SUCCESS;
}

/**************************************************************/

static Matrix compute_Nodal_Velocity(Particle MPM_Mesh, Mesh FEM_Mesh,
                                     Mask ActiveNodes, Matrix Lumped_Mass)
/*
  Compute the nodal velocity. The operation is linearized and
  all the dof split the velocity array in n components like :
  | M 0 |   |V.x|   | p.x |
  | 0 M | * |V.y| = | p.y |

*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Order = Ndim * Nnodes_mask;
  int Ap;
  int A_mask;
  int AB;
  Element Nodes_p;         /* Element for each particle */
  Matrix ShapeFunction_p;  /* Value of the shape-function */
  double ShapeFunction_pA; /* Evaluation of the particle in the node */
  double m_p;              /* Mass of the particle */

  /* Define and allocate the velocity vector */
  Matrix Velocity = allocZ__MatrixLib__(Nnodes_mask, Ndim);

  for (int p = 0; p < Np; p++) {

    /*
       Define element of the particle
    */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],
                                     MPM_Mesh.NumberNodes[p]);

    /*
      Evaluate the shape function in the coordinates of the particle
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Get the mass of the GP
    */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /*
      Get the nodal mommentum
    */
    for (int A = 0; A < Nodes_p.NumberNodes; A++) {
      /*
        Get the node with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Evaluate the GP function in the node
      */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      for (int i = 0; i < Ndim; i++) {
        Velocity.nM[A_mask][i] +=
            m_p * ShapeFunction_pA * MPM_Mesh.Phi.vel.nM[p][i];
      }
    }

    /*
      Free the value of the shape functions
    */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Compute the nodal velocities
  */
  for (int A = 0; A < Nnodes_mask; A++) {
    for (int i = 0; i < Ndim; i++) {
      AB = A * Ndim + i;
      Velocity.nM[A][i] = Velocity.nM[A][i] / Lumped_Mass.nM[AB][AB];
    }
  }


  return Velocity;
}

/**************************************************************/

static void impose_Dirichlet_Boundary_Conditions(
  double * D_Displacement,
  Mesh FEM_Mesh, 
  Mask ActiveNodes,
  int TimeStep, 
  int NumTimeStep)
/*
  Apply the boundary conditions over the nodes
*/
{

  /* 1º Define auxilar variables */
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int NumNodesBound;           /* Number of nodes of the bound */
  int Ndim = NumberDimensions; /* Number of dimensions */
  int Ndof = NumberDOF;        /* Number of degree of freedom */
  int Id_BCC;                  /* Index of the node where we apply the BCC */
  int Id_BCC_mask;

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
            Initialise increments using newmark and the value of the boundary
            condition
          */
          D_Displacement[Id_BCC_mask*Ndim + k] =
              FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep];
        }
      }
    }
  }
}

/**************************************************************/

static int __update_Local_State(
  const double * D_Displacement,
  Mask ActiveNodes,
  Particle MPM_Mesh, 
  Mesh FEM_Mesh) {

  /*
    Auxiliar variables
  */
  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Nnodes_p;
  int MatIndx_p;
  unsigned NumberNodes_p;
  unsigned Order_p;  
  int Idx_Element_p;
  int Idx_Patch_p;
  double rho_n_p;
  double Delta_J_p;
  double Vn_patch;
  double Vn1_patch;
  double J_patch;
  Element Nodes_p;
  Material MatProp_p;
  Matrix gradient_p;
  double * D_Displacement_Ap;
  double * F_n_p;
  double * F_n1_p;
  double * DF_p;

  /*
    Loop in the material point set to update kinematics
  */
  for (unsigned p = 0; p < Np; p++) {

    //  Define tributary nodes of the particle
    NumberNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],NumberNodes_p);
    Order_p = NumberNodes_p*Ndim;

    //  Get the nodal increment of displacement using the mask
    D_Displacement_Ap = (double *)calloc(Order_p, __SIZEOF_DOUBLE__);
      if(D_Displacement_Ap == NULL){
        fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
        return EXIT_FAILURE;
      } 
    get_set_field__MeshTools__(D_Displacement_Ap, D_Displacement, Nodes_p, ActiveNodes);

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

    /*
      Compute the increment of the deformation gradient
    */
      compute_strains_ctx ctx_shape_fun;
      ctx_shape_fun.gradient_p = gradient_p.nV;
      ctx_shape_fun.Nnodes_p = NumberNodes_p;

      update_increment_Deformation_Gradient__Particles__(
          DF_p, D_Displacement_Ap, &ctx_shape_fun);

    /*
      Update the deformation gradient in t = n + 1 with the information
      from t = n and the increment of deformation gradient.
    */
    update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, DF_p);

    /*
      Compute Jacobian of the deformation gradient
    */
    MPM_Mesh.Phi.J_n1.nV[p] = I3__TensorLib__(F_n1_p);

    if (MPM_Mesh.Phi.J_n1.nV[p] <= 0.0) {
      fprintf(stderr, ""RED"%s : %s %i"RESET" \n", 
      "Error in I3__TensorLib__(F_n1_p)",
      "Negative jacobian in particle", p);
      return EXIT_FAILURE;
    }

    /*
      Update patch
    */
    if (FEM_Mesh.Locking_Control_Fbar) {
      Idx_Element_p = MPM_Mesh.Element_p[p];
      Idx_Patch_p = FEM_Mesh.Idx_Patch[Idx_Element_p];
      FEM_Mesh.Vol_Patch_n[Idx_Patch_p] +=
          MPM_Mesh.Phi.J_n.nV[p] * MPM_Mesh.Phi.Vol_0.nV[p];
      FEM_Mesh.Vol_Patch_n1[Idx_Patch_p] +=
          MPM_Mesh.Phi.J_n1.nV[p] * MPM_Mesh.Phi.Vol_0.nV[p];
    }

    /*
      Update density with the jacobian of the increment deformation gradient
    */
    Delta_J_p = I3__TensorLib__(DF_p);
    rho_n_p = MPM_Mesh.Phi.rho.nV[p];
    MPM_Mesh.Phi.rho.nV[p] = rho_n_p / Delta_J_p;

    /*
      Free memory
    */
    free(D_Displacement_Ap);
    free__MatrixLib__(gradient_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Loop in the material point set to update stress
  */
  for (int p = 0; p < Np; p++) {
    MatIndx_p = MPM_Mesh.MatIdx[p];
    MatProp_p = MPM_Mesh.Mat[MatIndx_p];

    if (FEM_Mesh.Locking_Control_Fbar) {
      Idx_Element_p = MPM_Mesh.Element_p[p];
      Idx_Patch_p = FEM_Mesh.Idx_Patch[Idx_Element_p];

      Vn_patch = FEM_Mesh.Vol_Patch_n[Idx_Patch_p];
      Vn1_patch = FEM_Mesh.Vol_Patch_n1[Idx_Patch_p];
      J_patch = Vn1_patch / Vn_patch;

      get_locking_free_Deformation_Gradient_n1__Particles__(p, J_patch,
                                                            MPM_Mesh);

      MPM_Mesh.Phi.Jbar.nV[p] *= J_patch;
    }

    /*
      Update the first Piola-Kirchhoff stress tensor with an apropiate
      integration rule.
    */
    STATUS = Stress_integration__Constitutive__(p, MPM_Mesh, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, ""RED"Error in Stress_integration__Constitutive__(%i,,,)"RESET" \n",p);
      return EXIT_FAILURE;
    }

  }

  return STATUS;
}

/**************************************************************/

static Matrix compute_Nodal_Forces(Mask ActiveNodes, Particle MPM_Mesh,
                                   Mesh FEM_Mesh, int TimeStep,
                                   int NumTimeStep) {
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  Matrix Forces = allocZ__MatrixLib__(Nnodes_mask, Ndim);

  /*
    Add internal forces contribution
  */
  compute_Nodal_Internal_Forces(Forces, ActiveNodes, MPM_Mesh, FEM_Mesh);

  /*
    Add contact forces contribution
  */
  compute_Nodal_Nominal_traction_Forces(Forces, ActiveNodes, MPM_Mesh, FEM_Mesh,
                                        TimeStep, NumTimeStep);

  return Forces;
}

/**************************************************************/

static void compute_Nodal_Internal_Forces(Matrix Forces, Mask ActiveNodes,
                                          Particle MPM_Mesh, Mesh FEM_Mesh) {

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int NumNodes_p;
  int idx_A_mask_i;

  Tensor P_p; /* First Piola-Kirchhoff Stress tensor */
  Tensor InternalForcesDensity_Ap;

  Element Nodes_p;   /* List of nodes for particle */
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pA;
  Tensor GRADIENT_pA;
  Tensor F_n_p;
  Tensor transpose_F_n_p;

  double V0_p; /* Volume of the Gauss-Point */

  /*
    Loop in the particles
  */
  for (int p = 0; p < Np; p++) {
    /*
      Get the volume of the particle in the reference configuration
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /*
      Define nodes for each particle
    */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /*
      Compute gradient of the shape function in each node
    */
    gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Take the values of the deformation gradient ant t = n and t = n + 1.
      Later compute the midpoint deformation gradient and
      the transpose of the deformation gradient.
    */
    F_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p], 2);
    transpose_F_n_p = transpose__TensorLib__(F_n_p);

    /*
      Get the first Piola-Kirchhoff stress tensor.
    */
    Tensor P_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);

    for (int A = 0; A < NumNodes_p; A++) {

      /*
        Compute the gradient in the reference configuration
      */
      gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);
      GRADIENT_pA =
          vector_linear_mapping__TensorLib__(transpose_F_n_p, gradient_pA);

      /*
    Compute the nodal forces of the particle
  */
      InternalForcesDensity_Ap =
          vector_linear_mapping__TensorLib__(P_p, GRADIENT_pA);

      /*
    Get the node of the mesh for the contribution
  */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
    Asign the nodal forces contribution to the node
  */
      for (int i = 0; i < Ndim; i++) {
        Forces.nM[A_mask][i] -= InternalForcesDensity_Ap.n[i] * V0_p;
      }

      /*
    Free memory
  */
      free__TensorLib__(InternalForcesDensity_Ap);
      free__TensorLib__(GRADIENT_pA);
    }

    /*
      Free memory
    */
    free__TensorLib__(transpose_F_n_p);
    free__MatrixLib__(gradient_p);
    free(Nodes_p.Connectivity);
  }
}

/**************************************************************/

static void compute_Nodal_Nominal_traction_Forces(Matrix Forces,
                                                  Mask ActiveNodes,
                                                  Particle MPM_Mesh,
                                                  Mesh FEM_Mesh, int TimeStep,
                                                  int NumTimeStep) {

  int Ndim = NumberDimensions;
  Load Load_i;
  Element Nodes_p;        /* Element for each Gauss-Point */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pA;
  Tensor T = alloc__TensorLib__(1); // Nominal traction
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
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
        Fill vector of contact forces
      */
      for (int k = 0; k < Ndim; k++) {
        if (Load_i.Dir[k * NumTimeStep + TimeStep] == 1) {
          T.n[k] = Load_i.Value[k].Fx[TimeStep];
        }
      }

      /*
        Get the node of the mesh for the contribution
      */
      for (int A = 0; A < NumNodes_p; A++) {

        /*
          Pass the value of the nodal shape function to a scalar
        */
        ShapeFunction_pA = ShapeFunction_p.nV[A];

        /*
          Node for the contribution
        */
        Ap = Nodes_p.Connectivity[A];
        A_mask = ActiveNodes.Nodes2Mask[Ap];

        /*
          Compute Contact forces
        */
        for (int k = 0; k < Ndim; k++) {
          Forces.nM[A_mask][k] += ShapeFunction_pA * T.n[k] * A0_p;
        }
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }
  }

  free__TensorLib__(T);
}

/**************************************************************/

static Matrix solve_Nodal_Equilibrium(Matrix Lumped_Mass, Matrix Gravity_field,
                                      Matrix Total_Forces,
                                      Matrix D_Displacement, Particle MPM_Mesh,
                                      Mesh FEM_Mesh, Mask ActiveNodes,
                                      Mask Free_and_Restricted_Dofs)
/*
  Call the LAPACK solver to compute the accelerations and velocities
  Solve equilibrium equation to get the nodal values of the aceleration
  at t = n + 1
*/
{

  /*
    General varibles
  */
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Np = MPM_Mesh.NumGP;
  int Nnodes = ActiveNodes.Nactivenodes;
  int NumNodes_p;
  int Order = Nnodes * Ndim;
  int AB;
  int Ap;
  int A_mask;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix ShapeFunction_p;
  double ShapeFunction_pA;

  /*
    Solution nodal variable
  */
  Matrix Acceleration = allocZ__MatrixLib__(Nnodes, Ndim);

  /*
    Output
  */
  Matrix Reactions = allocZ__MatrixLib__(Nnodes, Ndim);

  /*
    The solution is now stored in the internal forces vector
  */
  for (int A = 0; A < Nnodes; A++) {
    for (int i = 0; i < Ndim; i++) {
      if (Free_and_Restricted_Dofs.Nodes2Mask[A * Ndof + i] != -1) {
        AB = A * Ndim + i;
        Acceleration.nM[A][i] = Gravity_field.nM[A][i] +
                                Total_Forces.nM[A][i] / Lumped_Mass.nM[AB][AB];
      } else {
        Acceleration.nM[A][i] = 0.0;
        Reactions.nM[A][i] = Total_Forces.nM[A][i];
      }
    }
  }

  /*
    Update particle acceleration
  */
  for (int p = 0; p < Np; p++) {

    /*
      Set to zero particle acceleration for interpolation
    */
    for (int i = 0; i < Ndim; i++) {
      MPM_Mesh.Phi.acc.nM[p][i] = 0.0;
      MPM_Mesh.Phi.D_dis.nM[p][i] = 0.0;
    }

    /*
      Define element of the particle
    */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /*
      Evaluate the shape function in the coordinates of the particle
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Interpolate the new value of the acceleration
    */
    for (int A = 0; A < NumNodes_p; A++) {
      /*
        Get the node in the nodal momentum with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Evaluate the GP function in the node
      */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      for (int i = 0; i < Ndim; i++) {
        MPM_Mesh.Phi.acc.nM[p][i] +=
            ShapeFunction_pA * Acceleration.nM[A_mask][i];
        MPM_Mesh.Phi.D_dis.nM[p][i] +=
            ShapeFunction_pA * D_Displacement.nM[A_mask][i];
      }
    }

    /*
      Free auxiliar data
    */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Free acceleration vector
  */
  free__MatrixLib__(Acceleration);

  return Reactions;
}

/**************************************************************/

static void compute_Explicit_Newmark_Corrector(Particle MPM_Mesh,
                                               double gamma) {

  unsigned Np = MPM_Mesh.NumGP;

#if NumberDimensions == 2
unsigned Size_vector = 2;
unsigned Size_tensor = 5;
#else
unsigned Size_vector = 3;
unsigned Size_tensor = 9;
#endif
  

  for (unsigned p = 0; p < Np; p++) {

    /*
      Replace the determinant of the deformation gradient
    */
    MPM_Mesh.Phi.J_n.nV[p] = MPM_Mesh.Phi.J_n1.nV[p];

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

    /*
      Update/correct vector variables
    */
    for (unsigned i = 0; i < Size_vector; i++) {
      /*
        Correct particle velocity
      */
      MPM_Mesh.Phi.vel.nM[p][i] +=
          gamma * DeltaTimeStep * MPM_Mesh.Phi.acc.nM[p][i];

      /*
        Update the particles position and displacement
      */
      MPM_Mesh.Phi.x_GC.nM[p][i] += MPM_Mesh.Phi.D_dis.nM[p][i];
      MPM_Mesh.Phi.dis.nM[p][i] += MPM_Mesh.Phi.D_dis.nM[p][i];
  
    }

    /*
      Update/correct tensor variables
    */
    for (unsigned i = 0; i < Size_tensor; i++) {
      MPM_Mesh.Phi.F_n.nM[p][i] = MPM_Mesh.Phi.F_n1.nM[p][i];
    }

  }
}

/**************************************************************/

static void output_selector(Particle MPM_Mesh, Mesh FEM_Mesh, Mask ActiveNodes,
                            Matrix Velocity, Matrix D_Displacement,
                            Matrix Forces, Matrix Reactions,
                            double DeltaTimeStep, int TimeStep,
                            int ResultsTimeStep) {

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

    particle_results_vtk__InOutFun__(MPM_Mesh, TimeStep, ResultsTimeStep);

    nodal_results_vtk__InOutFun__(FEM_Mesh, ActiveNodes, Reactions, TimeStep,
                                  ResultsTimeStep);


#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  if (Out_Partition_Unity) {
      free__MatrixLib__(MPM_Mesh.Phi.PU);
    }

  #endif
#endif
    


  }

  /*
    csv results
  */
  for (int i = 0; i < Number_Out_nodal_path_csv; i++) {

    if (Out_nodal_path_csv[i].Out_csv_nodes_path_Velocity) {
      path_nodes_analysis_csv__InOutFun__(
          Velocity, FEM_Mesh.Coordinates, "Nodal_path_velocity_csv",
          ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if (Out_nodal_path_csv[i].Out_csv_nodes_path_D_Displacement) {
      path_nodes_analysis_csv__InOutFun__(
          D_Displacement, FEM_Mesh.Coordinates, "Nodal_path_displacement_csv",
          ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if (Out_nodal_path_csv[i].Out_csv_nodes_path_Forces) {
      path_nodes_analysis_csv__InOutFun__(
          Forces, FEM_Mesh.Coordinates, "Nodal_path_forces_csv", ActiveNodes,
          Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if (Out_nodal_path_csv[i].Out_csv_nodes_path_Reactions) {
      path_nodes_analysis_csv__InOutFun__(
          Reactions, FEM_Mesh.Coordinates, "Nodal_path_reactions_csv",
          ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }
  }


  for (int i = 0; i < Number_Out_particles_path_csv; i++) {

/*
    if (Out_particles_path_csv[i].Out_csv_particles_path_Damage) {
      path_particles_analysis_csv__InOutFun__(
          MPM_Mesh.Phi.chi, MPM_Mesh.Phi.x_GC, "Particles_path_damage_csv",
          Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }
*/

    if (Out_particles_path_csv[i].Out_csv_particles_path_Velocity) {
      path_particles_analysis_csv__InOutFun__(
          MPM_Mesh.Phi.vel, MPM_Mesh.Phi.x_GC, "Particles_path_velocity_csv",
          Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if (Out_particles_path_csv[i].Out_csv_particles_path_Acceleration) {
      path_particles_analysis_csv__InOutFun__(
          MPM_Mesh.Phi.acc, MPM_Mesh.Phi.x_GC,
          "Particles_path_acceleration_csv", Out_particles_path_csv[i], i,
          TimeStep, DeltaTimeStep);
    }

    if (Out_particles_path_csv[i].Out_csv_particles_path_Displacement) {
      path_particles_analysis_csv__InOutFun__(
          MPM_Mesh.Phi.dis, MPM_Mesh.Phi.x_GC,
          "Particles_path_displacement_csv", Out_particles_path_csv[i], i,
          TimeStep, DeltaTimeStep);
    }

    if (Out_particles_path_csv[i].Out_csv_particles_path_Stress) {
      path_particles_analysis_csv__InOutFun__(
          MPM_Mesh.Phi.Stress, MPM_Mesh.Phi.x_GC, "Particles_path_stress_csv",
          Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if (Out_particles_path_csv[i].Out_csv_particles_path_Strain) {
      path_particles_analysis_csv__InOutFun__(
          MPM_Mesh.Phi.Strain, MPM_Mesh.Phi.x_GC, "Particles_path_strain_csv",
          Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if (Out_particles_path_csv[i].Out_csv_particles_path_Deformation_gradient) {
      path_particles_analysis_csv__InOutFun__(
          MPM_Mesh.Phi.F_n, MPM_Mesh.Phi.x_GC,
          "Particles_path_deformation_gradient_csv", Out_particles_path_csv[i],
          i, TimeStep, DeltaTimeStep);
    }
  }
}

/**************************************************************/