#include <math.h>
#include <string.h>
#include "nl-partsol.h"

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

typedef struct {

  Matrix value;
  Matrix d_value_dt;
  Matrix d2_value_dt2;

} Nodal_Field;

typedef struct {
  double alpha_1;
  double alpha_2;
  double alpha_3;
  double alpha_4;
  double alpha_5;
  double alpha_6;
} Newmark_parameters;

/*
  Auxiliar functions
*/
static Newmark_parameters compute_Newmark_parameters(double, double, double);
static Matrix compute_Nodal_Effective_Mass(Particle, Mesh, Mask, double);
static Nodal_Field compute_Nodal_Field(Particle, Mesh, Mask);
static Nodal_Field initialise_Nodal_Increments(Nodal_Field, Mesh, Mask,
                                               Newmark_parameters, int, int);

static int __local_deformation(
  Nodal_Field D_U, 
  Mask ActiveNodes,
  Particle MPM_Mesh, 
  Mesh FEM_Mesh,
  double TimeStep);

static Matrix compute_Nodal_Forces(Mask, Particle, Mesh, int, int);

static void compute_Nodal_Internal_Forces(Matrix, Mask, Particle, Mesh);

static void compute_Nodal_Nominal_traction_Forces(Matrix, Mask, Particle, Mesh,
                                                  int, int);

static void compute_Nodal_Body_Forces(Matrix, Mask, Particle, Mesh, int, int);

static Matrix compute_Nodal_Reactions(Mesh, Matrix, Mask, int, int);

static Matrix compute_Nodal_Residual(
  Nodal_Field U_n /**< */, 
  Nodal_Field D_U /**< */,
  Mask ActiveNodes /**< */, 
  Matrix Forces /**< */,
  Matrix Mass /**< */, 
  Newmark_parameters Params /**< */);

static double __error_residual(
  double * Residual /**< */,
  int Total_dof /**< */);

static int __assemble_tangent_stiffness(
  double * Tangent_Stiffness /**< */,
  const double * Effective_Mass /**< */,
  Mask ActiveNodes /**< */,
  Particle MPM_Mesh /**< */, 
  Mesh FEM_Mesh /**< */,
  Newmark_parameters Params /**< */);

static void system_reduction(
  double * Tangent_Stiffness /**< */,
  Mask ActiveNodes /**< */, 
  Mesh FEM_Mesh /**< */, 
  int TimeStep /**< */,
  int NumTimeStep /**< */);

static int __solve_equilibrium(
  Nodal_Field D_U /**< */, 
  double * Tangent_Stiffness /**< */,
  double * Residual /**< */,
  int Nnodes_mask /**< */);

static void __update_Nodal_Increments(
  Nodal_Field D_U /**< */,
  Nodal_Field U_n /**< */,
  Mask Free_and_Restricted_Dofs /**< */,
  Newmark_parameters Params /**< */);

static void __update_Particles(
  Nodal_Field D_U /**< */,
  Particle MPM_Mesh /**< */, 
  Mesh FEM_Mesh /**< */,
  Mask ActiveNodes /**< */);

static void output_selector(Particle, Mesh, Mask, Nodal_Field, Nodal_Field,
                            Matrix, Matrix, Matrix, int, int);
/**************************************************************/

int U_Newmark_beta_Finite_Strains(Mesh FEM_Mesh, Particle MPM_Mesh,
                                   Time_Int_Params Parameters_Solver) {

  int STATUS = EXIT_SUCCESS;

  /*
    Auxiliar variables for the solver
  */
  unsigned Ndim = NumberDimensions;
  unsigned Nactivenodes;
  unsigned Nactivedofs;
  unsigned InitialStep = Parameters_Solver.InitialTimeStep;
  unsigned NumTimeStep = Parameters_Solver.NumTimeStep;
  unsigned MaxIter = Parameters_Solver.MaxIter;
  unsigned Iter;
  unsigned Order;

  unsigned TimeStep = InitialStep;

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

  Matrix Effective_Mass;
  double * Tangent_Stiffness;
  Matrix Forces;
  Matrix Reactions;
  Nodal_Field U_n;
  Nodal_Field D_U;
  Matrix D_Displacement;
  Matrix Residual;

  Mask ActiveNodes;
  Mask Free_and_Restricted_Dofs;

  Newmark_parameters Params;

  /*
    Time step is defined at the init of the simulation throught the
    CFL condition. Notice that for this kind of solver, CFL confition is
    not required to be satisfied. The only purpose of it is to use the existing
    software interfase.
  */
  DeltaTimeStep = U_DeltaT__SolversLib__(MPM_Mesh, DeltaX, Parameters_Solver);

  /*
    Compute alpha parameters
  */
  Params = compute_Newmark_parameters(beta, gamma, DeltaTimeStep);

  while (TimeStep < NumTimeStep) {
    print_Status("*************************************************", TimeStep);
    print_step(TimeStep, DeltaTimeStep);

    local_search__MeshTools__(MPM_Mesh, FEM_Mesh);
    ActiveNodes = generate_NodalMask__MeshTools__(FEM_Mesh);
    Nactivenodes = ActiveNodes.Nactivenodes;
    Free_and_Restricted_Dofs =
        generate_Mask_for_static_condensation__MeshTools__(
            ActiveNodes, FEM_Mesh, TimeStep, NumTimeStep);
    Nactivedofs = Free_and_Restricted_Dofs.Nactivenodes;

    Order = NumberDimensions * Nactivenodes;

    // Compute kinematic nodal values
    U_n = compute_Nodal_Field(MPM_Mesh, FEM_Mesh, ActiveNodes);
    D_U = initialise_Nodal_Increments(U_n, FEM_Mesh, ActiveNodes, Params, TimeStep, NumTimeStep);
    Effective_Mass = compute_Nodal_Effective_Mass(MPM_Mesh, FEM_Mesh, ActiveNodes, epsilon);

    // First trial    
    Forces = compute_Nodal_Forces(ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep, NumTimeStep);
    Reactions = compute_Nodal_Reactions(FEM_Mesh, Forces, ActiveNodes, TimeStep, NumTimeStep);
    Residual = compute_Nodal_Residual(U_n, D_U, ActiveNodes, Forces, Effective_Mass, Params);

    // Compute error
    Error_0 = Error_i = __error_residual(Residual.nV,Order);  
    Error_relative = Error_i/Error_0;
    Iter = 0;

    while (Error_relative > TOL) {

      if ((Error_i < TOL * 100) 
      || (Error_relative < TOL) 
      || (Iter > MaxIter)) {
        break;
      }

      Tangent_Stiffness = (double *)calloc(Order * Order, __SIZEOF_DOUBLE__);
      if(STATUS == EXIT_FAILURE){
        fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
        return EXIT_FAILURE;
      } 

      STATUS = __assemble_tangent_stiffness(
        Tangent_Stiffness, Effective_Mass.nV,
        ActiveNodes, MPM_Mesh, FEM_Mesh, Params);
      if(STATUS == EXIT_FAILURE){
          fprintf(stderr, ""RED"Error in __assemble_tangent_stiffness()"RESET" \n");
          return EXIT_FAILURE;
      } 

      system_reduction(Tangent_Stiffness, ActiveNodes, FEM_Mesh, TimeStep, NumTimeStep);

      STATUS = __solve_equilibrium(D_U, Tangent_Stiffness, Residual.nV, Nactivenodes);
      if(STATUS == EXIT_FAILURE){
        fprintf(stderr, ""RED"Error in __solve_equilibrium()"RESET" \n");
        return EXIT_FAILURE;
      }         

      __update_Nodal_Increments(D_U, U_n, Free_and_Restricted_Dofs, Params);

      STATUS = __local_deformation(D_U, ActiveNodes, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
      if(STATUS == EXIT_FAILURE){
        fprintf(stderr, ""RED"Error in __local_deformation()"RESET" \n");
        return EXIT_FAILURE;
      } 

      free__MatrixLib__(Forces);
      free__MatrixLib__(Reactions);
      free__MatrixLib__(Residual);
      free(Tangent_Stiffness);

      Forces = compute_Nodal_Forces(ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep, NumTimeStep);
      Reactions = compute_Nodal_Reactions(FEM_Mesh, Forces, ActiveNodes, TimeStep, NumTimeStep);
      Residual = compute_Nodal_Residual(U_n, D_U, ActiveNodes, Forces, Effective_Mass, Params);

      Error_i = __error_residual(Residual.nV,Order);
      Error_relative = Error_i/Error_0;
      Iter++;
      printf("Iter: [%i/%i]. Total Error: %e, Relative Error: %e \n", 
      Iter, MaxIter, Error_i, Error_relative);

    }
    
    print_convergence_stats(TimeStep, Iter, Error_0, Error_i, Error_relative);

    if(Iter > MaxIter)
    {
      fprintf(stderr, ""RED"Convergence not reached in the maximum number of iterations"RESET" \n");
      return EXIT_FAILURE;
    }

    __update_Particles(D_U, MPM_Mesh, FEM_Mesh, ActiveNodes);

    output_selector(MPM_Mesh, FEM_Mesh, ActiveNodes, D_U, U_n, Forces,
                    Reactions, Residual, TimeStep, ResultsTimeStep);

    TimeStep++;

    free__MatrixLib__(Effective_Mass);
    free__MatrixLib__(U_n.value);
    free__MatrixLib__(U_n.d_value_dt);
    free__MatrixLib__(U_n.d2_value_dt2);
    free__MatrixLib__(D_U.value);
    free__MatrixLib__(D_U.d_value_dt);
    free__MatrixLib__(D_U.d2_value_dt2);
    free__MatrixLib__(Forces);
    free__MatrixLib__(Reactions);
    free__MatrixLib__(Residual);
    free(ActiveNodes.Nodes2Mask);
    free(Free_and_Restricted_Dofs.Nodes2Mask);

//    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static Newmark_parameters compute_Newmark_parameters(double beta, double gamma,
                                                     double DeltaTimeStep) {
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

static Matrix compute_Nodal_Effective_Mass(Particle MPM_Mesh, Mesh FEM_Mesh,
                                           Mask ActiveNodes, double epsilon)
/*
  This function computes the effective mass matrix as a convex combination
  of the lumped mass matrix and the consistent mass matrix. Later assemble
  a total mass matrix with the contribution of each degree of freedom.

  | M_eff |   0   |              | M_cons |   0    |          | M_lump |   0 |
  -----------------  = (1-eps) * -------------------  + eps *
  ------------------- |    0  | M_eff |	         |   0    | M_cons |	      |
  0    | M_lump |
*/
{

  int Nnodes_mask = ActiveNodes.Nactivenodes;
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

  /* Define and allocate the effective mass matrix */
  Matrix Effective_MassMatrix = allocZ__MatrixLib__(Order, Order);

  /* Define and allocate the lumped mass matrix */
  Matrix Lumped_MassMatrix = allocZ__MatrixLib__(Order, 1);

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
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

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

      /* Get the value of the shape function */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /*
        Compute the nodal A contribution of the particle p
      */
      m_A_p = m_p * ShapeFunction_pA;

      /*
         Fill the Lumped mass matrix considering the number of dofs
      */
      for (int i = 0; i < Ndof; i++) {
        Lumped_MassMatrix.nV[A_mask * Ndof + i] += m_A_p;
      }

      for (int B = 0; B < Nodes_p.NumberNodes; B++) {
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
        for (int i = 0; i < Ndof; i++) {
          /*
            Compute the vectorized index
          */
          Effective_MassMatrix.nM[A_mask * Ndof + i][B_mask * Ndof + i] +=
              m_AB_p;
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
      Effective_MassMatrix.nM[A][B] =
          (1 - epsilon) * Effective_MassMatrix.nM[A][B] +
          (A == B) * epsilon * Lumped_MassMatrix.nV[A];
    }
  }

  /*
    Free lumped mass matrix.
  */
  free__MatrixLib__(Lumped_MassMatrix);


  return Effective_MassMatrix;
}

/**************************************************************/

static Nodal_Field compute_Nodal_Field(Particle MPM_Mesh,
                                       Mesh FEM_Mesh, Mask ActiveNodes)
/*
  Call the LAPACK solver to compute the nodal velocity. The operation is
  linearized and all the dof split the velocity array in n components like : | M
  0 |   |V.x|   | p.x | | 0 M | * |V.y| = | p.y |

*/
{
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
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

  /* Define and allocate the output vector */
  Nodal_Field U_n;
  U_n.value = allocZ__MatrixLib__(Nnodes_mask, Ndof);
  U_n.d_value_dt = allocZ__MatrixLib__(Nnodes_mask, Ndof);
  U_n.d2_value_dt2 = allocZ__MatrixLib__(Nnodes_mask, Ndof);

  /* Iterate over the particles to get the nodal values */
  for (int p = 0; p < Np; p++) {

    /* Define element of the particle */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],
                                     MPM_Mesh.NumberNodes[p]);

    /* Evaluate the shape function in the coordinates of the particle */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Get the mass of the GP */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* Get the nodal mommentum */
    for (int A = 0; A < Nodes_p.NumberNodes; A++) {

      /*
        Get the node in the nodal momentum with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /* Evaluate the GP function in the node */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /* Nodal velocity and acceleration  */
      for (int i = 0; i < Ndim; i++) {
        U_n.value.nM[A_mask][i] +=
            m_p * ShapeFunction_pA * MPM_Mesh.Phi.dis.nM[p][i];
        U_n.d_value_dt.nM[A_mask][i] +=
            m_p * ShapeFunction_pA * MPM_Mesh.Phi.vel.nM[p][i];
        U_n.d2_value_dt2.nM[A_mask][i] +=
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
  int Order = Nnodes_mask * Ndim;
  int LDA = Order;
  int LDB = Order;
  char TRANS = 'N'; /* (Transpose) */
  int INFO = 3;
  int *IPIV = (int *)Allocate_Array(Order, sizeof(int));
  int NRHS = 1;

  /*
    Generate auxiliar copy of the mass matrix to avoid destructive operations
  */
 Matrix Effective_Mass =
        compute_Nodal_Effective_Mass(MPM_Mesh, FEM_Mesh, ActiveNodes, epsilon);

  /*
    Compute the LU factorization for the mass matrix
  */
  dgetrf_(&Order, &Order, Effective_Mass.nV, &LDA, IPIV, &INFO);

  if (INFO != 0) {
    if (INFO < 0) {
      printf("%s : \n", "Error in compute_Nodal_Field()");
      printf("the %i-th argument had an illegal value", abs(INFO));
    } else if (INFO > 0) {
      printf("%s :\n", "Error in compute_Nodal_Field()");
      printf(" M(%i,%i) %s \n %s \n %s \n %s \n", INFO, INFO,
             "is exactly zero. The factorization",
             "has been completed, but the factor M is exactly",
             "singular, and division by zero will occur if it is used",
             "to solve a system of equations.");
    }
    exit(EXIT_FAILURE);
  }

  /*
    Solve for the velocity
  */
  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass.nV, &LDA, IPIV, U_n.value.nV, &LDB,
          &INFO);
  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass.nV, &LDA, IPIV, U_n.d_value_dt.nV,
          &LDB, &INFO);
  dgetrs_(&TRANS, &Order, &NRHS, Effective_Mass.nV, &LDA, IPIV, U_n.d2_value_dt2.nV,
          &LDB, &INFO);

  /*
    Free auxiliar memory
  */
  free__MatrixLib__(Effective_Mass);
  free(IPIV);

  return U_n;
}

/**************************************************************/

static Nodal_Field initialise_Nodal_Increments(Nodal_Field U_n, Mesh FEM_Mesh,
                                               Mask ActiveNodes,
                                               Newmark_parameters Params,
                                               int TimeStep, int NumTimeStep)
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

  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

  double D_U_value_It;
  Nodal_Field D_U;
  D_U.value = allocZ__MatrixLib__(Nnodes_mask, Ndof);
  D_U.d_value_dt = allocZ__MatrixLib__(Nnodes_mask, Ndof);
  D_U.d2_value_dt2 = allocZ__MatrixLib__(Nnodes_mask, Ndof);

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
          U_n.value.nM[Id_BCC_mask][k] = 0.0;
          U_n.d_value_dt.nM[Id_BCC_mask][k] = 0.0;
          U_n.d2_value_dt2.nM[Id_BCC_mask][k] = 0.0;

          for (int t = 0; t < TimeStep; t++) {
            D_U_value_It = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[t];
            U_n.value.nM[Id_BCC_mask][k] += D_U_value_It;
            U_n.d2_value_dt2.nM[Id_BCC_mask][k] +=
                alpha_1 * D_U_value_It -
                alpha_2 * U_n.d_value_dt.nM[Id_BCC_mask][k] -
                (alpha_3 + 1) * U_n.d2_value_dt2.nM[Id_BCC_mask][k];
            U_n.d_value_dt.nM[Id_BCC_mask][k] +=
                alpha_4 * D_U_value_It +
                (alpha_5 - 1) * U_n.d_value_dt.nM[Id_BCC_mask][k] +
                alpha_6 * U_n.d2_value_dt2.nM[Id_BCC_mask][k];
          }

          /*
            Initialise increments using newmark and the value of the boundary
            condition
          */
          D_U.value.nM[Id_BCC_mask][k] =
              FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep];
          D_U.d2_value_dt2.nM[Id_BCC_mask][k] =
              alpha_1 * D_U.value.nM[Id_BCC_mask][k] -
              alpha_2 * U_n.d_value_dt.nM[Id_BCC_mask][k] -
              (alpha_3 + 1) * U_n.d2_value_dt2.nM[Id_BCC_mask][k];
          D_U.d_value_dt.nM[Id_BCC_mask][k] =
              alpha_4 * D_U.value.nM[Id_BCC_mask][k] +
              (alpha_5 - 1) * U_n.d_value_dt.nM[Id_BCC_mask][k] +
              alpha_6 * U_n.d2_value_dt2.nM[Id_BCC_mask][k];
        }
      }
    }
  }

  return D_U;
}

/**************************************************************/

static int __local_deformation(
  Nodal_Field D_U, 
  Mask ActiveNodes,
  Particle MPM_Mesh, 
  Mesh FEM_Mesh,
  double TimeStep) {

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int MatIndx_p;
  int Nnodes_p;
  int Idx_Element_p;
  int Idx_Patch_p;
  double Vn_patch;
  double Vn1_patch;
  double J_patch;
  Element Nodes_p;
  Material MatProp_p;
  Matrix gradient_p;
  Matrix D_Displacement_Ap;
  Matrix D_Velocity_Ap;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor DF_p;
  Tensor dFdt_n_p;
  Tensor dFdt_n1_p;
  Tensor dt_DF_p;

  /*
    Loop in the material point set
  */
  for (int p = 0; p < Np; p++) {
    /*
      Define tributary nodes of the particle
    */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],
                                     MPM_Mesh.NumberNodes[p]);

    /*
      Get the nodal increment of displacement using the mask
    */
    D_Displacement_Ap =
        get_set_field__MeshTools__(D_U.value, Nodes_p, ActiveNodes);
    D_Velocity_Ap =
        get_set_field__MeshTools__(D_U.d_value_dt, Nodes_p, ActiveNodes);

    /*
      Evaluate the shape function gradient in the coordinates of the particle
    */
    gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Take the values of the deformation gradient from the previous step
    */
    F_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p], 2);
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p], 2);
    DF_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.DF.nM[p], 2);
    dFdt_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n.nM[p], 2);
    dFdt_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n1.nM[p], 2);
    dt_DF_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_DF.nM[p], 2);

    /*
      Compute the increment of the deformation gradient
    */
    update_increment_Deformation_Gradient__Particles__(DF_p, D_Displacement_Ap,
                                                       gradient_p);
    update_rate_increment_Deformation_Gradient__Particles__(
        dt_DF_p, D_Velocity_Ap, gradient_p);

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
      fprintf(stderr, "%s %i\n", "Negative jacobian in particle", p);
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
             Free memory
    */
    free__MatrixLib__(D_Displacement_Ap);
    free__MatrixLib__(D_Velocity_Ap);
    free__MatrixLib__(gradient_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Loop in the material point set to update stress
  */
  for (int p = 0; p < Np; p++) {

    if (FEM_Mesh.Locking_Control_Fbar) {
      Idx_Element_p = MPM_Mesh.Element_p[p];
      Idx_Patch_p = FEM_Mesh.Idx_Patch[Idx_Element_p];

      Vn_patch = FEM_Mesh.Vol_Patch_n[Idx_Patch_p];
      Vn1_patch = FEM_Mesh.Vol_Patch_n1[Idx_Patch_p];
      J_patch = Vn1_patch / Vn_patch;

      get_locking_free_Deformation_Gradient_n1__Particles__(p, J_patch,MPM_Mesh);

      MPM_Mesh.Phi.Jbar.nV[p] *= J_patch;
    }
  }


  return EXIT_SUCCESS;
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
    Add nominal traction contribution
  */
  compute_Nodal_Nominal_traction_Forces(Forces, ActiveNodes, MPM_Mesh, FEM_Mesh,
                                        TimeStep, NumTimeStep);

  /*
    Add body forces contribution
  */
  compute_Nodal_Body_Forces(Forces, ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep,
                            NumTimeStep);

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
  int MatIndx_p;

  Tensor P_p; /* First Piola-Kirchhoff Stress tensor */
  Tensor InternalForcesDensity_Ap;

  Element Nodes_p;   /* List of nodes for particle */
  Material MatProp_p;  
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pA;
  Tensor GRADIENT_pA;
  Tensor F_n_p;
  Tensor transpose_F_n_p;
  double V0_p; /* Volume of the particle in the reference configuration */

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
      Take the value of the deformation gradient at t = n
    */
    F_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p], 2);
    transpose_F_n_p = transpose__TensorLib__(F_n_p);

    /*
      Update the first Piola-Kirchhoff stress tensor with an apropiate
      integration rule.
    */
    MatIndx_p = MPM_Mesh.MatIdx[p];
    MatProp_p = MPM_Mesh.Mat[MatIndx_p];   
    Stress_integration__Particles__(p, MPM_Mesh, FEM_Mesh, MatProp_p);

    /*
      Get the first Piola-Kirchhoff stress tensor
    */
    P_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);

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
        Forces.nM[A_mask][i] += InternalForcesDensity_Ap.n[i] * V0_p;
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
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_p;      /* Nodal values of the sahpe function */
  double N_pa;
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
      N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

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
        N_pa = N_p.nV[A];

        /*
          Node for the contribution
        */
        Ap = Nodes_p.Connectivity[A];
        A_mask = ActiveNodes.Nodes2Mask[Ap];

        /*
          Compute Contact forces
        */
        for (int k = 0; k < Ndim; k++) {
          Forces.nM[A_mask][k] -= N_pa * T.n[k] * A0_p;
        }
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(N_p);
      free(Nodes_p.Connectivity);
    }
  }

  free__TensorLib__(T);
}

/**************************************************************/

static void compute_Nodal_Body_Forces(Matrix Forces, Mask ActiveNodes,
                                      Particle MPM_Mesh, Mesh FEM_Mesh,
                                      int TimeStep, int NumTimeStep) {
  /* Define auxilar variables */
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int NumBodyForces = MPM_Mesh.NumberBodyForces;
  int NumParticles_i; /* Number of particles with the */
  int NumNodes_p;     /* Number of tributary nodes of p */
  int A_mask;         /* Index of the node where we apply the body force */
  int idx_A_mask_k;   /* Index of the node where we apply the body force */
  int Ap;             /* Tributary node A of particle p */
                      //  int p; /* Particle index */

  double m_p;              /* Mass of the particle */
  Load *B = MPM_Mesh.B;    /* List with the load cases */
  Element Nodes_p;         /* Element for each particle */
  Matrix ShapeFunction_p;  /* Nodal values of the sahpe function */
  double ShapeFunction_pA; /* Evaluation in the node I for the particle p */

  Tensor b = alloc__TensorLib__(1); /* Body forces vector */

  b.n[0] = 0.0;
  b.n[1] = -9.81;

  for (int p = 0; p < MPM_Mesh.NumGP; p++) {

    /* Get the value of the mass */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* Define tributary nodes of the particle */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /* Compute shape functions */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Get the node of the mesh for the contribution */
    for (int A = 0; A < NumNodes_p; A++) {

      /* Pass the value of the nodal shape function to a scalar */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /* Get the node of the mesh for the contribution */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /* Compute body forces */
      for (int k = 0; k < Ndim; k++) {
        Forces.nM[A_mask][k] -= ShapeFunction_pA * b.n[k] * m_p;
      }
    }

    /* Free the matrix with the nodal gradient of the element */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Free auxiliar tensor
  */
  free__TensorLib__(b);
}

/**********************************************************************/

static Matrix compute_Nodal_Reactions(Mesh FEM_Mesh, Matrix Forces,
                                      Mask ActiveNodes, int TimeStep,
                                      int NumTimeStep)
/*
  Compute the nodal reactions
*/
{
  /* 1º Define auxilar variables */
  int Ndim = NumberDimensions;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound;   /* Number of dimensions */
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

  Matrix Reactions = allocZ__MatrixLib__(Nnodes_mask, Ndim);

  /*
    Loop over the the boundaries
  */
  for (int i = 0; i < FEM_Mesh.Bounds.NumBounds; i++) {
    /*
       Get the number of nodes of this boundarie
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
        The boundary condition is not affecting any active node,
        continue interating
      */
      if (Id_BCC_mask == -1) {
        continue;
      }

      /*
         Loop over the dimensions of the boundary condition
      */
      for (int k = 0; k < NumDimBound; k++) {

        /*
           Apply only if the direction is active (1)
        */
        if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + TimeStep] == 1) {
          /*
             Set to zero the forces in the nodes where velocity is fixed
          */
          Reactions.nM[Id_BCC_mask][k] = Forces.nM[Id_BCC_mask][k];
          Forces.nM[Id_BCC_mask][k] = 0;
        }
      }
    }
  }

  return Reactions;
}

/**************************************************************/

static Matrix compute_Nodal_Residual(
  Nodal_Field U_n, 
  Nodal_Field D_U,
  Mask ActiveNodes, 
  Matrix Forces,
  Matrix Mass, 
  Newmark_parameters Params) {
  
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Ndim * Nnodes_mask;
  Matrix Acceleration_n1 = allocZ__MatrixLib__(Nnodes_mask, Ndim);
  Matrix Inertial_Forces = allocZ__MatrixLib__(Nnodes_mask, Ndim);
  Matrix Residual = allocZ__MatrixLib__(Nnodes_mask, Ndim);
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;

  /*
    Compute nodal acceleration (Vectorized)
  */
  for (int idx_B = 0; idx_B < Order; idx_B++) {
    Acceleration_n1.nV[idx_B] = alpha_1 * D_U.value.nV[idx_B] -
                                alpha_2 * U_n.d_value_dt.nV[idx_B] -
                                alpha_3 * U_n.d2_value_dt2.nV[idx_B];
  }

  /*
    Compute inertial forces (Vectorized)
  */
  for (int idx_A = 0; idx_A < Order; idx_A++) {
    for (int idx_B = 0; idx_B < Order; idx_B++) {
      Inertial_Forces.nV[idx_A] +=
          Mass.nM[idx_A][idx_B] * Acceleration_n1.nV[idx_B];
    }
  }

  /*
    Compute (-) residual (Vectorized). The minus symbol is due to
    solver purposes. See compute_D_Displacement
  */
  for (int idx_A = 0; idx_A < Order; idx_A++) {
    Residual.nV[idx_A] = Inertial_Forces.nV[idx_A] + Forces.nV[idx_A];
  }

  /*
    Free Memory
  */
  free__MatrixLib__(Acceleration_n1);
  free__MatrixLib__(Inertial_Forces);

  return Residual;
}

/**************************************************************/
static double __error_residual(double * Residual, int Total_dof) {
  
  double Error = 0;
  
  for (unsigned A = 0; A < Total_dof; A++) {
      Error += DSQR(Residual[A]);
  }
  Error = pow(Error, 0.5);

  return Error;
}

/**************************************************************/

static int __assemble_tangent_stiffness(
  double * Tangent_Stiffness,
  const double * Effective_Mass,
  Mask ActiveNodes,
  Particle MPM_Mesh, 
  Mesh FEM_Mesh,
  Newmark_parameters Params)
{
  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Ndof = NumberDOF;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned Order = Ndof * Nnodes_mask;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned Ap;
  unsigned A_mask;
  unsigned Bp;
  unsigned B_mask;
  unsigned NumNodes_p;
  unsigned MatIndx_p;
  unsigned idx_AB_mask_ij;

  Element Nodes_p;   /* List of nodes for particle */
  Matrix d_shapefunction_p; /* Shape functions gradients */
  double * d_shapefunction_pA;
  double * d_shapefunction_pB;    

#if NumberDimensions == 2
  double Stiffness_density_p[4];
#else
  double Stiffness_density_p[9];
#endif

  Material MatProp_p;
  State_Parameters IO_State_p;
  double m_p; /* Mass of the particle */
  double V0_p; /* Volume of the particle in the reference configuration */
  double J_p;  /* Jacobian of the deformation gradient */
  double alpha_1 = Params.alpha_1;
  double alpha_4 = Params.alpha_4;

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
    d_shapefunction_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);


    for (unsigned A = 0; A < NumNodes_p; A++) {
      
      // Get the gradient evaluation in node A
      // and the masked index of the node A
      d_shapefunction_pA = d_shapefunction_p.nM[A];
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      for (unsigned B = 0; B < NumNodes_p; B++) {

        // Get the gradient evaluation in node B
        // and the masked index of the node B
        d_shapefunction_pB = d_shapefunction_p.nM[B];
        Bp = Nodes_p.Connectivity[B];
        B_mask = ActiveNodes.Nodes2Mask[Bp];

        // Compute particle evaluation of the stiffness matrix for each node
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
        else {
          fprintf(stderr, "%s : %s %s %s \n",
          "Error in __assemble_tangent_stiffness()", "The material",
          MatProp_p.Type, "has not been yet implemnented");
          return EXIT_FAILURE;
        }

        //  Assembling process
        for (unsigned i = 0; i < Ndim; i++) {
          for (unsigned j = 0; j < Ndim; j++) {
            Tangent_Stiffness[(A_mask * Ndim + i)*Order + (B_mask * Ndim + j)] += 
            Stiffness_density_p[i*Ndim + j] * V0_p;
          }
        }

      }
    }

    free__MatrixLib__(d_shapefunction_p);
    free(Nodes_p.Connectivity);
  }

  for (unsigned A = 0; A < Order*Order; A++)
  {
    Tangent_Stiffness[A] += alpha_1 * Effective_Mass[A];
  }

  return EXIT_SUCCESS;
}

/**************************************************************/

static void system_reduction(
  double * Tangent_Stiffness,
  Mask ActiveNodes, 
  Mesh FEM_Mesh, 
  int TimeStep,
  int NumTimeStep) {

  unsigned Ndim = NumberDimensions;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned Order = Nnodes_mask * Ndim;
  unsigned Number_of_BCC = FEM_Mesh.Bounds.NumBounds;
  unsigned NumNodesBound; /* Number of nodes of the bound */
  unsigned NumDimBound;   /* Number of dimensions */
  unsigned Id_BCC;        /* Index of the node where we apply the BCC */
  unsigned Id_BCC_mask;
  unsigned Id_BCC_mask_k;

  /*
    Loop over the the boundaries to find the constrained dofs
  */
  for (unsigned i = 0; i < Number_of_BCC; i++) {

    /*
      Get the number of nodes of this boundary
    */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

    /*
      Get the number of dimensions where the BCC it is applied
    */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;

    for (unsigned j = 0; j < NumNodesBound; j++) {
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
        for (unsigned k = 0; k < NumDimBound; k++) {
          if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + TimeStep] == 1) {

            for (unsigned A_mask = 0; A_mask < Nnodes_mask; A_mask++) {
              for (unsigned l = 0; l < Ndim; l++) {
                Tangent_Stiffness[(A_mask * Ndim + l)*Order + (Id_BCC_mask * Ndim + k)] = 0.0;
                Tangent_Stiffness[(Id_BCC_mask * Ndim + k)*Order +(A_mask * Ndim + l)] = 0.0;
              }
            }

            Tangent_Stiffness[(Id_BCC_mask * Ndim + k)*Order + (Id_BCC_mask * Ndim + k)] = 1.0;
          }
        }
      }
    }
  }
}

/**************************************************************/

static int __solve_equilibrium(
  Nodal_Field D_U, 
  double * Tangent_Stiffness,
  double * Residual,
  int Nnodes_mask)
{
  int Ndof = NumberDOF;
  int Order = Nnodes_mask * Ndof;
  int LDA = Nnodes_mask * Ndof;
  int LDB = Nnodes_mask * Ndof;
  char TRANS = 'T'; /* (Transpose) */
  int INFO = 3;
  int *IPIV = (int *)Allocate_Array(Order, sizeof(int));
  int NRHS = 1;

  //  Compute the LU factorization
  dgetrf_(&Order, &Order, Tangent_Stiffness, &LDA, IPIV, &INFO);
  if (INFO) {
    free(IPIV);
    fprintf(stderr, "%s : %s %s %s \n", "Error in dgetrf_", "The function",
            "dgetrf_", "returned an error message !!!");
    return EXIT_FAILURE;
  }

  /*
    Solve the system
  */
  dgetrs_(&TRANS, &Order, &NRHS, Tangent_Stiffness, &LDA, IPIV, Residual, &LDB, &INFO);
  if (INFO) {
    free(IPIV);
    fprintf(stderr, "%s : %s %s %s \n", "Error in dgetrs_", "The function",
            "dgetrs_", "returned an error message !!!");
    exit(EXIT_FAILURE);
  }

  /*
    Update
  */
  for (int idx_A_i = 0; idx_A_i < Order; idx_A_i++) {
    D_U.value.nV[idx_A_i] -= Residual[idx_A_i];
  }

  // Free memory
  free(IPIV);
}
/**************************************************************/

static void __update_Nodal_Increments(
  Nodal_Field D_U,
  Nodal_Field U_n,
  Mask Free_and_Restricted_Dofs,
  Newmark_parameters Params) {

  unsigned Nnodes = U_n.value.N_rows;
  unsigned Ndim = NumberDimensions;
  unsigned Ndof = NumberDOF;
  unsigned Total_dof = Nnodes * NumberDOF;
  unsigned Mask_idx_A_i;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

  /*
    Update nodal variables
  */
  for (unsigned A = 0; A < Nnodes; A++) {
    for (unsigned i = 0; i < Ndim; i++) {

      Mask_idx_A_i = Free_and_Restricted_Dofs.Nodes2Mask[A * Ndim + i];

      if (Mask_idx_A_i == -1) {
        D_U.d2_value_dt2.nM[A][i] = 0.0;
        D_U.d_value_dt.nM[A][i] = alpha_4 * D_U.value.nM[A][i] +
                                  (alpha_5 - 1) * U_n.d_value_dt.nM[A][i] +
                                  alpha_6 * U_n.d2_value_dt2.nM[A][i];
      } else {
        D_U.d2_value_dt2.nM[A][i] = alpha_1 * D_U.value.nM[A][i] -
                                    alpha_2 * U_n.d_value_dt.nM[A][i] -
                                    (alpha_3 + 1) * U_n.d2_value_dt2.nM[A][i];
        D_U.d_value_dt.nM[A][i] = alpha_4 * D_U.value.nM[A][i] +
                                  (alpha_5 - 1) * U_n.d_value_dt.nM[A][i] +
                                  alpha_6 * U_n.d2_value_dt2.nM[A][i];
      }
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
  unsigned Ap;
  unsigned A_mask;
  unsigned idx_A_mask_i;
  unsigned idx_ij;

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
        D_U_pI = ShapeFunction_pI * D_U.value.nM[A_mask][i];
        D_V_pI = ShapeFunction_pI * D_U.d_value_dt.nM[A_mask][i];
        D_A_pI = ShapeFunction_pI * D_U.d2_value_dt2.nM[A_mask][i];

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

static void output_selector(Particle MPM_Mesh, Mesh FEM_Mesh, Mask ActiveNodes,
                            Nodal_Field D_U, Nodal_Field U_n, Matrix Forces,
                            Matrix Reactions, Matrix Residual, int TimeStep,
                            int ResultsTimeStep) {

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
