#include "Formulations/Displacements-WaterPressure/U-pw-Verlet.h"

/*
  Auxiliar functions and variables
*/

static char Error_message[MAXW];
static void standard_error();

/* Step 1 */
static Matrix compute_Mass_Matrix_Mixture(Particle, Mesh, Mask);
static Matrix compute_Compressibility_Matrix_Fluid(Particle, Mesh, Mask);
/* Step 2 */
static void compute_Explicit_Newmark_Predictor(Particle, double, double);
/* Step 3 */
static Matrix compute_Nodal_D_Displacement(Particle, Mesh, Mask, Matrix);
static Matrix compute_Nodal_Velocity(Particle, Mesh, Mask, Matrix);
static Matrix compute_Nodal_Pore_water_pressure(Particle, Mesh, Mask, Matrix);
static void impose_Dirichlet_Boundary_Conditions(Mesh, Matrix, Matrix, Mask,
                                                 int, int);
/* Step 4 */
static void update_Local_State(Matrix, Matrix, Mask, Particle, Mesh);
/* Step 5 */
static Matrix compute_Total_Forces_Mixture(Mask, Particle, Mesh, int, int);
static void compute_Internal_Forces_Mixture(Matrix, Mask, Particle, Mesh);
static void compute_Contact_Forces_Mixture(Matrix, Mask, Particle, Mesh, int,
                                           int);
static Tensor compute_total_first_Piola_Kirchhoff_stress(Tensor, double,
                                                         Tensor);
static Matrix solve_Nodal_Equilibrium_Mixture(Matrix, Matrix, Matrix, Particle,
                                              Mesh, Mask, Mask);
/* Step 6 */
static Matrix compute_Mass_exchanges_Source_Terms(Matrix, Mask, Particle, Mesh,
                                                  int);
static void compute_Jacobian_Rate_Mass_Balance(Matrix, Mask, Particle, Mesh);
static void compute_Permeability_Mass_Balance(Matrix, Mask, Particle, Mesh,
                                              Matrix);
static Tensor compute_Pore_water_pressure_gradient(Matrix, Matrix);
static void compute_Permeability_Inertial_Forces_Fluid(Matrix, Mask, Particle,
                                                       Mesh);
// static  void  compute_Fluid_Mass_Pump(Matrix,Mesh,Matrix,Mask,int);
static void solve_Nodal_Mass_Balance(Matrix, Matrix, Particle, Mesh, Mask,
                                     Mask);
/* Step 7 */
static void compute_Explicit_Newmark_Corrector(Particle, double, double);
/* Step 8 */
static void output_selector(Particle, Mesh, Mask, Matrix, Matrix, Matrix,
                            Matrix, double, int, int);

/**************************************************************/

void upw_Verlet(
    Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver) {

  /*
    Auxiliar variables for the solver
  */
  int Ndim = NumberDimensions;
  int Nactivenodes;
  int InitialStep = Parameters_Solver.InitialTimeStep;
  int NumTimeStep = Parameters_Solver.NumTimeStep;

  double gamma = 0.5;
  double CFL = Parameters_Solver.CFL;
  double DeltaTimeStep;
  double DeltaX = FEM_Mesh.DeltaX;

  Matrix Mass_Matrix_Mixture;
  Matrix Compressibility_Matrix_Fluid;
  Matrix Velocity;
  Matrix Pore_water_pressure;
  Matrix Rate_Pore_water_pressure;
  Matrix Gravity_field;
  Matrix D_Displacement;
  Matrix Total_Forces_Mixture;
  Matrix Reactions_Mixture;
  Matrix Mass_Exchanges_Source_Terms;
  Matrix Reactions_Fluid;

  Mask ActiveNodes;
  Mask Free_and_Restricted_Dofs;

  for (unsigned TimeStep = InitialStep; TimeStep < NumTimeStep; TimeStep++) {

    DeltaTimeStep = DeltaT_Coussy__SolversLib__(MPM_Mesh, DeltaX, 1.0, CFL);
    local_search__MeshTools__(MPM_Mesh, FEM_Mesh);
    ActiveNodes = get_active_nodes__MeshTools__(FEM_Mesh);
    Nactivenodes = ActiveNodes.Nactivenodes;
    Free_and_Restricted_Dofs =
        get_active_dofs__MeshTools__(
            ActiveNodes, FEM_Mesh, TimeStep, NumTimeStep);
    print_step(TimeStep, NumTimeStep, DeltaTimeStep);

    Mass_Matrix_Mixture =
        compute_Mass_Matrix_Mixture(MPM_Mesh, FEM_Mesh, ActiveNodes);

    Compressibility_Matrix_Fluid =
        compute_Compressibility_Matrix_Fluid(MPM_Mesh, FEM_Mesh, ActiveNodes);

    compute_Explicit_Newmark_Predictor(MPM_Mesh, gamma, DeltaTimeStep);

    D_Displacement = compute_Nodal_D_Displacement(
        MPM_Mesh, FEM_Mesh, ActiveNodes, Mass_Matrix_Mixture);

    Velocity = compute_Nodal_Velocity(MPM_Mesh, FEM_Mesh, ActiveNodes,
                                      Mass_Matrix_Mixture);

    Pore_water_pressure = compute_Nodal_Pore_water_pressure(
        MPM_Mesh, FEM_Mesh, ActiveNodes, Compressibility_Matrix_Fluid);

    impose_Dirichlet_Boundary_Conditions(FEM_Mesh, Velocity,
                                         Pore_water_pressure, ActiveNodes,
                                         TimeStep, NumTimeStep);

    update_Local_State(D_Displacement, Velocity, ActiveNodes, MPM_Mesh,
                       FEM_Mesh);

    Total_Forces_Mixture = compute_Total_Forces_Mixture(
        ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep, NumTimeStep);

    Reactions_Mixture = solve_Nodal_Equilibrium_Mixture(
        Mass_Matrix_Mixture, Gravity_field, Total_Forces_Mixture, MPM_Mesh,
        FEM_Mesh, ActiveNodes, Free_and_Restricted_Dofs);

    Mass_Exchanges_Source_Terms = compute_Mass_exchanges_Source_Terms(
        Pore_water_pressure, ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep);

    solve_Nodal_Mass_Balance(Compressibility_Matrix_Fluid,
                             Mass_Exchanges_Source_Terms, MPM_Mesh, FEM_Mesh,
                             ActiveNodes, Free_and_Restricted_Dofs);

    compute_Explicit_Newmark_Corrector(MPM_Mesh, gamma, DeltaTimeStep);

    output_selector(MPM_Mesh, FEM_Mesh, ActiveNodes, Velocity, D_Displacement,
                    Total_Forces_Mixture, Reactions_Mixture, DeltaTimeStep,
                    TimeStep, ResultsTimeStep);

    free__MatrixLib__(Mass_Matrix_Mixture);
    free__MatrixLib__(Compressibility_Matrix_Fluid);
    free__MatrixLib__(Gravity_field);
    free__MatrixLib__(Velocity);
    free__MatrixLib__(D_Displacement);
    free__MatrixLib__(Total_Forces_Mixture);
    free__MatrixLib__(Reactions_Mixture);
    free(ActiveNodes.Nodes2Mask);

    print_Status("DONE !!!", TimeStep);
  }
}

/**************************************************************/

static Matrix compute_Mass_Matrix_Mixture(
    Particle MPM_Mesh, // Variable with information of the particles
    Mesh FEM_Mesh,     // Variable with information of the nodes
    Mask ActiveNodes)  // Variable with information of the active nodes
/*
  This function computes the lumped mass matrix with the displacements degree of
  freedom
*/
{

  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Order = Ndim * Nnodes_mask;
  int Ap;
  int A_mask;
  int idx_A_mask_i;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;

  /* Evaluation of the particle in the node */
  double ShapeFunction_pA;
  /* Mass of the particle */
  double m_p;
  /* Nodal contribution A of the particle p */
  double m_A_p;
  /* Element for each particle */
  Element Nodes_p;

  /* Define and allocate the lumped mass matrix */
  Matrix Lumped_MassMatrix = allocZ__MatrixLib__(Order, Order);

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

      /*
        Get the value of the shape function
      */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /*
        Compute the nodal A contribution of the particle p
      */
      m_A_p = m_p * ShapeFunction_pA;

      /*
         Fill the Lumped mass matrix considering the number of dofs
      */
      for (int i = 0; i < Ndim; i++) {
        Lumped_MassMatrix.nM[A_mask * Ndim + i][A_mask * Ndim + i] += m_A_p;
      }
    }

    /*
      Free the value of the shape functions
    */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }


  return Lumped_MassMatrix;
}

/**************************************************************/

static Matrix compute_Compressibility_Matrix_Fluid(
    Particle MPM_Mesh, // Variable with information of the particles
    Mesh FEM_Mesh,     // Variable with information of the nodes
    Mask ActiveNodes)  // Variable with information of the active nodes
/*
  Compute the lumped compresibility matrix
*/
{
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Order = 1 * Nnodes_mask;
  int Mixture_idx;
  int Material_Water_idx;
  int Ap;
  int A_mask;
  Element Nodes_p;          /* Element for each particle */
  Matrix ShapeFunction_p;   /* Value of the shape-function */
  Material MatProp_Water_p; /* Variable with the material properties of the
                               fluid phase */
  double ShapeFunction_pA;  /* Evaluation of the particle in the node A */
  double V0_p;              /* Mass of the particle (mixture) */
  double rho_f_p;           /* Material density of the fluid */
  double relative_rho_f_p;  /* Relative density of the fluid */
  double phi_f_p;           /* Volume fractions of the fluid */
  double K_f;               /* Compressibility (fluid) */
  double compressibility_density_f_p; /* Compressibilidy density fo the fluid */
  double
      compressibility_density_A_p; /* Nodal contribution A of the particle p */
  /*
    Define and allocate the lumped Compressibility matrix
  */
  Matrix Lumped_Compressibility_Matrix_Fluid =
      allocZ__MatrixLib__(Order, Order);

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
      Get the initial volume of the particle
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /*
      Get the current material density, volume fraction
      and compressibility for each material point (fluid)
    */
    rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];
    phi_f_p = MPM_Mesh.Phi.phi_f.nV[p];

    /*
      Load intrinsic properties for the fluid phase to get the compressibility
    */
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Water_idx = Soil_Water_Mixtures[Mixture_idx].Water_Idx;
    MatProp_Water_p = MPM_Mesh.Mat[Material_Water_idx];
    K_f = MatProp_Water_p.Compressibility;

    /*
      Load intrinsic properties for the fluid phase to get the compressibility
    */
    K_f = MatProp_Water_p.Compressibility;

    /*
      Compute relative density
    */
    relative_rho_f_p = phi_f_p * rho_f_p;

    /*
      Compute the compressibility density
    */
    compressibility_density_f_p = (relative_rho_f_p / K_f) * V0_p;

    for (int A = 0; A < Nodes_p.NumberNodes; A++) {
      /*
        Get the node in the mass matrix with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Get the value of the shape function
       */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /*
         Compute the nodal A contribution of the particle p
       */
      compressibility_density_A_p =
          compressibility_density_f_p * ShapeFunction_pA;

      /*
        Fill the Lumped Compressibility matrix
      */
      Lumped_Compressibility_Matrix_Fluid.nM[A_mask][A_mask] +=
          compressibility_density_A_p;
    }

    /* Free the value of the shape functions */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }


  return Lumped_Compressibility_Matrix_Fluid;
}

/**************************************************************/

static void compute_Explicit_Newmark_Predictor(
    Particle MPM_Mesh, // Information related with particles
    double gamma,      // Newmark integration parameter
    double Dt)         // Time step
/*
  The predictor stage is computed in the particles
*/
{

  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;

  for (int p = 0; p < Np; p++) {

    /*
      Compute pore water pressure predictor
    */
    MPM_Mesh.Phi.Pw.nV[p] += (1 - gamma) * Dt * MPM_Mesh.Phi.d_Pw_dt_n1.nV[p];

    /*
      Compute velocity predictor and increment of displacements
    */
    for (int i = 0; i < Ndim; i++) {

      MPM_Mesh.Phi.D_dis.nM[p][i] = Dt * MPM_Mesh.Phi.vel.nM[p][i] +
                                    0.5 * DSQR(Dt) * MPM_Mesh.Phi.acc.nM[p][i];

      MPM_Mesh.Phi.dis.nM[p][i] += MPM_Mesh.Phi.D_dis.nM[p][i];

      MPM_Mesh.Phi.vel.nM[p][i] += (1 - gamma) * Dt * MPM_Mesh.Phi.acc.nM[p][i];
    }
  }
}


/**************************************************************/

static Matrix compute_Nodal_D_Displacement(Particle MPM_Mesh, Mesh FEM_Mesh,
                                           Mask ActiveNodes,
                                           Matrix Mass_Matrix_Mixture)
/*
  Compute the nodal increment of displacement. The operation is linearized and
  all the dof split the increment of displacement array in n components like :
  | M 0 |   |D_u.x|   | M*D_u.x |
  | 0 M | * |D_u.y| = | M*D_u.y |
*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Order = Nnodes_mask * Ndim;
  int Ap;
  int A_mask;
  int AB;
  Element Nodes_p;         /* Element for each particle */
  Matrix ShapeFunction_p;  /* Value of the shape-function */
  double ShapeFunction_pA; /* Evaluation of the particle in the node */
  double m_p;              /* Mass of the particle */

  /*
    Define and allocate the nodal increment of displacement vector
  */
  Matrix D_Displacement = allocZ__MatrixLib__(Nnodes_mask, Ndim);

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

    for (int A = 0; A < Nodes_p.NumberNodes; A++) {

      /*
        Get the node in the nodal momentum with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
         Evaluate the GP function in the node
      */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /*
        Nodal velocity and D_Displacement
      */
      for (int i = 0; i < Ndim; i++) {
        D_Displacement.nM[A_mask][i] +=
            m_p * ShapeFunction_pA * MPM_Mesh.Phi.D_dis.nM[p][i];
      }
    }

    /*
      Free the value of the shape functions
    */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Compute the D_Displacements
  */
  for (int A = 0; A < Nnodes_mask; A++) {
    for (int i = 0; i < Ndim; i++) {
      AB = A * Ndim + i;
      D_Displacement.nM[A][i] =
          D_Displacement.nM[A][i] / Mass_Matrix_Mixture.nM[AB][AB];
    }
  }


  return D_Displacement;
}

/**************************************************************/

static Matrix compute_Nodal_Velocity(Particle MPM_Mesh, Mesh FEM_Mesh,
                                     Mask ActiveNodes,
                                     Matrix Mass_Matrix_Mixture)
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
      Velocity.nM[A][i] = Velocity.nM[A][i] / Mass_Matrix_Mixture.nM[AB][AB];
    }
  }

  return Velocity;
}

/**************************************************************/

static Matrix
compute_Nodal_Pore_water_pressure(Particle MPM_Mesh, Mesh FEM_Mesh,
                                  Mask ActiveNodes,
                                  Matrix Compressibility_Matrix_Fluid)
/*

*/
{
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Mixture_idx;
  int Material_Water_idx;
  int Ap;
  int A_mask;
  Material MatProp_Water_p; /* Variable with the material properties of the
                               fluid phase */
  Element Nodes_p;          /* Element for each particle */
  Matrix ShapeFunction_p;   /* Value of the shape-function */
  double ShapeFunction_pA;  /* Evaluation of the particle in the node */
  double Pw_p;              /* Particle pore water pressure */
  double V0_p;              /* Volume of the particle */
  double rho_f_p;           /* Material density of the particle (fluid phase) */
  double phi_f_p;           /* Volume fraction of the particle (fluid phase) */
  double relative_rho_f_p;  /* Relative density of the particle (fluid phase) */
  double K_f;               /* Compressibility (fluid) */
  double compressibility_density_f_p; /* Compressibily density for the particle
                                         (fluid phase) */
  double compressibility_density_A_p; /* Contribution of the particle p to the
                                         compressibility in the node A */

  /*
    Define and allocate the Pore water pressure vector
  */
  Matrix Pore_water_pressure = allocZ__MatrixLib__(Nnodes_mask, 1);

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
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
       Get the pore water pressure
    */
    Pw_p = MPM_Mesh.Phi.Pw.nV[p];

    /*
      Get the initial volume of the particle
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /*
      Get the current material density, volume fraction
      and compressibility for each material point (fluid)
    */
    rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];
    phi_f_p = MPM_Mesh.Phi.phi_f.nV[p];

    /*
      Get intrinsic material properties for the fluid phase (Compressibility)
    */
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Water_idx = Soil_Water_Mixtures[Mixture_idx].Water_Idx;
    MatProp_Water_p = MPM_Mesh.Mat[Material_Water_idx];
    K_f = MatProp_Water_p.Compressibility;

    /*
      Compute relative density
    */
    relative_rho_f_p = phi_f_p * rho_f_p;

    /*
      Compute the compressibility density
    */
    compressibility_density_f_p = (relative_rho_f_p / K_f) * V0_p;

    /* Get the nodal mommentum */
    for (int A = 0; A < Nodes_p.NumberNodes; A++) {
      /*
        Get the node in the nodal momentum with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /* Evaluate the GP function in the node */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /*
        Compute the nodal A contribution of the particle p
      */
      compressibility_density_A_p =
          compressibility_density_f_p * ShapeFunction_pA;

      /* Nodal Pore water pressure */
      Pore_water_pressure.nV[A_mask] += compressibility_density_A_p * Pw_p;
    }

    /* Free the value of the shape functions */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Compute the nodal velocities
  */
  for (int A = 0; A < Nnodes_mask; A++) {
    Pore_water_pressure.nV[A] =
        Pore_water_pressure.nV[A] / Compressibility_Matrix_Fluid.nM[A][A];
  }


  return Pore_water_pressure;
}

/**************************************************************/

static void impose_Dirichlet_Boundary_Conditions(Mesh FEM_Mesh, Matrix Velocity,
                                                 Matrix Pore_water_pressure,
                                                 Mask ActiveNodes, int TimeStep,
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

  /* 2º Loop over the the boundaries */
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
            Assign the boundary condition
          */
          if (k < Ndim) {
            Velocity.nM[Id_BCC_mask][k] =
                FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep];
          } else {
            Pore_water_pressure.nV[Id_BCC_mask] =
                FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep];
          }
        }
      }
    }
  }
}

/**************************************************************/

static void update_Local_State(
    Matrix D_Displacement, // Nodal values of the increment of displacement
    Matrix Velocity,       // Nodal values of the velocity
    Mask ActiveNodes,      // Information with the active nodes
    Particle MPM_Mesh,     // Information related with the particles
    Mesh FEM_Mesh)         // Information related with the nodes
/*

*/
{

  /*
    Auxiliar variables
  */
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Mixture_idx;
  int Material_Soil_idx;
  int Material_Water_idx;
  int NumNodes_p;
  double J_n1_p;     /* Jacobian of the deformation gradient at t = n + 1 */
  double dJ_dt_n1_p; /* Rate of the Jacobian of the deformation gradient at t =
                        n + 1 */
  double K_f;        /* Compressibility (fluid) */
  double rho_f_0;    /* Initial density of the fluid */
  double phi_s_0;    /* Initial volume fraction (solid) */
  double phi_f_0;    /* Initial volume fraction (fluid) */
  State_Parameters
      Input_Plastic_Parameters; /* Input parameters for plasticity */
  State_Parameters
      Output_Plastic_Parameters; /* Output parameters for plasticity */
  Element Nodes_p;               /* Element for each particle */
  Material MatProp_Soil_p; /* Variable with the material properties of the solid
                              phase */
  Material MatProp_Water_p; /* Variable with the material properties of the
                               fluid phase */
  Matrix gradient_p;
  Matrix D_Displacement_Ap;
  Matrix Nodal_Velocity_p;
  Tensor grad_Nodal_Velocity_p;
  Tensor GRAD_Nodal_Velocity_p;
  Tensor F_n_p;  /* Deformation gradient of the soil skeleton (t = n) */
  Tensor F_n1_p; /* Deformation gradient of the soil skeleton (t = n + 1) */
  Tensor DF_p; /* Increment of the deformation gradient of the soil skeleton */
  Tensor FT_n1_p; /* Transpose of the deformation gradient of the soil skeleton
                     (t = n + 1) */
  Tensor F_plastic_p; /* Plastic deformation gradient of the soil skeleton */
  double Pw_0;        /* Cauchy pore water pressure at t = 0 */
  double Pw_n1;       /* Cauchy pore water pressure at t = n + 1 */

  for (int p = 0; p < Np; p++) {
    /*
      Load material properties for each phase
    */
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
    Material_Water_idx = Soil_Water_Mixtures[Mixture_idx].Water_Idx;
    MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];
    MatProp_Water_p = MPM_Mesh.Mat[Material_Water_idx];

    /*
      Read intrinsic material properties (fluid)
    */
    K_f = MatProp_Water_p.Compressibility;
    rho_f_0 = MatProp_Water_p.rho;

    /*
      Read reference volume fraction for each phase
    */
    phi_f_0 = Soil_Water_Mixtures[Mixture_idx].phi_f_0;
    phi_s_0 = Soil_Water_Mixtures[Mixture_idx].phi_s_0;

    /*
             Define tributary nodes of the particle
    */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /*
      Get the nodal increment of displacement using the mask
    */
//    D_Displacement_Ap =
//        get_set_field__MeshTools__(D_Displacement, Nodes_p, ActiveNodes);
//    Nodal_Velocity_p =
//        get_set_field__MeshTools__(Velocity, Nodes_p, ActiveNodes);

    /*
             Evaluate the shape function gradient in the coordinates of the
       particle
    */
    gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Take the values of the deformation gradient from the previous step
    */
    F_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p], 2);
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p], 2);
    DF_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.DF.nM[p], 2);
    FT_n1_p = transpose__TensorLib__(F_n1_p);

    /*
      Compute the increment of the deformation gradient to compute the
      deformation gradient
    */
//    update_increment_Deformation_Gradient__Particles__(DF_p, D_Displacement_Ap,
//                                                       gradient_p);
//    update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, DF_p);

    /*
      Compute the Jacobian of the deformation gradient
    */
//    J_n1_p = I3__TensorLib__(F_n1_p);

    /*
      Compute the rate of the jacobian
    */
    grad_Nodal_Velocity_p =
        interpolate_vectorial_magnitude_gradient__MeshTools__(Nodal_Velocity_p,
                                                              gradient_p);
    GRAD_Nodal_Velocity_p =
        vector_linear_mapping__TensorLib__(FT_n1_p, grad_Nodal_Velocity_p);
    dJ_dt_n1_p = I1__TensorLib__(GRAD_Nodal_Velocity_p);

    /*
      Update state parameters
    */
    Pw_0 = MPM_Mesh.Phi.Pw_0.nV[p]; /* Get the initial pressure */
    Pw_n1 =
        MPM_Mesh.Phi.Pw.nV[p] / J_n1_p; /* From the Kirchhoff pressure compure
                                           the cauchy pore water pressure */

    MPM_Mesh.Phi.rho_f.nV[p] =
        rho_f_0 * exp((Pw_n1 - Pw_0) / K_f); /* Update the fluid density */

    MPM_Mesh.Phi.phi_s.nV[p] =
        phi_s_0 / J_n1_p; /* Update the volume fraction of the solid phase */
    MPM_Mesh.Phi.phi_f.nV[p] =
        1 - (1 - phi_f_0) /
                J_n1_p; /* Update the volume fraction of the fluid phase */

    MPM_Mesh.Phi.rho.nV[p] =
        MPM_Mesh.Phi.rho_s.nV[p] * MPM_Mesh.Phi.phi_s.nV[p] +
        MPM_Mesh.Phi.rho_f.nV[p] *
            MPM_Mesh.Phi.phi_f.nV[p]; /* Update density of the mixture */

    MPM_Mesh.Phi.J_n1.nV[p] = J_n1_p; /* Update soil skeleton jacobian */
    MPM_Mesh.Phi.dJ_dt.nV[p] =
        dJ_dt_n1_p; /* Update soil skeleton rate of jacobian */

    /*
             Free memory
    */
    free__MatrixLib__(D_Displacement_Ap);
    free__MatrixLib__(Nodal_Velocity_p);
    free__MatrixLib__(gradient_p);
    free__TensorLib__(FT_n1_p);
    free__TensorLib__(grad_Nodal_Velocity_p);
    free__TensorLib__(GRAD_Nodal_Velocity_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Loop in the material point set to update stress
  */
  for (int p = 0; p < Np; p++) {
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
    MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];

    /*
      Update the first Piola-Kirchhoff stress tensor with an apropiate
      integration rule.
    */
    Stress_integration__Constitutive__(p, MPM_Mesh, MatProp_Soil_p);
  }
}

/**************************************************************/

static Matrix compute_Total_Forces_Mixture(Mask ActiveNodes, Particle MPM_Mesh,
                                           Mesh FEM_Mesh, int TimeStep,
                                           int NumTimeStep) {

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  Matrix Forces = allocZ__MatrixLib__(Nnodes_mask, Ndim);

  /*
    Add internal forces contribution
  */
  compute_Internal_Forces_Mixture(Forces, ActiveNodes, MPM_Mesh, FEM_Mesh);

  /*
    Add contact forces contribution
  */
  compute_Contact_Forces_Mixture(Forces, ActiveNodes, MPM_Mesh, FEM_Mesh,
                                 TimeStep, NumTimeStep);

  return Forces;
}

/**************************************************************/

static void compute_Internal_Forces_Mixture(Matrix Forces, Mask ActiveNodes,
                                            Particle MPM_Mesh, Mesh FEM_Mesh)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int NumNodes_p;

  Tensor P_p;           /* Total First Piola-Kirchhoff Stress tensor */
  Tensor P_effective_p; /* Effective First Piola-Kirchhoff Stress tensor */
  double theta_p;       /* Kirchhoff pore fluid pressure */

  Tensor InternalForcesDensity_Ap;

  Element Nodes_p;   /* List of nodes for particle */
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pA;
  Tensor GRADIENT_pA;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor transpose_F_n_p;
  double V0_p; /* Volume of the Gauss-Point */

  /*
    Loop in the particles
  */
  for (int p = 0; p < Np; p++) {

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
             Take the values of the deformation gradient at t = n. And transpose
       it
    */
    F_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p], 2);
    transpose_F_n_p = transpose__TensorLib__(F_n_p);
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p], 2);

    /*
      Get the volume of the particle in the reference configuration
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /*
             Get the first Piola-Kirchhoff stress tensor
    */
    P_effective_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);

    /*
      Get the Kirchhoff pressure
    */
    theta_p = MPM_Mesh.Phi.Pw.nV[p];

    /*
      Following Terzaghi's idea, the effective stress tensor in the reference
      configuration is computed:
    */
    P_p = compute_total_first_Piola_Kirchhoff_stress(P_effective_p, theta_p,
                                                     F_n1_p);

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
    free__TensorLib__(P_effective_p);
    free__MatrixLib__(gradient_p);
    free(Nodes_p.Connectivity);
  }
}

/**************************************************************/

static Tensor compute_total_first_Piola_Kirchhoff_stress(Tensor P_effective_p,
                                                         double theta_p,
                                                         Tensor F_n1_p)
/*
  This function returns : P = P' - theta*F^{-T}
*/
{
  int Ndim = NumberDimensions;

  Tensor P_p = alloc__TensorLib__(2);
  Tensor inverse_F_n1_p = Inverse__TensorLib__(F_n1_p);

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      P_p.N[i][j] = P_effective_p.N[i][j] - theta_p * inverse_F_n1_p.N[j][i];
    }
  }

  free__TensorLib__(inverse_F_n1_p);

  return P_p;
}

/*********************************************************************/

static void compute_Contact_Forces_Mixture(Matrix Forces, Mask ActiveNodes,
                                           Particle MPM_Mesh, Mesh FEM_Mesh,
                                           int TimeStep, int NumTimeStep) {
  int Ndim = NumberDimensions;
  Load T_i;
  Element Nodes_p;        /* Element for each Gauss-Point */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pA;
  Tensor t = alloc__TensorLib__(1); /* Body forces vector */
  double V0_p; /* Volumen of the particle in the reference configuration */
  double thickness_p; /* Thickness of the particle */
  double A0_p;        /* Area of the particle in the reference configuration */

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
    T_i = MPM_Mesh.Neumann_Contours.BCC_i[i];

    NumNodesLoad = T_i.NumNodes;

    for (int j = 0; j < NumNodesLoad; j++) {

      /*
        Get the index of the particle
      */
      p = T_i.Nodes[j];

      /*
        Get the volume of the particle in the reference configuration
      */
      V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

      /*
        Get the thickness of each particle
      */
      thickness_p = Thickness_Plain_Stress;

      /*
        Get the area of each particle
      */
      A0_p = V0_p / thickness_p;

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
        if (T_i.Dir[k * NumTimeStep + TimeStep]) {
          t.n[k] = T_i.Value[k].Fx[TimeStep];
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
          Forces.nM[A_mask][k] += ShapeFunction_pA * t.n[k] * A0_p;
        }
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }
  }

  free__TensorLib__(t);
}

/**************************************************************/

static Matrix solve_Nodal_Equilibrium_Mixture(Matrix Mass_Matrix_Mixture,
                                              Matrix Gravity_field,
                                              Matrix Total_Forces_Mixture,
                                              Particle MPM_Mesh, Mesh FEM_Mesh,
                                              Mask ActiveNodes,
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
        Acceleration.nM[A][i] =
            Gravity_field.nM[A][i] +
            Total_Forces_Mixture.nM[A][i] / Mass_Matrix_Mixture.nM[AB][AB];
      } else {
        Acceleration.nM[A][i] = 0.0;
        Reactions.nM[A][i] = Total_Forces_Mixture.nM[A][i];
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

static Matrix compute_Mass_exchanges_Source_Terms(Matrix Pore_water_pressure,
                                                  Mask ActiveNodes,
                                                  Particle MPM_Mesh,
                                                  Mesh FEM_Mesh, int TimeStep) {

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  Matrix Mass_Exchanges_Source_Terms = allocZ__MatrixLib__(Nnodes_mask, 1);

  /*
    Compute the diferrent contributions to the mass balance
  */

  compute_Jacobian_Rate_Mass_Balance(Mass_Exchanges_Source_Terms, ActiveNodes,
                                     MPM_Mesh, FEM_Mesh);

  //  compute_Fluid_Mass_Pump(Mass_Exchanges_Source_Terms,FEM_Mesh,Forces,ActiveNodes,TimeStep);

  compute_Permeability_Mass_Balance(Mass_Exchanges_Source_Terms, ActiveNodes,
                                    MPM_Mesh, FEM_Mesh, Pore_water_pressure);

  compute_Permeability_Inertial_Forces_Fluid(Mass_Exchanges_Source_Terms,
                                             ActiveNodes, MPM_Mesh, FEM_Mesh);

  return Mass_Exchanges_Source_Terms;
}

/**************************************************************/

static void
compute_Jacobian_Rate_Mass_Balance(Matrix Mass_Exchanges_Source_Terms,
                                   Mask ActiveNodes, Particle MPM_Mesh,
                                   Mesh FEM_Mesh)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int NumNodes_p; /* Number of tributary nodes of p */
  int A_mask;     /* Index of the node where we apply the body force */
  int Ap;         /* Tributary node A of particle p */

  Element Nodes_p;         /* Element for each particle p */
  Matrix ShapeFunction_p;  /* Matrix with the value of the shape function in the
                              particle p */
  double ShapeFunction_pA; /* Value of the shape funtion in node A for the
                              particle p */
  double rho_f_p;    /* Material density of the fluid phase for particle p */
  double V0_p;       /* Volume of the particle */
  double dJ_dt_n1_p; /* Rate of the jacobian */

  for (int p = 0; p < Np; p++) {

    /*
      Get the current material density for each material point (fluid)
    */
    rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];

    /*
      Get the reference volume for each material point (mixture)
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /*
      Define nodes for each particle
    */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /*
      Evaluate the shape function in the coordinates of the particle
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Get the rate of the jacobian for each material point
    */
    dJ_dt_n1_p = MPM_Mesh.Phi.dJ_dt.nV[p];

    for (int A = 0; A < NumNodes_p; A++) {

      /*
        Compute the nodal contribution of the shape function
      */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /*
        Get the node of the mesh for the contribution
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Add the contribution of the jacobian rate to the mass conservation
      */
      Mass_Exchanges_Source_Terms.nV[A_mask] -=
          ShapeFunction_pA * rho_f_p * dJ_dt_n1_p * V0_p;
    }

    /*
 Free memory
    */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }
}

/**************************************************************/

static void
compute_Permeability_Mass_Balance(Matrix Mass_Exchanges_Source_Terms,
                                  Mask ActiveNodes, Particle MPM_Mesh,
                                  Mesh FEM_Mesh, Matrix Pore_water_pressure)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int NumNodes_p;
  int Mixture_idx;
  int Material_Soil_idx;
  int Ap;
  int A_mask;

  Material MatProp_Soil_p; /* Variable with the material propertie of solid
                              phase for each particle */
  Element Nodes_p;         /* Element for each particle */
  Matrix Nodal_Pore_water_pressure_p;
  Matrix gradient_p;      /* Shape functions gradients */
  Tensor gradient_pA;     /* Shape functions gradients (Node A), def config */
  Tensor GRADIENT_pA;     /* Shape functions gradients (Node A), ref config */
  Tensor F_n_p;           /* Deformation gradient t = n */
  Tensor F_n1_p;          /* Deformation gradient t = n + 1 */
  Tensor transpose_F_n_p; /* Transpose of the deformation gradient t = n */
  Tensor inverse_F_n1_p;  /* */
  Tensor k_p;             /* Spatial permebility tensor */
  Tensor Fk_p; /* Product of the inverse of the defomration gradient and the
                  permeability tensor */
  Tensor gradPw;
  Tensor Fk__x__gradPw;
  double GRADIENT_pA__x__Fk__x__gradPw;
  double g = -9.81;
  double V0_p; /* Volume of the particle at the reference configuration */

  /*
    Iterate over the particles to get the nodal values
  */
  for (int p = 0; p < Np; p++) {

    /*
      Define tributary nodes of the particle
    */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /*
      Compute gradient of the shape function in each node
    */
    gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Take the values of the deformation gradient at t = n and t = n + 1.
      Later compute the transpose of the deformation gradient and the
      determinant.
    */
    F_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p], 2);
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p], 2);
    transpose_F_n_p = transpose__TensorLib__(F_n_p);
    inverse_F_n1_p = Inverse__TensorLib__(F_n1_p);

    /*
      Get the reference volume for each material point (mixture)
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /*
      Load intrinsic properties for the solid phase to compute the permeability
    */
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
    MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];
    k_p = Soil_Water_Mixtures[Mixture_idx].Permeability;

    /*
      Intermediate result 1
    */
    Fk_p = matrix_product_old__TensorLib__(inverse_F_n1_p, k_p);

    /*
      Compute particle pore water pressure gradient
    */
//    Nodal_Pore_water_pressure_p =
//        get_set_field__MeshTools__(Pore_water_pressure, Nodes_p, ActiveNodes);
//    gradPw = compute_Pore_water_pressure_gradient(Nodal_Pore_water_pressure_p,
//                                                  gradient_p);

    /*
      Intermediate result 1
    */
    Fk__x__gradPw = vector_linear_mapping__TensorLib__(Fk_p, gradPw);

    for (int A = 0; A < NumNodes_p; A++) {

      /*
        Get the node in the mass matrix with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Compute the gradient in the reference configuration for the node A
      */
      gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);
      GRADIENT_pA =
          vector_linear_mapping__TensorLib__(transpose_F_n_p, gradient_pA);

      /*
        Intermediate result 2
      */
      GRADIENT_pA__x__Fk__x__gradPw =
          inner_product__TensorLib__(GRADIENT_pA, Fk__x__gradPw);

      /*
        Compute nodal contribution to the mass conservation
      */
      Mass_Exchanges_Source_Terms.nV[A_mask] +=
          (1 / g) * GRADIENT_pA__x__Fk__x__gradPw * V0_p;

      free__TensorLib__(GRADIENT_pA);
    }

    /*
      Free the value of the shape functions
    */
    free__MatrixLib__(gradient_p);
    free__MatrixLib__(Nodal_Pore_water_pressure_p);
    free__TensorLib__(transpose_F_n_p);
    free__TensorLib__(inverse_F_n1_p);
    free__TensorLib__(Fk_p);
    free__TensorLib__(gradPw);
    free__TensorLib__(Fk__x__gradPw);
    free(Nodes_p.Connectivity);
  }
}

/**************************************************************/

static Tensor
compute_Pore_water_pressure_gradient(Matrix Nodal_Pore_water_pressure_p,
                                     Matrix gradient_p)
/*

*/
{

  /* Variable definition */
  int Ndim = NumberDimensions;
  int Nnodes_p = Nodal_Pore_water_pressure_p.N_rows;
  double Pore_water_pressure_pA;
  Tensor gradient_A;
  Tensor gradPw = alloc__TensorLib__(1);

  for (int A = 0; A < Nnodes_p; A++) {

    /*
      Assign from matrix
    */
    Pore_water_pressure_pA = Nodal_Pore_water_pressure_p.nV[A];
    gradient_A = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);

    /*
      Compute nodal contribution
    */
    for (int i = 0; i < Ndim; i++) {
      gradPw.n[i] += Pore_water_pressure_pA * gradient_A.n[i];
    }
  }

  return gradPw;
}

/**************************************************************/

static void
compute_Permeability_Inertial_Forces_Fluid(Matrix Mass_Exchanges_Source_Terms,
                                           Mask ActiveNodes, Particle MPM_Mesh,
                                           Mesh FEM_Mesh)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int NumNodes_p;  /* Number of tributary nodes of p */
  int A_mask;      /* Index of the node where we apply the body force */
  int Ap;          /* Tributary node A of particle p */
  int B_mask;      /* Index of the node where we apply the body force */
  int Bp;          /* Tributary node B of particle p */
  int Mixture_idx; /* Index for the material point mixture parameters */

  Element Nodes_p;        /* Element for each particle */
  Matrix gradient_p;      /* Shape functions gradients */
  Tensor gradient_pA;     /* Shape functions gradients (Node A), def config */
  Tensor GRADIENT_pA;     /* Shape functions gradients (Node A), ref config */
  Tensor F_n_p;           /* Deformation gradient t = n */
  Tensor F_n1_p;          /* Deformation gradient t = n + 1 */
  Tensor transpose_F_n_p; /* Transpose of the deformation gradient t = n */
  Tensor inverse_F_n1_p;  /* Inverse of the deformation gradient t = n + 1 */
  Tensor k_p;             /* Spatial permebility tensor */
  Tensor Fk_p;            /* One side pull-back of the permeability tensor */
  Tensor a_p;             /* Particle acceleration */
  Tensor b_p;             /* External acceleration */
  Tensor dyn_p;           /* Total acceleration of the particle */
  Tensor Fk_dyn_p;        /* Axiliar tensor for intermediate result */
  double GRADIENT_pA__x__Fk_dyn_p; /* Auxiliar scalar for intermediate result */
  double rho_f_p;   /* Intrinsic or material density (fluid phase) */
  double g = -9.81; /* Gravity constant */
  double V0_p;      /* Inital volume of the particle p */
  double J_n1_p; /* Determianant of the soil skeleton deformation gradient at t
                    = n + 1 */

  for (int p = 0; p < Np; p++) {

    /*
      Get the current material density for each material point (fluid)
    */
    rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];

    /*
      Get the reference volume for each material point (mixture)
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
      Take the value of the deformation gradient at t = n + 1, and t = n.
      Get intermediate results
    */
    F_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p], 2);
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p], 2);
    transpose_F_n_p = transpose__TensorLib__(F_n_p);
    inverse_F_n1_p = Inverse__TensorLib__(F_n1_p);
    J_n1_p = MPM_Mesh.Phi.J_n1.nV[p];

    /*
      Get the permeability tensor
    */
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    k_p = Soil_Water_Mixtures[Mixture_idx].Permeability;

    /*
      Compute intermediate result
    */
    Fk_p = matrix_product_old__TensorLib__(inverse_F_n1_p, k_p);

    /*
      Compute total particle acceleration
    */
    a_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.acc.nM[p], 1);
//    b_p = MPM_Mesh.b;
    dyn_p = subtraction__TensorLib__(a_p, b_p);

    /*
      Compute intermediate result
    */
    Fk_dyn_p = vector_linear_mapping__TensorLib__(Fk_p, dyn_p);

    for (int A = 0; A < NumNodes_p; A++) {

      /*
        Get the gradient of the shape function and multiply it by the
        permeability tensor
      */
      gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);
      GRADIENT_pA =
          vector_linear_mapping__TensorLib__(transpose_F_n_p, gradient_pA);

      /*
        Get the node of the mesh for the contribution
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Compute intermediate result
      */
      GRADIENT_pA__x__Fk_dyn_p =
          inner_product__TensorLib__(GRADIENT_pA, Fk_dyn_p);

      /*
        Add nodal contributions
      */
      Mass_Exchanges_Source_Terms.nV[A_mask] +=
          (J_n1_p * rho_f_p / g) * GRADIENT_pA__x__Fk_dyn_p * V0_p;

      /*
        Free some auxiliar resutls
      */
      free__TensorLib__(GRADIENT_pA);
    }

    /*
     Free memory
    */
    free__TensorLib__(transpose_F_n_p);
    free__TensorLib__(inverse_F_n1_p);
    free__TensorLib__(Fk_p);
    free__TensorLib__(dyn_p);
    free__TensorLib__(Fk_dyn_p);
    free__MatrixLib__(gradient_p);
    free(Nodes_p.Connectivity);
  }
}

/**************************************************************/

static void solve_Nodal_Mass_Balance(Matrix Compressibility_Matrix_Fluid,
                                     Matrix Mass_Exchanges_Source_Terms,
                                     Particle MPM_Mesh, Mesh FEM_Mesh,
                                     Mask ActiveNodes,
                                     Mask Free_and_Restricted_Dofs)
/*

*/
{
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nnodes_mask = Mass_Exchanges_Source_Terms.N_rows;
  int AB_indx;
  int Np = MPM_Mesh.NumGP;
  int NumNodes_p;
  int Ap;
  int A_mask;
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p;
  double ShapeFunction_pA;

  Matrix Rate_Pore_water_pressure = allocZ__MatrixLib__(Nnodes_mask, 1);

  /*
    The solution is now stored in the fluid forces vector
  */
  for (int A = 0; A < Nnodes_mask; A++) {
    if (Free_and_Restricted_Dofs.Nodes2Mask[A * Ndof + Ndim] != -1) {
      Rate_Pore_water_pressure.nV[A] = Mass_Exchanges_Source_Terms.nV[A] /
                                       Compressibility_Matrix_Fluid.nM[A][A];
    } else {
      Rate_Pore_water_pressure.nV[A] = 0.0;
    }
  }

  /*
    Compute the particle rate of pore water pressure
  */
  for (int p = 0; p < Np; p++) {
    /*
      Set to zero the rate of pore water pressure
    */
    MPM_Mesh.Phi.d_Pw_dt_n1.nV[p] = 0.0;

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
      Iterate over the nodes of the particle
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

      /*
        Update the particle rate of pore water pressure
      */
      MPM_Mesh.Phi.d_Pw_dt_n1.nV[p] +=
          ShapeFunction_pA * Rate_Pore_water_pressure.nV[A_mask];
    }

    /*
      Free shape function and particle connectivity
    */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

  /*
    Free auxiliar rate of pore water pressure
  */
  free__MatrixLib__(Rate_Pore_water_pressure);
}

/**************************************************************/

static void compute_Explicit_Newmark_Corrector(Particle MPM_Mesh, double gamma,
                                               double Dt) {
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  Tensor F_n_p;
  Tensor F_n1_p;

  for (int p = 0; p < Np; p++) {

    /*
      Replace the deformation gradient at t = n with the new one
    */
    F_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p], 2);
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p], 2);

    /*
      Update/correct tensor and vector variables
    */
    for (int i = 0; i < Ndim; i++) {

      /*
        Correct particle velocity
      */
      MPM_Mesh.Phi.vel.nM[p][i] += gamma * Dt * MPM_Mesh.Phi.acc.nM[p][i];

      /*
        Update the particles position
      */
      MPM_Mesh.Phi.x_GC.nM[p][i] += MPM_Mesh.Phi.D_dis.nM[p][i];

      /* Update deformation gradient tensor */
      for (int j = 0; j < Ndim; j++) {
        F_n_p.N[i][j] = F_n1_p.N[i][j];
      }
    }

    /*
      Replace the determinant of the deformation gradient
    */
    MPM_Mesh.Phi.J_n.nV[p] = MPM_Mesh.Phi.J_n1.nV[p];

    /*
      Correct pressure field
    */
    MPM_Mesh.Phi.Pw.nV[p] += gamma * Dt * MPM_Mesh.Phi.d_Pw_dt_n1.nV[p];
  }
}

/**************************************************************/

static void output_selector(Particle MPM_Mesh, Mesh FEM_Mesh, Mask ActiveNodes,
                            Matrix Velocity, Matrix D_Displacement,
                            Matrix Forces, Matrix Reactions,
                            double DeltaTimeStep, int TimeStep,
                            int ResultsTimeStep)
/*

*/
{
  /*
    vtk results
  */
  if (TimeStep % ResultsTimeStep == 0) {
    particle_results_vtk__InOutFun__(MPM_Mesh, TimeStep, ResultsTimeStep);

    nodal_results_vtk__InOutFun__(FEM_Mesh, ActiveNodes, Reactions, TimeStep,
                                  ResultsTimeStep);
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

static void standard_error() {
  fprintf(stderr, "%s !!! \n", Error_message);
  exit(EXIT_FAILURE);
}

/***************************************************************************/
