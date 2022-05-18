#include "Formulations/Displacements/U-Forward-Euler.h"


/*
  Auxiliar functions
*/
static void update_Particles(Particle, Mesh, Matrix, Matrix, double);
static Matrix compute_NodalMomentumMass(Particle, Mesh);
static void imposed_Momentum(Mesh, Matrix, int, int);
static Matrix compute_Nodal_Velocity(Mesh, Matrix);
static void update_Nodal_Momentum(Mesh, Matrix, Matrix, double);
static void update_LocalState(Matrix, Particle, Mesh, double);
static Matrix compute_InternalForces(Matrix, Particle, Mesh);
static Matrix compute_ContacForces(Matrix, Particle, Mesh, int, int);
static Matrix compute_Reactions(Mesh, Matrix, int, int);

/**************************************************************/

void U_Forward_Euler(Mesh FEM_Mesh, Particle MPM_Mesh,
                     Time_Int_Params Parameters_Solver) {

  /*
    Auxiliar variables
  */
  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int InitialTimeStep = Parameters_Solver.InitialTimeStep;
  int NumTimeStep = Parameters_Solver.NumTimeStep;
  double CFL = Parameters_Solver.CFL;
  double DeltaTimeStep;

  /*
    Auxiliar variable for the mass and momentum
  */
  Matrix Phi_I;
  Matrix V_I;
  Matrix F_I;
  Matrix R_I;

  for (int TimeStep = InitialTimeStep; TimeStep < NumTimeStep; TimeStep++) {

    DeltaTimeStep =
        U_DeltaT__SolversLib__(MPM_Mesh, FEM_Mesh.DeltaX, Parameters_Solver);
    print_step(TimeStep, NumTimeStep, DeltaTimeStep);
    local_search__MeshTools__(MPM_Mesh, FEM_Mesh);

    print_Status("*************************************************", TimeStep);
    print_Status(" First step : Get the nodal fields ... WORKING", TimeStep);
    Phi_I = compute_NodalMomentumMass(MPM_Mesh, FEM_Mesh);
    imposed_Momentum(FEM_Mesh, Phi_I, TimeStep, NumTimeStep);
    V_I = compute_Nodal_Velocity(FEM_Mesh, Phi_I);
    print_Status("DONE !!!", TimeStep);

    print_Status("*************************************************", TimeStep);
    print_Status("Second step : Compute equilibrium ... WORKING", TimeStep);
    update_LocalState(V_I, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
    F_I = allocZ__MatrixLib__(Nnodes, Ndim);
    F_I = compute_InternalForces(F_I, MPM_Mesh, FEM_Mesh);
    F_I = compute_ContacForces(F_I, MPM_Mesh, FEM_Mesh, TimeStep, NumTimeStep);
    R_I = compute_Reactions(FEM_Mesh, F_I, TimeStep, NumTimeStep);
    print_Status("DONE !!!", TimeStep);

    print_Status("*************************************************", TimeStep);
    print_Status(" Third step : Update nodal momentum ... WORKING", TimeStep);
    update_Nodal_Momentum(FEM_Mesh, Phi_I, F_I, DeltaTimeStep);
    print_Status("DONE !!!", TimeStep);

    print_Status("*************************************************", TimeStep);
    print_Status(" Four step : Update lagrangian ... WORKING", TimeStep);
    update_Particles(MPM_Mesh, FEM_Mesh, Phi_I, F_I, DeltaTimeStep);
    print_Status("DONE !!!", TimeStep);

    if (TimeStep % ResultsTimeStep == 0) {
      /*!
        Print Nodal values after appling the BCCs
      */
      //       nodal_results_vtk__InOutFun__("Mesh",FEM_Mesh,R_I,TimeStep,(int)TimeStep/ResultsTimeStep);
      /*!
        Print particle results
      */
      particle_results_vtk__InOutFun__(
          MPM_Mesh, (int)TimeStep / ResultsTimeStep, ResultsTimeStep);
    }

    print_Status("*************************************************", TimeStep);
    print_Status(" Five step : Reset nodal values ... WORKING", TimeStep);
    free__MatrixLib__(Phi_I);
    free__MatrixLib__(V_I);
    free__MatrixLib__(F_I);
    free__MatrixLib__(R_I);
    print_Status("DONE !!!", TimeStep);
  }
}

/*******************************************************/

static Matrix compute_NodalMomentumMass(Particle MPM_Mesh, Mesh FEM_Mesh)
/*
  This function performs a information trasference from the lagrangian particles
  to the nodes of the eulerian mesh. Output : Phi_I = {P_I | M_I}
*/
{
  /* Define some dimensions */
  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Np = MPM_Mesh.NumGP;
  int Ip;

  /* Auxiliar variables */
  Matrix ShapeFunction_p;  /* Value of the shape-function */
  double ShapeFunction_pI; /* Evaluation of the GP in the node */
  double m_p;              /* Mass of the GP */
  double M_Ip;
  Element Nodes_p; /* Element for each Gauss-Point */

  /* Output */
  Matrix Phi_I;

  /* Allocate the output list of fields */
  Phi_I = allocZ__MatrixLib__(Nnodes, Ndim + 1);

  /* Iterate over the GP to get the nodal values */
  for (int p = 0; p < Np; p++) {

    /* Define element of the GP */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],
                                     MPM_Mesh.NumberNodes[p]);

    /* Evaluate the shape function in the coordinates of the GP */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Get the mass of the GP */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* Get the nodal mass and mommentum */
    for (int A = 0; A < Nodes_p.NumberNodes; A++) {

      /* Get the node for the GP */
      Ip = Nodes_p.Connectivity[A];

      /* Evaluate the GP function in the node */
      ShapeFunction_pI = ShapeFunction_p.nV[A];

      /* If this node has a null Value of the SHF continue */
      if (ShapeFunction_pI == 0) {
        continue;
      }

      /* Gauss-point contribution to node I */
      M_Ip = m_p * ShapeFunction_pI;

      /* Nodal mass */
      Phi_I.nM[Ip][Ndim] += M_Ip;

      /* Nodal momentum */
      for (int i = 0; i < Ndim; i++) {
        Phi_I.nM[Ip][i] += M_Ip * MPM_Mesh.Phi.vel.nM[p][i];
      }
    }

    /* Free the value of the shape functions */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

  return Phi_I;
}

/**********************************************************************/

static void imposed_Momentum(Mesh FEM_Mesh, Matrix Phi_I, int TimeStep,
                             int NumTimeStep)
/*
  Apply the boundary conditions over the nodes
*/
{

  /* 1º Define auxilar variables */
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound;   /* Number of dimensions */
  int Id_BCC;        /* Index of the node where we apply the BCC */

  /* 2º Loop over the the boundaries */
  for (int i = 0; i < FEM_Mesh.Bounds.NumBounds; i++) {
    /* 3º Get the number of nodes of this boundarie */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;
    /* 4º Get the number of dimensions where the BCC it is applied */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;
    for (int j = 0; j < NumNodesBound; j++) {
      /* 5º Get the index of the node */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      /* 6º Loop over the dimensions of the boundary condition */
      for (int k = 0; k < NumDimBound; k++) {
        /* 7º Apply only if the direction is active (1) */
        if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + TimeStep] == 1) {
          /* 9º Assign the boundary condition */
          Phi_I.nM[Id_BCC][k] = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep];
        }
      }
    }
  }
}

/**************************************************************/

static void update_Particles(Particle MPM_Mesh, Mesh FEM_Mesh, Matrix Phi_I,
                             Matrix F_I, double Dt) {

  int Ndim = NumberDimensions;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_p;      /* Value of the shape-function in the GP */
  double N_pI;     /* Nodal value for the GP */
  double M_I;      /* Value of the nodal mass */
  double D_U_pI;
  int Np = MPM_Mesh.NumGP;
  int Nnodes;
  int Ip; /* Index of each tributary node for the GP */

  /* 1º iterate over the Gauss-Points */
  for (int p = 0; p < Np; p++) {

    /* 2º Define element of the GP */
    Nnodes = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], Nnodes);

    /* 3º Evaluate shape function in the GP i */
    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* 4º Iterate over the nodes of the element */
    for (int A = 0; A < Nnodes; A++) {
      /* Node of the GP */
      Ip = Nodes_p.Connectivity[A];
      /* Evaluate the GP function in the node */
      N_pI = N_p.nV[A];
      /* If this node has a null Value of the SHF continue */
      if (fabs(N_pI) <= TOL_zero) {
        continue;
      }
      /* Get the nodal mass */
      M_I = Phi_I.nM[Ip][Ndim];
      /* Update GP cuantities with nodal values */
      for (int i = 0; i < Ndim; i++) {
        /* Update the GP velocities */
        MPM_Mesh.Phi.vel.nM[p][i] += Dt * N_pI * F_I.nM[Ip][i] / M_I;

        D_U_pI = Dt * N_pI * Phi_I.nM[Ip][i] / M_I;

        /* Update the GP displacemetn */
        MPM_Mesh.Phi.dis.nM[p][i] += D_U_pI;

        /* Update the GP position */
        MPM_Mesh.Phi.x_GC.nM[p][i] += D_U_pI;
      }
    }

    /* 5º Free memory */
    free(Nodes_p.Connectivity);
    free__MatrixLib__(N_p);
  }
}

/*******************************************************/

static Matrix compute_Nodal_Velocity(Mesh FEM_Mesh, Matrix Phi_I) {
  /*
    Get the nodal velocity using :
    v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
    Initialize nodal velocities
  */

  /* Define some dimensions */
  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  /* Define output */
  Matrix V_I = allocZ__MatrixLib__(Nnodes, Ndim);

  /* Value of the nodal mass */
  double M_I;

  /* 1º Get nodal values of the velocity */
  for (int i = 0; i < Nnodes; i++) {
    M_I = Phi_I.nM[i][Ndim];
    if (M_I > 0) {
      for (int j = 0; j < Ndim; j++) {
        V_I.nM[i][j] = Phi_I.nM[i][j] / M_I;
      }
    }
  }

  return V_I;
}

/*******************************************************/

static void update_Nodal_Momentum(Mesh FEM_Mesh, Matrix Phi_I, Matrix F_I,
                                  double DeltaTimeStep)
/*!
 * \brief Brief description of UpdateGridNodalMomentum.
 *        Compute the nodal contribution of each GP to the total forces.
 *
 *  The parameters for this functions are  :
 *  @param MPM_Mesh : Mesh with the material points.
 *  @param Phi_I {P_I | M_I}
 *  @param F_I : Nodal value of the total forces.
 *
 */
{
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Ndim = NumberDimensions;

  /* Update the grid nodal momentum */
  for (int i = 0; i < Nnodes; i++) {
    for (int j = 0; j < Ndim; j++) {
      if (FEM_Mesh.ActiveNode[i]) {
        Phi_I.nM[i][j] += DeltaTimeStep * F_I.nM[i][j];
      }
    }
  }
}

/**************************************************************/

static void update_LocalState(Matrix V_I, Particle MPM_Mesh, Mesh FEM_Mesh,
                              double TimeStep) {
  Element Nodes_p;         /* Element for each Gauss-Point */
  Matrix Gradient_p;       /* Shape functions gradients */
  Matrix Nodal_Velocity_p; /* Velocity of the element nodes */
  Material Material_p;     /* Properties of the Gauss-Point material */
  Tensor Rate_Strain_p;    /* Increment of strain tensor */
  Tensor Strain_p;         /*  Strain tensor */
  Tensor Stress_p;         /* Stress tensor */
  double rho_p;            /* Density of the Gauss-Point */
  int Np = MPM_Mesh.NumGP;
  int Nn;
  int Idx_Mat_p;

  /* Loop in the GPs */
  for (int p = 0; p < Np; p++) {

    /* Get the value of the density */
    rho_p = MPM_Mesh.Phi.rho.nV[p];

    /* Asign memory to tensors */
    Strain_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Strain.nM[p], 2);
    Stress_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);

    /* Define element for each GP */
    Nn = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], Nn);

    /* Get the velocity of the nodes of the element */
    Nodal_Velocity_p = get_set_field_old__MeshTools__(V_I, Nodes_p);

    /* Compute gradient of the shape function in each node */
    Gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Get the material properties */
    Idx_Mat_p = MPM_Mesh.MatIdx[p];
    Material_p = MPM_Mesh.Mat[Idx_Mat_p];


    /* Update density field */
//    MPM_Mesh.Phi.rho.nV[p] =
//        update_density__Particles__(rho_p, TimeStep, Rate_Strain_p);
    free__TensorLib__(Rate_Strain_p);

    /* Compute stress tensor */
//int Stress_integration__Constitutive__(
//    int p, 
//    Particle MPM_Mesh, 
//    Material MatProp_p);

    /* Compute deformation energy */
    MPM_Mesh.Phi.W[p] = 0.5 * inner_product__TensorLib__(Strain_p, Stress_p);

    /* Free the matrix with the nodal velocity of the element */
    free__MatrixLib__(Nodal_Velocity_p);

    /* Free the matrix with the nodal gradient of the element */
    free__MatrixLib__(Gradient_p);
    free(Nodes_p.Connectivity);
  }

}

/*************************************************************/

static Matrix compute_InternalForces(Matrix F_I, Particle MPM_Mesh,
                                     Mesh FEM_Mesh) {
  int Ndim = NumberDimensions;
  Element Nodes_p;   /* Element for each Gauss-Point */
  Matrix Gradient_p; /* Shape functions gradients */
  Tensor Stress_p;   /* Stress tensor */
  Tensor Gradient_pI;
  Tensor InternalForcesDensity_Ip;
  double m_p;   /* Mass of the Gauss-Point */
  double rho_p; /* Density of the Gauss-Point */
  double V_p;   /* Volume of the Gauss-Point */
  int Np = MPM_Mesh.NumGP;
  int Ip;
  int Nn;

  /* Loop in the GPs */
  for (int p = 0; p < Np; p++) {

    /* Get the value of the density */
    rho_p = MPM_Mesh.Phi.rho.nV[p];

    /* Get the value of the mass */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /* Asign memory to tensors */
    Stress_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);

    /* Define element for each GP */
    Nn = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], Nn);

    /* Compute gradient of the shape function in each node */
    Gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Compute the volume of the Gauss-Point */
    V_p = m_p / rho_p;

    /* Compute nodal forces */
    for (int A = 0; A < Nn; A++) {
      /* Pass by reference the nodal gradient to the tensor */
      Gradient_pI = memory_to_tensor__TensorLib__(Gradient_p.nM[A], 1);

      /* Compute the nodal forces of the Gauss-Point */
      InternalForcesDensity_Ip =
          vector_linear_mapping__TensorLib__(Stress_p, Gradient_pI);

      /* Get the node of the mesh for the contribution */
      Ip = Nodes_p.Connectivity[A];

      /* Asign the nodal forces contribution to the node */
      for (int i = 0; i < Ndim; i++) {
        F_I.nM[Ip][i] -= InternalForcesDensity_Ip.n[i] * V_p;
      }

      /* Free the internal forces density */
      free__TensorLib__(InternalForcesDensity_Ip);
    }

    /* Free the matrix with the nodal gradient of the element */
    free__MatrixLib__(Gradient_p);
    free(Nodes_p.Connectivity);
  }

  return F_I;
}

/*********************************************************************/

static Matrix compute_ContacForces(Matrix F_I, Particle MPM_Mesh, Mesh FEM_Mesh,
                                   int TimeStep, int NumTimeStep) {
  int Ndim = NumberDimensions;
  Load Load_i;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_p;      /* Nodal values of the sahpe function */
  double N_pa;
  Tensor t = alloc__TensorLib__(1); /* Body forces vector */
  double m_p;                       /* Mass of the Gauss-Point */
  double rho_p;                     /* Density of the Gauss-Point */
  double V_p;                       /* Volumen of the Gauss-Point */
  double A_p;                       /* Area of the Gauss-Point */
  int NumContactForces = MPM_Mesh.Neumann_Contours.NumBounds;
  int NumNodesLoad;
  int p;
  int Ip;
  int Nn; /* Number of nodes of each Gauss-Point */

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
      m_p = MPM_Mesh.Phi.mass.nV[p];
      rho_p = MPM_Mesh.Phi.mass.nV[p];
      V_p = m_p / rho_p;
      A_p = V_p / Thickness_Plain_Stress;

      /* Define element for each GP */
      Nn = MPM_Mesh.NumberNodes[p];
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], Nn);

      /* Compute shape functions */
      N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Fill vector of body forces */
      for (int k = 0; k < Ndim; k++) {
        if (Load_i.Dir[k * NumTimeStep + TimeStep]) {
          t.n[k] = Load_i.Value[k].Fx[TimeStep];
        }
      }

      /* Get the node of the mesh for the contribution */
      for (int A = 0; A < Nn; A++) {
        /* Node for the contribution */
        Ip = Nodes_p.Connectivity[A];

        /* Pass the value of the nodal shape function to a scalar */
        N_pa = N_p.nV[A];

        /* Compute Contact forces */
        for (int k = 0; k < Ndim; k++) {
          F_I.nM[Ip][k] += N_pa * t.n[k] * A_p;
        }
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(N_p);
      free(Nodes_p.Connectivity);
    }
  }

  free__TensorLib__(t);

  return F_I;
}

/**********************************************************************/

static Matrix compute_Reactions(Mesh FEM_Mesh, Matrix F_I, int TimeStep,
                                int NumTimeStep)
/*
  Compute the nodal reactions
*/
{
  /* 1º Define auxilar variables */
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound;   /* Number of dimensions */
  int Id_BCC;        /* Index of the node where we apply the BCC */
  int Ndim = NumberDimensions;

  Matrix R_I = allocZ__MatrixLib__(FEM_Mesh.NumNodesMesh, Ndim);


  /* 2º Loop over the the boundaries */
  for (int i = 0; i < FEM_Mesh.Bounds.NumBounds; i++) {
    /* 3º Get the number of nodes of this boundarie */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;
    /* 4º Get the number of dimensions where the BCC it is applied */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;
    for (int j = 0; j < NumNodesBound; j++) {
      /* 5º Get the index of the node */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      /* 6º Loop over the dimensions of the boundary condition */
      for (int k = 0; k < NumDimBound; k++) {
        /* 7º Apply only if the direction is active (1) */
        if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + TimeStep] == 1) {
          /* 8º Set to zero the forces in the nodes where velocity is fixed */
          R_I.nM[Id_BCC][k] = F_I.nM[Id_BCC][k];
          F_I.nM[Id_BCC][k] = 0;
        }
      }
    }
  }

  return R_I;
}

/**********************************************************************/
