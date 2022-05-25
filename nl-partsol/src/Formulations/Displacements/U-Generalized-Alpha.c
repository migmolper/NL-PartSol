#include "Formulations/Displacements/U-Generalized-Alpha.h"
#include "Globals.h"

/*
  Auxiliar functions
*/
static Matrix compute_Nodal_Kinetics(Particle, Mesh);
static void update_Kinetics(Mesh, Matrix, Matrix, double, double);
static Matrix GetNodalVelocityDisplacement(Particle, Mesh);
static void update_Particles(Particle, Mesh, Matrix, double, double);
static void update_LocalState(Matrix, Particle, Mesh, double);
static Matrix compute_InternalForces(Matrix, Particle, Mesh);
static Matrix compute_ContacForces(Matrix, Particle, Mesh, int, int);
static Matrix compute_Reactions(Mesh, Matrix, int, int);

/**************************************************************/

void U_Generalized_alpha(Mesh FEM_Mesh, Particle MPM_Mesh,
                         Time_Int_Params Parameters_Solver) {

  /*!
    Time step
  */
  int TimeStep;

  /*!
    Auxiliar variable for the nodal kinetics
    Nodal_Kinetics = {mass, a0, a1, v}
  */
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Ndim = NumberDimensions;
  int InitialStep = Parameters_Solver.InitialTimeStep;
  int NumTimeStep = Parameters_Solver.NumTimeStep;

  double rb = Parameters_Solver.rb_Generalized_alpha;
  double CFL = Parameters_Solver.CFL;
  double DeltaTimeStep;
  double DeltaX = FEM_Mesh.DeltaX;

  Matrix V_I;
  Matrix Nodal_Kinetics;
  Matrix F_I = memory_to_matrix__MatrixLib__(Ndim, Nnodes, NULL);
  Matrix R_I = memory_to_matrix__MatrixLib__(Ndim, Nnodes, NULL);

  puts("*************************************************");
  puts(" First step : Get the nodal kinetics");
  puts(" \t WORKING ...");
  Nodal_Kinetics = compute_Nodal_Kinetics(MPM_Mesh, FEM_Mesh);
  /* V_I =  MatAssign(Nnodes,Ndim,NAN,NULL, */
  /* 		   (double**)malloc(Nnodes*sizeof(double *))); */
  for (int i = 0; i < Ndim; i++) {
    V_I.nM[i] = Nodal_Kinetics.nM[1 + 2 * Ndim + i];
  }
  puts(" \t DONE !!! \n");

  for (TimeStep = InitialStep; TimeStep < NumTimeStep; TimeStep++) {

    DeltaTimeStep = U_DeltaT__SolversLib__(MPM_Mesh, DeltaX, Parameters_Solver);
    print_step(TimeStep,NumTimeStep, DeltaTimeStep);
    local_search__MeshTools__(MPM_Mesh, FEM_Mesh);
    print_Status("*************************************************", TimeStep);
    print_Status("First step : Compute nodal kinetics ... WORKING", TimeStep);
    /* BCC_Nod_VALUE(FEM_Mesh, V_I, TimeStep); */
    print_Status("DONE !!!", TimeStep);

    print_Status("*************************************************", TimeStep);
    print_Status("Second step : Compute equilibrium ... WORKING", TimeStep);
    F_I = allocZ__MatrixLib__(Nnodes, Ndim);
    F_I = compute_InternalForces(F_I, MPM_Mesh, FEM_Mesh);
    F_I = compute_ContacForces(F_I, MPM_Mesh, FEM_Mesh, TimeStep, NumTimeStep);
    R_I = compute_Reactions(FEM_Mesh, F_I, TimeStep, NumTimeStep);
    print_Status("DONE !!!", TimeStep);

    print_Status("*************************************************", TimeStep);
    print_Status(" Third step : Update kinetics ... WORKING", TimeStep);
    update_Kinetics(FEM_Mesh, Nodal_Kinetics, F_I, rb, DeltaTimeStep);
    print_Status("DONE !!!", TimeStep);

    print_Status("*************************************************", TimeStep);
    print_Status("Four step : Update particles lagrangian ... WORKING",
                 TimeStep);
    update_Particles(MPM_Mesh, FEM_Mesh, Nodal_Kinetics, rb, DeltaTimeStep);
    print_Status("DONE !!!", TimeStep);

    if (TimeStep % ResultsTimeStep == 0) {
      /*!
       Print Nodal values after appling the BCCs
       */
      //        nodal_results_vtk__InOutFun__("Mesh",FEM_Mesh,R_I,TimeStep,(int)TimeStep/ResultsTimeStep);
      /*!
      Print GPs results
      */
      particle_results_vtk__InOutFun__(
          MPM_Mesh, (int)TimeStep / ResultsTimeStep, ResultsTimeStep);
    }

    print_Status("*************************************************", TimeStep);
    print_Status("Five step : Reset nodal values ... WORKING", TimeStep);
    free__MatrixLib__(F_I);
    free__MatrixLib__(R_I);
    print_Status("DONE !!!", TimeStep);
  }
}

/*******************************************************/

static void update_Kinetics(Mesh FEM_Mesh, Matrix Nodal_Kinetics,
                            Matrix Nodal_Forces, double rb,
                            double DeltaTimeStep)
/*!
 * \brief Brief description of Update_Nodal_Acceleration_Velocity.
 *
 *  The parameters for this functions are  :
 *  @param FEM_Mesh
 *  @param Nodal_Kinetics = {m, a0, a1, v}
 *  @param Nodal_Forces
 *  @param Params
 */
{
  int N_Nodes = FEM_Mesh.NumNodesMesh;
  int Ndim = NumberDimensions;
  double SizeTable = Ndim * sizeof(double *);
  double Mass_I;

  /* Time integration parameters */
  double alpha = (2 * rb - 1) / (1 + rb);
  double beta = (5 - 3 * rb) / (pow((1 + rb), 2) * (2 - rb));
  double gamma = 3 / 2 - alpha;

  /*!
  Control parameters of the generalized-alpha algorithm
  all the parameters are controled by a simple parameter :
  SpectralRadius
*/
  Time_Int_Params Params;

  /* Asign forces to an auxiliar variable */
  Matrix F = Nodal_Forces;

  /* Nodal values the fields */
  Matrix Nodal_Acceleration_t0;
  /* MatAssign(Ndim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */
  Matrix Nodal_Acceleration_t1;
  /* MatAssign(Ndim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */
  Matrix Nodal_Velocity;
  /* MatAssign(Ndim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */

  /* 1º Asign Kinetics values to matricial tables */
  for (int i = 0; i < Ndim; i++) {
    Nodal_Acceleration_t0.nM[i] = Nodal_Kinetics.nM[1 + i];
    Nodal_Acceleration_t1.nM[i] = Nodal_Kinetics.nM[1 + Ndim + i];
    Nodal_Velocity.nM[i] = Nodal_Kinetics.nM[1 + 2 * Ndim + i];
  }
  /* 2º Update the grid nodal variales */
  for (int i = 0; i < N_Nodes; i++) {
    Mass_I = Nodal_Kinetics.nM[0][i];
    if (Mass_I > 0) {
      for (int j = 0; j < Ndim; j++) {
        /* Get the nodal acceleration t1 */
        Nodal_Acceleration_t1.nM[j][i] =
            (F.nM[j][i] / Mass_I - alpha * Nodal_Acceleration_t0.nM[j][i]) /
            (1 - alpha);
        /* Update nodal velocity */
        Nodal_Velocity.nM[j][i] +=
            ((1 - gamma) * Nodal_Acceleration_t0.nM[j][i] +
             gamma * Nodal_Acceleration_t1.nM[j][i]) *
            DeltaTimeStep;
      }
    }
  }

  /* 3º Free tables */
  free(Nodal_Acceleration_t0.nM);
  free(Nodal_Acceleration_t1.nM);
  free(Nodal_Velocity.nM);
}

/*********************************************************************/

static Matrix compute_Nodal_Kinetics(Particle MPM_Mesh, Mesh FEM_Mesh)
/*!
 *  Nodal_Kinetics = {mass, a0, a1, v}
 */
{
  /* */
  int N_Nodes_Mesh = FEM_Mesh.NumNodesMesh;
  int N_GPs = MPM_Mesh.NumGP;

  /* Sizes */
  int Ndim = NumberDimensions;
  int N_Acceleration_dim = Ndim;
  int N_Velocity_dim = Ndim;

  /* Nodal values the fields */
  Matrix Nodal_Mass = allocZ__MatrixLib__(1, N_Nodes_Mesh);
  Matrix Nodal_Acceleration_t0 =
      allocZ__MatrixLib__(N_Acceleration_dim, N_Nodes_Mesh);
  Matrix Nodal_Acceleration_t1 =
      allocZ__MatrixLib__(N_Acceleration_dim, N_Nodes_Mesh);
  Matrix Nodal_Velocity = allocZ__MatrixLib__(N_Velocity_dim, N_Nodes_Mesh);

  /* GPs values of the fields */
  double Mass_GP; /* Mass of the GP */
  Matrix Vel_GP;  /* Velocity of the GP */
  Matrix Acc_GP;  /* Stress of the GP */

  /* Shape function auxiliar variables */
  Element Nodes_p;  /* Element for each Gauss-Point */
  Matrix N_GP;      /* Value of the shape-function */
  double N_GP_I;    /* Evaluation of the GP in the node */
  double Mass_GP_I; /* Nodal contribution of each GP */
  int GP_I;

  /* 1º Iterate over the GP to get the nodal values */
  for (int i = 0; i < N_GPs; i++) {

    /* 2º Define element of the GP */
    Nodes_p = nodal_set__Particles__(i, MPM_Mesh.ListNodes[i],
                                     MPM_Mesh.NumberNodes[i]);

    /* 3º Evaluate the shape function in the coordinates of the GP */
    N_GP = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* 4º Get the properties of the GP */
    Mass_GP = MPM_Mesh.Phi.mass.nV[i];
    Acc_GP.nV = MPM_Mesh.Phi.acc.nM[i];
    Vel_GP.nV = MPM_Mesh.Phi.vel.nM[i];

    /* 5º Get the nodal kinetics (I) */
    for (int k = 0; k < Nodes_p.NumberNodes; k++) {
      /* Get the node for the GP */
      GP_I = Nodes_p.Connectivity[k];
      /* Evaluate the GP function in the node */
      N_GP_I = N_GP.nV[k];
      /* If this node has a null Value of the SHF continue */
      if (N_GP_I == 0) {
        continue;
      }
      /* Nodal constribution */
      Mass_GP_I = Mass_GP * N_GP_I;
      /* Nodal mass */
      Nodal_Mass.nV[GP_I] += Mass_GP_I;
      /* Nodal acceleration t0 */
      for (int l = 0; l < N_Acceleration_dim; l++) {
        Nodal_Acceleration_t0.nM[l][GP_I] += Acc_GP.nV[l] * Mass_GP_I;
      }
      /* Nodal velocity */
      for (int l = 0; l < N_Velocity_dim; l++) {
        Nodal_Velocity.nM[l][GP_I] += Vel_GP.nV[l] * Mass_GP_I;
      }
    }
    /* 6º Free the value of the shape functions */
    free__MatrixLib__(N_GP);
    free(Nodes_p.Connectivity);
  }

  /* 7º Get the nodal kinetics (II) */
  for (int i = 0; i < N_Nodes_Mesh; i++) {
    if (Nodal_Mass.nV[i] > 0) {
      /* Get the nodal acceleration */
      for (int j = 0; j < N_Acceleration_dim; j++) {
        Nodal_Acceleration_t0.nM[j][i] =
            Nodal_Acceleration_t0.nM[j][i] / Nodal_Mass.nV[i];
      }
      /* Get the nodal velocity */
      for (int j = 0; j < N_Velocity_dim; j++) {
        Nodal_Velocity.nM[j][i] = Nodal_Velocity.nM[j][i] / Nodal_Mass.nV[i];
      }
    }
  }

  /* 8º Asign Kinetics values to matricial tables */

  int N_Kinetics_dim = 1 + 2 * N_Acceleration_dim + N_Velocity_dim;
  double SizeKinetics = N_Kinetics_dim * sizeof(double *);

  Matrix Nodal_Kinetics;
  
  Nodal_Kinetics.nM[0] = Nodal_Mass.nV;
  for (int i = 0; i < Ndim; i++) {
    Nodal_Kinetics.nM[1 + i] = Nodal_Acceleration_t0.nM[i];
    Nodal_Kinetics.nM[1 + Ndim + i] = Nodal_Acceleration_t1.nM[i];
    Nodal_Kinetics.nM[1 + 2 * Ndim + i] = Nodal_Velocity.nM[i];
  }

  /* Free table of pointers */
  free(Nodal_Acceleration_t0.nM);
  free(Nodal_Acceleration_t1.nM);
  free(Nodal_Velocity.nM);

  /* Return the kinetics variables in the nodes */
  return Nodal_Kinetics;
}

/*******************************************************/

static Matrix GetNodalVelocityDisplacement(Particle MPM_Mesh, Mesh FEM_Mesh)
/*!
 *  Nodal_Kinetics = {m, d, v}
 */
{
  /* */
  int N_Nodes_Mesh = FEM_Mesh.NumNodesMesh;
  int N_GPs = MPM_Mesh.NumGP;

  /* Sizes */
  int Ndim = NumberDimensions;

  /* Nodal values the fields */
  Matrix Nodal_Mass = allocZ__MatrixLib__(1, N_Nodes_Mesh);
  Matrix Nodal_Displacement = allocZ__MatrixLib__(Ndim, N_Nodes_Mesh);
  Matrix Nodal_Velocity = allocZ__MatrixLib__(Ndim, N_Nodes_Mesh);
  Matrix Nodal_Acceleration = allocZ__MatrixLib__(Ndim, N_Nodes_Mesh);

  /* GPs values of the fields */
  double Mass_GP; /* Mass of the GP */
  Matrix Disp_GP; /* Stress of the GP */
  Matrix Vel_GP;  /* Velocity of the GP */
  Matrix Acel_GP; /* Velocity of the GP */

  /* Shape function auxiliar variables */
  Element Nodes_p;  /* Element for each Gauss-Point */
  Matrix N_GP;      /* Value of the shape-function */
  double N_GP_I;    /* Evaluation of the GP in the node */
  double Mass_GP_I; /* Nodal contribution of each GP */
  int GP_I;

  /* 1º Iterate over the GP to get the nodal values */
  for (int i = 0; i < N_GPs; i++) {

    /* 2º Define element of the GP */
    Nodes_p = nodal_set__Particles__(i, MPM_Mesh.ListNodes[i],
                                     MPM_Mesh.NumberNodes[i]);

    /* 3º Evaluate the shape function in the coordinates of the GP */
    N_GP = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* 4º Get the properties of the GP */
    Mass_GP = MPM_Mesh.Phi.mass.nV[i];
    Vel_GP.nV = MPM_Mesh.Phi.vel.nM[i];
    Disp_GP.nV = MPM_Mesh.Phi.dis.nM[i];
    Acel_GP.nV = MPM_Mesh.Phi.acc.nM[i];

    /* 5º Get the nodal kinetics (I) */
    for (int k = 0; k < Nodes_p.NumberNodes; k++) {

      /* Get the node for the GP */
      GP_I = Nodes_p.Connectivity[k];

      /* Evaluate the GP function in the node */
      N_GP_I = N_GP.nV[k];

      /* If this node has a null Value of the SHF continue */
      if (N_GP_I == 0) {
        continue;
      }

      /* Nodal constribution */
      Mass_GP_I = Mass_GP * N_GP_I;

      /* Nodal mass */
      Nodal_Mass.nV[GP_I] += Mass_GP_I;

      for (int l = 0; l < Ndim; l++) {
        /* Nodal displacement */
        Nodal_Displacement.nM[l][GP_I] += Disp_GP.nV[l] * Mass_GP_I;
        /* Nodal velocity */
        Nodal_Velocity.nM[l][GP_I] += Vel_GP.nV[l] * Mass_GP_I;
        /* Nodal acceleration */
        Nodal_Acceleration.nM[l][GP_I] += Acel_GP.nV[l] * Mass_GP_I;
      }
    }
    /* 6º Free the value of the shape functions */
    free__MatrixLib__(N_GP);
    free(Nodes_p.Connectivity);
  }

  /* 7º Get the nodal kinetics (II) */
  for (int i = 0; i < N_Nodes_Mesh; i++) {
    if (Nodal_Mass.nV[i] > 0) {
      /* Get the nodal displacement */
      for (int j = 0; j < Ndim; j++) {
        Nodal_Displacement.nM[j][i] =
            Nodal_Displacement.nM[j][i] / Nodal_Mass.nV[i];
      }
      /* Get the nodal velocity */
      for (int j = 0; j < Ndim; j++) {
        Nodal_Velocity.nM[j][i] = Nodal_Velocity.nM[j][i] / Nodal_Mass.nV[i];
      }
      /* Get the nodal acceleration */
      for (int j = 0; j < Ndim; j++) {
        Nodal_Acceleration.nM[j][i] =
            Nodal_Acceleration.nM[j][i] / Nodal_Mass.nV[i];
      }
    }
  }

  /* 8º Asign Kinetics values to matricial tables */

  int N_Kinetics_dim = 1 + 3 * Ndim;
  double SizeKinetics = N_Kinetics_dim * sizeof(double *);

  Matrix Nodal_Kinetics;

  Nodal_Kinetics.nM[0] = Nodal_Mass.nV;
  for (int i = 0; i < Ndim; i++) {
    Nodal_Kinetics.nM[1 + i] = Nodal_Displacement.nM[i];
    Nodal_Kinetics.nM[1 + Ndim + i] = Nodal_Velocity.nM[i];
    Nodal_Kinetics.nM[1 + 2 * Ndim + i] = Nodal_Acceleration.nM[i];
  }

  /* Free table of pointers */
  free(Nodal_Displacement.nM);
  free(Nodal_Velocity.nM);
  free(Nodal_Acceleration.nM);

  /* Return the kinetics variables in the nodes */
  return Nodal_Kinetics;
}

/*******************************************************/

static void update_Particles(Particle MPM_Mesh, Mesh FEM_Mesh,
                             Matrix Nodal_Kinetics, double rb,
                             double DeltaTimeStep) {

  int N_Nodes = FEM_Mesh.NumNodesMesh;
  int N_dim = NumberDimensions;
  double SizeTable = N_dim * sizeof(double *);

  /* Shape function nodal parameters */
  Matrix N_GP;        /* Value of the shape-function in the GP */
  Element GP_Element; /* Element for each Gauss-Point */
  int GP_I;           /* Index of each tributary node for the GP */
  double N_I_GP;      /* Nodal value for the GP */

  /* Time integration parameters */
  double alpha = (2 * rb - 1) / (1 + rb);
  double beta = (5 - 3 * rb) / (pow((1 + rb), 2) * (2 - rb));
  double gamma = 3 / 2 - alpha;

  /* Asign Kinetics values to matricial tables */
  Matrix Nodal_Acceleration_t0;
  /* = MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */
  Matrix Nodal_Acceleration_t1;
  /* = MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */
  Matrix Nodal_Velocity;
  /* = MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable)); */
  for (int i = 0; i < N_dim; i++) {
    Nodal_Acceleration_t0.nM[i] = Nodal_Kinetics.nM[1 + i];
    Nodal_Acceleration_t1.nM[i] = Nodal_Kinetics.nM[1 + N_dim + i];
    Nodal_Velocity.nM[i] = Nodal_Kinetics.nM[1 + 2 * N_dim + i];
  }

  /* 1º iterate over the Gauss-Points */
  for (int i = 0; i < MPM_Mesh.NumGP; i++) {

    /* 2º Define element of the GP */
    GP_Element = nodal_set__Particles__(i, MPM_Mesh.ListNodes[i],
                                        MPM_Mesh.NumberNodes[i]);

    /* 3º Evaluate shape function in the GP i */
    N_GP = compute_N__MeshTools__(GP_Element, MPM_Mesh, FEM_Mesh);

    /* 4º Set to zero the GPs acceleration and velocity of the previous step */
    for (int k = 0; k < N_dim; k++) {
      MPM_Mesh.Phi.acc.nM[i][k] = 0.0;
    }

    /* 5º Iterate over the nodes of the element */
    for (int j = 0; j < GP_Element.NumberNodes; j++) {
      /* Node of the GP */
      GP_I = GP_Element.Connectivity[j];
      /* Evaluate the GP function in the node */
      N_I_GP = N_GP.nV[j];
      /* If this node has a null Value of the SHF continue */
      if (fabs(N_I_GP) <= TOL_zero) {
        continue;
      }
      /* Update GP cuantities with nodal values */
      for (int k = 0; k < N_dim; k++) {
        /* Get nodal values
           Nodal_Kinetics = {m, a0, a1, v}
        */
        /* Get the GP accelerations */
        MPM_Mesh.Phi.acc.nM[i][k] += N_I_GP * Nodal_Acceleration_t1.nM[k][GP_I];
      }
    }

    /* 5º Iterate over the nodes of the element */
    for (int j = 0; j < GP_Element.NumberNodes; j++) {
      /* Node of the GP */
      GP_I = GP_Element.Connectivity[j];
      /* Evaluate the GP function in the node */
      N_I_GP = N_GP.nV[j];
      /* If this node has a null Value of the SHF continue */
      if (fabs(N_I_GP) <= TOL_zero) {
        continue;
      }
      /* Update GP cuantities with nodal values */
      for (int k = 0; k < N_dim; k++) {
        /* Get nodal values
           Nodal_Kinetics = {m, a0, a1, v}
         */
        /* Update the GP velocities */
        MPM_Mesh.Phi.vel.nM[i][k] +=
            N_I_GP *
            ((1 - gamma) * Nodal_Acceleration_t0.nM[k][GP_I] +
             gamma * Nodal_Acceleration_t1.nM[k][GP_I]) *
            DeltaTimeStep;
      }
    }

    /* 5º Iterate over the nodes of the element */
    for (int j = 0; j < GP_Element.NumberNodes; j++) {
      /* Node of the GP */
      GP_I = GP_Element.Connectivity[j];
      /* Evaluate the GP function in the node */
      N_I_GP = N_GP.nV[j];
      /* If this node has a null Value of the SHF continue */
      if (fabs(N_I_GP) <= TOL_zero) {
        continue;
      }
      /* Update GP cuantities with nodal values */
      for (int k = 0; k < N_dim; k++) {
        /* Get nodal values
           Nodal_Kinetics = {m, a0, a1, v}
         */
        /* Update the GP position I */
        MPM_Mesh.Phi.x_GC.nM[i][k] +=
            N_I_GP * Nodal_Velocity.nM[k][GP_I] * DeltaTimeStep +
            N_I_GP *
                ((0.5 - beta) * Nodal_Acceleration_t0.nM[k][GP_I] +
                 beta * Nodal_Acceleration_t1.nM[k][GP_I]) *
                DeltaTimeStep * DeltaTimeStep;
      }
    }

    /* 6º Free memory */
    free(GP_Element.Connectivity);
    free__MatrixLib__(N_GP);
  }

  /* Update the grid nodal variales */
  for (int i = 0; i < N_Nodes; i++) {
    for (int j = 0; j < N_dim; j++) {
      Nodal_Acceleration_t0.nM[j][i] = Nodal_Acceleration_t1.nM[j][i];
      Nodal_Acceleration_t1.nM[j][i] = 0.0;
    }
  }

  /* Free tables */
  free(Nodal_Acceleration_t0.nM);
  free(Nodal_Acceleration_t1.nM);
  free(Nodal_Velocity.nM);
}

/*************************************************************/

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
//    free__TensorLib__(Rate_Strain_p);

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
  Load *F = MPM_Mesh.F;
  Element Nodes_p;        /* Element for each Gauss-Point */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pI;
  Tensor t = alloc__TensorLib__(1); /* Body forces vector */
  double m_p;                       /* Mass of the Gauss-Point */
  double rho_p;                     /* Density of the Gauss-Point */
  double V_p;                       /* Volumen of the Gauss-Point */
  double thickness_p;               /* Thickness of the Gauss-Point */
  int NumContactForces = MPM_Mesh.NumNeumannBC;
  int NumNodesLoad;
  int p;
  int Ip;
  int Nn; /* Number of nodes of each Gauss-Point */

  for (int i = 0; i < NumContactForces; i++) {

    NumNodesLoad = F[i].NumNodes;

    for (int j = 0; j < NumNodesLoad; j++) {

      /* Get the index of the Gauss-Point */
      p = F[i].Nodes[j];

      /* Get the value of the mass */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      /* Get the value of the density */
      rho_p = MPM_Mesh.Phi.mass.nV[p];

      /* Get the value of the volum */
      V_p = m_p / rho_p;

      /* Get the thickness of the material point */
      thickness_p = Thickness_Plain_Stress;

      /* Define element for each GP */
      Nn = MPM_Mesh.NumberNodes[p];
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], Nn);

      /* Compute shape functions */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Fill vector of body forces */
      for (int k = 0; k < Ndim; k++) {
        if (F[i].Dir[k * NumTimeStep + TimeStep]) {
          t.n[k] = F[i].Value[k].Fx[TimeStep];
        }
      }

      /* Get the node of the mesh for the contribution */
      for (int A = 0; A < Nn; A++) {

        /* Node for the contribution */
        Ip = Nodes_p.Connectivity[A];

        /* Pass the value of the nodal shape function to a scalar */
        ShapeFunction_pI = ShapeFunction_p.nV[A];

        /* Compute Contact forces */
        for (int k = 0; k < Ndim; k++) {
          F_I.nM[Ip][k] += ShapeFunction_pI * (t.n[k] / thickness_p) * V_p;
        }
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(ShapeFunction_p);
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
