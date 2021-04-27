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
double epsilon_Mass_Matrix; 
double beta_Newmark_beta;   
double gamma_Newmark_beta;
double TOL_Newmark_beta;
Event * Out_nodal_path_csv;
Event * Out_particles_path_csv;
int Number_Out_nodal_path_csv;
int Number_Out_particles_path_csv;

/*
  Define local global variable for the relative error
*/
double Error0;

typedef struct
{
  double alpha_1;
  double alpha_2;
  double alpha_3;
  double alpha_4;
  double alpha_5;
  double alpha_6;
} Newmark_parameters;

typedef struct
{

  Matrix value;
  Matrix d_value_dt;
  Matrix d2_value_dt2;

} Nodal_Field;
  
/*
  Auxiliar functions and variables
*/

static Newmark_parameters compute_Newmark_parameters(double, double, double);
static Matrix compute_Nodal_Effective_Mass(Particle, Mesh, Mask, double);
static  void  compute_Gravity_field(Mask, Particle, int);
static Nodal_Field compute_Nodal_Field(Matrix, Particle, Mesh, Mask);
static Nodal_Field initialise_Nodal_Increments(Mask, Mesh, int);
static  void  update_Local_State(Nodal_Field, Mask, Particle, Mesh);
static Matrix compute_Residual(Nodal_Field,Nodal_Field,Mask,Particle,Mesh,int);
static  void  compute_Inertial_Forces_Mixture(Nodal_Field,Matrix,Mask,Particle,Mesh);
static  void  compute_Internal_Forces_Mixture(Matrix,Mask,Particle,Mesh);
static Tensor compute_total_first_Piola_Kirchhoff_stress(Tensor,double,Tensor);
static  void  compute_Contact_Forces_Mixture(Matrix,Mask,Particle,Mesh,int);
static  void  compute_Compresibility_Mass_Balance(Nodal_Field, Matrix, Mask, Particle, Mesh);
static  void  compute_Jacobian_Rate_Mass_Balance(Matrix,Mask,Particle,Mesh);
static  void  compute_Permeability_Mass_Balance(Nodal_Field, Nodal_Field,Matrix,Mask,Particle,Mesh);
static Tensor compute_Pore_water_pressure_gradient_n1(Matrix,Matrix,Matrix);
static  void  compute_Permeability_Inertial_Forces_Fluid(Nodal_Field,Matrix,Mask,Particle,Mesh);
static  bool  check_convergence(Matrix,double,int,int,int);
static Matrix assemble_Nodal_Tangent_Stiffness(Nodal_Field,Nodal_Field,Mask,Particle,Mesh,Newmark_parameters);
static Tensor compute_stiffness_density(Tensor, Tensor, Tensor, double, Material);
static Tensor compute_H_AB(Tensor,Tensor);
static Tensor compute_L_AB(double,Tensor,Tensor,Tensor);
static  void  solve_non_reducted_system(Nodal_Field,Matrix,Matrix);
static  void  solve_reducted_system(Nodal_Field,Matrix,Matrix,Mask);
static  void  update_Newmark_Nodal_Increments(Nodal_Field,Nodal_Field,Newmark_parameters);
static  void  update_Particles(Nodal_Field,Particle,Mesh,Mask);
static  void  output_selector(Nodal_Field,Nodal_Field,Particle,Mesh,Mask,int,int);
static  char  Error_message[MAXW];
static  void  standard_error();

/**************************************************************/

void upw_Newmark_beta_Finite_Strains(Mesh FEM_Mesh, Particle MPM_Mesh, int InitialStep)
{

  /*
    Integer variables 
  */
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nactivenodes;
  /*
    Auxiliar variables for the solver
  */
  Matrix Effective_Mass;
  Nodal_Field D_upw;
  Nodal_Field upw_n;
  Matrix Residual;
  Matrix Tangent_Stiffness;
  Mask ActiveNodes;
  Mask Free_and_Restricted_Dofs;
  double TOL = TOL_Newmark_beta;
  double epsilon = epsilon_Mass_Matrix;
  double beta = beta_Newmark_beta;
  double gamma = gamma_Newmark_beta;

  /*
    Alpha parameters for the Newmark-beta
  */
  Newmark_parameters Params;
  
  double DeltaTimeStep;
  bool Convergence;
  int Iter = 0;
  int MaxIter = 100;

  /*
    Time step is defined at the init of the simulation throught the
    CFL condition. Notice that for this kind of solver, CFL confition is
    not required to be satisfied. The only purpose of it is to use the existing
    software interfase.
  */
  DeltaTimeStep = DeltaT_Coussy__SolversLib__(MPM_Mesh, FEM_Mesh.DeltaX, 1.0); 
  Params = compute_Newmark_parameters(beta, gamma, DeltaTimeStep);

  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
  {

    print_Status("*************************************************",TimeStep);
    print_step(TimeStep,DeltaTimeStep);

    print_Status("*************************************************",TimeStep);
    print_Status("First step : Generate Mask ... WORKING",TimeStep);

    ActiveNodes = generate_NodalMask__MeshTools__(FEM_Mesh);
    Nactivenodes = ActiveNodes.Nactivenodes;
    Free_and_Restricted_Dofs = generate_Mask_for_static_condensation__MeshTools__(ActiveNodes,FEM_Mesh);

    print_Status("DONE !!!",TimeStep);

    print_Status("*************************************************",TimeStep);
    print_Status("Second step : Compute effective mass ... WORKING",TimeStep);

    Effective_Mass = compute_Nodal_Effective_Mass(MPM_Mesh,FEM_Mesh,ActiveNodes,epsilon);

    print_Status("DONE !!!",TimeStep);

    print_Status("*************************************************",TimeStep);
    print_Status("Third step : Compute nodal kinetics ... WORKING",TimeStep);
    
    compute_Gravity_field(ActiveNodes, MPM_Mesh, TimeStep);
    upw_n = compute_Nodal_Field(Effective_Mass,MPM_Mesh,FEM_Mesh,ActiveNodes);
    D_upw = initialise_Nodal_Increments(ActiveNodes, FEM_Mesh, TimeStep);

    print_Status("DONE !!!",TimeStep);  

    print_Status("*************************************************",TimeStep);
    print_Status("Four step : Compute equilibrium ... WORKING",TimeStep);
    

    Convergence = false;
    Iter = 0;
    while(Convergence == false)
    {

      update_Local_State(D_upw,ActiveNodes,MPM_Mesh,FEM_Mesh);

      Residual = compute_Residual(upw_n,D_upw,ActiveNodes,MPM_Mesh,FEM_Mesh,TimeStep);

      Convergence = check_convergence(Residual,TOL,Iter,MaxIter,TimeStep);

      if(Convergence == false)
      {

        Tangent_Stiffness = assemble_Nodal_Tangent_Stiffness(upw_n,D_upw,ActiveNodes,MPM_Mesh,FEM_Mesh,Params);

        if((Free_and_Restricted_Dofs.Nactivenodes - Ndim*Nactivenodes) == 0)
        {
          solve_non_reducted_system(D_upw,Tangent_Stiffness,Residual);
        }
        else
        {
          solve_reducted_system(D_upw,Tangent_Stiffness,Residual,Free_and_Restricted_Dofs);
        }

        update_Newmark_Nodal_Increments(D_upw,upw_n,Params);

        Iter++;

        free__MatrixLib__(Residual);
        free__MatrixLib__(Tangent_Stiffness);
      }
    }

    print_Status("DONE !!!",TimeStep);

    print_Status("*************************************************",TimeStep);
    print_Status("Seven step : Update particles lagrangian ... WORKING",TimeStep);

    update_Particles(D_upw,MPM_Mesh,FEM_Mesh,ActiveNodes);

    local_search__Particles__(MPM_Mesh,FEM_Mesh);

    print_Status("DONE !!!",TimeStep);
      

    /*
      Outputs
    */
    output_selector(upw_n, D_upw, MPM_Mesh, FEM_Mesh, ActiveNodes, TimeStep, ResultsTimeStep);
  
    print_Status("*************************************************",TimeStep);
    print_Status("Eight step : Reset nodal values ... WORKING",TimeStep);
    
    /*
      Free memory.
    */
    free__MatrixLib__(Effective_Mass); 
    free__MatrixLib__(upw_n.value);
    free__MatrixLib__(upw_n.d_value_dt);
    free__MatrixLib__(upw_n.d2_value_dt2);
    free__MatrixLib__(D_upw.value);
    free__MatrixLib__(D_upw.d_value_dt);
    free__MatrixLib__(D_upw.d2_value_dt2);
    free__MatrixLib__(Residual);
    free(ActiveNodes.Nodes2Mask);
    free(Free_and_Restricted_Dofs.Nodes2Mask);

    print_Status("DONE !!!",TimeStep);

  }

}

/**************************************************************/

static Newmark_parameters compute_Newmark_parameters(
  double beta,
  double gamma,
  double DeltaTimeStep)
{
  Newmark_parameters Params;
  
  Params.alpha_1 = 1/(beta*DSQR(DeltaTimeStep));
  Params.alpha_2 = 1/(beta*DeltaTimeStep);
  Params.alpha_3 = (1-2*beta)/(2*beta);
  Params.alpha_4 = gamma/(beta*DeltaTimeStep);
  Params.alpha_5 = 1-gamma/beta;
  Params.alpha_6 = (1-gamma/(2*beta))*DeltaTimeStep;

  return Params;
}

/**************************************************************/

static Matrix compute_Nodal_Effective_Mass(
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  double epsilon)
/*
  This function computes the effective mass matrix as a convex combination
  of the lumped mass matrix and the consistent mass matrix. Later assemble
  a total mass matrix with the contribution of each degree of freedom.

  | M_eff |   0   |              | M_cons |   0    |          | M_lump |   0    |
  -----------------  = (1-eps) * -------------------  + eps * -------------------
  |    0  | M_eff |              |   0    | M_cons |        |   0    | M_lump |
*/
{

  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ndof = NumberDOF;
  int Np = MPM_Mesh.NumGP;
  int Order = Ndof*Nnodes_mask;
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
  Matrix Lumped_MassMatrix = allocZ__MatrixLib__(Order, Order);

  /*
    Iterate over the particles to get the nodal values 
  */
  for(int p = 0 ; p<Np ; p++) 
  {

    /*
      Define tributary nodes of the particle 
    */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    /* 
      Evaluate the shape function in the coordinates of the particle
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*  
      Get the mass of the particle 
    */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
    {

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
      m_A_p = m_p*ShapeFunction_pA;

      /* 
        Fill the Lumped mass matrix considering the number of dofs
      */
      for(int i = 0 ; i<Ndof ; i++)
      {
        Lumped_MassMatrix.nM[A_mask*Ndof + i][A_mask*Ndof + i] += m_A_p;
      }

      for(int B = 0 ; B<Nodes_p.NumberNodes ; B++)
      {       
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
        m_AB_p = m_p*ShapeFunction_pA*ShapeFunction_pB;

        /* 
          Fill the effective mass matrix considering the number of dofs
        */
        for(int i = 0 ; i<Ndof ; i++)
        {
          Effective_MassMatrix.nM[A_mask*Ndof+i][B_mask*Ndof+i] += m_AB_p;
        }

      }
    }

    /* 
      Free the value of the shape functions
    */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);      
  }

  /*
    At this point the effective mass matrix coincides with the consistent mass
    matrix. We can tune it by a convecx combination with the lumped mass matrix
  */
  for(int A = 0 ; A<Order ; A++)
  {
    for(int B = 0 ; B<Order ; B++)
    {    
      Effective_MassMatrix.nM[A][B] = (1-epsilon)*Effective_MassMatrix.nM[A][B] + (A == B)*epsilon*Lumped_MassMatrix.nM[A][B];
    }
  }

  /*
    Free lumped mass matrix.
  */
  free__MatrixLib__(Lumped_MassMatrix);

  /* 
    Add some usefulll info 
  */
  strcpy(Effective_MassMatrix.Info,"Effective-Mass-Matrix");

  return Effective_MassMatrix; 
}

/**************************************************************/

static void compute_Gravity_field(
  Mask ActiveNodes,
  Particle MPM_Mesh,
  int TimeStep)
/*

*/
{
  /* Define auxilar variables */
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int NumBodyForces = MPM_Mesh.NumberBodyForces;
  
  double m_p; /* Mass of the particle */
  Load * B = MPM_Mesh.B; /* List with the load cases */

  /*
    Initialise distance accelerations
  */
  for(int k = 0 ; k<Ndim ; k++)
  {
    MPM_Mesh.b.n[k] = 0.0;
  } 
  
  for(int i = 0 ; i<NumBodyForces ; i++)
  {

    /* Fill vector b of body acclerations */
    for(int k = 0 ; k<Ndim ; k++)
    {
      if(B[i].Dir[k])
      {
        if( (TimeStep < 0) || (TimeStep > B[i].Value[k].Num))
        {
          printf("%s : %s\n", "Error in compute_Gravity_field()","The time step is out of the curve !!");
          exit(EXIT_FAILURE);
        }
         
        MPM_Mesh.b.n[k] += B[i].Value[k].Fx[TimeStep];
      }
    }
      
  }

}

/**************************************************************/

static Nodal_Field compute_Nodal_Field(
  Matrix Effective_Mass,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes)
/*
  Call the LAPACK solver to compute simultanesly :

  -> The nodal Acceleration/Second time derivative of the pore water pressure.
  The operation is linearized and all the dof split the solution array in n components like :
  | M 0 0 |   |d2 (u.x) dt2|   |m * d2 (u.x) dt2|
  | 0 M 0 | * |d2 (u.y) dt2| = |m * d2 (u.y) dt2|
  | 0 0 M |   |d2 (p.w) dt2|   |m * d2 (p.w) dt2|

  -> The nodal Velocity/First time derivative of the pore water pressure.
  The operation is linearized and all the dof split the solution array in n components like :
  | M 0 0 |   |d (u.x) dt|   |m * d (u.x) dt|
  | 0 M 0 | * |d (u.y) dt| = |m * d (u.y) dt|
  | 0 0 M |   |d (p.w) dt|   |m * d (p.w) dt|

  -> The nodal displacement/pore water pressure.
  The operation is linearized and all the dof split the solution array in n components like :
  | M 0 0 |   |u.x|   |m * u.x|
  | 0 M 0 | * |u.y| = |m * u.y|
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
  Matrix ShapeFunction_p;

  /* Evaluation of the particle in the node */
  double ShapeFunction_pA;
  /* Mass of the particle */
  double m_p;
  /* Element for each particle */
  Element Nodes_p;

  /* Define and allocate the output vector */
  Nodal_Field upw;
  upw.value        = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  upw.d_value_dt   = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  upw.d2_value_dt2 = allocZ__MatrixLib__(Nnodes_mask,Ndof);


  /* Iterate over the particles to get the nodal values */
  for(int p = 0 ; p<Np ; p++)
  {
    /* 
      Define element of the particle
    */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);
  
    /* 
      Evaluate the shape function in the coordinates of the particle 
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* 
      Get the mass of the GP
    */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
    {

      /*
        Get the node with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Evaluate the GP function in the node
      */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /*
        Nodal displacement and pressure (and its rates)
      */
      for(int i = 0 ; i<Ndof ; i++)
      {
        if(i<Ndim)
        {
          upw.value.nM[A_mask][i]        += m_p*ShapeFunction_pA*MPM_Mesh.Phi.dis.nM[p][i];
          upw.d_value_dt.nM[A_mask][i]   += m_p*ShapeFunction_pA*MPM_Mesh.Phi.vel.nM[p][i];
          upw.d2_value_dt2.nM[A_mask][i] += m_p*ShapeFunction_pA*MPM_Mesh.Phi.acc.nM[p][i];
        } 
        else
        {
          upw.value.nM[A_mask][i]        += m_p*ShapeFunction_pA*MPM_Mesh.Phi.Pw.nV[p]; 
          upw.d_value_dt.nM[A_mask][i]   += m_p*ShapeFunction_pA*MPM_Mesh.Phi.d_Pw_dt.nV[p];
          upw.d2_value_dt2.nM[A_mask][i] += m_p*ShapeFunction_pA*MPM_Mesh.Phi.d2_Pw_dt2.nV[p];          
        } 
        
      }
    }

    /*
      Free the value of the shape functions
    */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);

  }

  /*
    Call the LAPACK solver to compute the accelerations and second derivative of the pore water pressure
  */
  int Order = Nnodes_mask*Ndof;
  int LDA   = Order;
  int LDB = Order;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO= 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;

  /*
    Compute the LU factorization for the mass matrix
  */
  dgetrf_(&Order,&Order,Effective_Mass.nV,&LDA,IPIV,&INFO);

  /*
    Solve for the acceleration and second derivative of the pore water pressure
  */
  dgetrs_(&TRANS,&Order,&NRHS,Effective_Mass.nV,&LDA,IPIV,upw.value.nV,&LDB,&INFO);
  dgetrs_(&TRANS,&Order,&NRHS,Effective_Mass.nV,&LDA,IPIV,upw.d_value_dt.nV,&LDB,&INFO);
  dgetrs_(&TRANS,&Order,&NRHS,Effective_Mass.nV,&LDA,IPIV,upw.d2_value_dt2.nV,&LDB,&INFO);
  free(IPIV);
 
  /*
    Add some usefulll info
  */
  strcpy(upw.value.Info,"Nodal-fields");
  strcpy(upw.d_value_dt.Info,"First-time-derivative-nodal-fields");
  strcpy(upw.d2_value_dt2.Info,"Second-time-derivative-nodal-fields");

  return upw;
}

/**************************************************************/

static Nodal_Field initialise_Nodal_Increments(
  Mask ActiveNodes,
  Mesh FEM_Mesh,
  int TimeStep)
/*
  Allocate the increments nodal values for the time integration scheme,
  and apply the boundary conditions over the nodes.
*/
{

  /* 1º Define auxilar variables */
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int NumNodesBound; /* Number of nodes of the bound */
  int Ndim = NumberDimensions; /* Number of dimensions */
  int Ndof = NumberDOF; /* Number of degree of freedom */
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Id_BCC_mask;

  /* 
    Allocate memory
  */
  Nodal_Field D_upw;
  D_upw.value        = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  D_upw.d_value_dt   = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  D_upw.d2_value_dt2 = allocZ__MatrixLib__(Nnodes_mask,Ndof);

  /* 
    Loop over the the boundaries to set boundary conditions
  */
  for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++)
    {

    /* 
      Get the number of nodes of this boundarie 
    */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

    /* 
      Get the number of dimensions where the BCC it is applied 
    */
    Ndof = FEM_Mesh.Bounds.BCC_i[i].Dim;

    for(int j = 0 ; j<NumNodesBound ; j++)
    {
      /* 
        Get the index of the node 
      */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

      /*
        The boundary condition is not affecting any active node,
        continue interating
      */
      if(Id_BCC_mask == -1)
      {
        continue;
      }

      /* 
        Loop over the dimensions of the boundary condition 
      */
      for(int k = 0 ; k<Ndof ; k++)
      {

        /* 
          Apply only if the direction is active (1) 
        */
        if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1)
        {
    
          /* 
            Check if the curve it is on time 
          */
          if( (TimeStep < 0) || (TimeStep > FEM_Mesh.Bounds.BCC_i[i].Value[k].Num))
          {
            printf("%s : %s \n","Error in initialise_Nodal_Increments()","The time step is out of the curve !!");
            exit(EXIT_FAILURE);
          }

          /* 
            Assign the boundary condition 
          */
          D_upw.value.nM[Id_BCC_mask][k] = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*(double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];                    
               
        }
      }
    }    
  }

  /*
    Add some usefull info
  */
  strcpy(D_upw.value.Info,"Nodal-fields");
  strcpy(D_upw.d_value_dt.Info,"First-time-derivative-nodal-fields");
  strcpy(D_upw.d2_value_dt2.Info,"Second-time-derivative-nodal-fields");

  return D_upw;
}

/**************************************************************/

static void update_Local_State(
  Nodal_Field D_upw, // Structure with the nodal values of the increment of displacementm and pore water pressure
  Mask ActiveNodes, // Information with the active nodes
  Particle MPM_Mesh, // Information related with the particles
  Mesh FEM_Mesh) // Information related with the nodes
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
  double J_n1_p; /* Jacobian of the deformation gradient at t = n + 1 */
  double dJ_dt_n1_p; /* Rate of the Jacobian of the deformation gradient at t = n + 1 */
  double K_f; /* Compressibility (fluid) */
  double rho_f_0; /* Initial density of the fluid */
  double phi_s_0; /* Initial volume fraction (solid) */
  double phi_f_0; /* Initial volume fraction (fluid) */
  Plastic_status Input_Plastic_Parameters; /* Input parameters for plasticity */
  Plastic_status Output_Plastic_Parameters; /* Output parameters for plasticity */
  Element Nodes_p; /* Element for each particle */
  Material MatProp_Soil_p; /* Variable with the material properties of the solid phase */
  Material MatProp_Water_p; /* Variable with the material properties of the fluid phase */
  Matrix gradient_p;
  Matrix ShapeFunction_p;
  Matrix Nodal_D_Displacement_p;
  Matrix Nodal_D_Velocity_p;
  Matrix Nodal_D_Pw_p;
  Tensor F_n_p; /* Deformation gradient of the soil skeleton (t = n) */
  Tensor F_n1_p; /* Deformation gradient of the soil skeleton (t = n + 1) */
  Tensor DF_p; /* Increment of the deformation gradient of the soil skeleton */
  Tensor dFdt_n_p; /* Rate of the deformation gradient of the soil skeleton (t = n) */
  Tensor dFdt_n1_p; /* Rate of the deformation gradient of the soil skeleton (t = n + 1) */
  Tensor dt_DF_p; /* Rate of the increment of the deformation gradient of the soil skeleton */
  Tensor F_plastic_p; /* Plastic deformation gradient of the soil skeleton */
  Tensor P_n1_p; /* First Piola-Kirchhoff stress tensor (t = n + 1) */
  double Pw_0; /* Cauchy pore water pressure at t = 0 */
  double Pw_n1; /* Cauchy pore water pressure at t = n + 1 */

  for(int p = 0 ; p<Np ; p++)
  {
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
      Get the nodal increment of displacement and velocity using the mask
    */
    Nodal_D_Displacement_p = get_U_set_field_upw__MeshTools__(D_upw.value, Nodes_p, ActiveNodes);
    Nodal_D_Velocity_p     = get_U_set_field_upw__MeshTools__(D_upw.d_value_dt, Nodes_p, ActiveNodes);

    /*
      Evaluate the shape function gradient in the coordinates of the particle
    */
    gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    
    /*
      Take the values of the deformation gradient and its rates from the previous step
    */
    F_n_p     = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    F_n1_p    = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    DF_p      = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.DF.nM[p],2);
    dFdt_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n.nM[p],2);
    dFdt_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n1.nM[p],2);
    dt_DF_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_DF.nM[p],2);
    
    /*
      Compute the increment of the deformation gradient and its rate
    */
    update_increment_Deformation_Gradient__Particles__(DF_p, Nodal_D_Displacement_p, gradient_p);
    update_rate_increment_Deformation_Gradient__Particles__(dt_DF_p, Nodal_D_Velocity_p, gradient_p);

    /*
      Update the deformation gradient in t = n + 1 with the information
      from t = n and the increment of deformation gradient.
    */  
    update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, DF_p);
    update_rate_Deformation_Gradient_n1__Particles__(dFdt_n1_p, dt_DF_p, F_n_p, DF_p, dFdt_n_p);

    /*
      Compute the Jacobian of the deformation gradient and its rate
    */
    J_n1_p     = I3__TensorLib__(F_n1_p);
    dJ_dt_n1_p = compute_Jacobian_Rate__Particles__(J_n1_p,F_n1_p,dFdt_n1_p);

    /*
      Compute the first Piola-Kirchhoff stress tensor (P).
    */
    P_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2);      

    
    if(strcmp(MatProp_Soil_p.Type,"Neo-Hookean-Wriggers") == 0)
    {
      P_n1_p = compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(P_n1_p, F_n1_p, J_n1_p, MatProp_Soil_p);
    }
    else
    {
      fprintf(stderr,"%s : %s %s %s \n","Error in update_Local_State()","The material",MatProp_Soil_p.Type,"has not been yet implemnented");
      exit(EXIT_FAILURE);
    }

    /*
      Get the Kirchhoff pore water pressure at t = n + 1
    */
    Nodal_D_Pw_p = get_Pw_set_field_upw__MeshTools__(D_upw.value, Nodes_p, ActiveNodes);
    MPM_Mesh.Phi.Pw_n1.nV[p] = MPM_Mesh.Phi.Pw.nV[p] + interpolate_scalar_magnitude__MeshTools__(Nodal_D_Pw_p, ShapeFunction_p);

    /*
      Update state parameters
    */
    Pw_0 = MPM_Mesh.Phi.Pw_0.nV[p]; /* Get the initial pressure */
    Pw_n1 = MPM_Mesh.Phi.Pw_n1.nV[p]/J_n1_p; /* From the Kirchhoff pressure compure the cauchy pore water pressure */

    MPM_Mesh.Phi.rho_f.nV[p] = rho_f_0*exp((Pw_n1-Pw_0)/K_f); /* Update the fluid density */

    MPM_Mesh.Phi.phi_s.nV[p] = phi_s_0/J_n1_p; /* Update the volume fraction of the solid phase */
    MPM_Mesh.Phi.phi_f.nV[p] = 1 - (1 - phi_f_0)/J_n1_p; /* Update the volume fraction of the fluid phase */

    MPM_Mesh.Phi.rho.nV[p] = MPM_Mesh.Phi.rho_s.nV[p]*MPM_Mesh.Phi.phi_s.nV[p] +
                             MPM_Mesh.Phi.rho_f.nV[p]*MPM_Mesh.Phi.phi_f.nV[p];  /* Update density of the mixture */

    MPM_Mesh.Phi.J.nV[p] = J_n1_p; /* Update soil skeleton jacobian */
    MPM_Mesh.Phi.dJ_dt.nV[p] = dJ_dt_n1_p; /* Update soil skeleton rate of jacobian */

    /*
      Free memory 
    */
    free__MatrixLib__(Nodal_D_Displacement_p);
    free__MatrixLib__(Nodal_D_Velocity_p);
    free__MatrixLib__(Nodal_D_Pw_p);
    free__MatrixLib__(gradient_p);
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);

  }
  
}

/**************************************************************/

static Matrix compute_Residual(
  Nodal_Field upw_n,
  Nodal_Field D_upw,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  int TimeStep)
{

  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;

  Matrix Residual = allocZ__MatrixLib__(Nnodes_mask,Ndof);

  compute_Inertial_Forces_Mixture(D_upw,Residual,ActiveNodes,MPM_Mesh,FEM_Mesh);

  compute_Internal_Forces_Mixture(Residual,ActiveNodes,MPM_Mesh,FEM_Mesh);

  compute_Contact_Forces_Mixture(Residual,ActiveNodes,MPM_Mesh,FEM_Mesh,TimeStep);

  compute_Compresibility_Mass_Balance(D_upw, Residual, ActiveNodes, MPM_Mesh, FEM_Mesh);

  compute_Jacobian_Rate_Mass_Balance(Residual, ActiveNodes, MPM_Mesh, FEM_Mesh);

  compute_Permeability_Mass_Balance(upw_n, D_upw, Residual, ActiveNodes, MPM_Mesh, FEM_Mesh);

  compute_Permeability_Inertial_Forces_Fluid(D_upw, Residual, ActiveNodes, MPM_Mesh, FEM_Mesh);
  
  return Residual;
}

/**************************************************************/

static void compute_Inertial_Forces_Mixture(
  Nodal_Field D_upw,
  Matrix Residual,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int idx_A_mask_i;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;
  /* Evaluation of the particle in the node */
  double ShapeFunction_pA;
  /* Increment of nodal acceleration */
  Matrix Nodal_D_Acceleration_p;
  /* Mass of the particle */
  double m_p;
  /* Acceleration vectors */
  Tensor a_n1_p;
  Tensor a_n_p;
  Tensor D_a_p;
  Tensor b_p;
  Tensor dyn_p;
  /* Element for each particle */
  Element Nodes_p;


  /* Iterate over the particles to get the nodal values */
  for(int p = 0 ; p<Np ; p++)
  {
    /* 
      Define element of the particle
    */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);
  
    /* 
      Evaluate the shape function in the coordinates of the particle 
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* 
      Get the mass of the GP
    */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /*
      Compute dynamic terms
    */
    Nodal_D_Acceleration_p = get_U_set_field_upw__MeshTools__(D_upw.d2_value_dt2, Nodes_p, ActiveNodes);
    D_a_p = interpolate_vectorial_magnitude__MeshTools__(Nodal_D_Acceleration_p, ShapeFunction_p);
    a_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.acc.nM[p],1);
    a_n1_p = addition__TensorLib__(a_n_p, D_a_p);
    b_p = MPM_Mesh.b;
    dyn_p = subtraction__TensorLib__(a_n1_p, b_p);

    for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
    {

      /*
        Get the node with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Evaluate the GP function in the node
      */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /*
        Nodal velocity and acceleration
      */
      for(int i = 0 ; i<Ndim ; i++)
      {
        Residual.nM[A_mask][i] += ShapeFunction_pA*m_p*dyn_p.n[i];
      }
    }

    /*
      Free memory
    */
    free__MatrixLib__(ShapeFunction_p);
    free__MatrixLib__(Nodal_D_Acceleration_p);
    free__TensorLib__(D_a_p);
    free__TensorLib__(a_n1_p);
    free__TensorLib__(dyn_p);
    free(Nodes_p.Connectivity);

  }

}

/**************************************************************/

static void compute_Internal_Forces_Mixture(
  Matrix Residual,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int NumNodes_p;

  Tensor P_p; /* Total First Piola-Kirchhoff Stress tensor */
  Tensor P_effective_p; /* Effective First Piola-Kirchhoff Stress tensor */
  double theta_n1_p; /* Kirchhoff pore fluid pressure t = n + 1 */
  Tensor InternalForcesDensity_Ap;
  Element Nodes_p; /* List of nodes for particle */
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pA;
  Tensor GRADIENT_pA;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor FT_n_p;
  double V0_p; /* Volume of the Gauss-Point */

  /*
    Loop in the particles 
  */
  for(int p = 0 ; p<Np ; p++)
  {

    /*
      Define nodes for each particle
    */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /*
      Compute shape function gradient in each node 
    */
    gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
        
    /*
      Take the values of the deformation gradient at t = n. And transpose it
    */
    F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    FT_n_p = transpose__TensorLib__(F_n_p);
    F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);

    /*
      Get the volume of the particle in the reference configuration 
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /*
      Get the first Piola-Kirchhoff stress tensor
    */
    P_effective_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);

    /*
      Get the Kirchhoff pore water pressure 
    */
    theta_n1_p = MPM_Mesh.Phi.Pw_n1.nV[p];

    /*
      Following Terzaghi's idea, the effective stress tensor in the reference configuration
      is computed:
    */
    P_p = compute_total_first_Piola_Kirchhoff_stress(P_effective_p,theta_n1_p,F_n1_p);

    for(int A = 0 ; A<NumNodes_p ; A++)
    {
      
      /*
        Compute the gradient in the reference configuration 
      */
      gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);
      GRADIENT_pA = vector_linear_mapping__TensorLib__(FT_n_p,gradient_pA);
      
      /*
        Compute the nodal forces of the particle 
      */
      InternalForcesDensity_Ap = vector_linear_mapping__TensorLib__(P_p, GRADIENT_pA);
      
      /*
        Get the node of the mesh for the contribution 
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];
      
      /*
        Asign the nodal forces contribution to the node 
      */
      for(int i = 0 ; i<Ndim ; i++)
      {
        Residual.nM[A_mask][i] += InternalForcesDensity_Ap.n[i]*V0_p;
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
    free__TensorLib__(FT_n_p);
    free__TensorLib__(P_p);
    free__MatrixLib__(gradient_p);
    free(Nodes_p.Connectivity);
  }
  
}

/**************************************************************/

static Tensor compute_total_first_Piola_Kirchhoff_stress(
  Tensor P_effective_p,
  double theta_p,
  Tensor F_n1_p)
/*
  This function returns : P = P' - theta*F^{-T}
*/
{
  int Ndim = NumberDimensions;

  Tensor P_p = alloc__TensorLib__(2);
  Tensor Fm1_n1_p = Inverse__TensorLib__(F_n1_p);

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      P_p.N[i][j] = P_effective_p.N[i][j] - theta_p*Fm1_n1_p.N[j][i];
    }
  }

  free__TensorLib__(Fm1_n1_p);

  return P_p;
}

/*********************************************************************/

static void compute_Contact_Forces_Mixture(
  Matrix Residual,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  int TimeStep)
{
  int Ndim = NumberDimensions;
  Load T_i;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pA;
  Tensor t = alloc__TensorLib__(1); /* Body forces vector */
  double V0_p; /* Volumen of the particle in the reference configuration */
  double thickness_p; /* Thickness of the particle */
  double A0_p; /* Area of the particle in the reference configuration */ 

  int NumContactForces = MPM_Mesh.Neumann_Contours.NumBounds;
  int NumNodesLoad;
  int p;
  int Ap;
  int A_mask;
  int NumNodes_p; /* Number of nodes of each particle */

  for(int i = 0 ; i<NumContactForces; i++)
  {

    /*
      Read load i
    */
    T_i = MPM_Mesh.Neumann_Contours.BCC_i[i];

    NumNodesLoad = T_i.NumNodes;
      
    for(int j = 0 ; j<NumNodesLoad ; j++)
    {

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
      A0_p = V0_p/thickness_p;

      /*
        Define tributary nodes of the particle 
      */
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],NumNodes_p);

      /* 
        Evaluate the shape function in the coordinates of the particle
      */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
        Fill vector of contact forces
      */
      for(int k = 0 ; k<Ndim ; k++)
      {
        if(T_i.Dir[k])
        {
          if( (TimeStep < 0) || (TimeStep > T_i.Value[k].Num))
          {
            sprintf(Error_message,"%s : %s",
              "Error in compute_Contact_Forces_Mixture()",
              "The time step is out of the curve !!");
            standard_error();
          }
          t.n[k] = T_i.Value[k].Fx[TimeStep];
        }
      }

      /*
        Get the node of the mesh for the contribution
      */
      for(int A = 0 ; A<NumNodes_p ; A++)
      {

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
        for(int k = 0 ; k<Ndim ; k++)
        {
          Residual.nM[A_mask][k] -= ShapeFunction_pA*t.n[k]*A0_p;
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

static void compute_Compresibility_Mass_Balance(
  Nodal_Field D_upw,
  Matrix Residual,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
{
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ndim = NumberDimensions;
  int Mixture_idx;
  int Material_Water_idx;
  int Ap;
  int A_mask;
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Value of the shape-function */
  Material MatProp_Water_p; /* Variable with the material properties of the fluid phase */
  Matrix D_d_Pw_dt_n_I;
  double ShapeFunction_pA; /* Evaluation of the particle in the node A */
  double V0_p; /* Mass of the particle (mixture) */
  double rho_f_p; /* Material density of the fluid */
  double phi_f_p; /* Volume fractions of the fluid */
  double relative_rho_f_p; /* Relative density of the fluid */
  double K_f; /* Compressibility (fluid) */
  double compressibility_density_f_p; /* Compressibilidy density fo the fluid */
  double compressibility_density_f_A_p; /* Nodal contribution A of the particle p */
  double d_Pw_dt_n_p; /* Rate of pore water pressure (t = n) in particle p */
  double d_Pw_dt_n1_p; /* Rate of pore water pressure (t = n + 1) in particle p */

  for(int p = 0 ; p<Np ; p++)
  {

    /* 
      Define tributary nodes of the particle 
    */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    /* 
      Evaluate the shape function in the coordinates of the particle 
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
   
    /* 
      Get the initial volume of the particle 
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /*
      Get the current intrinsic density, volume fraction 
      and compressibility for each material point (fluid).
      Compute relative density.
    */
    rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];
    phi_f_p = MPM_Mesh.Phi.phi_f.nV[p];
    relative_rho_f_p = phi_f_p*rho_f_p;

    /*
      Load intrinsic properties for the fluid phase to 
      get the fluid compressibility
    */
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Water_idx = Soil_Water_Mixtures[Mixture_idx].Water_Idx;
    MatProp_Water_p = MPM_Mesh.Mat[Material_Water_idx];
    K_f = MatProp_Water_p.Compressibility;

    /*
      Compute the value of the rate of pore water pressure
      at t = n + 1, in each particle
    */
    d_Pw_dt_n_p = MPM_Mesh.Phi.d_Pw_dt.nV[p];
    D_d_Pw_dt_n_I = get_Pw_set_field_upw__MeshTools__(D_upw.d_value_dt, Nodes_p, ActiveNodes);
    d_Pw_dt_n1_p = d_Pw_dt_n_p + interpolate_scalar_magnitude__MeshTools__(D_d_Pw_dt_n_I, ShapeFunction_p);


    /* 
      Compute the compressibility density 
    */
    compressibility_density_f_p = relative_rho_f_p/K_f;

      
    for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
    {
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
      compressibility_density_f_A_p = compressibility_density_f_p*ShapeFunction_pA;

      /* 
        Add the compressubility contribution
      */
      Residual.nM[A_mask][Ndim] += compressibility_density_f_A_p*d_Pw_dt_n1_p*V0_p;      
    }

    /* Free the value of the shape functions */
    free__MatrixLib__(ShapeFunction_p);
    free__MatrixLib__(D_d_Pw_dt_n_I);
    free(Nodes_p.Connectivity);      

  }

}

/**************************************************************/

static void compute_Jacobian_Rate_Mass_Balance(
  Matrix Residual,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int NumNodes_p; /* Number of tributary nodes of p */
  int A_mask; /* Index of the node where we apply the body force */
  int Ap; /* Tributary node A of particle p */

  Element Nodes_p; /* Element for each particle p */
  Matrix ShapeFunction_p; /* Matrix with the value of the shape function in the particle p */ 
  double ShapeFunction_pA; /* Value of the shape funtion in node A for the particle p */
  double rho_f_p; /* Material density of the fluid phase for particle p */
  double V0_p; /* Volume of the particle */
  double dJ_dt_n1_p; /* Rate of the jacobian */

  for(int p = 0 ; p<Np ; p++)
  {
      
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

    
    for(int A = 0 ; A<NumNodes_p ; A++)
    {
      
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
      Residual.nM[A_mask][Ndim] += ShapeFunction_pA*rho_f_p*dJ_dt_n1_p*V0_p;

    }
        
    /* 
      Free memory 
    */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

}

/**************************************************************/

static void compute_Permeability_Mass_Balance(
  Nodal_Field upw_n,
  Nodal_Field D_upw,
  Matrix Residual,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
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

  Material MatProp_Soil_p; /* Variable with the material propertie of solid phase for each particle */
  Element Nodes_p; /* Element for each particle */
  Matrix Nodal_Pw_n_p;
  Matrix Nodal_D_Pw_p;
  Matrix gradient_p;  /* Shape functions gradients */
  Tensor gradient_pA; /* Shape functions gradients (Node A), def config */
  Tensor GRADIENT_pA; /* Shape functions gradients (Node A), ref config */
  Tensor F_n_p; /* Deformation gradient t = n */
  Tensor F_n1_p; /* Deformation gradient t = n + 1 */
  Tensor FT_n_p; /* Transpose of the deformation gradient t = n */
  Tensor Fm1_n1_p; /* */
  Tensor k_p; /* Spatial permebility tensor */
  Tensor Fk_p; /* Product of the inverse of the defomration gradient and the permeability tensor */
  Tensor gradPw_n1; 
  Tensor Fk__x__gradPw_n1;
  double GRADIENT_pA__x__Fk__x__gradPw_n1;
  double g = -9.81;
  double V0_p; /* Volume of the particle at the reference configuration */

  /*
    Iterate over the particles to get the nodal values 
  */
  for(int p = 0 ; p<Np ; p++)
    {

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
        Later compute the transpose of the deformation gradient and the determinant.
      */
      F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
      FT_n_p = transpose__TensorLib__(F_n_p);
      Fm1_n1_p = Inverse__TensorLib__(F_n1_p);

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
      Fk_p = matrix_product__TensorLib__(Fm1_n1_p,k_p);

      /* 
        Compute particle pore water pressure gradient
      */
      Nodal_Pw_n_p = get_Pw_set_field_upw__MeshTools__(upw_n.value, Nodes_p, ActiveNodes);
      Nodal_D_Pw_p = get_Pw_set_field_upw__MeshTools__(D_upw.value, Nodes_p, ActiveNodes);
      gradPw_n1 = compute_Pore_water_pressure_gradient_n1(Nodal_Pw_n_p,Nodal_D_Pw_p,gradient_p);

      /*
        Intermediate result 1
      */
      Fk__x__gradPw_n1 = vector_linear_mapping__TensorLib__(Fk_p,gradPw_n1);


      for(int A = 0 ; A<NumNodes_p ; A++)
      {

        /* 
          Get the node in the mass matrix with the mask
        */
        Ap = Nodes_p.Connectivity[A];
        A_mask = ActiveNodes.Nodes2Mask[Ap];
    
        /*
          Compute the gradient in the reference configuration for the node A
        */
        gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);
        GRADIENT_pA = vector_linear_mapping__TensorLib__(FT_n_p,gradient_pA);
  
        /*
          Intermediate result 2
        */
        GRADIENT_pA__x__Fk__x__gradPw_n1 = inner_product__TensorLib__(GRADIENT_pA,Fk__x__gradPw_n1);

        /* 
          Compute nodal contribution to the mass conservation
        */
        Residual.nM[A_mask][Ndim] -= (1/g)*GRADIENT_pA__x__Fk__x__gradPw_n1*V0_p; 

        free__TensorLib__(GRADIENT_pA);
      }

      /* 
        Free the value of the shape functions 
      */
      free__MatrixLib__(gradient_p);
      free__MatrixLib__(Nodal_Pw_n_p);
      free__MatrixLib__(Nodal_D_Pw_p);
      free__TensorLib__(gradPw_n1);
      free__TensorLib__(FT_n_p);
      free__TensorLib__(Fm1_n1_p);
      free__TensorLib__(Fk_p);
      free__TensorLib__(Fk__x__gradPw_n1);
      free(Nodes_p.Connectivity);      
    }

}

/**************************************************************/

static Tensor compute_Pore_water_pressure_gradient_n1(
  Matrix Nodal_Pore_water_pressure_p,
  Matrix Nodal_D_Pore_water_pressure_p,
  Matrix gradient_p)
/*

*/
 {

  /* Variable definition */
  int Ndim = NumberDimensions;
  int Nnodes_p = Nodal_Pore_water_pressure_p.N_rows; 
  double Pore_water_pressure_pI;
  double D_Pore_water_pressure_pI;
  Tensor gradient_I;
  Tensor gradPw_n1 = alloc__TensorLib__(1);

  for(int I = 0 ; I<Nnodes_p ; I++)
  {

    /*
      Assign from matrix
    */
    Pore_water_pressure_pI = Nodal_Pore_water_pressure_p.nV[I];
    D_Pore_water_pressure_pI = Nodal_D_Pore_water_pressure_p.nV[I];
    gradient_I = memory_to_tensor__TensorLib__(gradient_p.nM[I], 1);

    /*
      Compute nodal contribution
    */
    for(int i = 0 ; i<Ndim ; i++)
     {
      gradPw_n1.n[i] += (Pore_water_pressure_pI + D_Pore_water_pressure_pI)*gradient_I.n[i];
     } 

  }

  return gradPw_n1;
 } 

/**************************************************************/

static void compute_Permeability_Inertial_Forces_Fluid(
  Nodal_Field D_upw,
  Matrix Residual,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
/*

*/
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int NumNodes_p; /* Number of tributary nodes of p */
  int A_mask; /* Index of the node where we apply the body force */
  int Ap; /* Tributary node A of particle p */
  int B_mask; /* Index of the node where we apply the body force */
  int Bp; /* Tributary node B of particle p */
  int Mixture_idx; /* Index for the material point mixture parameters */

  Element Nodes_p; /* Element for each particle */
  Matrix Nodal_D_Acceleration_p;
  Matrix ShapeFunction_p;
  Matrix gradient_p;  /* Shape functions gradients */
  Tensor gradient_pA; /* Shape functions gradients (Node A), def config */
  Tensor GRADIENT_pA; /* Shape functions gradients (Node A), ref config */
  Tensor F_n_p; /* Deformation gradient t = n */
  Tensor F_n1_p; /* Deformation gradient t = n + 1 */
  Tensor FT_n_p; /* Transpose of the deformation gradient t = n */
  Tensor Fm1_n1_p; /* Inverse of the deformation gradient t = n + 1 */
  Tensor k_p; /* Spatial permebility tensor */
  Tensor Fk_p; /* One side pull-back of the permeability tensor */
  Tensor a_n_p; /* Particle acceleration */
  Tensor a_n1_p; /* Particle acceleration */
  Tensor D_a_p; /* Particle acceleration */
  Tensor b_p; /* External acceleration */
  Tensor dyn_p; /* Total acceleration of the particle */
  Tensor Fk_dyn_p; /* Axiliar tensor for intermediate result */
  double GRADIENT_pA__x__Fk_dyn_p; /* Auxiliar scalar for intermediate result */
  double rho_f_p; /* Intrinsic or material density (fluid phase) */
  double g = -9.81; /* Gravity constant */
  double V0_p; /* Inital volume of the particle p */
  double J_n1_p; /* Determianant of the soil skeleton deformation gradient at t = n + 1 */

  for(int p = 0 ; p<Np ; p++)
  {
      
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
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
        
    /*
      Take the value of the deformation gradient at t = n + 1, and t = n. 
      Get intermediate results
    */
    F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    FT_n_p = transpose__TensorLib__(F_n_p);
    Fm1_n1_p = Inverse__TensorLib__(F_n1_p);
    J_n1_p = MPM_Mesh.Phi.J.nV[p];

    /*
      Get the permeability tensor
    */
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    k_p = Soil_Water_Mixtures[Mixture_idx].Permeability;

    /*
      Compute intermediate result
    */
    Fk_p = matrix_product__TensorLib__(Fm1_n1_p,k_p);

    /*
      Compute total particle acceleration at n+1 using nodal variables
    */
    Nodal_D_Acceleration_p = get_U_set_field_upw__MeshTools__(D_upw.d2_value_dt2, Nodes_p, ActiveNodes);
    D_a_p = interpolate_vectorial_magnitude__MeshTools__(Nodal_D_Acceleration_p, ShapeFunction_p);
    a_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.acc.nM[p],1);
    a_n1_p = addition__TensorLib__(a_n_p, D_a_p);
    b_p = MPM_Mesh.b;
    dyn_p = subtraction__TensorLib__(a_n1_p, b_p);

    /*
      Compute intermediate result
    */
    Fk_dyn_p = vector_linear_mapping__TensorLib__(Fk_p,dyn_p);

    for(int A = 0 ; A<NumNodes_p ; A++)
    {

      /*
        Get the gradient of the shape function and multiply it by the permeability tensor
      */
      gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A],1);
      GRADIENT_pA = vector_linear_mapping__TensorLib__(FT_n_p,gradient_pA);
     
      /*
        Get the node of the mesh for the contribution 
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Compute intermediate result
      */
      GRADIENT_pA__x__Fk_dyn_p = inner_product__TensorLib__(GRADIENT_pA, Fk_dyn_p);

      /*
        Add nodal contributions
      */
      Residual.nM[A_mask][Ndim] -= (J_n1_p*rho_f_p/g)*GRADIENT_pA__x__Fk_dyn_p*V0_p;

      /*
        Free some auxiliar resutls
      */  
      free__TensorLib__(GRADIENT_pA);

    }
        
    /* 
      Free memory 
    */
    free__MatrixLib__(Nodal_D_Acceleration_p);
    free__TensorLib__(FT_n_p);
    free__TensorLib__(Fm1_n1_p);
    free__TensorLib__(Fk_p);
    free__TensorLib__(D_a_p);
    free__TensorLib__(a_n1_p);
    free__TensorLib__(dyn_p);
    free__TensorLib__(Fk_dyn_p);
    free__MatrixLib__(gradient_p);
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);
  }

}

/**************************************************************/

static bool check_convergence(
  Matrix Residual,
  double TOL,
  int Iter,
  int MaxIter,
  int Step)
{
  bool convergence;
  int Ndof = NumberDOF;
  int Nnodes_mask = Residual.N_rows;
  int Total_dof = Ndof*Nnodes_mask;
  double Error = 0;
  double Error_relative = 0;

  if(Iter > MaxIter)
  {
      fprintf(stderr,"%s : %s !!! \n","Error in upw_Newmark_beta_Finite_Strains()",
        "Convergence not reached in the maximum number of iterations");
      exit(EXIT_FAILURE);
    }
  else
    {
      /*
        Compute absolute error 
      */
      for(int A = 0 ; A<Total_dof ; A++)
      {
        Error += DSQR(Residual.nV[A]);
      }
      Error = pow(Error,0.5);

      /*
        Compute relative error
      */
      if(Iter == 0)
      {
        Error0 = Error;
        Error_relative = Error/Error0;
        printf("Error iter %i : %1.4e ; %1.4e \n",Iter,Error,Error_relative);
      }
      else
      {
        Error_relative = Error/Error0;
        printf("Error iter %i : %1.4e ; %1.4e \n",Iter,Error,Error_relative);
      }
      
      /*
        Check convergence using the relative error
      */
      if(Error_relative > TOL)
      {
        return false;
      }
      else
      {
        print_convergence_stats(Step, Iter, Error, Error_relative);
        return true;
      }
    }
}


/**************************************************************/

static Matrix assemble_Nodal_Tangent_Stiffness(
  Nodal_Field upw_n,
  Nodal_Field D_upw,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Newmark_parameters Params)
/*
  This function computes the tangent stiffness matrix as a combination
  of the geometrical stiffness matrix and the material stiffness matrix. 

  | K_xx | K_xy |   | K_geom_xx |    0      |     | K_mat_xx | K_mat_xy |
  ---------------- = -------------------------  +  ----------------------- + ...
  | K_yx | K_yy |   |     0     | K_geom_yy |     | K_mat_yx | K_mat_y  |
*/
{

  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Ndof*Nnodes_mask;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int Bp;
  int B_mask;
  int NumNodes_p;
  int Mixture_idx;
  int Material_Soil_idx;
  int Material_Water_idx;

  Element Nodes_p; /* List of nodes for particle */

  Matrix ShapeFunction_p; /* Shape functions */
  Matrix gradient_p; /* Shape functions gradients */
  Matrix Nodal_D_Acceleration_p; // Nodal values of the acceleration increments
  Matrix Nodal_Pw_n_p; // Nodal values of the pore-water pressure
  Matrix Nodal_D_Pw_p; // Nodal values of the pore-water pressure increments

  Material MatProp_Soil_p; // structure with the material properties of the soil
  Material MatProp_Water_p; // structure with the material properties of the water

  Tensor k_p; // Particle permeability
  Tensor D_a_p; // Particle increment of acceleration
  Tensor a_n1_p; // Particle acceleration at t = n + 1 
  Tensor a_n_p; // Particle acceleration at t = n 
  Tensor gradPw_n1; // Particle pore water pressure
  Tensor b_p; // External accelerations
  Tensor dyn_p; // Particle acceleration minus externa accelerations 
  Tensor kdyn_p; // Product of the permeability tensor and the total acceleration
  Tensor kgradPw_p; // Product of the permeability tensor and the gradient of the pore water pressure
  Tensor F_n_p; // Particle deformation gradient at t = n 
  Tensor FT_n_p; // Transpose of the particle deformation gradient at t = n 
  Tensor F_n1_p; // Particle deformation gradient at t = n + 1
  Tensor FT_n1_p; // Transpose of the particle deformation gradient at t = n + 1
  Tensor Fm1_n1_p; // Inverse of the particle deformation gradient at t = n + 1
  Tensor FmT_n1_p; // Inverse of the transposed particle deformation gradient at t = n + 1
  Tensor dFdt_n1_p; // Rate of the particle deformation gradient at t = n + 1
  double FmT__dd__dFdt_n1_p; // FmT : dFdt_n1
  Tensor dFdt__x__Fm1_n1_p; // dFdt x Fm1_n1
  Tensor gradient_pA, GRADIENT_pA; // Shape function gradient in the deformed/reference configurations (Node A)
  Tensor gradient_pB, GRADIENT_pB; // Shape function gradient in the deformed/reference configurations (Node B)
  Tensor FmTGRADIENT_pA; // Covariant push-forward of the GRADIENT_pA vector
  Tensor FmTGRADIENT_pB; // Covariant push-forward of the GRADIENT_pB vector
  Tensor dyn__o__FmTGRADIENT_pB;
  Tensor FmTGRADIENT__o__FmTGRADIENT_pAB; // Diadic product between FmTGRADIENT_A and FmTGRADIENT_B
  Tensor FmTGRADIENT__o__FmTGRADIENT_pBA; 
  Tensor FmTGRADIENT__o__FmTGRADIENT_AB__d__kdyn_p; // Auxiliar tensors
  Tensor FmTGRADIENT__o__FmTGRADIENT_BA__d__kdyn_p;
  Tensor FmTGRADIENT__o__FmTGRADIENT_AB__d__kgradPw_p; // Auxiliar tensors
  Tensor FmTGRADIENT__o__FmTGRADIENT_BA__d__kgradPw_p;
  Tensor dFdt__x__Fm1_n1_dot_FmTGRADIENT_pB;
  double FmTGRADIENT__o__FmTGRADIENT_AB__dd__k_p;

  Tensor Stiffness_density_pAB; // Stiffness density contribution to the nodes A and B

  double V0_p; // Volume of the particle in the reference configuration
  double m_p; // Mass of the particle
  double J_p; // Jacobian of the deformation gradient
  double dJ_dt_n1_p;// Rate of the Jacobian of the deformation gradient
  double K_f_p; // Compresibility of the fluid
  double phi_f_0; // Reference volume fraction of the fluid phase
  double phi_s_0; // Reference volume fraction of the solid phase
  double phi_f_p; // Current volume fraction of the fluid phase
  double rho_f_p; // Current intrinsic density of the fluid phase
  double phi_s_p; // Current volume fraction of the solid phase
  double rho_s_p; // Current intrinsic density of the solid phase
  double relative_rho_f_p; // Current relative density of the fluid phase 
  double relative_rho_s_p; // Current relative density of the solid phase 
  double relative_rho_p; // Current relative density of the mixture
  double theta_n1_p; // Kirchhoff pore water pressure
  double d_theta_n1_p_dt; // Rate of the Kirchhoff pore water pressure
  double alpha_1 = Params.alpha_1; // Time integration parameter
  double alpha_4 = Params.alpha_4; // Time integration parameter
  double ShapeFunction_pA; // Nodal value of the shape function in node A
  double ShapeFunction_pB; // Nodal value of the shape function in node B
  double g = -9.81; // Gravity constant

  /*
    Define and allocate the tangent stiffness matrix
  */
  Matrix Tangent_Stiffness = allocZ__MatrixLib__(Order, Order);

  /*
    Loop in the particles for the assembling process
  */
  for(int p = 0 ; p<Np ; p++)
  {

    /*
      Define tributary nodes of the particle 
    */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /* 
      Evaluate the shape function in the coordinates of the particle
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    gradient_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
      
    /* 
      Load material properties for each phase
    */
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
    Material_Water_idx = Soil_Water_Mixtures[Mixture_idx].Water_Idx;
    MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];
    MatProp_Water_p = MPM_Mesh.Mat[Material_Water_idx];

    /*
      Get some material properties
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p]; 
    m_p = MPM_Mesh.Phi.mass.nV[p]; 
    K_f_p = MatProp_Water_p.Compressibility; 
    phi_f_0 = Soil_Water_Mixtures[Mixture_idx].phi_f_0;
    phi_s_0 = Soil_Water_Mixtures[Mixture_idx].phi_s_0;
    phi_f_p = MPM_Mesh.Phi.phi_f.nV[p]; // Volume fraction of fluid for particle p
    phi_s_p = MPM_Mesh.Phi.phi_s.nV[p]; // Volume fraction of solid for particle p
    rho_f_p = MPM_Mesh.Phi.rho_f.nV[p]; // Intrinsic density of fluid for particle p
    rho_s_p = MPM_Mesh.Phi.rho_s.nV[p]; // Intrinsic density of solid for particle p
    relative_rho_f_p = phi_f_p*rho_f_p; // Relative density of fluid for particle p
    relative_rho_s_p = phi_s_p*rho_s_p; // Relative density of solid for particle p
    relative_rho_p = MPM_Mesh.Phi.rho.nV[p]; // Relative density of mixture for particle p
    k_p = Soil_Water_Mixtures[Mixture_idx].Permeability; // Particle permeability

    /*
      Get some nodal values
    */
    Nodal_D_Acceleration_p = get_U_set_field_upw__MeshTools__(D_upw.d2_value_dt2, Nodes_p, ActiveNodes);
    Nodal_Pw_n_p = get_Pw_set_field_upw__MeshTools__(upw_n.value, Nodes_p, ActiveNodes);
    Nodal_D_Pw_p = get_Pw_set_field_upw__MeshTools__(D_upw.value, Nodes_p, ActiveNodes);
    
    /*
      Compute auxiliar particle velocity in the next time step
    */
    D_a_p = interpolate_vectorial_magnitude__MeshTools__(Nodal_D_Acceleration_p, ShapeFunction_p);
    a_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.acc.nM[p],1);
    a_n1_p = addition__TensorLib__(a_n_p, D_a_p);

    /*
      Compute current particle gradient of pore water pressure
    */
    gradPw_n1 = compute_Pore_water_pressure_gradient_n1(Nodal_Pw_n_p,Nodal_D_Pw_p,gradient_p);

    /*
      Compute auxiliar terms
    */
    b_p = MPM_Mesh.b;
    dyn_p = subtraction__TensorLib__(a_n1_p, b_p);
    kdyn_p = vector_linear_mapping__TensorLib__(k_p,dyn_p);
    kgradPw_p = vector_linear_mapping__TensorLib__(k_p,gradPw_n1);
      
    /*
      Take the values of the deformation gradient ant t = n and t = n + 1. 
      Later compute the midpoint deformation gradient and 
      the transpose of the deformation gradient at the midpoint.
    */
    F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    dFdt_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n1.nM[p],2);
    FT_n_p = transpose__TensorLib__(F_n_p);
    FT_n1_p = transpose__TensorLib__(F_n1_p);
    Fm1_n1_p = Inverse__TensorLib__(F_n1_p);
    FmT_n1_p = transpose__TensorLib__(Fm1_n1_p);
    FmT__dd__dFdt_n1_p = inner_product__TensorLib__(FmT_n1_p, dFdt_n1_p);
    dFdt__x__Fm1_n1_p = matrix_product__TensorLib__(dFdt_n1_p,Fm1_n1_p);

    /*
      Get the jacobian of the deformation gradient
      and its rate
    */
    J_p = MPM_Mesh.Phi.J.nV[p];
    dJ_dt_n1_p = MPM_Mesh.Phi.dJ_dt.nV[p];

    /*
      Get Kirchhoff pore water pressure and its rate
    */
    theta_n1_p = MPM_Mesh.Phi.Pw_n1.nV[p];
    d_theta_n1_p_dt = MPM_Mesh.Phi.d_Pw_dt.nV[p];


    for(int A = 0 ; A<NumNodes_p ; A++)
    {
      /* 
        Get the value of the shape function in node A
      */
      ShapeFunction_pA = ShapeFunction_p.nV[A];

      /*
        Compute the gradient in the reference configuration and the covariant pushforward
      */
      gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);
      GRADIENT_pA = vector_linear_mapping__TensorLib__(FT_n_p,gradient_pA);
      FmTGRADIENT_pA = vector_linear_mapping__TensorLib__(FmT_n1_p,GRADIENT_pA);
      /*
        Get the node of the mesh for the contribution 
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];


      for(int B = 0 ; B<NumNodes_p ; B++)
      {
        /* 
          Get the value of the shape function in node B
        */
        ShapeFunction_pB = ShapeFunction_p.nV[B];

        /*
          Compute the gradient in the reference configuration and the covariant pushforward
        */
        gradient_pB = memory_to_tensor__TensorLib__(gradient_p.nM[B], 1);
        GRADIENT_pB = vector_linear_mapping__TensorLib__(FT_n_p,gradient_pB);
        FmTGRADIENT_pB = vector_linear_mapping__TensorLib__(FmT_n1_p,GRADIENT_pB);

        /*
          Get the node of the mesh for the contribution 
        */
        Bp = Nodes_p.Connectivity[B];
        B_mask = ActiveNodes.Nodes2Mask[Bp];
        
        /*
          Compute operators
        */
        dyn__o__FmTGRADIENT_pB = dyadic_Product__TensorLib__(dyn_p, FmTGRADIENT_pB); 
        FmTGRADIENT__o__FmTGRADIENT_pAB = dyadic_Product__TensorLib__(FmTGRADIENT_pA, FmTGRADIENT_pB); 
        FmTGRADIENT__o__FmTGRADIENT_pBA = dyadic_Product__TensorLib__(FmTGRADIENT_pB, FmTGRADIENT_pA); 
        dFdt__x__Fm1_n1_dot_FmTGRADIENT_pB = vector_linear_mapping__TensorLib__(dFdt__x__Fm1_n1_p,FmTGRADIENT_pB);
        FmTGRADIENT__o__FmTGRADIENT_AB__d__kdyn_p = vector_linear_mapping__TensorLib__(FmTGRADIENT__o__FmTGRADIENT_pAB,kdyn_p);
        FmTGRADIENT__o__FmTGRADIENT_BA__d__kdyn_p = vector_linear_mapping__TensorLib__(FmTGRADIENT__o__FmTGRADIENT_pBA,kdyn_p);
        FmTGRADIENT__o__FmTGRADIENT_AB__d__kgradPw_p = vector_linear_mapping__TensorLib__(FmTGRADIENT__o__FmTGRADIENT_pAB,kgradPw_p);
        FmTGRADIENT__o__FmTGRADIENT_BA__d__kgradPw_p = vector_linear_mapping__TensorLib__(FmTGRADIENT__o__FmTGRADIENT_pBA,kgradPw_p);
        FmTGRADIENT__o__FmTGRADIENT_AB__dd__k_p = inner_product__TensorLib__(FmTGRADIENT__o__FmTGRADIENT_pAB, k_p);
        

        if(strcmp(MatProp_Soil_p.Type,"Neo-Hookean-Wriggers") == 0)
        {
          Stiffness_density_pAB = compute_stiffness_density_Neo_Hookean_Wriggers(GRADIENT_pA,GRADIENT_pB,F_n1_p,J_p,MatProp_Soil_p);
        }
        else
        {
          fprintf(stderr,"%s : %s %s %s \n","Error in assemble_Nodal_Tangent_Stiffness()",
            "The material",MatProp_Soil_p.Type,"has not been yet implemnented");
          exit(EXIT_FAILURE);
        }
        
        /*
          Add the geometric and material contribution to each dof for the assembling process
        */
        for(int i = 0 ; i<Ndim ; i++)
        {
          /*
            Compute the contribution of the linearised terms to the D_phi-D_theta direcction. 
          */
          Tangent_Stiffness.nM[A_mask*Ndof+i][B_mask*Ndof+Ndim] +=  
          + (relative_rho_f_p/(K_f_p*J_p))*ShapeFunction_pA*ShapeFunction_pB*dyn_p.n[i]*V0_p
          - FmTGRADIENT_pA.n[i]*ShapeFunction_pB*V0_p;

          for(int j = 0 ; j<Ndim ; j++)
          {
            /*
              Compute the contribution of the linearised terms to the D_phi-D_phi direcction. 
            */
            Tangent_Stiffness.nM[A_mask*Ndof+i][B_mask*Ndof+j] += 
            + (phi_f_0*rho_f_p/J_p -
              phi_s_0*rho_s_p/J_p -
              (relative_rho_f_p*theta_n1_p)/(K_f_p*J_p))*ShapeFunction_pA*dyn__o__FmTGRADIENT_pB.N[i][j]*V0_p
            + alpha_1*relative_rho_p*ShapeFunction_pA*ShapeFunction_pB*(i==j)*V0_p
            + theta_n1_p*FmTGRADIENT__o__FmTGRADIENT_pAB.N[j][i]*V0_p
            + Stiffness_density_pAB.N[i][j]*V0_p;
            /*
              Compute the contribution of the linearised terms to the D_theta-D_phi direcction. 
            */
            Tangent_Stiffness.nM[A_mask*Ndof+Ndim][B_mask*Ndof+j] += 
            + ((d_theta_n1_p_dt*phi_f_0*rho_f_p)/(K_f_p*J_p) - 
              (theta_n1_p*d_theta_n1_p_dt*relative_rho_f_p)/(K_f_p*K_f_p*J_p) - 
              (rho_f_p*theta_n1_p*dJ_dt_n1_p)/(K_f_p*J_p) + 
              rho_f_p*J_p*alpha_4 + 
              rho_f_p*J_p*FmT__dd__dFdt_n1_p)*ShapeFunction_pA*FmTGRADIENT_pB.n[j]*V0_p
            - rho_f_p*J_p*ShapeFunction_pA*dFdt__x__Fm1_n1_dot_FmTGRADIENT_pB.n[j]*V0_p 
            + (1/g)*FmTGRADIENT__o__FmTGRADIENT_AB__d__kgradPw_p.n[j]*V0_p
            + (1/g)*FmTGRADIENT__o__FmTGRADIENT_AB__dd__k_p*gradPw_n1.n[j]*V0_p 
            + (rho_f_p/(g*K_f_p))*FmTGRADIENT__o__FmTGRADIENT_BA__d__kdyn_p.n[j]*V0_p
            - (rho_f_p*J_p/g)*FmTGRADIENT__o__FmTGRADIENT_BA__d__kdyn_p.n[j]*V0_p
            + (rho_f_p*J_p/g)*FmTGRADIENT__o__FmTGRADIENT_AB__d__kdyn_p.n[j]*V0_p;
          }
        }

        /*
          Compute the contribution of the linearised terms to the D_theta-D_theta direcction.
        */
        Tangent_Stiffness.nM[A_mask*Ndof+Ndim][B_mask*Ndof+Ndim] =+ 
        ((d_theta_n1_p_dt*relative_rho_f_p)/(K_f_p*K_f_p*J_p) + 
          relative_rho_f_p*alpha_4/K_f_p +
          (dJ_dt_n1_p*rho_f_p)/(K_f_p*J_p))*ShapeFunction_pA*ShapeFunction_pB*V0_p
        - (1/g)*FmTGRADIENT__o__FmTGRADIENT_AB__dd__k_p*V0_p;

        /*
          Free memory 
        */
        free__TensorLib__(GRADIENT_pB);
        free__TensorLib__(FmTGRADIENT_pB);
        free__TensorLib__(dyn__o__FmTGRADIENT_pB);
        free__TensorLib__(FmTGRADIENT__o__FmTGRADIENT_pAB); 
        free__TensorLib__(FmTGRADIENT__o__FmTGRADIENT_pBA); 
        free__TensorLib__(FmTGRADIENT__o__FmTGRADIENT_AB__d__kdyn_p);
        free__TensorLib__(FmTGRADIENT__o__FmTGRADIENT_BA__d__kdyn_p);
        free__TensorLib__(FmTGRADIENT__o__FmTGRADIENT_AB__d__kgradPw_p);
        free__TensorLib__(FmTGRADIENT__o__FmTGRADIENT_BA__d__kgradPw_p);
        free__TensorLib__(dFdt__x__Fm1_n1_dot_FmTGRADIENT_pB);
        free__TensorLib__(Stiffness_density_pAB);
      }

      /*
        Free memory 
      */
      free__TensorLib__(GRADIENT_pA);
      free__TensorLib__(FmTGRADIENT_pA);
    }
      

    /* 
      Free memory 
    */
    free__TensorLib__(FT_n_p);
    free__TensorLib__(FT_n1_p);
    free__TensorLib__(Fm1_n1_p);
    free__TensorLib__(FmT_n1_p);
    free__TensorLib__(dFdt__x__Fm1_n1_p);
    free__TensorLib__(D_a_p);
    free__TensorLib__(a_n1_p);
    free__TensorLib__(dyn_p);
    free__TensorLib__(gradPw_n1);
    free__TensorLib__(kdyn_p);
    free__TensorLib__(kgradPw_p);
    free__MatrixLib__(ShapeFunction_p);
    free__MatrixLib__(gradient_p);
    free__MatrixLib__(Nodal_D_Acceleration_p);
    free__MatrixLib__(Nodal_Pw_n_p);
    free__MatrixLib__(Nodal_D_Pw_p);
    free(Nodes_p.Connectivity);

  }

  /*
    Return tangent matrix
  */
  return Tangent_Stiffness;
}

/**************************************************************/

static void solve_non_reducted_system(
  Nodal_Field D_upw,
  Matrix Tangent_Stiffness,
  Matrix Residual)
/*


*/
{
  int Nnodes_mask = Residual.N_rows;
  int Ndof = NumberDOF;
  int Order = Nnodes_mask*Ndof;
  int LDA   = Nnodes_mask*Ndof;
  int LDB   = Nnodes_mask*Ndof;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO = 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;

  /*
    Compute the LU factorization 
  */
  dgetrf_(&Order,&Order,Tangent_Stiffness.nV,&LDA,IPIV,&INFO);

  /*
    Check error messages in the LAPACK LU descompistion  
  */
  if(INFO)
  {
    fprintf(stderr,"%s : %s %s %s \n",
      "Error in solve_non_reducted_system",
      "The function","dgetrf_","returned an error message !!!" );
    exit(EXIT_FAILURE);
  }

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order,&NRHS,Tangent_Stiffness.nV,&LDA,IPIV,Residual.nV,&LDB,&INFO);
  free(IPIV);

  /*
    Check error messages in the LAPACK solver  
  */
  if(INFO)
  {
    fprintf(stderr,"%s : %s %s %s \n","Error in solve_non_reducted_system",
      "The function","dgetrs_","returned an error message !!!" );
    exit(EXIT_FAILURE);
  }

  /*
    Update 
  */
  for(int idx_A_i = 0 ; idx_A_i < Order ; idx_A_i++)
  {
    D_upw.value.nV[idx_A_i] -= Residual.nV[idx_A_i];
  }
  
}


/**************************************************************/

static void solve_reducted_system(
  Nodal_Field D_upw,
  Matrix Tangent_Stiffness,
  Matrix Residual,
  Mask Free_and_Restricted_Dofs)
/*


*/
{
  int Nnodes_mask = Residual.N_rows;
  int Ndof = NumberDOF;
  int Order = Nnodes_mask*Ndof;
  int Num_Free_dofs = Free_and_Restricted_Dofs.Nactivenodes;
  int Free_A_ij;
  int idx_A_ij, idx_B_ij;
  int Mask_idx_A_ij, Mask_idx_B_ij;

  /*
    Guyan reduction : static condensation
    Here, we will work directly with vectorized matrix
  */
  Matrix Tangent_Stiffness_FF = allocZ__MatrixLib__(Num_Free_dofs,Num_Free_dofs);
  Matrix Residual_F  = allocZ__MatrixLib__(Num_Free_dofs,1);

  for(idx_A_ij = 0 ; idx_A_ij < Order ; idx_A_ij++)
  {

    /* Get the index mask of the dof */
    Mask_idx_A_ij = Free_and_Restricted_Dofs.Nodes2Mask[idx_A_ij];
    
    /*
      Get the Residual with the Free dofs
    */
    if(Mask_idx_A_ij != - 1)
    {
      Residual_F.nV[Mask_idx_A_ij] = Residual.nV[idx_A_ij];
    } 

    for(idx_B_ij = 0 ; idx_B_ij < Order ; idx_B_ij++)
    { 

      /* Get the index mask of the dof */
      Mask_idx_B_ij = Free_and_Restricted_Dofs.Nodes2Mask[idx_B_ij];

      /*
          Get the K matrix with the Free-Free dofs
      */
      if((Mask_idx_A_ij != - 1) && (Mask_idx_B_ij != - 1))
      {
        Tangent_Stiffness_FF.nM[Mask_idx_A_ij][Mask_idx_B_ij] = Tangent_Stiffness.nM[idx_A_ij][idx_B_ij];
      }

    }

  }


/*
  Parameters for the solver
 */
 int LDA   = Num_Free_dofs;
 int LDB   = Num_Free_dofs;
 int Order_FF = Num_Free_dofs;
 char  TRANS = 'T'; /* (Transpose) */
 int   INFO = 3;
 int * IPIV = (int *)Allocate_Array(Num_Free_dofs,sizeof(int));
 int NRHS = 1;

  /*
    Compute the LU factorization 
  */
  dgetrf_(&Order_FF,&Order_FF,Tangent_Stiffness_FF.nV,&LDA,IPIV,&INFO);

  /*
    Check error messages in the LAPACK LU descompistion  
  */
  if(INFO)
  {
    fprintf(stderr,"%s : %s %s %s \n","Error in solve_reducted_system",
      "The function","dgetrf_","returned an error message !!!" );
    exit(EXIT_FAILURE);
  }

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order_FF,&NRHS,Tangent_Stiffness_FF.nV,&LDA,IPIV,Residual_F.nV,&LDB,&INFO);
  free(IPIV);
  
  /*
    Check error messages in the LAPACK solver  
  */
  if(INFO)
    {
      fprintf(stderr,"%s : %s %s %s \n","Error in solve_reducted_system",
        "The function","dgetrs_","returned an error message !!!" );
      exit(EXIT_FAILURE);
    }

  /*
    Update 
  */
  for(idx_A_ij =  0 ; idx_A_ij<Order ; idx_A_ij ++)
  { 

    /* Get the index mask of the dof */
    Mask_idx_A_ij = Free_and_Restricted_Dofs.Nodes2Mask[idx_A_ij];

    /*
      Get the Residual with the Free dofs
    */
    if(Mask_idx_A_ij != - 1)
    {
      D_upw.value.nV[idx_A_ij] -= Residual_F.nV[Mask_idx_A_ij];
    } 
      
  }

  /*
    Free auxiliar
  */
  free__MatrixLib__(Tangent_Stiffness_FF);
  free__MatrixLib__(Residual_F);

}

/**************************************************************/

static void update_Newmark_Nodal_Increments(
  Nodal_Field D_upw,
  Nodal_Field upw_n,
  Newmark_parameters Params)
{

  int Nnodes = upw_n.value.N_rows;
  int Ndof = NumberDOF;
  int Total_dof = Nnodes*NumberDOF;
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

  /*
    Update nodal variables
  */
  for(int A = 0 ; A<Total_dof ; A++)
  {  
    D_upw.d2_value_dt2.nV[A] = alpha_1*D_upw.value.nV[A] - (alpha_2 + 1)*upw_n.d_value_dt.nV[A] - alpha_3*upw_n.d2_value_dt2.nV[A];
    D_upw.d_value_dt.nV[A]   = alpha_4*D_upw.value.nV[A] + alpha_5*upw_n.d_value_dt.nV[A] + (alpha_6 - 1)*upw_n.d2_value_dt2.nV[A];
  }
}


/**************************************************************/

static void update_Particles(
  Nodal_Field D_upw,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes)
{
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ap;
  int A_mask;
  int idx_A_mask_i;
  int idx_ij;
  int Mixture_idx;
  int Material_Soil_idx;
  Element Nodes_p; /* Element for each particle */
  Material MatProp_Soil_p;
  Matrix ShapeFunction_p; /* Value of the shape-function in the particle */
  Matrix gradient_p;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor dFdt_n_p;
  Tensor dFdt_n1_p;
  double ShapeFunction_pI; /* Nodal value for the particle */
  double D_upw_pI;
  double D_dt_upw_pI;
  double D_dt2_upw_pI;
  double Vol_0_p;

  /* iterate over the particles */
  for(int p = 0 ; p<Np ; p++)
  {
      
    /* Define element of the particle */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);
      
    /*
      Evaluate the shape function and gradient in the coordinates of the particle 
    */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    gradient_p      = compute_dN__MeshTools__(Nodes_p,MPM_Mesh,FEM_Mesh);

    /*
      Take the values of the deformation gradient and its rates from the previous step
    */
    F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    dFdt_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n.nM[p],2);
    dFdt_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n1.nM[p],2);

    /*
      Replace the deformation gradient and it rate at t = n with the converged deformation gradient
    */
    for(int i = 0 ; i<Ndim  ; i++)
    {
      for(int j = 0 ; j<Ndim  ; j++)
      {
        F_n_p.N[i][j] = F_n1_p.N[i][j];
        dFdt_n_p.N[i][j] = dFdt_n1_p.N[i][j];
      }
    }

    /*
      Compute the deformation energy (reference volume + material properties (solid phase))
    */
    Vol_0_p = MPM_Mesh.Phi.Vol_0.nV[p];
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
    MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];
    MPM_Mesh.Phi.W.nV[p]= finite_strains_internal_energy__Particles__(F_n_p, MatProp_Soil_p,Vol_0_p);
      
    /*
      Iterate over the nodes of the particle
    */
    for(int A = 0; A<Nodes_p.NumberNodes; A++)
    {

      /*
        Get the node in the nodal momentum with the mask
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];
    
      /*
        Evaluate the GP function in the node 
      */
      ShapeFunction_pI = ShapeFunction_p.nV[A];
    
      /*
        Update the particle primitive fields using nodal values of the increments
      */
      for(int i = 0 ; i<Ndim ; i++)
      {
        D_upw_pI = ShapeFunction_pI*D_upw.value.nM[A_mask][i];
        D_dt_upw_pI = ShapeFunction_pI*D_upw.d_value_dt.nM[A_mask][i];
        D_dt2_upw_pI = ShapeFunction_pI*D_upw.d2_value_dt2.nM[A_mask][i];
                
        if(i<Ndim)
        {
          MPM_Mesh.Phi.acc.nM[p][i]  += D_dt2_upw_pI;
          MPM_Mesh.Phi.vel.nM[p][i]  += D_dt_upw_pI;
          MPM_Mesh.Phi.dis.nM[p][i]  += D_upw_pI;
          MPM_Mesh.Phi.x_GC.nM[p][i] += D_upw_pI;
        }
        else
        {
          MPM_Mesh.Phi.d2_Pw_dt2.nV[p] = D_dt2_upw_pI;
          MPM_Mesh.Phi.d_Pw_dt.nV[p] = D_dt_upw_pI;
          MPM_Mesh.Phi.Pw.nV[p] = D_upw_pI;
        }
      } 

    }

    /*
      Free memory
    */
    free(Nodes_p.Connectivity);
    free__MatrixLib__(ShapeFunction_p);
    free__MatrixLib__(gradient_p);
  }  
}

/**************************************************************/

static void output_selector(
  Nodal_Field D_upw,
  Nodal_Field upw_n,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  int TimeStep,
  int ResultsTimeStep)
{
  int Ndof = NumberDOF;
//  int Nnodes_mask = ;

  /*
    Define and allocate the output vector
  */
//  Nodal_Field upw_n1;
//  upw_n1.value        = allocZ__MatrixLib__(Nnodes_mask,Ndof);
//  upw_n1.d_value_dt   = allocZ__MatrixLib__(Nnodes_mask,Ndof);
//  upw_n1.d2_value_dt2 = allocZ__MatrixLib__(Nnodes_mask,Ndof);

//  for(int I = 0 ; I<Nnodes_mask*Ndof ; I++)
//  {
//    upw_n1.value.nV[I] = D_upw.value.nV[I] + upw_n.value.nV[I];
//    upw_n1.d_value_dt.nV[I] = D_upw.d_value_dt.nV[I] + upw_n.d_value_dt.nV[I];
//    upw_n1.d2_value_dt2.nV[I] = D_upw.d2_value_dt2.nV[I] + upw_n.d2_value_dt2.nV[I];
//  }

  /*
    vtk results
  */
  if(TimeStep % ResultsTimeStep == 0)
  {
    particle_results_vtk__InOutFun__(MPM_Mesh,TimeStep,ResultsTimeStep);

//    nodal_results_upw_vtk__InOutFun__(upw_n1.value, upw_n1.d_value_dt, FEM_Mesh, ActiveNodes, TimeStep, ResultsTimeStep);
  }

  /* 
    csv results 
  */
  for(int i = 0 ; i<Number_Out_nodal_path_csv ; i++)
  {

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Velocity)
    {
//      path_nodes_analysis_csv__InOutFun__(Velocity, FEM_Mesh.Coordinates,"Nodal_path_velocity_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Acceleration)
    {
 //     path_nodes_analysis_csv__InOutFun__(Acceleration, FEM_Mesh.Coordinates,"Nodal_path_acceleration_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_D_Displacement)
    {
//      path_nodes_analysis_csv__InOutFun__(D_Displacement, FEM_Mesh.Coordinates,"Nodal_path_displacement_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Forces)
    {
//      path_nodes_analysis_csv__InOutFun__(Forces, FEM_Mesh.Coordinates,"Nodal_path_forces_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Reactions)
    {
//      path_nodes_analysis_csv__InOutFun__(Reactions, FEM_Mesh.Coordinates,"Nodal_path_reactions_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Residual)
    {
//      path_nodes_analysis_csv__InOutFun__(Residual, FEM_Mesh.Coordinates,"Nodal_path_residual_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }
  }


  for(int i = 0 ; i<Number_Out_particles_path_csv ; i++)
  {
    if(Out_particles_path_csv[i].Out_csv_particles_path_Damage)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.chi, MPM_Mesh.Phi.x_GC, "Particles_path_damage_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Velocity)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.vel, MPM_Mesh.Phi.x_GC, "Particles_path_velocity_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Acceleration)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.acc, MPM_Mesh.Phi.x_GC, "Particles_path_acceleration_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Displacement)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.dis, MPM_Mesh.Phi.x_GC, "Particles_path_displacement_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Stress)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.Stress, MPM_Mesh.Phi.x_GC, "Particles_path_stress_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Strain)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.Strain, MPM_Mesh.Phi.x_GC, "Particles_path_strain_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_particles_path_csv[i].Out_csv_particles_path_Deformation_gradient)
    {
      path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.F_n, MPM_Mesh.Phi.x_GC, "Particles_path_deformation_gradient_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

  }

}

/**************************************************************/

static void standard_error()
{
  fprintf(stderr,"%s !!! \n",Error_message);
  exit(EXIT_FAILURE);
}

/**************************************************************/
