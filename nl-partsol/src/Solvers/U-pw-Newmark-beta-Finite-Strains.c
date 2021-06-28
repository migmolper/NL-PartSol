#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>

#endif

#include <omp.h>

/*
  Call global variables
*/
double Thickness_Plain_Stress;
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
static Matrix compute_Residual(Nodal_Field,Nodal_Field,Matrix,Mask,Particle,Mesh,Newmark_parameters,int);
static  void  compute_Inertial_Forces_Mixture(Nodal_Field,Nodal_Field,Matrix,Matrix,Mask,Particle,Newmark_parameters);
static  void  compute_Internal_Forces_Mixture(Matrix,Mask,Particle,Mesh);
static Tensor compute_total_first_Piola_Kirchhoff_stress(Tensor,double,Tensor);
static  void  compute_Rate_Mass_Fluid(Matrix,Mask,Particle,Mesh);
static  void  compute_Flow_contribution_Fluid(Nodal_Field, Nodal_Field,Matrix,Mask,Particle,Mesh);
static Tensor compute_Kirchoff_Pore_water_pressure_gradient_n1(Matrix,Matrix,Matrix);
static  void  compute_nominal_traction_and_fluid_flux(Matrix,Mask,Particle,Mesh,int);
static Matrix compute_Nodal_Reactions(Mesh,Matrix,Mask);
static  bool  check_convergence(Matrix,double,int,int,int);

static Matrix assemble_Tangent_Stiffness(Nodal_Field,Nodal_Field,Matrix,Mask,Particle,Mesh,double,Newmark_parameters);

static Matrix compute_mixture_stiffness_density(Tensor,Tensor,Tensor,double,double);
static Matrix compute_mixture_inertial_density(Tensor,Tensor,double,double,double,double,double,double,double,double,double,double,double);
static Matrix compute_water_flux_density(Tensor,Tensor,Tensor,Tensor,Tensor,double,double,double,double,double,double);
static Matrix compute_water_inertial_density(Tensor,Tensor,double,double,double,double,double,double,double,double,double,double,double);

static Tensor compute_stiffness_density(Tensor, Tensor, Tensor, double, Material);
static  void  system_reduction(Matrix,Matrix,Mask,Mesh);
static  void  solve_system(Nodal_Field,Matrix,Matrix);

static  void  update_Newmark_Nodal_Increments(Nodal_Field,Nodal_Field,Newmark_parameters);
static  void  update_Particles(Nodal_Field,Particle,Mesh,Mask);
static  void  output_selector(Nodal_Field,Nodal_Field,Particle,Mesh,Mask,int,int);
static  char  Error_message[MAXW];
static  void  standard_error();

/**************************************************************/

void upw_Newmark_beta_Finite_Strains(
  Mesh FEM_Mesh,
  Particle MPM_Mesh,
  Time_Int_Params Parameters_Solver)
{

  /*
    Auxiliar variables for the solver
  */
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nactivenodes;
  int InitialStep = Parameters_Solver.InitialTimeStep;
  int NumTimeStep = Parameters_Solver.NumTimeStep;  
  int MaxIter = Parameters_Solver.MaxIter;
  int Iter;

  double TOL = Parameters_Solver.TOL_Newmark_beta;
  double epsilon = Parameters_Solver.epsilon_Mass_Matrix;
  double beta = Parameters_Solver.beta_Newmark_beta;
  double gamma = Parameters_Solver.gamma_Newmark_beta;
  double CFL = Parameters_Solver.CFL;
  double DeltaTimeStep;
  double DeltaX = FEM_Mesh.DeltaX;

  bool Convergence;

  Matrix Effective_Mass;
  Nodal_Field D_upw;
  Nodal_Field upw_n;
  Matrix Residual;
  Matrix Reactions;
  Matrix Tangent_Stiffness;

  Mask ActiveNodes;

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
  DeltaTimeStep = 0.01;// DeltaT_Coussy__SolversLib__(MPM_Mesh, DeltaX, 1.0, CFL); 
  Params = compute_Newmark_parameters(beta, gamma, DeltaTimeStep);

  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
  {

    print_Status("*************************************************",TimeStep);
    print_step(TimeStep,DeltaTimeStep);

    print_Status("*************************************************",TimeStep);
    print_Status("First step : Generate Mask ... WORKING",TimeStep);

    ActiveNodes = generate_NodalMask__MeshTools__(FEM_Mesh);
    Nactivenodes = ActiveNodes.Nactivenodes;

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

      Residual = compute_Residual(upw_n,D_upw,Effective_Mass,ActiveNodes,MPM_Mesh,FEM_Mesh,Params,TimeStep);

      Reactions = compute_Nodal_Reactions(FEM_Mesh,Residual,ActiveNodes);

      Convergence = check_convergence(Residual,TOL,Iter,MaxIter,TimeStep);

      if(Convergence == false)
      {

        Tangent_Stiffness = assemble_Tangent_Stiffness(upw_n,D_upw,Effective_Mass,ActiveNodes,MPM_Mesh,FEM_Mesh,epsilon,Params);

        system_reduction(Tangent_Stiffness,Residual,ActiveNodes,FEM_Mesh);

        solve_system(D_upw,Tangent_Stiffness,Residual);

        update_Newmark_Nodal_Increments(D_upw,upw_n,Params);
        
        Iter++;

        free__MatrixLib__(Residual);
        free__MatrixLib__(Reactions);
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
    free__MatrixLib__(Reactions);
    free(ActiveNodes.Nodes2Mask);

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
  of the lumped mass matrix and the consistent mass matrix. 
*/
{

  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ndim = NumberDimensions;
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
  for(int p = 0 ; p<Np ; p++) 
  {

    /*
      Define tributary nodes of the particle 
    */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    /* 
      Evaluate the shape function in the coordinates of the particle
    */
    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

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
      Na_p = N_p.nV[A];


      for(int B = 0 ; B<Nodes_p.NumberNodes ; B++)
      {       
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
        for(int i = 0 ; i<Ndof ; i++)
        {
          Effective_MassMatrix.nM[A_mask*Ndof+i][B_mask*Ndof+i] += (1-epsilon)*m_p*Na_p*Nb_p + (A_mask==B_mask)*epsilon*m_p*Na_p;
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
  Matrix N_p;

  /* Evaluation of the particle in the node */
  double Na_p;
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
    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

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
      Na_p = N_p.nV[A];

      /*
        Nodal displacement and pressure (and its rates)
      */
      for(int i = 0 ; i<Ndof ; i++)
      {
        if(i<Ndim)
        {
          upw.value.nM[A_mask][i]        += m_p*Na_p*MPM_Mesh.Phi.dis.nM[p][i];
          upw.d_value_dt.nM[A_mask][i]   += m_p*Na_p*MPM_Mesh.Phi.vel.nM[p][i];
          upw.d2_value_dt2.nM[A_mask][i] += m_p*Na_p*MPM_Mesh.Phi.acc.nM[p][i];
        } 
        else
        {
          upw.value.nM[A_mask][i]        += m_p*Na_p*MPM_Mesh.Phi.Pw.nV[p]; 
          upw.d_value_dt.nM[A_mask][i]   += m_p*Na_p*MPM_Mesh.Phi.d_Pw_dt_n.nV[p];
          upw.d2_value_dt2.nM[A_mask][i] += m_p*Na_p*MPM_Mesh.Phi.d2_Pw_dt2.nV[p];          
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
    Call the LAPACK solver to compute the accelerations and second derivative of the pore water pressure
  */
  int Order = Nnodes_mask*Ndof;
  int LDA   = Order;
  int LDB   = Order;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO= 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;
  double * AUX_MEMORY = (double *)calloc(Order*Order,sizeof(double));

  /*
    Generate auxiliar copy of the mass matrix to avoid destructive operations
  */
  memcpy(AUX_MEMORY, Effective_Mass.nV, Order*Order*sizeof(double));

  /*
    Compute the LU factorization for the mass matrix
  */
  dgetrf_(&Order,&Order,AUX_MEMORY,&LDA,IPIV,&INFO);

  if(INFO != 0)
  {
    if(INFO < 0)
    {
      printf("%s : \n","Error in compute_Nodal_Field()");
      printf("the %i-th argument had an illegal value",abs(INFO));
    }
    else if(INFO > 0)
    {
      printf("%s :\n","Error in compute_Nodal_Field()");
      printf(" M(%i,%i) %s \n %s \n %s \n %s \n",INFO,INFO,"is exactly zero. The factorization",
        "has been completed, but the factor M is exactly",
        "singular, and division by zero will occur if it is used",
        "to solve a system of equations.");
    }    
    exit(EXIT_FAILURE);
  }

  /*
    Solve for the acceleration and second derivative of the pore water pressure
  */
  dgetrs_(&TRANS,&Order,&NRHS,AUX_MEMORY,&LDA,IPIV,upw.value.nV,&LDB,&INFO);
  dgetrs_(&TRANS,&Order,&NRHS,AUX_MEMORY,&LDA,IPIV,upw.d_value_dt.nV,&LDB,&INFO);
  dgetrs_(&TRANS,&Order,&NRHS,AUX_MEMORY,&LDA,IPIV,upw.d2_value_dt2.nV,&LDB,&INFO);

  /*
    Free auxiliar memory
  */
  free(AUX_MEMORY);
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
          D_upw.d_value_dt.nM[Id_BCC_mask][k] = 0.0;
          D_upw.d2_value_dt2.nM[Id_BCC_mask][k] = 0.0;
               
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
  Element Nodes_p; /* Element for each particle */
  Material MatProp_Soil_p; /* Variable with the material properties of the solid phase */
  Material MatProp_Water_p; /* Variable with the material properties of the fluid phase */
  Matrix gradient_p;
  Matrix N_p;
  Matrix Nodal_D_Displacement_p;
  Matrix Nodal_D_Velocity_p;
  Matrix Nodal_D_theta_p;
  Matrix Nodal_D_theta_dt;
  Tensor F_n_p; /* Deformation gradient of the soil skeleton (t = n) */
  Tensor F_n1_p; /* Deformation gradient of the soil skeleton (t = n + 1) */
  Tensor DF_p; /* Increment of the deformation gradient of the soil skeleton */
  Tensor dFdt_n_p; /* Rate of the deformation gradient of the soil skeleton (t = n) */
  Tensor dFdt_n1_p; /* Rate of the deformation gradient of the soil skeleton (t = n + 1) */
  Tensor dt_DF_p; /* Rate of the increment of the deformation gradient of the soil skeleton */
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
    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

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
      Get the Kirchhoff pore water pressure at t = n + 1
    */
    Nodal_D_theta_p = get_Pw_set_field_upw__MeshTools__(D_upw.value, Nodes_p, ActiveNodes);
    Nodal_D_theta_dt = get_Pw_set_field_upw__MeshTools__(D_upw.d_value_dt, Nodes_p, ActiveNodes);

    MPM_Mesh.Phi.Pw_n1.nV[p] = MPM_Mesh.Phi.Pw.nV[p] + interpolate_scalar_magnitude__MeshTools__(Nodal_D_theta_p, N_p);
    MPM_Mesh.Phi.d_Pw_dt_n1.nV[p] = MPM_Mesh.Phi.d_Pw_dt_n.nV[p] + interpolate_scalar_magnitude__MeshTools__(Nodal_D_theta_dt, N_p);

    /*
      Update state parameters
    */
    Pw_0 = MPM_Mesh.Phi.Pw_0.nV[p]; /* Get the initial pressure */
    Pw_n1 = MPM_Mesh.Phi.Pw_n1.nV[p]/J_n1_p; /* From the Kirchhoff pressure compure the cauchy pore water pressure */

    MPM_Mesh.Phi.rho_f.nV[p] = rho_f_0*exp((Pw_n1-Pw_0)/K_f); /* Update the fluid density */

    MPM_Mesh.Phi.phi_s.nV[p] = phi_s_0/J_n1_p; /* Update the volume fraction of the solid phase */
    MPM_Mesh.Phi.phi_f.nV[p] = 1.0 - (1.0 - phi_f_0)/J_n1_p; /* Update the volume fraction of the fluid phase */

    MPM_Mesh.Phi.rho.nV[p] = MPM_Mesh.Phi.rho_s.nV[p]*MPM_Mesh.Phi.phi_s.nV[p] +
                             MPM_Mesh.Phi.rho_f.nV[p]*MPM_Mesh.Phi.phi_f.nV[p];  /* Update density of the mixture */

    MPM_Mesh.Phi.J.nV[p] = J_n1_p; /* Update soil skeleton jacobian */
    MPM_Mesh.Phi.dJ_dt.nV[p] = dJ_dt_n1_p; /* Update soil skeleton rate of jacobian */

    /*
      Free memory 
    */
    free__MatrixLib__(Nodal_D_Displacement_p);
    free__MatrixLib__(Nodal_D_Velocity_p);
    free__MatrixLib__(Nodal_D_theta_p);
    free__MatrixLib__(Nodal_D_theta_dt);
    free__MatrixLib__(gradient_p);
    free__MatrixLib__(N_p);
    free(Nodes_p.Connectivity);

  }
  
  /*
    Loop in the material point set to update stress
  */
  for(int p = 0 ; p<Np ; p++)
  {

    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
    MatProp_Soil_p = MPM_Mesh.Mat[Material_Soil_idx];

    /*
      Update the first Piola-Kirchhoff stress tensor with an apropiate
      integration rule.
    */
    Stress_integration__Particles__(p,MPM_Mesh,FEM_Mesh,MatProp_Soil_p); 
  }

}

/**************************************************************/

static Matrix compute_Residual(
  Nodal_Field upw_n,
  Nodal_Field D_upw,
  Matrix Effective_Mass,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Newmark_parameters Params,
  int TimeStep)
{

  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;

  Matrix Residual = allocZ__MatrixLib__(Nnodes_mask,Ndof);

  compute_Inertial_Forces_Mixture(D_upw,upw_n,Effective_Mass,Residual,ActiveNodes,MPM_Mesh,Params);

  compute_Internal_Forces_Mixture(Residual,ActiveNodes,MPM_Mesh,FEM_Mesh);

  compute_Rate_Mass_Fluid(Residual, ActiveNodes, MPM_Mesh, FEM_Mesh);

  compute_Flow_contribution_Fluid(upw_n, D_upw, Residual, ActiveNodes, MPM_Mesh, FEM_Mesh);

  compute_nominal_traction_and_fluid_flux(Residual,ActiveNodes,MPM_Mesh,FEM_Mesh,TimeStep);

  return Residual;
}

/**************************************************************/

static void compute_Inertial_Forces_Mixture(
  Nodal_Field D_upw,
  Nodal_Field upw_n,
  Matrix Effective_Mass,
  Matrix Residual,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Newmark_parameters Params)
{
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Ndof*Nnodes_mask;
  Matrix Acceleration_n1 = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;

  /*
    Compute nodal acceleration
  */
  for(int A = 0 ; A<Nnodes_mask ; A++)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      Acceleration_n1.nM[A][i] = 
      alpha_1*D_upw.value.nM[A][i] - 
      alpha_2*upw_n.d_value_dt.nM[A][i] - 
      alpha_3*upw_n.d2_value_dt2.nM[A][i] - 
      MPM_Mesh.b.n[i];
    }
  }

  /*
    Compute inertial forces
  */
  for(int idx_A = 0 ; idx_A<Order ; idx_A++)
  {
    for(int idx_B = 0 ; idx_B<Order ; idx_B++)
    {
      Residual.nV[idx_A] += Effective_Mass.nM[idx_A][idx_B]*Acceleration_n1.nV[idx_B];
    }
  }


  free__MatrixLib__(Acceleration_n1);

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

  Tensor P_total_p; /* Total First Piola-Kirchhoff Stress tensor */
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
      Get the first Piola-Kirchhoff stress tensor.
    */
    P_effective_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2); 

    /*
      Get the Kirchhoff pore water pressure 
    */
    theta_n1_p = MPM_Mesh.Phi.Pw_n1.nV[p];

    /*
      Following Terzaghi's idea, the effective stress tensor in the reference configuration
      is computed:
    */
    P_total_p = compute_total_first_Piola_Kirchhoff_stress(P_effective_p,theta_n1_p,F_n1_p);

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
      InternalForcesDensity_Ap = vector_linear_mapping__TensorLib__(P_total_p, GRADIENT_pA);
      
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
    free__TensorLib__(P_total_p);
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

/**************************************************************/

static void compute_Rate_Mass_Fluid(
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
  int Mixture_idx;
  int Material_Water_idx;
  int NumNodes_p; /* Number of tributary nodes of p */
  int A_mask; /* Index of the node where we apply the body force */
  int Ap; /* Tributary node A of particle p */
  Material MatProp_Water_p; /* Variable with the material properties of the fluid phase */
  Element Nodes_p; /* Element for each particle p */
  Matrix N_p; /* Matrix with the value of the shape function in the particle p */ 
  double Na_p; /* Value of the shape funtion in node A for the particle p */
  double phi_f_p; /* Volume fractions of the fluid */
  double intrinsic_rho_f_p; /* Initrinsic density of the fluid */
  double relative_rho_f_p; /* Relative density of the fluid */
  double K_f; /* Compressibility (fluid) */
  double theta_dt_n1_p; /* Rate of pore water pressure (t = n + 1) in particle p */
  double V0_p; /* Volume of the particle */
  double dJ_dt_n1_p; /* Rate of the jacobian */

  for(int p = 0 ; p<Np ; p++)
  {

    /*
      Get the reference volume for each material point (mixture) 
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /*
      Get the rate of the jacobian for each material point
    */
    dJ_dt_n1_p = MPM_Mesh.Phi.dJ_dt.nV[p];
      
    /*
      Load intrinsic properties for the fluid phase to 
      get the fluid compressibility
    */
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Water_idx = Soil_Water_Mixtures[Mixture_idx].Water_Idx;
    MatProp_Water_p = MPM_Mesh.Mat[Material_Water_idx];
    K_f = MatProp_Water_p.Compressibility;

    /*
      Get the current intrinsic density, volume fraction 
      and compressibility for each material point (fluid).
      Compute relative density.
    */
    phi_f_p = MPM_Mesh.Phi.phi_f.nV[p];
    intrinsic_rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];
    relative_rho_f_p = phi_f_p*intrinsic_rho_f_p;

    /*
      Define nodes for each particle
    */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /*
      Evaluate the shape function in the coordinates of the particle 
    */
    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
        
    /*
      Get the value of the rate of pore water pressure
      at t = n + 1, in each particle
    */
    theta_dt_n1_p = MPM_Mesh.Phi.d_Pw_dt_n1.nV[p];

    for(int A = 0 ; A<NumNodes_p ; A++)
    {
      
      /*
        Compute the nodal contribution of the shape function
      */
      Na_p = N_p.nV[A];

      /*
        Get the node of the mesh for the contribution 
      */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];

      /*
        Add the contribution of the jacobian rate to the mass conservation
      */
      Residual.nM[A_mask][Ndim] += Na_p*(intrinsic_rho_f_p*dJ_dt_n1_p + theta_dt_n1_p*relative_rho_f_p/K_f)*V0_p;

    }
        
    /* 
      Free memory 
    */
    free__MatrixLib__(N_p);
    free(Nodes_p.Connectivity);
  }

}

/**************************************************************/

static void compute_Flow_contribution_Fluid(
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
  Matrix Nodal_theta_n;
  Matrix Nodal_D_theta; 
  Matrix Nodal_D_acceleration_p;

  /* Shape function variables */
  Matrix N_p; /* Shape function */
  Matrix gradient_N_p;  /* Shape function gradient */
  Tensor gradient_Na_n_p; /* Shape function gradient evaluated in Node A, def config */
  Tensor GRADIENT_Na_n_p; /* Shape function gradient evaluated in Node A, ref config */

  /* Deformation gradient tensor */
  Tensor F_n_p; /* Deformation gradient t = n */
  Tensor F_n1_p; /* Deformation gradient t = n + 1 */
  Tensor DF_p;
  Tensor FT_n_p; /* Transpose of the deformation gradient t = n */
  Tensor Fm1_n1_p; /* Inverse of the deformation agradient a t = n + 1 */
  Tensor FmT_n1_p; /* Transpose and inverse of the deformation agradient a t = n + 1 */
  Tensor k_p; /* Spatial permebility tensor */

  /* Variables to compute the relative flux */
  double g = - 10.0;
  double intrinsic_rho_f_p;
  double J_n1_p; /* Determianant of the soil skeleton deformation gradient at t = n + 1 */
  Tensor gradient_theta_n1_p;
  Tensor GRADIENT_theta_n1_p;
  Tensor a_n_p; /* Particle acceleration */
  Tensor D_a_p; /* Increment of the particle acceleration */
  Tensor a_n1_p; /* Particle acceleration */
  Tensor b_p; /* External acceleration */
  Tensor dyn_p;
  double V0_p; /* Volume of the particle at the reference configuration */

  /* Auxiliar values */
  Tensor Fm1__x__k_p;
  Tensor Fm1__x__k__x__FmT_p;
  Tensor Fm1K__x__dyn_p;
  Tensor Fm1KFmT__x__GRAD_theta;
  double GRADIENT_Na__x__Fm1KF__x__dyn;
  double GRADIENT_Na__x__Fm1KFmT__x__GRAD_theta;


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
      Compute shape function and its gradient in each node 
    */
    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    gradient_N_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Take the values of the deformation gradient at t = n and t = n + 1. 
      Later compute the transpose of the deformation gradient and the determinant.
    */
    F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    DF_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.DF.nM[p],2);
    FT_n_p = transpose__TensorLib__(F_n_p);
    Fm1_n1_p = Inverse__TensorLib__(F_n1_p);
    FmT_n1_p = transpose__TensorLib__(Fm1_n1_p);
    J_n1_p = MPM_Mesh.Phi.J.nV[p];

    /*
      Get the reference volume for each material point (mixture) 
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /* 
      Load intrinsic density of the fluid phase and the eulerian permeability
    */
    intrinsic_rho_f_p = MPM_Mesh.Phi.rho_f.nV[p];
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    k_p = Soil_Water_Mixtures[Mixture_idx].Permeability;

    /* 
      Compute particle pore water pressure gradient
    */
    Nodal_theta_n = get_Pw_set_field_upw__MeshTools__(upw_n.value, Nodes_p, ActiveNodes);
    Nodal_D_theta = get_Pw_set_field_upw__MeshTools__(D_upw.value, Nodes_p, ActiveNodes);
    gradient_theta_n1_p = compute_Kirchoff_Pore_water_pressure_gradient_n1(Nodal_theta_n,Nodal_D_theta,gradient_N_p);
    GRADIENT_theta_n1_p = vector_linear_mapping__TensorLib__(FT_n_p,gradient_theta_n1_p);

    /*
      Compute total particle acceleration at n+1 using nodal variables
    */
    Nodal_D_acceleration_p = get_U_set_field_upw__MeshTools__(D_upw.d2_value_dt2, Nodes_p, ActiveNodes);
    D_a_p = interpolate_vectorial_magnitude__MeshTools__(Nodal_D_acceleration_p, N_p);
    a_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.acc.nM[p],1);
    a_n1_p = addition__TensorLib__(D_a_p,a_n_p);
    b_p = MPM_Mesh.b;
    dyn_p = subtraction__TensorLib__(a_n1_p,b_p);


    /*
      Compute intermediate values
    */
    Fm1__x__k_p = matrix_product__TensorLib__(Fm1_n1_p,k_p);
    Fm1__x__k__x__FmT_p = matrix_product__TensorLib__(Fm1__x__k_p,FmT_n1_p);
    Fm1K__x__dyn_p = vector_linear_mapping__TensorLib__(Fm1__x__k_p,dyn_p);
    Fm1KFmT__x__GRAD_theta = vector_linear_mapping__TensorLib__(Fm1__x__k__x__FmT_p,GRADIENT_theta_n1_p);
    

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
      gradient_Na_n_p = memory_to_tensor__TensorLib__(gradient_N_p.nM[A], 1);
      GRADIENT_Na_n_p = vector_linear_mapping__TensorLib__(FT_n_p,gradient_Na_n_p);
  
      /*
        Intermediate results
      */
      GRADIENT_Na__x__Fm1KF__x__dyn = inner_product__TensorLib__(GRADIENT_Na_n_p,Fm1K__x__dyn_p);
      GRADIENT_Na__x__Fm1KFmT__x__GRAD_theta = inner_product__TensorLib__(GRADIENT_Na_n_p,Fm1KFmT__x__GRAD_theta);

      /* 
        Compute nodal contribution to the mass conservation
      */
      Residual.nM[A_mask][Ndim] -= (1.0/g)*(GRADIENT_Na__x__Fm1KFmT__x__GRAD_theta + J_n1_p*intrinsic_rho_f_p*GRADIENT_Na__x__Fm1KF__x__dyn)*V0_p;

      free__TensorLib__(GRADIENT_Na_n_p);
    }

    /* 
      Free the value of the shape functions 
    */
    free__MatrixLib__(N_p);
    free__MatrixLib__(gradient_N_p);
    free__MatrixLib__(Nodal_theta_n);
    free__MatrixLib__(Nodal_D_theta);
    free__MatrixLib__(Nodal_D_acceleration_p);

    free__TensorLib__(FT_n_p);
    free__TensorLib__(Fm1_n1_p);
    free__TensorLib__(FmT_n1_p);

    free__TensorLib__(gradient_theta_n1_p);
    free__TensorLib__(GRADIENT_theta_n1_p);
    free__TensorLib__(D_a_p);
    free__TensorLib__(a_n1_p);
    free__TensorLib__(dyn_p);

    free__TensorLib__(Fm1__x__k_p);
    free__TensorLib__(Fm1__x__k__x__FmT_p);
    free__TensorLib__(Fm1K__x__dyn_p);
    free__TensorLib__(Fm1KFmT__x__GRAD_theta);


    free(Nodes_p.Connectivity);      
  }

}

/**************************************************************/

static Tensor compute_Kirchoff_Pore_water_pressure_gradient_n1(
  Matrix Nodal_theta_n,
  Matrix Nodal_D_theta,
  Matrix gradient_p)
/*

*/
 {

  /* Variable definition */
  int Ndim = NumberDimensions;
  int Nnodes_p = Nodal_theta_n.N_rows; 
  double theta_a_n;
  double D_theta_a;
  double theta_a_n1;
  Tensor gradient_Na_n_p;
  Tensor gradient_theta_n1_p = alloc__TensorLib__(1);

  for(int a = 0 ; a<Nnodes_p ; a++)
  {

    /*
      Get nodal value (a) of the kirchhoff pore water pressure at n+1
    */
    theta_a_n = Nodal_theta_n.nV[a];
    D_theta_a = Nodal_D_theta.nV[a];
    theta_a_n1 = theta_a_n + D_theta_a;

    /*
      Get nodal value (a) of the gradient evaluated in the particle position (p) at n+1
    */
    gradient_Na_n_p = memory_to_tensor__TensorLib__(gradient_p.nM[a], 1);

    /*
      Compute nodal contribution
    */
    for(int i = 0 ; i<Ndim ; i++)
    {
     gradient_theta_n1_p.n[i] += theta_a_n1*gradient_Na_n_p.n[i];
    } 

  }

  return gradient_theta_n1_p;

 } 


/**************************************************************/


static void compute_nominal_traction_and_fluid_flux(
  Matrix Residual,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  int TimeStep)
{
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  Load Load_i;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_p; /* Nodal values of the sahpe function */
  double Na_p;
  Tensor T = alloc__TensorLib__(1); // Nominal traction
  double Q; // Nominal flux
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
    Load_i = MPM_Mesh.Neumann_Contours.BCC_i[i];

    NumNodesLoad = Load_i.NumNodes;
      
    for(int j = 0 ; j<NumNodesLoad ; j++)
    {

      /*
        Get the index of the particle
      */
      p = Load_i.Nodes[j];

      /*
        Get the area of each particle
      */      
      if(Ndim == 2)
      {
        A0_p = MPM_Mesh.Phi.Vol_0.nV[p]/Thickness_Plain_Stress;
      }
      else if(Ndim == 3)
      {
        A0_p = MPM_Mesh.Phi.Area_0.nV[p];
      }

      /*
        Define tributary nodes of the particle 
      */
      NumNodes_p = MPM_Mesh.NumberNodes[p];
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p],NumNodes_p);

      /* 
        Evaluate the shape function in the coordinates of the particle
      */
      N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
        Fill vector of contact forces
      */
      for(int k = 0 ; k<Ndof ; k++)
      {
        if(Load_i.Dir[k] == 1)
        {
          if((TimeStep < 0) || (TimeStep > Load_i.Value[k].Num))
          {
            sprintf(Error_message,"%s : %s",
              "Error in compute_Contact_Forces_Mixture()",
              "The time step is out of the curve !!");
            standard_error();
          }

          if(k<Ndim)
          {
            T.n[k] = Load_i.Value[k].Fx[TimeStep];
          }
          else if(k == Ndim)
          {
            Q = Load_i.Value[k].Fx[TimeStep];
          }
          
        }
        else
        {
          if(k<Ndim)
          {
            T.n[k] = 0.0;
          }
          else if(k == Ndim)
          {
            Q = 0.0;
          }
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
        Na_p = N_p.nV[A];

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
          Residual.nM[A_mask][k] -= Na_p*T.n[k]*A0_p;
        }

        Residual.nM[A_mask][Ndim] -= Na_p*Q*A0_p;  
      }

      /* Free the matrix with the nodal gradient of the element */
      free__MatrixLib__(N_p);
      free(Nodes_p.Connectivity);
  
    }

      
  }

  free__TensorLib__(T);

}


/**************************************************************/

static Matrix compute_Nodal_Reactions(
  Mesh FEM_Mesh,
  Matrix Residual,
  Mask ActiveNodes)
/*
  Compute the nodal reactions
*/
{
  /* 1º Define auxilar variables */
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

  Matrix Reactions = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  strcpy(Reactions.Info,"REACTIONS");

  /*
    Loop over the the boundaries 
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
      NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;
    
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
    for(int k = 0 ; k<NumDimBound ; k++)
      {

        /* 
     Apply only if the direction is active (1) 
        */
        if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1)
    {
      /* 
         Set to zero the Residual in the nodes where velocity is fixed 
      */
      Reactions.nM[Id_BCC_mask][k] = Residual.nM[Id_BCC_mask][k];
      Residual.nM[Id_BCC_mask][k] = 0;
    }
      }
  }    
    }

  return Reactions;
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
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nnodes_mask = Residual.N_cols;
  double Error = 0.0;
  double Error_relative = 0.0;

  if(Iter > MaxIter)
  {
      fprintf(stderr,"%s : %s !!! \n","Error in check_convergence()",
        "Convergence not reached in the maximum number of iterations");
      return true;
      //exit(EXIT_FAILURE);
    }
  else
    {
      /*
        Compute absolute error 
      */
      for(int A = 0 ; A<Nnodes_mask ; A++)
      {
        for(int i = 0 ; i<Ndim ; i++)
        {
          Error += DSQR(Residual.nM[A][i]);
        }
      }

      Error = pow(Error,0.5);


      /*
        Compute relative error
      */
      if(Iter == 0)
      {
        Error0 = Error;
        Error_relative = Error/Error0;      
      }
      else
      {
        Error_relative = Error/Error0;
      }
      
      /*
        Check convergence using the relative error
      */
      if(Error_relative > TOL)
      {
        printf("%e -> %e\n",Error,Error_relative);
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

static Matrix assemble_Tangent_Stiffness(
  Nodal_Field upw_n,
  Nodal_Field D_upw,
  Matrix Effective_Mass,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  double epsilon,
  Newmark_parameters Params)
/*
  This function computes the tangent stiffness matrix as a combination
  of the geometrical stiffness matrix and the material stiffness matrix.
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

  /*
    Nodal values
  */
  Matrix Nodal_D_Acceleration_p; // Nodal values of the acceleration increments
  Matrix Nodal_Pw_n_p; // Nodal values of the pore-water pressure
  Matrix Nodal_D_Pw_p; // Nodal values of the pore-water pressure increments

  /*
    Particle properties
  */
  Material MatProp_Soil_p; // structure with the material properties of the soil
  Material MatProp_Water_p; // structure with the material properties of the water
  Tensor k_p; // Particle permeability and its transpose
  Tensor D_a_p; // Particle increment of acceleration
  Tensor a_n1_p; // Particle acceleration at t = n + 1 
  Tensor a_n_p; // Particle acceleration at t = n 
  Tensor b_p; // External accelerations
  Tensor dyn_p;
  Tensor gradPw_n1; // Particle pore water pressure
  Tensor gradient_theta_n1;
  Tensor F_n_p; // Particle deformation gradient at t = n 
  Tensor DF_p;
  Tensor FT_n_p; // Transpose of the particle deformation gradient at t = n 
  Tensor F_n1_p; // Particle deformation gradient at t = n + 1
  Tensor Fm1_n1_p; // Inverse of the particle deformation gradient at t = n + 1
  Tensor FmT_n1_p; // Inverse of the transposed particle deformation gradient at t = n + 1
  Tensor dFdt_n1_p; // Rate of the particle deformation gradient at t = n + 1
  Tensor Stiffness_density_pAB;
  /*
    Auxiliar variables
  */
  double div_v_p; // FmT : dFdt_n1
  Tensor grad_v_p; // dFdt x Fm1_n1
  Tensor transpose_grad_v_p; // (dFdt x Fm1_n1)^T
  

  /* Shape functions variables */
  Matrix N_p; /* Shape functions */
  Matrix gradient_N_p; /* Shape functions gradients */
  Tensor gradient_Na_p, GRADIENT_Na_p; // Shape function gradient in the deformed/reference configurations (Node A)
  Tensor gradient_Nb_p, GRADIENT_Nb_p; // Shape function gradient in the deformed/reference configurations (Node B)
  Tensor FmTGRADIENT_Na_p; // Covariant push-forward of the GRADIENT_pA vector
  Tensor FmTGRADIENT_Nb_p; // Covariant push-forward of the GRADIENT_pB vector
  double Na_p; // Nodal value of the shape function in node A
  double Nb_p; // Nodal value of the shape function in node B

  double V0_p; // Volume of the particle in the reference configuration
  double m_p; // Particle mass
  double J_p; // Jacobian of the deformation gradient
  double K_f_p; // Compresibility of the fluid
  double phi_f_0; // Reference volume fraction of the fluid phase
  double phi_s_0; // Reference volume fraction of the solid phase
  double phi_f_p; // Volume fraction of the fluid phase
  double phi_s_p; // Volume fraction of the solid phase
  double intrinsic_rho_f_p; // Intrinsic density of the fluid phase
  double intrinsic_rho_s_p; // Intrinsic density of the solid phase
  double theta_n1_p; // Kirchhoff pore water pressure
  double d_theta_n1_p_dt; // Rate of the Kirchhoff pore water pressure
  double alpha_1 = Params.alpha_1; // Time integration parameter
  double alpha_4 = Params.alpha_4; // Time integration parameter

  /*
    Define and allocate the tangent stiffness matrix
  */
  Matrix Tangent_Stiffness = allocZ__MatrixLib__(Order, Order);

  Matrix mixture_inertial_density;
  Matrix mixture_stiffness_density;
  Matrix water_inertial_density;
  Matrix water_flux_density;

  /*
    Loop in the particles for the assembling process
  */
  for(int p = 0 ; p<Np ; p++)
  {

    /*
      Particle mass
     */
    m_p = MPM_Mesh.Phi.mass.nV[p];

    /*
      Define tributary nodes of the particle 
    */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /* 
      Evaluate the shape function in the coordinates of the particle
    */
    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    gradient_N_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
      
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
    K_f_p = MatProp_Water_p.Compressibility; 
    phi_f_0 = Soil_Water_Mixtures[Mixture_idx].phi_f_0;
    phi_s_0 = Soil_Water_Mixtures[Mixture_idx].phi_s_0;
    phi_f_p = MPM_Mesh.Phi.phi_f.nV[p]; // Volume fraction of fluid for particle p
    phi_s_p = MPM_Mesh.Phi.phi_s.nV[p]; // Volume fraction of solid for particle p
    intrinsic_rho_f_p = MPM_Mesh.Phi.rho_f.nV[p]; // Intrinsic density of fluid for particle p
    intrinsic_rho_s_p = MPM_Mesh.Phi.rho_s.nV[p]; // Intrinsic density of solid for particle p
    k_p = Soil_Water_Mixtures[Mixture_idx].Permeability; // Particle permeability
    
    /*
      Get some nodal values
    */
    Nodal_D_Acceleration_p = get_U_set_field_upw__MeshTools__(D_upw.d2_value_dt2, Nodes_p, ActiveNodes);
    
    /*
      Take the values of the deformation gradient ant t = n and t = n + 1. 
      Later compute the midpoint deformation gradient and 
      the transpose of the deformation gradient at the midpoint.
    */
    F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    DF_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.DF.nM[p],2);
    dFdt_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n1.nM[p],2);
    FT_n_p = transpose__TensorLib__(F_n_p);
    Fm1_n1_p = Inverse__TensorLib__(F_n1_p);
    FmT_n1_p = transpose__TensorLib__(Fm1_n1_p);

    /*
      Get the jacobian of the deformation gradient
    */
    J_p = MPM_Mesh.Phi.J.nV[p];

    /*
      Get some usefull intermediate results
    */
    div_v_p = inner_product__TensorLib__(FmT_n1_p, dFdt_n1_p);
    grad_v_p = matrix_product__TensorLib__(dFdt_n1_p,Fm1_n1_p);
    transpose_grad_v_p = transpose__TensorLib__(grad_v_p);

    /*
      Get Kirchhoff pore water pressure and its rate
    */
    theta_n1_p = MPM_Mesh.Phi.Pw_n1.nV[p];
    d_theta_n1_p_dt = MPM_Mesh.Phi.d_Pw_dt_n1.nV[p];

    /* 
      Compute particle pore water pressure gradient
    */
    Matrix Nodal_theta_n = get_Pw_set_field_upw__MeshTools__(upw_n.value, Nodes_p, ActiveNodes);
    Matrix Nodal_D_theta = get_Pw_set_field_upw__MeshTools__(D_upw.value, Nodes_p, ActiveNodes);
    Tensor gradient_theta_n1_p = compute_Kirchoff_Pore_water_pressure_gradient_n1(Nodal_theta_n,Nodal_D_theta,gradient_N_p);
    Tensor GRADIENT_theta_n1_p = vector_linear_mapping__TensorLib__(FT_n_p,gradient_theta_n1_p);
    Tensor FmTGRADIENT_theta_p = vector_linear_mapping__TensorLib__(FmT_n1_p,GRADIENT_theta_n1_p);
    
    /*
      Compute mixture/fluid dynamics in the next time step
    */
    D_a_p = interpolate_vectorial_magnitude__MeshTools__(Nodal_D_Acceleration_p, N_p);
    a_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.acc.nM[p],1);
    a_n1_p = addition__TensorLib__(D_a_p,a_n_p);
    b_p = MPM_Mesh.b;
    dyn_p = subtraction__TensorLib__(a_n1_p,b_p);    


    for(int A = 0 ; A<NumNodes_p ; A++)
    {
      /* 
        Get the value of the shape function in node A
      */
      Na_p = N_p.nV[A];

      /*
        Compute the gradient in the reference configuration and the covariant pushforward
      */
      gradient_Na_p = memory_to_tensor__TensorLib__(gradient_N_p.nM[A], 1);
      GRADIENT_Na_p = vector_linear_mapping__TensorLib__(FT_n_p,gradient_Na_p);
      FmTGRADIENT_Na_p = vector_linear_mapping__TensorLib__(FmT_n1_p,GRADIENT_Na_p);

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
        Nb_p = N_p.nV[B];

        /*
          Compute the gradient in the reference configuration and the covariant pushforward
        */
        gradient_Nb_p = memory_to_tensor__TensorLib__(gradient_N_p.nM[B], 1);
        GRADIENT_Nb_p = vector_linear_mapping__TensorLib__(FT_n_p,gradient_Nb_p);
        FmTGRADIENT_Nb_p = vector_linear_mapping__TensorLib__(FmT_n1_p,GRADIENT_Nb_p);

        /*
          Get the node of the mesh for the contribution 
        */
        Bp = Nodes_p.Connectivity[B];
        B_mask = ActiveNodes.Nodes2Mask[Bp];

        if(strcmp(MatProp_Soil_p.Type,"Neo-Hookean-Wriggers") == 0)
        {
          Stiffness_density_pAB = compute_stiffness_density_Neo_Hookean_Wriggers(GRADIENT_Na_p,GRADIENT_Nb_p,F_n1_p,J_p,MatProp_Soil_p);
        }
        else
        {
          fprintf(stderr,"%s : %s %s %s \n","Error in assemble_Tangent_Stiffness()",
            "The material",MatProp_Soil_p.Type,"has not been yet implemnented");
          exit(EXIT_FAILURE);
        }        

        /*
          Compute the contribution of the mixture to the tangent matrix
        */
        mixture_inertial_density = compute_mixture_inertial_density(
          FmTGRADIENT_Nb_p,dyn_p,Na_p,Nb_p,J_p,
          intrinsic_rho_f_p,intrinsic_rho_s_p,
          phi_f_p,phi_s_p,phi_f_0,phi_s_0,K_f_p,
          theta_n1_p);

        mixture_stiffness_density = compute_mixture_stiffness_density(
          FmTGRADIENT_Na_p,FmTGRADIENT_Nb_p,
          Stiffness_density_pAB,Nb_p,theta_n1_p);

        /*
          Compute the contribution of the fluid phase to the tangent matrix
        */
        water_inertial_density = compute_water_inertial_density(
          FmTGRADIENT_Nb_p,transpose_grad_v_p,Na_p,
          Nb_p,div_v_p,J_p,intrinsic_rho_f_p,
          phi_f_p,phi_f_0,theta_n1_p,d_theta_n1_p_dt,
          K_f_p,alpha_4);

        water_flux_density = compute_water_flux_density(
          FmTGRADIENT_Na_p,FmTGRADIENT_Nb_p,
          FmTGRADIENT_theta_p,k_p,dyn_p,Nb_p,
          J_p,intrinsic_rho_f_p,theta_n1_p,
          K_f_p,alpha_1);

        for(int i = 0 ; i<Ndof ; i++)
        { 

          for(int j = 0 ; j<Ndof ; j++)
          {
            if(i<Ndim)
            {
              Tangent_Stiffness.nM[A_mask*Ndof+i][B_mask*Ndof+j] += 
              (i==j)*((1-epsilon)*Na_p*Nb_p + (A_mask==B_mask)*epsilon*Na_p)*m_p +
              (mixture_inertial_density.nM[i][j] + mixture_stiffness_density.nM[i][j])*V0_p;
            }
            else
            {
              Tangent_Stiffness.nM[A_mask*Ndof+i][B_mask*Ndof+j] += 
              (water_inertial_density.nV[j] - water_flux_density.nV[j])*V0_p;
            }        
          }
        }

        /*
          Free memory 
        */
        free__TensorLib__(GRADIENT_Nb_p);
        free__TensorLib__(FmTGRADIENT_Nb_p);
        free__TensorLib__(Stiffness_density_pAB);
        free__MatrixLib__(mixture_inertial_density);
        free__MatrixLib__(mixture_stiffness_density);
        free__MatrixLib__(water_inertial_density);
        free__MatrixLib__(water_flux_density);
      }

      /*
        Free memory 
      */
      free__TensorLib__(GRADIENT_Na_p);
      free__TensorLib__(FmTGRADIENT_Na_p);
    }
      

    /* 
      Free memory 
    */
    free(Nodes_p.Connectivity);
    free__TensorLib__(gradient_theta_n1_p);
    free__MatrixLib__(N_p);
    free__MatrixLib__(gradient_N_p);
  }

  /*
    Return tangent matrix
  */
  return Tangent_Stiffness;

}

/**************************************************************/

static Matrix compute_mixture_inertial_density(
  Tensor FmTGRADIENT_Nb_p,
  Tensor dyn_p,
  double Na_p,
  double Nb_p,
  double J_p,
  double intrinsic_rho_f_p,
  double intrinsic_rho_s_p,
  double phi_f_p,
  double phi_s_p,
  double phi_f_0,
  double phi_s_0,
  double K_f_p,
  double theta_n1_p)
{
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;

  double rho_0 = J_p*(intrinsic_rho_f_p*phi_f_p + intrinsic_rho_s_p*phi_s_p);
  double C1 = phi_f_p*intrinsic_rho_f_p/(K_f_p*J_p);
  double C2 = (1/J_p)*(intrinsic_rho_f_p*(1 - phi_f_0) - phi_f_p*intrinsic_rho_f_p*theta_n1_p/K_f_p);
  double C3 = C2 - intrinsic_rho_s_p*phi_s_0/J_p;

  Matrix mixture_inertial_density = allocZ__MatrixLib__(Ndim,Ndof);

  Tensor dyn__o__FmTGRADIENT_Nb_p = dyadic_Product__TensorLib__(dyn_p,FmTGRADIENT_Nb_p);

  for(int i = 0; i<Ndim ; i++)
  {
    for(int j = 0; j<Ndim ; j++)
    {
      mixture_inertial_density.nM[i][j] = Na_p*(J_p*C3 + rho_0)*dyn__o__FmTGRADIENT_Nb_p.N[i][j];
    }

    mixture_inertial_density.nM[i][Ndim] = Na_p*J_p*C1*dyn_p.n[i]*Nb_p;

  }

  free__TensorLib__(dyn__o__FmTGRADIENT_Nb_p);

  return mixture_inertial_density;
}

/**************************************************************/

static Matrix compute_mixture_stiffness_density(
  Tensor FmTGRADIENT_Na_p,
  Tensor FmTGRADIENT_Nb_p,
  Tensor Stiffness_density_pAB,
  double Nb_p,
  double theta_n1_p)
{  
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;

  Matrix mixture_stiffness_density = allocZ__MatrixLib__(Ndim,Ndof);

  Tensor FmTGRADIENT_Na_p__o__FmTGRADIENT_Nb_p = dyadic_Product__TensorLib__(FmTGRADIENT_Na_p,FmTGRADIENT_Nb_p);

  for(int i = 0; i<Ndim ; i++)
  {
    for(int j = 0; j<Ndim ; j++)
    {
      mixture_stiffness_density.nM[i][j] = Stiffness_density_pAB.N[i][j] + theta_n1_p*FmTGRADIENT_Na_p__o__FmTGRADIENT_Nb_p.N[i][j];
    }

    mixture_stiffness_density.nM[i][Ndim] = - FmTGRADIENT_Na_p.n[i]*Nb_p;

  }

  free__TensorLib__(FmTGRADIENT_Na_p__o__FmTGRADIENT_Nb_p);

  return mixture_stiffness_density;
}

/**************************************************************/

static Matrix compute_water_inertial_density(
  Tensor FmTGRADIENT_Nb_p,
  Tensor transpose_grad_v_p,
  double Na_p,
  double Nb_p,
  double div_v_p,
  double J_p,
  double intrinsic_rho_f_p,
  double phi_f_p,
  double phi_f_0,
  double theta_n1_p,
  double d_theta_n1_p_dt,
  double K_f_p,
  double alpha_4)
{
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;

  Matrix water_flux_density = allocZ__MatrixLib__(1,Ndof);

  Tensor T_grad_v__x__FmTGRADIENT_Nb_p = vector_linear_mapping__TensorLib__(transpose_grad_v_p,FmTGRADIENT_Nb_p);

  for(int i = 0; i<Ndim ; i++)
  {
    water_flux_density.nV[i] = 
    + (d_theta_n1_p_dt*(1-phi_f_0)*intrinsic_rho_f_p/(K_f_p*J_p))*(Na_p*FmTGRADIENT_Nb_p.n[i])
    - (d_theta_n1_p_dt*phi_f_p*intrinsic_rho_f_p*theta_n1_p/(K_f_p*K_f_p*J_p))*(Na_p*FmTGRADIENT_Nb_p.n[i])
    - (intrinsic_rho_f_p*div_v_p*theta_n1_p/(K_f_p))*(Na_p*FmTGRADIENT_Nb_p.n[i])
    + (intrinsic_rho_f_p*J_p*div_v_p)*(Na_p*FmTGRADIENT_Nb_p.n[i])
    + (intrinsic_rho_f_p*alpha_4*J_p)*(Na_p*FmTGRADIENT_Nb_p.n[i])
    - (intrinsic_rho_f_p*J_p)*(Na_p*T_grad_v__x__FmTGRADIENT_Nb_p.n[i]);
  }

  water_flux_density.nV[Ndim] = 
  + (d_theta_n1_p_dt*phi_f_p*intrinsic_rho_f_p/(J_p*K_f_p*K_f_p))*(Na_p*Nb_p)
  + (phi_f_p*intrinsic_rho_f_p*alpha_4/K_f_p)*(Na_p*Nb_p)
  + (intrinsic_rho_f_p*div_v_p/K_f_p)*(Na_p*Nb_p);


  free__TensorLib__(T_grad_v__x__FmTGRADIENT_Nb_p);

  return water_flux_density; 
}

/**************************************************************/

static Matrix compute_water_flux_density(
  Tensor FmTGRADIENT_Na_p,
  Tensor FmTGRADIENT_Nb_p,
  Tensor FmTGRADIENT_theta_p,
  Tensor k_p,
  Tensor dyn_p,
  double Nb_p,
  double J_p,
  double intrinsic_rho_f_p,
  double theta_n1_p,
  double K_f_p,
  double alpha_1)
{
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  double g = - 10.0; // Gravity constant

  Matrix water_flux_density = allocZ__MatrixLib__(1,Ndof);

  Tensor FmTGRADIENT_Na_p__o__FmTGRADIENT_Nb_p = dyadic_Product__TensorLib__(FmTGRADIENT_Na_p,FmTGRADIENT_Nb_p);
  Tensor FmTGRADIENT_Nb_p__o__FmTGRADIENT_Na_p = transpose__TensorLib__(FmTGRADIENT_Na_p__o__FmTGRADIENT_Nb_p);
  Tensor K__x__FmTGRADIENT_Na_p = vector_linear_mapping__TensorLib__(k_p,FmTGRADIENT_Na_p);
  Tensor K__x__FmTGRADIENT_Nb_p = vector_linear_mapping__TensorLib__(k_p,FmTGRADIENT_Nb_p);
  Tensor K__x__FmTGRADIENT_theta_p = vector_linear_mapping__TensorLib__(k_p,FmTGRADIENT_theta_p);
  Tensor K__x__dyn_p = vector_linear_mapping__TensorLib__(k_p,dyn_p);
  Tensor FmTGRADIENT_a__o__FmTGRADIENT_b__x__KFmTGRADIENT_theta_p = vector_linear_mapping__TensorLib__(FmTGRADIENT_Na_p__o__FmTGRADIENT_Nb_p,K__x__FmTGRADIENT_theta_p);

  Tensor FmTGRADIENT_a__o__FmTGRADIENT_b__x__Kdyn_p = vector_linear_mapping__TensorLib__(FmTGRADIENT_Na_p__o__FmTGRADIENT_Nb_p,K__x__dyn_p);
  Tensor FmTGRADIENT_b__o__FmTGRADIENT_a__x__Kdyn_p = vector_linear_mapping__TensorLib__(FmTGRADIENT_Nb_p__o__FmTGRADIENT_Na_p,K__x__dyn_p);

  double FmTGRADIENT__x__K__x__FmTGRADIENT_p = inner_product__TensorLib__(FmTGRADIENT_Na_p,K__x__FmTGRADIENT_Nb_p);
  double FmTGRADIENT_a__x__K__x__dyn_p = inner_product__TensorLib__(FmTGRADIENT_Na_p,K__x__dyn_p);

  for(int i = 0; i<Ndim ; i++)
  { 
    water_flux_density.nV[i] = 
    - (1.0/g)*FmTGRADIENT_a__o__FmTGRADIENT_b__x__KFmTGRADIENT_theta_p.n[i]
    - (1.0/g)*FmTGRADIENT__x__K__x__FmTGRADIENT_p*FmTGRADIENT_theta_p.n[i]
    - (J_p*intrinsic_rho_f_p/g)*FmTGRADIENT_a__o__FmTGRADIENT_b__x__Kdyn_p.n[i]
    + (J_p*intrinsic_rho_f_p/g)*FmTGRADIENT_b__o__FmTGRADIENT_a__x__Kdyn_p.n[i]
    - (intrinsic_rho_f_p*theta_n1_p/(g*K_f_p))*FmTGRADIENT_b__o__FmTGRADIENT_a__x__Kdyn_p.n[i]
    + (J_p*intrinsic_rho_f_p*alpha_1/g)*K__x__FmTGRADIENT_Na_p.n[i]*Nb_p;
  }

  water_flux_density.nV[Ndim] = 
  + (1.0/g)*FmTGRADIENT__x__K__x__FmTGRADIENT_p
  + (intrinsic_rho_f_p/(K_f_p*g)*FmTGRADIENT_a__x__K__x__dyn_p*Nb_p);

  free__TensorLib__(FmTGRADIENT_Na_p__o__FmTGRADIENT_Nb_p);
  free__TensorLib__(FmTGRADIENT_Nb_p__o__FmTGRADIENT_Na_p);
  free__TensorLib__(K__x__FmTGRADIENT_Na_p);
  free__TensorLib__(K__x__FmTGRADIENT_Nb_p);
  free__TensorLib__(K__x__FmTGRADIENT_theta_p);
  free__TensorLib__(K__x__dyn_p);
  free__TensorLib__(FmTGRADIENT_a__o__FmTGRADIENT_b__x__KFmTGRADIENT_theta_p);
  free__TensorLib__(FmTGRADIENT_a__o__FmTGRADIENT_b__x__Kdyn_p);
  free__TensorLib__(FmTGRADIENT_b__o__FmTGRADIENT_a__x__Kdyn_p);

  return water_flux_density;
}

/**************************************************************/

static void system_reduction(
  Matrix Tangent_Stiffness,
  Matrix Residual,
  Mask ActiveNodes,
  Mesh FEM_Mesh)
{

    /* 
    Define auxilar variables 
  */
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Nnodes_mask*Ndof;
  int Number_of_BCC = FEM_Mesh.Bounds.NumBounds;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

  /* 
    Loop over the the boundaries to find the constrained dofs
  */
  for(int i = 0 ; i<Number_of_BCC ; i++)
    {

      /* 
        Get the number of nodes of this boundary
      */
      NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

      /* 
        Get the number of dimensions where the BCC it is applied 
      */
      NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;

      for(int j = 0 ; j<NumNodesBound ; j++)
        {
          /* 
            Get the index of the node 
          */
          Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
          Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

          /*
            If the boundary condition is under an active node 
          */
          if(Id_BCC_mask != -1)
          {
            /* 
              Loop over the dimensions of the boundary condition 
            */
            for(int k = 0 ; k<NumDimBound ; k++)
            {
              if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1)
              {

                for(int A_mask = 0 ; A_mask < Nnodes_mask ; A_mask++)
                {
                  for(int l = 0 ; l<Ndof ; l++)
                  {
                    Tangent_Stiffness.nM[A_mask*Ndof+l][Id_BCC_mask*Ndof+k] = 0.0;
                    Tangent_Stiffness.nM[Id_BCC_mask*Ndof+k][A_mask*Ndof+l] = 0.0;
                  }
                }

                Tangent_Stiffness.nM[Id_BCC_mask*Ndof+k][Id_BCC_mask*Ndof+k] = 1.0;

              }
            }
          }

        }    
    }

}

/**************************************************************/

static void solve_system(
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
    Check error messages in the LAPACK LU descomposition  
  */
  if(INFO)
  {
    fprintf(stderr,"%s : %s %s %s \n",
      "Error in solve_system",
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
    fprintf(stderr,"%s : %s %s %s \n","Error in solve_system",
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
  for(int A = 0 ; A<Nnodes ; A++)
  {  
    for(int i = 0 ; i<Ndof ; i++)
    {
      D_upw.d2_value_dt2.nM[A][i] = alpha_1*D_upw.value.nM[A][i] - alpha_2*upw_n.d_value_dt.nM[A][i] - (alpha_3+1)*upw_n.d2_value_dt2.nM[A][i];
      D_upw.d_value_dt.nM[A][i]   = alpha_4*D_upw.value.nM[A][i] + (alpha_5-1)*upw_n.d_value_dt.nM[A][i] + alpha_6*upw_n.d2_value_dt2.nM[A][i];
    }
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
  int Ndof = NumberDOF;
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
  for(int p = 0 ; p<Np ; p++)
  {
      
    /* Define element of the particle */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);
      
    /*
      Evaluate the shape function and gradient in the coordinates of the particle 
    */
    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    gradient_p = compute_dN__MeshTools__(Nodes_p,MPM_Mesh,FEM_Mesh);

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
      N_pI = N_p.nV[A];
    
      /*
        Update the particle primitive fields using nodal values of the increments
      */
      for(int i = 0 ; i<Ndof ; i++)
      {
        D_upw_pI = N_pI*D_upw.value.nM[A_mask][i];
        D_dt_upw_pI = N_pI*D_upw.d_value_dt.nM[A_mask][i];
        D_dt2_upw_pI = N_pI*D_upw.d2_value_dt2.nM[A_mask][i];
                
        if(i<Ndim)
        {
          MPM_Mesh.Phi.acc.nM[p][i]  += D_dt2_upw_pI;
          MPM_Mesh.Phi.vel.nM[p][i]  += D_dt_upw_pI;
          MPM_Mesh.Phi.dis.nM[p][i]  += D_upw_pI;
          MPM_Mesh.Phi.x_GC.nM[p][i] += D_upw_pI;
        }
        else
        {
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
/*  for(int i = 0 ; i<Number_Out_nodal_path_csv ; i++)
  {

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Velocity)
    {
      path_nodes_analysis_csv__InOutFun__(Velocity, FEM_Mesh.Coordinates,"Nodal_path_velocity_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Acceleration)
    {
      path_nodes_analysis_csv__InOutFun__(Acceleration, FEM_Mesh.Coordinates,"Nodal_path_acceleration_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_D_Displacement)
    {
      path_nodes_analysis_csv__InOutFun__(D_Displacement, FEM_Mesh.Coordinates,"Nodal_path_displacement_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Forces)
    {
      path_nodes_analysis_csv__InOutFun__(Forces, FEM_Mesh.Coordinates,"Nodal_path_forces_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Reactions)
    {
      path_nodes_analysis_csv__InOutFun__(Reactions, FEM_Mesh.Coordinates,"Nodal_path_reactions_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Residual)
    {
      path_nodes_analysis_csv__InOutFun__(Residual, FEM_Mesh.Coordinates,"Nodal_path_residual_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
    }
  }*/


  // for(int i = 0 ; i<Number_Out_particles_path_csv ; i++)
  // {
  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Damage)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.chi, MPM_Mesh.Phi.x_GC, "Particles_path_damage_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Velocity)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.vel, MPM_Mesh.Phi.x_GC, "Particles_path_velocity_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Acceleration)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.acc, MPM_Mesh.Phi.x_GC, "Particles_path_acceleration_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Displacement)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.dis, MPM_Mesh.Phi.x_GC, "Particles_path_displacement_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Stress)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.Stress, MPM_Mesh.Phi.x_GC, "Particles_path_stress_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Strain)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.Strain, MPM_Mesh.Phi.x_GC, "Particles_path_strain_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_particles_path_csv[i].Out_csv_particles_path_Deformation_gradient)
  //   {
  //     path_particles_analysis_csv__InOutFun__(MPM_Mesh.Phi.F_n, MPM_Mesh.Phi.x_GC, "Particles_path_deformation_gradient_csv", Out_particles_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  // }

}

/**************************************************************/

static void standard_error()
{
  fprintf(stderr,"%s !!! \n",Error_message);
  exit(EXIT_FAILURE);
}

/**************************************************************/
