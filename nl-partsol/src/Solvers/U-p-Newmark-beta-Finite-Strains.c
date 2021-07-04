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

static Matrix compute_Residual(Nodal_Field,Nodal_Field,Matrix,Mask,Particle,Mesh,Newmark_parameters,int);
static  void  compute_Inertial_Forces(Nodal_Field,Nodal_Field,Matrix,Matrix,Mask,Particle,Newmark_parameters);
static  void  compute_Internal_Forces(Matrix,Mask,Particle,Mesh);
static  void  compute_Volumetric_Constrain_Forces(Matrix,Mask,Particle,Mesh);
static  void  compute_Nodal_Nominal_traction_Forces(Matrix,Mask,Particle,Mesh,int);
static Matrix compute_Nodal_Reactions(Mesh,Matrix,Mask);
static  bool  check_convergence(Matrix,double,int,int,int);

static Matrix assemble_Tangent_Stiffness(Mask,Particle,Mesh,Newmark_parameters);
static  void  system_reduction(Matrix,Matrix,Mask,Mesh);
static  void  solve_system(Nodal_Field,Matrix,Matrix);

static  void  update_Newmark_Nodal_Increments(Nodal_Field,Nodal_Field,Newmark_parameters);
static  void  update_Local_State(Nodal_Field, Mask, Particle, Mesh);
static Matrix get_U_set_field_Up(Matrix, Element, Mask);
static Matrix get_P_set_field_Up(Matrix, Element, Mask);

static  void  update_Particles(Nodal_Field,Particle,Mesh,Mask);
static  void  output_selector(Nodal_Field,Nodal_Field,Particle,Mesh,Mask,int,int);
static  char  Error_message[MAXW];
static  void  standard_error();

/**************************************************************/

void Up_Newmark_beta_Finite_Strains(
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
  Nodal_Field D_up;
  Nodal_Field up_n;
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
  DeltaTimeStep = U_DeltaT__SolversLib__(MPM_Mesh, DeltaX, CFL);
  Params = compute_Newmark_parameters(beta, gamma, DeltaTimeStep);

  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
  {

    print_Status("*************************************************",TimeStep);
    print_step(TimeStep,DeltaTimeStep);

    print_Status("*************************************************",TimeStep);
    print_Status("First step : Generate Mask ... WORKING",TimeStep);
    local_search__MeshTools__(MPM_Mesh,FEM_Mesh);
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
    up_n = compute_Nodal_Field(Effective_Mass,MPM_Mesh,FEM_Mesh,ActiveNodes);
    D_up = initialise_Nodal_Increments(ActiveNodes, FEM_Mesh, TimeStep);

    print_Status("DONE !!!",TimeStep);  

    print_Status("*************************************************",TimeStep);
    print_Status("Four step : Compute equilibrium ... WORKING",TimeStep);
    

    Convergence = false;
    Iter = 0;
    while(Convergence == false)
    {

      Residual = compute_Residual(up_n,D_up,Effective_Mass,ActiveNodes,MPM_Mesh,FEM_Mesh,Params,TimeStep);

      Reactions = compute_Nodal_Reactions(FEM_Mesh,Residual,ActiveNodes);

      Convergence = check_convergence(Residual,TOL,Iter,MaxIter,TimeStep);

      if(Convergence == false)
      {

        Tangent_Stiffness = assemble_Tangent_Stiffness(ActiveNodes,MPM_Mesh,FEM_Mesh,Params);

        system_reduction(Tangent_Stiffness,Residual,ActiveNodes,FEM_Mesh);

        solve_system(D_up,Tangent_Stiffness,Residual);

        update_Newmark_Nodal_Increments(D_up,up_n,Params);

        update_Local_State(D_up,ActiveNodes,MPM_Mesh,FEM_Mesh);
        
        Iter++;

        free__MatrixLib__(Residual);
        free__MatrixLib__(Reactions);
        free__MatrixLib__(Tangent_Stiffness);
      }
    }

    print_Status("DONE !!!",TimeStep);

    print_Status("*************************************************",TimeStep);
    print_Status("Seven step : Update particles lagrangian ... WORKING",TimeStep);

    update_Particles(D_up,MPM_Mesh,FEM_Mesh,ActiveNodes);

    print_Status("DONE !!!",TimeStep);

    /*
      Outputs
    */
    output_selector(up_n, D_up, MPM_Mesh, FEM_Mesh, ActiveNodes, TimeStep, ResultsTimeStep);
  
    print_Status("*************************************************",TimeStep);
    print_Status("Eight step : Reset nodal values ... WORKING",TimeStep);
    
    /*
      Free memory.
    */
    free__MatrixLib__(Effective_Mass); 
    free__MatrixLib__(up_n.value);
    free__MatrixLib__(up_n.d_value_dt);
    free__MatrixLib__(up_n.d2_value_dt2);
    free__MatrixLib__(D_up.value);
    free__MatrixLib__(D_up.d_value_dt);
    free__MatrixLib__(D_up.d2_value_dt2);
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
  Nodal_Field D_up;
  D_up.value        = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  D_up.d_value_dt   = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  D_up.d2_value_dt2 = allocZ__MatrixLib__(Nnodes_mask,Ndof);

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
          D_up.value.nM[Id_BCC_mask][k] = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*(double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];  
          D_up.d_value_dt.nM[Id_BCC_mask][k] = 0.0;
          D_up.d2_value_dt2.nM[Id_BCC_mask][k] = 0.0;
               
        }
      }
    }    
  }

  /*
    Add some usefull info
  */
  strcpy(D_up.value.Info,"Nodal-fields");
  strcpy(D_up.d_value_dt.Info,"First-time-derivative-nodal-fields");
  strcpy(D_up.d2_value_dt2.Info,"Second-time-derivative-nodal-fields");

  return D_up;
}

/**************************************************************/

static void update_Local_State(
  Nodal_Field D_Up,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
{

  /*
    Auxiliar variables
  */
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int MatIndx_p;
  int Nnodes_p;
  double J_n1_p;  
  Element Nodes_p;
  Material MatProp_p;
  Matrix N_p;
  Matrix gradient_p;
  Matrix D_Displacement_Ap;
  Matrix D_Velocity_Ap;
  Matrix D_lambda_pressure_n1;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor DF_p;
  Tensor dFdt_n_p;
  Tensor dFdt_n1_p;
  Tensor dt_DF_p;
  double lambda_pressure_n1;
  double lambda_pressure_n;
  
  /*
    Loop in the material point set 
  */
  for(int p = 0 ; p<Np ; p++)
    {
      /*
  Define tributary nodes of the particle 
      */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);
      
      /*
        Get the nodal increment of displacement using the mask
      */
      D_Displacement_Ap = get_U_set_field_Up(D_Up.value, Nodes_p, ActiveNodes);
      D_Velocity_Ap = get_U_set_field_Up(D_Up.d_value_dt, Nodes_p, ActiveNodes);
      D_lambda_pressure_n1 = get_P_set_field_Up(D_Up.value, Nodes_p, ActiveNodes);

      /*
        Evaluate the shape function gradient in the coordinates of the particle
      */
      gradient_p = compute_dN__MeshTools__(Nodes_p,MPM_Mesh,FEM_Mesh);
      N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    
      /*
        Take the values of the deformation gradient from the previous step
      */
      F_n_p     = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p    = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
      DF_p      = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.DF.nM[p],2);
      dFdt_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n.nM[p],2);
      dFdt_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n1.nM[p],2);
      dt_DF_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_DF.nM[p],2);
      
      /*
        Compute the increment of the deformation gradient
      */
      update_increment_Deformation_Gradient__Particles__(DF_p,D_Displacement_Ap,gradient_p);
      update_rate_increment_Deformation_Gradient__Particles__(dt_DF_p,D_Velocity_Ap,gradient_p);

      /*
        Update the deformation gradient in t = n + 1 with the information
        from t = n and the increment of deformation gradient.
      */  
      update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, DF_p);
      update_rate_Deformation_Gradient_n1__Particles__(dFdt_n1_p, dt_DF_p, F_n_p, DF_p, dFdt_n_p);

      /*
        Update pressure lagrange multiplier
      */
      lambda_pressure_n = MPM_Mesh.Phi.lambda_pressure_n.nV[p];
      lambda_pressure_n1 = lambda_pressure_n + interpolate_scalar_magnitude__MeshTools__(D_lambda_pressure_n1, N_p);
      MPM_Mesh.Phi.lambda_pressure_n1.nV[p] = lambda_pressure_n1;

      /*
        Compute Jacobian of the deformation gradient
      */
      MPM_Mesh.Phi.J.nV[p] = I3__TensorLib__(F_n1_p);
            
      /*
         Free memory 
      */
      free__MatrixLib__(D_Displacement_Ap);
      free__MatrixLib__(D_Velocity_Ap);
      free__MatrixLib__(D_lambda_pressure_n1);
      free__MatrixLib__(gradient_p);
      free__MatrixLib__(N_p);
      free(Nodes_p.Connectivity);
    
    }

  /*
    Loop in the material point set to update stress
  */
  for(int p = 0 ; p<Np ; p++)
  {

    MatIndx_p = MPM_Mesh.MatIdx[p];
    MatProp_p = MPM_Mesh.Mat[MatIndx_p];

    /*
      Update the first Piola-Kirchhoff stress tensor with an apropiate
      integration rule.
    */
    Stress_integration__Particles__(p,MPM_Mesh,FEM_Mesh,MatProp_p); 
  }
  
}

/*********************************************************************/

static Matrix get_U_set_field_Up(
  Matrix Field_Up,
  Element Nodes_p, 
  Mask ActiveNodes)
{
  int Nnodes = Nodes_p.NumberNodes;
  int Ndim = NumberDimensions;
  Matrix Field_U_Ap = allocZ__MatrixLib__(Nnodes,Ndim);
  int Ap;
  int A_mask;


  for(int A = 0 ; A<Nnodes ; A++)
  {
    /* 
      Get the node in the mass matrix with the mask
    */
    Ap = Nodes_p.Connectivity[A];
    A_mask = ActiveNodes.Nodes2Mask[Ap];
     
    /*
      Get nodal field for particle p
    */
    for(int i = 0 ; i<Ndim ; i++)
    {
      Field_U_Ap.nM[A][i] = Field_Up.nM[A_mask][i];
    }
  }

  return Field_U_Ap;
}

/*********************************************************************/

static Matrix get_P_set_field_Up(
  Matrix Field_Up, 
  Element Nodes_p, 
  Mask ActiveNodes)
{
  int Nnodes = Nodes_p.NumberNodes;
  int Ndim = NumberDimensions;
  Matrix Field_p_Ap = allocZ__MatrixLib__(Nnodes,1);
  int Ap;
  int A_mask;

  for(int A = 0 ; A<Nnodes ; A++)
  {
    /* 
      Get the node in the mass matrix with the mask
    */
    Ap = Nodes_p.Connectivity[A];
    A_mask = ActiveNodes.Nodes2Mask[Ap];
    
    /*
      Get nodal field for particle p
    */
    Field_p_Ap.nV[A] = Field_Up.nM[A_mask][Ndim];
  }

  return Field_p_Ap;
}

/**************************************************************/

static Matrix compute_Residual(
  Nodal_Field up_n,
  Nodal_Field D_up,
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

  compute_Inertial_Forces(D_up,up_n,Effective_Mass,Residual,ActiveNodes,MPM_Mesh,Params);

  compute_Internal_Forces(Residual,ActiveNodes,MPM_Mesh,FEM_Mesh);

  compute_Volumetric_Constrain_Forces(Residual, ActiveNodes, MPM_Mesh, FEM_Mesh);

  compute_Nodal_Nominal_traction_Forces(Residual,ActiveNodes,MPM_Mesh,FEM_Mesh,TimeStep);

  return Residual;
}

/**************************************************************/

static void compute_Inertial_Forces(
  Nodal_Field D_up,
  Nodal_Field up_n,
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
      alpha_1*D_up.value.nM[A][i] - 
      alpha_2*up_n.d_value_dt.nM[A][i] - 
      alpha_3*up_n.d2_value_dt2.nM[A][i] - 
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

static void compute_Internal_Forces(
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

  Tensor P_constrained_p; /* Constrained First Piola-Kirchhoff Stress tensor */
  Tensor InternalForcesDensity_Ap;
  Element Nodes_p; /* List of nodes for particle */
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pA;
  Tensor GRADIENT_pA;
  Tensor F_n_p;
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

    /*
      Get the volume of the particle in the reference configuration 
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /*
      Get the first Piola-Kirchhoff stress tensor with the volumetric constrain.
    */
    P_constrained_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2); 

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
      InternalForcesDensity_Ap = vector_linear_mapping__TensorLib__(P_constrained_p, GRADIENT_pA);
      
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
    free__MatrixLib__(gradient_p);
    free(Nodes_p.Connectivity);
  }
  
}

/**************************************************************/

static void compute_Volumetric_Constrain_Forces(
  Matrix Residual,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
/*!
 * Weak form of the lagrange multiplier constrain
 * */
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int NumNodes_p;
  Element Nodes_p; /* List of nodes for particle */
  Matrix N_p; /* Shape function */
  double Na_p;
  double J_n1_p;
  double Volumetric_Constrain;
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
      Evaluate the shape function in the coordinates of the particle 
    */
    N_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /*
      Get the volume of the particle in the reference configuration 
    */
    V0_p = MPM_Mesh.Phi.Vol_0.nV[p];

    /*
      Get the Jacobian
    */
    J_n1_p = MPM_Mesh.Phi.J.nV[p];

    /*
      Impose the volumetric constrain
    */
    Volumetric_Constrain = (J_n1_p - 1.0);

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
        Asign the nodal forces contribution to the node 
      */
      Residual.nM[A_mask][Ndim] += Na_p*Volumetric_Constrain*V0_p;

    }
        
    /* 
      Free memory 
    */
    free__MatrixLib__(N_p);
    free(Nodes_p.Connectivity);
  }

}


/**************************************************************/

static void compute_Nodal_Nominal_traction_Forces(
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
      for(int k = 0 ; k<Ndim ; k++)
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
          
          T.n[k] = Load_i.Value[k].Fx[TimeStep];
          
        }
        else
        {
          T.n[k] = 0.0;
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
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Newmark_parameters Params)

/*
  This function computes the tangent stiffness matrix. 

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
  int MatIndx_p;
  int idx_AB_mask_ij;

  Element Nodes_p; /* List of nodes for particle */
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pA;
  Tensor GRADIENT_pA;
  Tensor gradient_pB;
  Tensor GRADIENT_pB;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor dFdt_n1_p;
  Tensor transpose_F_n_p;
  Matrix Stiffness_density_p;

  Material MatProp_p;
  double V0_p; /* Volume of the particle in the reference configuration */
  double J_p; /* Jacobian of the deformation gradient */
  double alpha_4 = Params.alpha_4; /* Newmark parameter (rate-dependent models) */

  Matrix Tangent_Stiffness = allocZ__MatrixLib__(Order, Order);

  /*
    Loop in the particles for the assembling process
  */
  for(int p = 0 ; p<Np ; p++)
    {

      /*
  Get the volume of the particle in the reference configuration
       */
      V0_p = MPM_Mesh.Phi.Vol_0.nV[p];
      
      /*
  Material properties of the particle
       */      
      MatIndx_p = MPM_Mesh.MatIdx[p];
      MatProp_p = MPM_Mesh.Mat[MatIndx_p];

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
  the transpose of the deformation gradient at the midpoint.
      */
      F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
      dFdt_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n1.nM[p],2);
      transpose_F_n_p = transpose__TensorLib__(F_n_p);

      /*
  Compute the jacobian of the deformation gradient in the deformed configuration
      */
      J_p = I3__TensorLib__(F_n1_p);

      for(int A = 0 ; A<NumNodes_p ; A++)
  {
    /*
      Compute the gradient in the reference configuration 
    */
    gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);
    GRADIENT_pA = vector_linear_mapping__TensorLib__(transpose_F_n_p,gradient_pA);

    /*
      Get the node of the mesh for the contribution 
    */
    Ap = Nodes_p.Connectivity[A];
    A_mask = ActiveNodes.Nodes2Mask[Ap];

    
    for(int B = 0 ; B<NumNodes_p ; B++)
      {

        /*
    Compute the gradient in the reference configuration 
        */
        gradient_pB = memory_to_tensor__TensorLib__(gradient_p.nM[B], 1);
        GRADIENT_pB = vector_linear_mapping__TensorLib__(transpose_F_n_p,gradient_pB);

        /*
    Get the node of the mesh for the contribution 
        */
        Bp = Nodes_p.Connectivity[B];
        B_mask = ActiveNodes.Nodes2Mask[Bp];

        /*
          Get the stiffness density of each particle
        */
        if(strcmp(MatProp_p.Type,"Newtonian-Fluid-Incompressible") == 0)
        {
          Stiffness_density_p = compute_stiffness_density_Newtonian_Fluid_Incompressible(GRADIENT_pA,GRADIENT_pB,F_n1_p,dFdt_n1_p,J_p,alpha_4,MatProp_p);
        }
        else
        {
          fprintf(stderr,"%s : %s %s %s \n",
            "Error in assemble_Tangent_Stiffness()",
            "The material",MatProp_p.Type,"has not been yet implemnented");
          exit(EXIT_FAILURE);
        }
        
        /*
    Add the geometric contribution to each dof for the assembling process
        */
        for(int i = 0 ; i<Ndim ; i++)
        {
          for(int j = 0 ; j<Ndim ; j++)
          {
            Tangent_Stiffness.nM[A_mask*Ndim+i][B_mask*Ndim+j] += Stiffness_density_p.nM[i][j]*V0_p;
          }
        }

        /*
    Free memory 
        */
        free__TensorLib__(GRADIENT_pB);
        free__MatrixLib__(Stiffness_density_p);
      }

    /*
      Free memory 
    */
    free__TensorLib__(GRADIENT_pA);   
  }
      

      /* 
   Free memory 
      */
      free__TensorLib__(transpose_F_n_p);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);

    }

  return Tangent_Stiffness;
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
  Nodal_Field D_up,
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
    D_up.value.nV[idx_A_i] -= Residual.nV[idx_A_i];
  }

}


/**************************************************************/

static void update_Newmark_Nodal_Increments(
  Nodal_Field D_up,
  Nodal_Field up_n,
  Newmark_parameters Params)
{

  int Nnodes = up_n.value.N_rows;
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
      D_up.d2_value_dt2.nM[A][i] = alpha_1*D_up.value.nM[A][i] - alpha_2*up_n.d_value_dt.nM[A][i] - (alpha_3+1)*up_n.d2_value_dt2.nM[A][i];
      D_up.d_value_dt.nM[A][i]   = alpha_4*D_up.value.nM[A][i] + (alpha_5-1)*up_n.d_value_dt.nM[A][i] + alpha_6*up_n.d2_value_dt2.nM[A][i];
    }
  }
}


/**************************************************************/

static void update_Particles(
  Nodal_Field D_up,
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
  double D_up_pI;
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
        D_up_pI = N_pI*D_up.value.nM[A_mask][i];
        D_dt_upw_pI = N_pI*D_up.d_value_dt.nM[A_mask][i];
        D_dt2_upw_pI = N_pI*D_up.d2_value_dt2.nM[A_mask][i];
                
        if(i<Ndim)
        {
          MPM_Mesh.Phi.acc.nM[p][i]  += D_dt2_upw_pI;
          MPM_Mesh.Phi.vel.nM[p][i]  += D_dt_upw_pI;
          MPM_Mesh.Phi.dis.nM[p][i]  += D_up_pI;
          MPM_Mesh.Phi.x_GC.nM[p][i] += D_up_pI;
        }
        else
        {
          MPM_Mesh.Phi.d2_Pw_dt2.nV[p] += D_dt2_upw_pI;
          MPM_Mesh.Phi.d_Pw_dt_n.nV[p] += D_dt_upw_pI;
          MPM_Mesh.Phi.Pw.nV[p] += D_up_pI;
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
  Nodal_Field D_up,
  Nodal_Field up_n,
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
//  Nodal_Field up_n1;
//  up_n1.value        = allocZ__MatrixLib__(Nnodes_mask,Ndof);
//  up_n1.d_value_dt   = allocZ__MatrixLib__(Nnodes_mask,Ndof);
//  up_n1.d2_value_dt2 = allocZ__MatrixLib__(Nnodes_mask,Ndof);

//  for(int I = 0 ; I<Nnodes_mask*Ndof ; I++)
//  {
//    up_n1.value.nV[I] = D_up.value.nV[I] + up_n.value.nV[I];
//    up_n1.d_value_dt.nV[I] = D_up.d_value_dt.nV[I] + up_n.d_value_dt.nV[I];
//    up_n1.d2_value_dt2.nV[I] = D_up.d2_value_dt2.nV[I] + up_n.d2_value_dt2.nV[I];
//  }

  /*
    vtk results
  */
  if(TimeStep % ResultsTimeStep == 0)
  {
    particle_results_vtk__InOutFun__(MPM_Mesh,TimeStep,ResultsTimeStep);

//    nodal_results_upw_vtk__InOutFun__(up_n1.value, up_n1.d_value_dt, FEM_Mesh, ActiveNodes, TimeStep, ResultsTimeStep);
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
