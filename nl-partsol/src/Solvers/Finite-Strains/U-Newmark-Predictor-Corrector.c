#include "nl-partsol.h"

#ifdef _OPENMP
    #include <omp.h>
    #define UD_Num_THREADS 10
#else
    #define omp_get_thread_num() 0
#endif

#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/*
  Call global variables
*/
double DeltaTimeStep;
double Thickness_Plain_Stress;
Event * Out_nodal_path_csv;
Event * Out_particles_path_csv;
int Number_Out_nodal_path_csv;
int Number_Out_particles_path_csv;

/*
  Auxiliar functions 
*/

/* Step 1 */
static int compute_Mass_Matrix(double *,Particle,Mesh,Mask);
/* Step 2 */
static int compute_Explicit_Newmark_Predictor(Particle,double);
/* Step 3 */
static int compute_Nodal_Gravity_field(double *, Mask, Particle, int);
static int compute_Nodal_D_Displacement(double *,double *,Particle,Mesh,Mask);
static int compute_Nodal_Velocity(double *,double *,Particle,Mesh,Mask);
static int impose_Dirichlet_Boundary_Conditions(double *, double *,Mesh,Mask,double,int);
/* Step 4 */
static int update_Local_State(double *,Mask,Particle,Mesh);
/* Step 5 */
static int compute_Nodal_Forces(double *,Mask,Particle,Mesh,int);
static void   compute_Nodal_Internal_Forces(Matrix,Mask,Particle,Mesh);
static void   compute_Nodal_Nominal_traction_Forces(Matrix,Mask,Particle,Mesh,int);
static int solve_Nodal_Equilibrium(double *,double *,double *,double *,double *,Particle,Mesh,Mask,Mask);
/* Step 6 */
static int compute_Explicit_Newmark_Corrector(Particle,double);
/* Step 7 */
static int output_selector(Particle, Mesh, Mask, double *, double *, double *, double *, double, int, int);

/**************************************************************/

int U_Newmark_Predictor_Corrector_Finite_Strains(
  Mesh FEM_Mesh,
  Particle MPM_Mesh,
  Time_Int_Params Parameters_Solver)
{

  int status = 0;

  /*
    Auxiliar variables for the solver
  */
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
  double * Velocity;
  double * Forces;
  double * Reactions;

  Mask ActiveNodes;
  Mask Free_and_Restricted_Dofs;


  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {

      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = U_DeltaT__SolversLib__(MPM_Mesh, DeltaX, Parameters_Solver);
      print_step(TimeStep,DeltaTimeStep);
      local_search__MeshTools__(MPM_Mesh,FEM_Mesh);
      ActiveNodes = generate_NodalMask__MeshTools__(FEM_Mesh);
      Nactivenodes = ActiveNodes.Nactivenodes;
      Free_and_Restricted_Dofs = generate_Mask_for_static_condensation__MeshTools__(ActiveNodes,FEM_Mesh);

      /* Define memory */
      Lumped_Mass = (double *)calloc(Nactivenodes*Ndim,sizeof(double));
      Gravity_field = (double *)calloc(Nactivenodes*Ndim,sizeof(double));
      D_Displacement = (double *)calloc(Nactivenodes*Ndim,sizeof(double));
      Velocity = (double *)calloc(Nactivenodes*Ndim,sizeof(double));
      Forces = (double *)calloc(Nactivenodes*Ndim,sizeof(double));
      Reactions = (double *)calloc(Nactivenodes*Ndim,sizeof(double));

      status = compute_Mass_Matrix(Lumped_Mass,MPM_Mesh,FEM_Mesh,ActiveNodes);
      if(status)
      {
        fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
        "Error in the line",__LINE__,"of the file",__FILE__, 
        "The function","compute_Mass_Matrix",
        "returned an error message !!!" );
        return EXIT_FAILURE;
      }
      
      status = compute_Explicit_Newmark_Predictor(MPM_Mesh, gamma);
      if(status)
      {
        fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
        "Error in the line",__LINE__,"of the file",__FILE__, 
        "The function","compute_Explicit_Newmark_Predictor",
        "returned an error message !!!" );
        return EXIT_FAILURE;
      }

      status = compute_Nodal_Gravity_field(Gravity_field, ActiveNodes, MPM_Mesh, TimeStep);
      if(status)
      {
        fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
        "Error in the line",__LINE__,"of the file",__FILE__, 
        "The function","compute_Nodal_Gravity_field",
        "returned an error message !!!" );
        return EXIT_FAILURE;
      }

      status = compute_Nodal_D_Displacement(D_Displacement,Lumped_Mass,MPM_Mesh,FEM_Mesh,ActiveNodes);
      if(status)
      {
        fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
        "Error in the line",__LINE__,"of the file",__FILE__, 
        "The function","compute_Nodal_D_Displacement",
        "returned an error message !!!" );
        return EXIT_FAILURE;
      }

      status = compute_Nodal_Velocity(Velocity,Lumped_Mass,MPM_Mesh,FEM_Mesh,ActiveNodes);
      if(status)
      {
        fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
        "Error in the line",__LINE__,"of the file",__FILE__, 
        "The function","compute_Nodal_Velocity",
        "returned an error message !!!" );
        return EXIT_FAILURE;
      }

      status = impose_Dirichlet_Boundary_Conditions(D_Displacement,Velocity,FEM_Mesh,ActiveNodes,DeltaTimeStep,TimeStep);
      if(status)
      {
        fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
        "Error in the line",__LINE__,"of the file",__FILE__, 
        "The function","impose_Dirichlet_Boundary_Conditions",
        "returned an error message !!!" );
        return EXIT_FAILURE;
      }

      status = update_Local_State(D_Displacement, ActiveNodes, MPM_Mesh, FEM_Mesh);
      if(status)
      {
        fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
        "Error in the line",__LINE__,"of the file",__FILE__, 
        "The function","update_Local_State",
        "returned an error message !!!" );
        return EXIT_FAILURE;
      }

      status = compute_Nodal_Forces(Forces,ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep);
      if(status)
      {
        fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
        "Error in the line",__LINE__,"of the file",__FILE__, 
        "The function","compute_Nodal_Forces",
        "returned an error message !!!" );
        return EXIT_FAILURE;
      }

      status = solve_Nodal_Equilibrium(Lumped_Mass,Gravity_field,Forces,Reactions,D_Displacement,MPM_Mesh,FEM_Mesh,ActiveNodes,Free_and_Restricted_Dofs);
      if(status)
      {
        fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
        "Error in the line",__LINE__,"of the file",__FILE__, 
        "The function","solve_Nodal_Equilibrium",
        "returned an error message !!!" );
        return EXIT_FAILURE;
      }
        
      status = compute_Explicit_Newmark_Corrector(MPM_Mesh,gamma);
      if(status)
      {
        fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
        "Error in the line",__LINE__,"of the file",__FILE__, 
        "The function","compute_Explicit_Newmark_Corrector",
        "returned an error message !!!" );
        return EXIT_FAILURE;
      }

      output_selector(MPM_Mesh, FEM_Mesh, ActiveNodes, Velocity, D_Displacement,Forces, Reactions, DeltaTimeStep, TimeStep, ResultsTimeStep);
      if(status)
      {
        fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
        "Error in the line",__LINE__,"of the file",__FILE__, 
        "The function","output_selector",
        "returned an error message !!!" );
        return EXIT_FAILURE; 
      }


      free(Lumped_Mass); 
      free(Gravity_field);
      free(D_Displacement);
      free(Velocity);
      free(Forces);
      free(Reactions);
      free(ActiveNodes.Nodes2Mask); 
      free(Free_and_Restricted_Dofs.Nodes2Mask);

    }
  
  return EXIT_SUCCESS;
}

/**************************************************************/

static int compute_Mass_Matrix(
  double * M_IJ,
  Particle MPM_Mesh, // Variable with information of the particles 
  Mesh FEM_Mesh, // Variable with information of the nodes
  Mask ActiveNodes) // Variable with information of the active nodes
/*
  This function computes the lumped mass matrix with the displacements degree of freedom
*/
{

unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
unsigned Ndim = NumberDimensions;
unsigned Np = MPM_Mesh.NumGP;

#pragma omp parallel 
{
  /* Particle variables */
  Element Nodes_p;
  Matrix ShapeFunction_p;
  double ShapeFunction_pA;
  double m_p;
  unsigned Idx_p;

#ifdef _OPENMP

  unsigned Machine_threads = omp_get_max_threads();
  unsigned Threads = IMIN(UD_Num_THREADS,Machine_threads);

  // Additional work to set the number of threads.
	// We hard-code to 4 for illustration purposes only.
	omp_set_num_threads(Threads);

	// determine how many elements each process will work on
	unsigned n_per_thread = Np/Threads;

#endif

  #pragma omp parallel for shared (M_IJ) schedule(static,n_per_thread)
  for(Idx_p = 0 ; Idx_p<Np ; Idx_p++)
    {

      /*
        Define tributary nodes of the particle 
      */
      Nodes_p = nodal_set__Particles__(Idx_p, MPM_Mesh.ListNodes[Idx_p], MPM_Mesh.NumberNodes[Idx_p]);

      /* 
        Evaluate the shape function in the coordinates of the particle
      */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
       Get the mass of the particle 
      */
      m_p = MPM_Mesh.Phi.mass.nV[Idx_p];


      for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
        {

          /* 
             Get the node in the mass matrix with the mask
          */
          unsigned Ap = Nodes_p.Connectivity[A];
          unsigned A_mask = ActiveNodes.Nodes2Mask[Ap];

          /* 
            Get the value of the shape function 
          */
          ShapeFunction_pA = ShapeFunction_p.nV[A];

          /*
            Compute the nodal A contribution of the particle p
          */
          double m_A_p = m_p*ShapeFunction_pA;

          /* 
             Fill the Lumped mass matrix considering the number of dofs
          */
          for(int i = 0 ; i<Ndim ; i++)
            {
              M_IJ[A_mask*Ndim + i] += m_A_p;       
            }

          }

        /* 
          Free the value of the shape functions 
        */
        free__MatrixLib__(ShapeFunction_p);
        free(Nodes_p.Connectivity);   

    }


}
  
  return EXIT_SUCCESS; 
}

/**************************************************************/

static int compute_Explicit_Newmark_Predictor(
  Particle MPM_Mesh,
  double gamma)
/*
  The predictor stage is computed in the particles
*/
{

  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Order = Ndim*Np;

#pragma omp parallel 
{
  /* Particle variables */
  unsigned idx_p;
  double * D_dis_p = MPM_Mesh.Phi.D_dis.nV;
  double * Vel_p = MPM_Mesh.Phi.vel.nV;
  double * Acc_p = MPM_Mesh.Phi.acc.nV;

#ifdef _OPENMP

  unsigned Machine_threads = omp_get_max_threads();
  unsigned Threads = IMIN(UD_Num_THREADS,Machine_threads);

  // Additional work to set the number of threads.
	// We hard-code to 4 for illustration purposes only.
	omp_set_num_threads(Threads);

	// determine how many elements each process will work on
	unsigned n_per_thread = Order/Threads;

#endif

  #pragma omp parallel for shared (D_dis_p,Vel_p,Acc_p) private(idx_p) schedule(static,n_per_thread)
  for(idx_p = 0 ; idx_p<Order ; idx_p++)
  {
    D_dis_p[idx_p] = DeltaTimeStep*Vel_p[idx_p] + 0.5*DSQR(DeltaTimeStep)*Acc_p[idx_p];

    Vel_p[idx_p] += (1-gamma)*DeltaTimeStep*Acc_p[idx_p];
  }

}

  return EXIT_SUCCESS; 
}

/**************************************************************/

static int compute_Nodal_Gravity_field(
  double * Gravity_field,
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
  Tensor b = alloc__TensorLib__(1);
  

  for(int i = 0 ; i<NumBodyForces ; i++)
    {

     /* Fill vector b of body acclerations */
     for(int k = 0 ; k<Ndim ; k++)
     {
       if(B[i].Dir[k])
       {
         if( (TimeStep < 0) || (TimeStep > B[i].Value[k].Num))
         {
            printf("%s : %s\n", "Error in compute_Nodal_Gravity_field()","The time step is out of the curve !!");
            exit(EXIT_FAILURE);
         }
         
         b.n[k] += B[i].Value[k].Fx[TimeStep];

       }
     }
      
    }

  /* Get the node of the mesh for the contribution */
     for(int A = 0 ; A<Nnodes_mask ; A++)
    {
      /* Compute body forces */
      for(int k = 0 ; k<Ndim ; k++)
        {
          Gravity_field[A*Ndim + k] = b.n[k];
        }  
    }

    free__TensorLib__(b);


  return EXIT_SUCCESS; 
}

/**************************************************************/

static int compute_Nodal_D_Displacement(
  double * D_dis_I,
  double * Mass_IJ,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes)
/*
  Compute the nodal increment of displacement. The operation is linearized and
  all the dof split the increment of displacement array in n components like :
  | M 0 |   |D_u.x|   | M*D_u.x |
  | 0 M | * |D_u.y| = | M*D_u.y |
*/
{
  unsigned Ndim = NumberDimensions;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned Order = Nnodes_mask*Ndim;
  unsigned Ap;
  unsigned A_mask;
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Value of the shape-function */
  double ShapeFunction_pA; /* Evaluation of the particle in the node */
  double m_p; /* Mass of the particle */

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
          for(int i = 0 ; i<Ndim ; i++)
            {
              D_dis_I[A_mask*Ndim + i] += m_p*ShapeFunction_pA*MPM_Mesh.Phi.D_dis.nM[p][i];
            }
        }

      /* 
        Free the value of the shape functions 
      */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }

#pragma omp parallel 
{

  /*
    Compute the D_Displacements
  */
#ifdef _OPENMP

  unsigned Machine_threads = omp_get_max_threads();
  unsigned Threads = IMIN(UD_Num_THREADS,Machine_threads);

  // Additional work to set the number of threads.
	// We hard-code to 4 for illustration purposes only.
	omp_set_num_threads(Threads);

	// determine how many elements each process will work on
	unsigned n_per_thread = Order/Threads;

#endif

  unsigned idx;

  #pragma omp parallel for shared (D_dis_I,Mass_IJ,Order) private(idx) schedule(static,n_per_thread)
  for(idx = 0 ; idx<Order ; idx++)
  {
    D_dis_I[idx] = D_dis_I[idx]/Mass_IJ[idx];
  }

}

 return EXIT_SUCCESS; 
}


/**************************************************************/

static int compute_Nodal_Velocity(
  double * Vel_I,
  double * Mass_IJ,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes)
/*
  Compute the nodal velocity. The operation is linearized and
  all the dof split the velocity array in n components like :
  | M 0 |   |V.x|   | p.x |
  | 0 M | * |V.y| = | p.y |
*/
{
  unsigned Ndim = NumberDimensions;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned Order = Ndim*Nnodes_mask;
  unsigned Ap;
  unsigned A_mask;
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Value of the shape-function */
  double ShapeFunction_pA; /* Evaluation of the particle in the node */
  double m_p; /* Mass of the particle */

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
        Get the nodal mommentum 
      */
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


          for(int i = 0 ; i<Ndim ; i++)
            {
              Vel_I[A_mask*Ndim + i] += m_p*ShapeFunction_pA*MPM_Mesh.Phi.vel.nM[p][i];
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
 #ifdef _OPENMP

  unsigned Machine_threads = omp_get_max_threads();
  unsigned Threads = IMIN(UD_Num_THREADS,Machine_threads);

  // Additional work to set the number of threads.
	// We hard-code to 4 for illustration purposes only.
	omp_set_num_threads(Threads);

	// determine how many elements each process will work on
	unsigned n_per_thread = Order/Threads;

#endif

  unsigned idx;

  #pragma omp parallel for shared (Vel_I,Mass_IJ,Order) private(idx) schedule(static,n_per_thread)
  for(idx = 0 ; idx<Order ; idx++)
  {
    Vel_I[idx] = Vel_I[idx]/Mass_IJ[idx];
  }

  return EXIT_SUCCESS;
}


/**************************************************************/

static int impose_Dirichlet_Boundary_Conditions(
  double * D_dis_I,
  double * Vel_I,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  double Time_Step,
  int TimeStep)
/*
  Apply the boundary conditions over the nodes 
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
            return EXIT_FAILURE;
          }

          /*
            Initialise increments using newmark and the value of the boundary condition
          */
          D_dis_I[Id_BCC_mask*Ndim + k] = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*(double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];
          Vel_I[Id_BCC_mask*Ndim + k] = D_dis_I[Id_BCC_mask*Ndim + k]/Time_Step;  

        }
      }
    }    
  }

  return EXIT_SUCCESS;

}

/**************************************************************/

static int update_Local_State(
  double * D_dis_I,
	Mask ActiveNodes,
	Particle MPM_Mesh,
	Mesh FEM_Mesh)
{

  /*
    Auxiliar variables
  */
  int status = 0;
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Nnodes_p;
  int MatIndx_p;
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
  Matrix D_Displacement_Ap;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor DF_p;

  Matrix D_Displacement = memory_to_matrix__MatrixLib__(Nnodes_mask,Ndim,D_dis_I);

  /*
    Loop in the material point set to update kinematics
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
    D_Displacement_Ap = get_set_field__MeshTools__(D_Displacement, Nodes_p, ActiveNodes);

    /*
      Evaluate the shape function gradient in the coordinates of the particle
    */
    gradient_p = compute_dN__MeshTools__(Nodes_p,MPM_Mesh,FEM_Mesh);
	  
    /*
      Take the values of the deformation gradient from the previous step
    */
    F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    DF_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.DF.nM[p],2);

    /*
      Compute the increment of the deformation gradient
    */
    update_increment_Deformation_Gradient__Particles__(DF_p,D_Displacement_Ap,gradient_p);

    /*
      Update the deformation gradient in t = n + 1 with the information
      from t = n and the increment of deformation gradient.
    */  
    update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, DF_p);

    /*
      Compute Jacobian of the deformation gradient
    */
    MPM_Mesh.Phi.J_n1.nV[p] = I3__TensorLib__(F_n1_p);

    /*
      Check non-pentrability condition
    */
    if(MPM_Mesh.Phi.J_n1.nV[p] <= 0.0)
    {
      fprintf(stderr,"%s : %s %i\n",
        "Error in update_Local_State()",
        "Negative jacobian in particle",p);
      return EXIT_FAILURE;
    }

    /*
      Update patch
    */
    if(FEM_Mesh.Locking_Control_Fbar)
    {
      Idx_Element_p = MPM_Mesh.Element_p[p];
      Idx_Patch_p = FEM_Mesh.Idx_Patch[Idx_Element_p];
      FEM_Mesh.Vol_Patch_n[Idx_Patch_p] += MPM_Mesh.Phi.J_n.nV[p]*MPM_Mesh.Phi.Vol_0.nV[p];
      FEM_Mesh.Vol_Patch_n1[Idx_Patch_p] += MPM_Mesh.Phi.J_n1.nV[p]*MPM_Mesh.Phi.Vol_0.nV[p];
    }

    /*
      Update density with the jacobian of the increment deformation gradient
    */
    Delta_J_p = I3__TensorLib__(DF_p);
    rho_n_p = MPM_Mesh.Phi.rho.nV[p];
    MPM_Mesh.Phi.rho.nV[p] = rho_n_p/Delta_J_p;

    /*
      Free memory 
    */
    free__MatrixLib__(D_Displacement_Ap);
    free__MatrixLib__(gradient_p);
    free(Nodes_p.Connectivity);
	  
  }

  /*
    Loop in the material point set to update stress
  */
  for(int p = 0 ; p<Np ; p++)
  {
    MatIndx_p = MPM_Mesh.MatIdx[p];
    MatProp_p = MPM_Mesh.Mat[MatIndx_p];
    
    if(FEM_Mesh.Locking_Control_Fbar)
    {
      Idx_Element_p = MPM_Mesh.Element_p[p];
      Idx_Patch_p = FEM_Mesh.Idx_Patch[Idx_Element_p];

      Vn_patch = FEM_Mesh.Vol_Patch_n[Idx_Patch_p];
      Vn1_patch = FEM_Mesh.Vol_Patch_n1[Idx_Patch_p];
      J_patch = Vn1_patch/Vn_patch;

      get_locking_free_Deformation_Gradient_n1__Particles__(p,J_patch,MPM_Mesh);

      MPM_Mesh.Phi.Jbar.nV[p] *= J_patch;
    }

    /*
      Update the first Piola-Kirchhoff stress tensor with an apropiate
      integration rule.
    */
    status = Stress_integration__Particles__(p,MPM_Mesh,FEM_Mesh,MatProp_p); 
    if(status)
    {
      fprintf(stderr,"%s %i %s %s : \n\t %s %s %s \n",
      "Error in the line",__LINE__,"of the file",__FILE__, 
      "The function","Stress_integration__Particles__",
      "returned an error message !!!" );
      return EXIT_FAILURE;
    }

  }

  return EXIT_SUCCESS;

}

/**************************************************************/

static int compute_Nodal_Forces(
  double * F_I,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  int TimeStep)
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;

  Matrix Forces = memory_to_matrix__MatrixLib__(Nnodes_mask,Ndim,F_I);

  /*
    Add internal forces contribution
  */
  compute_Nodal_Internal_Forces(Forces,ActiveNodes,MPM_Mesh,FEM_Mesh);

  /*
    Add contact forces contribution
  */
  compute_Nodal_Nominal_traction_Forces(Forces,ActiveNodes,MPM_Mesh,FEM_Mesh,TimeStep);

  
  return EXIT_SUCCESS;
}


/**************************************************************/

static void compute_Nodal_Internal_Forces(
  Matrix Forces,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int NumNodes_p;
  int idx_A_mask_i;

  Tensor P_p; /* First Piola-Kirchhoff Stress tensor */
  Tensor InternalForcesDensity_Ap;

  Element Nodes_p; /* List of nodes for particle */
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pA;
  Tensor GRADIENT_pA;
  Tensor F_n_p;
  Tensor transpose_F_n_p;
  
  double V0_p; /* Volume of the Gauss-Point */

  /*
    Loop in the particles 
  */
  for(int p = 0 ; p<Np ; p++)
  {
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
    F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    transpose_F_n_p = transpose__TensorLib__(F_n_p);

    /*
      Get the first Piola-Kirchhoff stress tensor.
    */
    Tensor P_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2); 

    for(int A = 0 ; A<NumNodes_p ; A++)
    {

      /*
        Compute the gradient in the reference configuration 
      */
      gradient_pA = memory_to_tensor__TensorLib__(gradient_p.nM[A], 1);
      GRADIENT_pA = vector_linear_mapping__TensorLib__(transpose_F_n_p,gradient_pA);
      
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
	      Forces.nM[A_mask][i] -= InternalForcesDensity_Ap.n[i]*V0_p;
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


static void compute_Nodal_Nominal_traction_Forces(
  Matrix Forces,
  Mask ActiveNodes,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  int TimeStep)
{

  int Ndim = NumberDimensions;
  Load Load_i;
  Element Nodes_p; /* Element for each Gauss-Point */
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
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /*
        Fill vector of contact forces
      */
      for(int k = 0 ; k<Ndim ; k++)
      {
        if(Load_i.Dir[k] == 1)
        {
          if((TimeStep < 0) || (TimeStep > Load_i.Value[k].Num))
          {
            printf("%s : %s",
              "Error in compute_Nodal_Nominal_traction_Forces()",
              "The time step is out of the curve !!");
            exit(EXIT_FAILURE);
          }

          T.n[k] = Load_i.Value[k].Fx[TimeStep];
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
          Forces.nM[A_mask][k] += ShapeFunction_pA*T.n[k]*A0_p;
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

static int solve_Nodal_Equilibrium(
  double * Mass_IJ,
  double * G_I,
  double * F_I,
  double * R_I,
  double * D_dis_I,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
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
  unsigned Ndim = NumberDimensions;
  unsigned Ndof = NumberDOF;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned Nnodes = ActiveNodes.Nactivenodes;
  unsigned NumNodes_p;
  unsigned Order = Nnodes*Ndim;
  unsigned AB;
  unsigned Ap;
  unsigned A_mask;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix ShapeFunction_p;
  double ShapeFunction_pA;

  /* 
    Solution nodal variable
  */
  double * Acc_I = (double *)calloc(Nnodes*Ndim,sizeof(double));

  /*
    The solution is now stored in the internal forces vector
  */
 #ifdef _OPENMP

  unsigned Machine_threads = omp_get_max_threads();
  unsigned Threads = IMIN(UD_Num_THREADS,Machine_threads);

  // Additional work to set the number of threads.
	// We hard-code to 4 for illustration purposes only.
	omp_set_num_threads(Threads);

	// determine how many elements each process will work on
	unsigned n_per_thread = Order/Threads;

#endif

  unsigned idx;
  int * Nodes2Mask = Free_and_Restricted_Dofs.Nodes2Mask;

  #pragma omp parallel for shared (R_I,F_I,Acc_I,G_I,Mass_IJ,Order) private(idx) schedule(static,n_per_thread)
  for(idx = 0 ; idx<Order ; idx++)
  {
    if(Nodes2Mask[idx] != -1)
    {
      Acc_I[idx] = G_I[idx] + F_I[idx]/Mass_IJ[idx];
    }
    else
    {
      Acc_I[idx] = 0.0;
      R_I[idx] = F_I[idx];
    }
  }


  /*
    Update particle acceleration
  */
 double * Acc_p = MPM_Mesh.Phi.acc.nV;
 double * D_dis_p = MPM_Mesh.Phi.D_dis.nV;

  for(int p = 0 ; p<Np ; p++)
  {

    /* Set to zero particle acceleration for interpolation */
    for(int i = 0 ; i<Ndim ; i++)
    {
      MPM_Mesh.Phi.acc.nM[p][i] = 0.0;
      MPM_Mesh.Phi.D_dis.nM[p][i] = 0.0;
    }

    /* Define element of the particle */
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

    /* Evaluate the shape function in the coordinates of the particle */
    ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    /* Interpolate the new value of the acceleration */
    for(int A = 0; A<NumNodes_p; A++)
    {
      /* Get the node in the nodal momentum with the mask */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];
    
      /* Evaluate the GP function in the node */
      ShapeFunction_pA = ShapeFunction_p.nV[A];
      
      for(int i = 0 ; i<Ndim ; i++)
       {
         Acc_p[p*Ndim + i]   += ShapeFunction_pA*Acc_I[A_mask*Ndim + i];
         D_dis_p[p*Ndim + i] += ShapeFunction_pA*D_dis_I[A_mask*Ndim + i];
       }
    }

    /* Free auxiliar data */
    free__MatrixLib__(ShapeFunction_p);
    free(Nodes_p.Connectivity);

  }

  /* Free acceleration vector */
  free(Acc_I);

  return EXIT_SUCCESS;
}

/**************************************************************/

static int compute_Explicit_Newmark_Corrector(
  Particle MPM_Mesh,
  double gamma)
{
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  Tensor F_n_p;
  Tensor F_n1_p;

  for(int p = 0 ; p<Np ; p++)
  {
      
    /*
      Replace the deformation gradient at t = n with the new one
    */
    F_n_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);

    if(MPM_Mesh.Mat[MPM_Mesh.MatIdx[p]].Locking_Control_Fbar)
    {
      F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Fbar.nM[p],2);
    }
    else
    {
      F_n1_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    }


    /*
      Replace the determinant of the deformation gradient
    */
    MPM_Mesh.Phi.J_n.nV[p] = MPM_Mesh.Phi.J_n1.nV[p];

    /* 
      Update/correct tensor and vector variables
    */
    for(int i = 0 ; i<Ndim  ; i++)
    {
      /* 
        Correct particle velocity
      */
      MPM_Mesh.Phi.vel.nM[p][i] += gamma*DeltaTimeStep*MPM_Mesh.Phi.acc.nM[p][i];

      /*
        Update the particles position and displacement
      */
      MPM_Mesh.Phi.x_GC.nM[p][i] += MPM_Mesh.Phi.D_dis.nM[p][i];
      MPM_Mesh.Phi.dis.nM[p][i]  += MPM_Mesh.Phi.D_dis.nM[p][i];

      /* Update deformation gradient tensor */
      for(int j = 0 ; j<Ndim  ; j++)
      {
        F_n_p.N[i][j] = F_n1_p.N[i][j];
      }
    }

  }

  return EXIT_SUCCESS;

}  

/**************************************************************/

static int output_selector(
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  double * Velocity,
  double * D_Displacement,
  double * F_I,
  double * R_I,
  double DeltaTimeStep,
  int TimeStep,
  int ResultsTimeStep)
{


  /*
    vtk results
  */
  if(TimeStep % ResultsTimeStep == 0)
  {

    int Nnodes = ActiveNodes.Nactivenodes;
    int p_idx = 2075;
    int NumNodes_p = MPM_Mesh.NumberNodes[p_idx];
    Element Nodes_p = nodal_set__Particles__(p_idx, MPM_Mesh.ListNodes[p_idx], NumNodes_p);
    Matrix ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    Matrix ShapeFunction_Ip = allocZ__MatrixLib__(Nnodes,1);

    for(int A = 0; A<NumNodes_p; A++)
    {
      ShapeFunction_Ip.nV[ActiveNodes.Nodes2Mask[Nodes_p.Connectivity[A]]] = ShapeFunction_p.nV[A];
    }

    Matrix Reactions = memory_to_matrix__MatrixLib__(Nnodes,NumberDimensions,R_I);

    particle_results_vtk__InOutFun__(MPM_Mesh,TimeStep,ResultsTimeStep);

    nodal_results_vtk__InOutFun__(FEM_Mesh, ActiveNodes, Reactions, ShapeFunction_Ip, TimeStep, ResultsTimeStep);

    free(Nodes_p.Connectivity);
    free__MatrixLib__(ShapeFunction_p);
    free__MatrixLib__(ShapeFunction_Ip);

  }

  /* 
    csv results 
  */

/*
  for(int i = 0 ; i<Number_Out_nodal_path_csv ; i++)
  {

    if(Out_nodal_path_csv[i].Out_csv_nodes_path_Velocity)
    {
      path_nodes_analysis_csv__InOutFun__(Velocity, FEM_Mesh.Coordinates,"Nodal_path_velocity_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
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

  }

*/

/*
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

  */


  return EXIT_SUCCESS;

}

/**************************************************************/
