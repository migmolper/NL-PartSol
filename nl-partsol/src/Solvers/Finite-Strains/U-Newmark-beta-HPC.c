#include "nl-partsol.h"

#include <petscsys.h>
#include <petscmat.h>
#include <petscviewer.h>

/*
  Call global variables
*/
int argc_copy;
char ** argv_copy;
double Thickness_Plain_Stress;
Event * Out_nodal_path_csv;
Event * Out_particles_path_csv;
int Number_Out_nodal_path_csv;
int Number_Out_particles_path_csv;

static char help[] = "Test VecConcatenate both in serial and parallel.\n";

/*
  Define local global variable for the relative error
*/
double Error0;

typedef struct
{

  Matrix value;
  Matrix d_value_dt;
  Matrix d2_value_dt2;

} Nodal_Field;

typedef struct
{
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
static PetscErrorCode Generate_List_Nonzero_Entries(PetscInt *,Particle,Mesh,Mask);
static PetscErrorCode compute_Nodal_Effective_Mass(Mat *,PetscInt *,Particle,Mesh,Mask,double);

/**************************************************************/

void U_Newmark_beta_Finite_Strains_HPC(
  Mesh FEM_Mesh,
  Particle MPM_Mesh,
  Time_Int_Params Parameters_Solver)
{

  /*
    Initializes the PETSc data base and MPI
  */
  PetscErrorCode ierr;
//  PetscInitialize(&argc_copy,&argv_copy,(char*)0,help);
  PetscInitializeNoArguments();

  if(ierr)
  {
    fprintf(stderr,"%s : %s %s \n \t -> %s \n",
    "Error in U_Newmark_beta_Finite_Strains_HPC",
    "The function","  PetscInitializeNoArguments() returned an error message !!!",help);
    exit(EXIT_FAILURE);
  }

  /*
    Auxiliar variables for the solver
  */
  int Ndim = NumberDimensions;
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

  PetscInt * nnz;
  Mat * Effective_MassMatrix;


  Matrix Tangent_Stiffness;
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


  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {
      print_Status("*************************************************",TimeStep);
      print_step(TimeStep,DeltaTimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("First step : Generate Mask ... WORKING",TimeStep);
      /*
	       With the active set of nodes generate a mask to help the algorithm to compute
	       the equilibrium only in the active nodes
      */
      local_search__MeshTools__(MPM_Mesh,FEM_Mesh);
      ActiveNodes = generate_NodalMask__MeshTools__(FEM_Mesh);
      Nactivenodes = ActiveNodes.Nactivenodes;
      Free_and_Restricted_Dofs = generate_Mask_for_static_condensation__MeshTools__(ActiveNodes,FEM_Mesh);
      Generate_List_Nonzero_Entries(nnz,MPM_Mesh,FEM_Mesh,ActiveNodes);
      print_Status("DONE !!!",TimeStep);

//exit(0);

      print_Status("*************************************************",TimeStep);
      print_Status("Second step : Compute effective mass ... WORKING",TimeStep);
      /*
	       Compute the effective mass matrix as a convex combination of the consistent mass
	       matrix and the lumped mass matrix.
      */
     compute_Nodal_Effective_Mass(Effective_MassMatrix,nnz,MPM_Mesh,FEM_Mesh,ActiveNodes,epsilon);
      print_Status("DONE !!!",TimeStep);

exit(0);

//      print_Status("*************************************************",TimeStep);
//      print_Status("Third step : Compute nodal kinetics ... WORKING",TimeStep);
      /*
        Compute the nodal values of velocity and acceleration
      */
//      U_n = compute_Nodal_Field(Effective_Mass,MPM_Mesh,FEM_Mesh,ActiveNodes);
//      D_U = initialise_Nodal_Increments(U_n,FEM_Mesh,ActiveNodes,Params,TimeStep);
//      print_Status("DONE !!!",TimeStep);

//      print_Status("*************************************************",TimeStep);
//      print_Status("Four step : Compute equilibrium ... WORKING",TimeStep);
           
//      Convergence = false;
//      Iter = 0;

//      while(Convergence == false)
//      {

//        update_Local_State(D_U,ActiveNodes,MPM_Mesh,FEM_Mesh,DeltaTimeStep);

//        Forces = compute_Nodal_Forces(ActiveNodes,MPM_Mesh,FEM_Mesh,TimeStep);

//        Reactions = compute_Nodal_Reactions(FEM_Mesh,Forces,ActiveNodes);

//        Residual = compute_Nodal_Residual(U_n,D_U,ActiveNodes,Forces,Effective_Mass,Params);

//        Convergence = check_convergence(Residual,TOL,Iter,MaxIter,TimeStep);

//        if(Convergence == false)
//        {

//	        Tangent_Stiffness = assemble_Nodal_Tangent_Stiffness(ActiveNodes,MPM_Mesh,FEM_Mesh,Params);

//	        if((Free_and_Restricted_Dofs.Nactivenodes - Ndim*Nactivenodes) == 0)
//		      {
//		        solve_non_reducted_system(D_U,Tangent_Stiffness,Effective_Mass,Residual,Params);
//		      }
//	        else
//		      {
//		        solve_reducted_system(D_U,Tangent_Stiffness,Effective_Mass,Residual,Free_and_Restricted_Dofs,Params);
//		      }

//          update_Newmark_Nodal_Increments(D_U,U_n,Params);

//	        Iter++;

//	        free__MatrixLib__(Forces);
//	        free__MatrixLib__(Reactions);
//	        free__MatrixLib__(Residual);
//	        free__MatrixLib__(Tangent_Stiffness);
//        }
//     }
      
//      print_Status("DONE !!!",TimeStep);


//      print_Status("*************************************************",TimeStep);
//      print_Status("Seven step : Update particles lagrangian ... WORKING",TimeStep);

//      update_Particles(D_U,MPM_Mesh,FEM_Mesh,ActiveNodes);
    
//      print_Status("DONE !!!",TimeStep);

      /*
	Outputs
      */
//      output_selector(MPM_Mesh, FEM_Mesh, ActiveNodes, D_U, U_n, Forces, Reactions, Residual,TimeStep, ResultsTimeStep);

    

//      print_Status("*************************************************",TimeStep);
//      print_Status("Eight step : Reset nodal values ... WORKING",TimeStep);
      /*
	Free memory.
      */
//      free__MatrixLib__(Effective_Mass); 
//      free__MatrixLib__(U_n.value);
//      free__MatrixLib__(U_n.d_value_dt);
//      free__MatrixLib__(U_n.d2_value_dt2);
//      free__MatrixLib__(D_U.value);
//      free__MatrixLib__(D_U.d_value_dt);
//      free__MatrixLib__(D_U.d2_value_dt2);
//      free__MatrixLib__(Forces);
//      free__MatrixLib__(Reactions);
//      free__MatrixLib__(Residual);
//      free(ActiveNodes.Nodes2Mask);
//      free(Free_and_Restricted_Dofs.Nodes2Mask);

//      print_Status("DONE !!!",TimeStep);
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

static PetscErrorCode Generate_List_Nonzero_Entries(
  PetscInt *nnz,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes)
{
  PetscErrorCode ierr;

  PetscInt Nnodes_mask = ActiveNodes.Nactivenodes;
  PetscInt Ndof = NumberDOF;
  PetscInt Np = MPM_Mesh.NumGP;
  PetscInt Order = Ndof*Nnodes_mask;
  PetscInt Ap;
  PetscInt Bp;
  PetscInt A_mask;
  PetscInt B_mask;
  PetscInt idx;

  /* Element for each particle */
  Element Nodes_p;

  int * Matrix_Layout = (int *)calloc(Order*Order,sizeof(int));

  PetscMalloc1(Order,&nnz);

  if(ierr)
  {
    fprintf(stderr,"%s : %s %s %s \n",
    "Error in Generate_List_Nonzero_Entries",
    "The function","PetscMalloc1",
    "returned an error message !!!" );
    return ierr;
  }

  /* Iterate over the particles to get the nodal values */
  for(int p = 0 ; p<Np ; p++)
  {

    /* Define tributary nodes of the particle */
    Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

    for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
    {

	    /* Get the node in the mass matrix with the mask */
      Ap = Nodes_p.Connectivity[A];
      A_mask = ActiveNodes.Nodes2Mask[Ap];
      
      for(int B = 0 ; B<Nodes_p.NumberNodes ; B++)
	    {	      
	      /* Get the node in the mass matrix with the mask */
	      Bp = Nodes_p.Connectivity[B];
	      B_mask = ActiveNodes.Nodes2Mask[Bp];

	      /* Fill the effective mass matrix considering the number of dofs */
	      for(int i = 0 ; i<Ndof ; i++)
        {

  	      for(int j = 0 ; j<Ndof ; j++)
           {
            /* Compute the vectorized index */
            idx = (A_mask*Ndof+i)*Order + (B_mask*Ndof+j);

            if(Matrix_Layout[idx] != 1)
            {
              Matrix_Layout[idx] = 1;
            }
          }
        }
	    }
      
    }

    free(Nodes_p.Connectivity);      
    
  }

  for(int i = 0 ; i<Order ; i++)
  {
    for(int j = 0 ; j<Order ; j++)
    {
      nnz[i] += Matrix_Layout[i*Order + j];
    }
  }

  return ierr;
}

/**************************************************************/

static PetscErrorCode compute_Nodal_Effective_Mass(
  Mat * Effective_MassMatrix,
  PetscInt * nnz,
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
  |    0  | M_eff |	         |   0    | M_cons |	      |   0    | M_lump |
*/
{

  PetscInt Nnodes_mask = ActiveNodes.Nactivenodes;
  PetscInt Ndof = NumberDOF;
  PetscInt Np = MPM_Mesh.NumGP;
  PetscInt Order = Ndof*Nnodes_mask;
  PetscInt Ap;
  PetscInt Bp;
  PetscInt A_mask;
  PetscInt B_mask;
  PetscInt idx_A_mask_i;
  PetscInt idx_B_mask_i;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;  

  /* Evaluation of the particle in the node */
  PetscScalar ShapeFunction_pA, ShapeFunction_pB;
  /* Mass of the particle */
  PetscScalar m_p;
  /* Nodal contribution A of the particle p */
  PetscScalar m_A_p;
  /* Nodal contribution A-B of the particle p */
  PetscScalar m_AB_p;
  /* Element for each particle */
  Element Nodes_p;

  /* Allocate the effective mass matrix */
  MatCreateSeqAIJ(PETSC_COMM_SELF,Order,Order,Order,nnz,Effective_MassMatrix);

  /* Iterate over the particles to get the nodal values */
  for(PetscInt p = 0 ; p<Np ; p++)
  {

      /* Define tributary nodes of the particle */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /* Evaluate the shape function in the coordinates of the particle */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Get the mass of the particle */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      for(PetscInt A = 0 ; A<Nodes_p.NumberNodes ; A++)
      {
        /* Get the node in the mass matrix with the mask */
        Ap = Nodes_p.Connectivity[A];
        A_mask = ActiveNodes.Nodes2Mask[Ap];
        
        /* Get the value of the shape function */
        ShapeFunction_pA = ShapeFunction_p.nV[A];
        
        /* Compute the nodal A contribution of the particle p */
        m_A_p = m_p*ShapeFunction_pA;
        
        for(PetscInt  B = 0 ; B<Nodes_p.NumberNodes ; B++)
        {	      
          /* Get the node in the mass matrix with the mask */
          Bp = Nodes_p.Connectivity[B];
          B_mask = ActiveNodes.Nodes2Mask[Bp];

	      /* Get the value of the shape function */
	      ShapeFunction_pB = ShapeFunction_p.nV[B];

	      /* Compute the nodal AB contribution of the particle p */
	      m_AB_p = m_p*ShapeFunction_pA*ShapeFunction_pB;

	      /* Fill the effective mass matrix considering the number of dofs */
	      for(PetscInt i = 0 ; i<Ndof ; i++)
        {
          idx_A_mask_i = A_mask*Ndof+i;
          idx_B_mask_i = B_mask*Ndof+i;
          
          MatSetValues(*Effective_MassMatrix,1,&idx_A_mask_i,1,&idx_B_mask_i,&m_AB_p,ADD_VALUES);

        }

	    }
	}

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);      

    }

    /*
      Start assembling process
    */
   MatAssemblyBegin(*Effective_MassMatrix,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(*Effective_MassMatrix,MAT_FINAL_ASSEMBLY);

   MatView(*Effective_MassMatrix,PETSC_VIEWER_STDOUT_SELF);

  /*
    At this point the effective mass matrix coincides with the consistent mass
    matrix. We can tune it by a convecx combination with the lumped mass matrix
  */
 /*
  for(PetscInt A = 0 ; A<Order ; A++)
    {
      for(PetscInt B = 0 ; B<Order ; B++)
    	{    
	     Effective_MassMatrix.nM[A][B] = (1-epsilon)*Effective_MassMatrix.nM[A][B] + (A == B)*epsilon*Lumped_MassMatrix.nV[A];
	     }
    }
*/

  /*
    Free lumped mass matrix.
  */
 /*
  free__MatrixLib__(Lumped_MassMatrix);
  */

}

/**************************************************************/


