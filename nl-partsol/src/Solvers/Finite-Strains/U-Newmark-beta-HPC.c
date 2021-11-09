#include "nl-partsol.h"

//#include <petscsys.h>

#include <petscmat.h>

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

/**************************************************************/

void U_Newmark_beta_Finite_Strains_HPC(
  Mesh FEM_Mesh,
  Particle MPM_Mesh,
  Time_Int_Params Parameters_Solver)
{

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

  Matrix Effective_Mass;
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
      print_Status("DONE !!!",TimeStep);

exit(0);

//      print_Status("*************************************************",TimeStep);
//      print_Status("Second step : Compute effective mass ... WORKING",TimeStep);
      /*
	       Compute the effective mass matrix as a convex combination of the consistent mass
	       matrix and the lumped mass matrix.
      */
//      Effective_Mass = compute_Nodal_Effective_Mass(MPM_Mesh,FEM_Mesh,ActiveNodes,epsilon);
//      print_Status("DONE !!!",TimeStep);

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
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ndof = NumberDOF;
  int Np = MPM_Mesh.NumGP;
  int Order = Ndof*Nnodes_mask;
  int Ap;
  int Bp;
  int A_mask;
  int B_mask;
  int idx;

  /* Element for each particle */
  Element Nodes_p;

  int * Matrix_Layout = (int *)calloc(Order*Order,sizeof(int));

  PetscMalloc1(Order,&nnz);

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
            idx = (A_mask*Ndof+i)*Nnodes_mask + B_mask*Ndof + j;

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

}



/**************************************************************/