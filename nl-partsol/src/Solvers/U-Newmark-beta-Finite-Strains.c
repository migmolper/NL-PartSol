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
static Matrix compute_Nodal_Effective_Mass(Particle, Mesh, Mask, double);
static Nodal_Field compute_Nodal_Field(Matrix,Particle, Mesh, Mask);
static Nodal_Field initialise_Nodal_Increments(Nodal_Field,Mesh,Mask,Newmark_parameters,int);
static void   update_Local_State(Nodal_Field,Mask,Particle,Mesh,double);
static Matrix compute_Nodal_Forces(Mask, Particle, Mesh, int);
static void   compute_Nodal_Internal_Forces(Matrix,Mask,Particle, Mesh);
static void   compute_Nodal_Nominal_traction_Forces(Matrix,Mask,Particle,Mesh,int);
static void   compute_Nodal_Body_Forces(Matrix,Mask,Particle,Mesh, int);
static Matrix compute_Nodal_Reactions(Mesh,Matrix,Mask);
static Matrix compute_Nodal_Residual(Nodal_Field,Nodal_Field,Mask,Matrix, Matrix, Newmark_parameters);
static bool   check_convergence(Matrix,double,int,int,int);
static Matrix assemble_Nodal_Tangent_Stiffness(Mask,Particle,Mesh,Newmark_parameters);
static void   solve_non_reducted_system(Nodal_Field,Matrix, Matrix, Matrix, Newmark_parameters);
static void   solve_reducted_system(Nodal_Field,Matrix,Matrix,Matrix,Mask,Newmark_parameters);
static void   update_Newmark_Nodal_Increments(Nodal_Field,Nodal_Field,Newmark_parameters);
static void   update_Particles(Nodal_Field,Particle,Mesh,Mask);
static void output_selector(Particle,Mesh,Mask,Nodal_Field,Nodal_Field,Matrix,Matrix,Matrix,int,int);
/**************************************************************/

void U_Newmark_beta_Finite_Strains(
  Mesh FEM_Mesh,
  Particle MPM_Mesh,
  int InitialStep)
{

  /*
    Integer variables 
  */
  int Ndim = NumberDimensions;
  int Nactivenodes;
  /*
    Auxiliar variables for the solver
  */
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
  DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
 
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
      ActiveNodes = generate_NodalMask__MeshTools__(FEM_Mesh);
      Nactivenodes = ActiveNodes.Nactivenodes;
      Free_and_Restricted_Dofs = generate_Mask_for_static_condensation__MeshTools__(ActiveNodes,FEM_Mesh);
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Second step : Compute effective mass ... WORKING",TimeStep);
      /*
	       Compute the effective mass matrix as a convex combination of the consistent mass
	       matrix and the lumped mass matrix.
      */
      Effective_Mass = compute_Nodal_Effective_Mass(MPM_Mesh,FEM_Mesh,ActiveNodes,epsilon);
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Third step : Compute nodal kinetics ... WORKING",TimeStep);
      /*
        Compute the nodal values of velocity and acceleration
      */
      U_n = compute_Nodal_Field(Effective_Mass,MPM_Mesh,FEM_Mesh,ActiveNodes);
      D_U = initialise_Nodal_Increments(U_n,FEM_Mesh,ActiveNodes,Params,TimeStep);
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Four step : Compute equilibrium ... WORKING",TimeStep);
           
      /*
	      Set the convergence false by default and start the iterations to compute
      	the incement of velocity and displacement in the nodes of the mesh
      */
      Convergence = false;
      Iter = 0;
      while(Convergence == false)
      {
        /*
          Compute the stress-strain state for each particle
        */
        update_Local_State(D_U,ActiveNodes,MPM_Mesh,FEM_Mesh,DeltaTimeStep);

        /*
          Compute the nodal forces
        */
        Forces = compute_Nodal_Forces(ActiveNodes,MPM_Mesh,FEM_Mesh,TimeStep);

        /*
          Compute nodal reactions and set to zero those DOF of the
          nodal forces with imposed displacements.
        */
        Reactions = compute_Nodal_Reactions(FEM_Mesh,Forces,ActiveNodes);

        /*
          Compute the numerical residual to check the equilibrium
        */
        Residual = compute_Nodal_Residual(U_n,D_U,ActiveNodes,Forces,Effective_Mass,Params);

        /*
          If the norm of the residual for each nodal value is below a tolerace skip
        */
        Convergence = check_convergence(Residual,TOL,Iter,MaxIter,TimeStep);

        /*
          If not, solve the linearized equilibrium to compute the next step
        */
        if(Convergence == false)
        {

          /*
		        Assemble the tangent stiffness matrix as a sum of the material tangent
		        stiffness and the geometrical tangent stiffness.
	        */
	        Tangent_Stiffness = assemble_Nodal_Tangent_Stiffness(ActiveNodes,MPM_Mesh,FEM_Mesh,Params);

          /*
        		Solve the resulting equation 
	        */
	        if((Free_and_Restricted_Dofs.Nactivenodes - Ndim*Nactivenodes) == 0)
		      {
		        solve_non_reducted_system(D_U,Tangent_Stiffness,Effective_Mass,Residual,Params);
		      }
	        else
		      {
		        solve_reducted_system(D_U,Tangent_Stiffness,Effective_Mass,Residual,Free_and_Restricted_Dofs,Params);
		      }

          update_Newmark_Nodal_Increments(D_U,U_n,Params);

	        Iter++;

	        free__MatrixLib__(Forces);
	        free__MatrixLib__(Reactions);
	        free__MatrixLib__(Residual);
	        free__MatrixLib__(Tangent_Stiffness);
        }
      }

      print_Status("DONE !!!",TimeStep);


      print_Status("*************************************************",TimeStep);
      print_Status("Seven step : Update particles lagrangian ... WORKING",TimeStep);

      update_Particles(D_U,MPM_Mesh,FEM_Mesh,ActiveNodes);
      local_search__Particles__(MPM_Mesh,FEM_Mesh);

      print_Status("DONE !!!",TimeStep);

      /*
	Outputs
      */
      output_selector(MPM_Mesh, FEM_Mesh, ActiveNodes, D_U, U_n, Forces, Reactions, Residual,TimeStep, ResultsTimeStep);

    

      print_Status("*************************************************",TimeStep);
      print_Status("Eight step : Reset nodal values ... WORKING",TimeStep);
      /*
	Free memory.
      */
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
  |    0  | M_eff |	         |   0    | M_cons |	      |   0    | M_lump |
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
  Matrix Lumped_MassMatrix = allocZ__MatrixLib__(Order, 1);

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

	  /* Get the value of the shape function */
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
	      Lumped_MassMatrix.nV[A_mask*Ndof + i] += m_A_p;	      
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
		  /*
		    Compute the vectorized index
		  */
		  Effective_MassMatrix.nM[A_mask*Ndof+i][B_mask*Ndof+i] += m_AB_p;
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
  for(int A = 0 ; A<Order ; A++)
    {
      for(int B = 0 ; B<Order ; B++)
    	{    
	     Effective_MassMatrix.nM[A][B] = (1-epsilon)*Effective_MassMatrix.nM[A][B] + (A == B)*epsilon*Lumped_MassMatrix.nV[A];
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

static Nodal_Field compute_Nodal_Field(
  Matrix Effective_Mass,
  Particle MPM_Mesh, 
  Mesh FEM_Mesh, 
  Mask ActiveNodes)
/*
  Call the LAPACK solver to compute the nodal velocity. The operation is linearized and
  all the dof split the velocity array in n components like :
  | M 0 |   |V.x|   | p.x |
  | 0 M | * |V.y| = | p.y |

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
  Nodal_Field U_n;
  U_n.value        = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  U_n.d_value_dt   = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  U_n.d2_value_dt2 = allocZ__MatrixLib__(Nnodes_mask,Ndof);

  /* Iterate over the particles to get the nodal values */
  for(int p = 0 ; p<Np ; p++)
    {

      /* Define element of the particle */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /* Evaluate the shape function in the coordinates of the particle */
      ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Get the mass of the GP */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      /* Get the nodal mommentum */
      for(int A = 0 ; A<Nodes_p.NumberNodes ; A++)
	{

	  /*
	    Get the node in the nodal momentum with the mask
	  */
	  Ap = Nodes_p.Connectivity[A];
	  A_mask = ActiveNodes.Nodes2Mask[Ap];

	  /* Evaluate the GP function in the node */
	  ShapeFunction_pA = ShapeFunction_p.nV[A];

	  /* Nodal velocity and acceleration  */
	  for(int i = 0 ; i<Ndim ; i++)
    {
        U_n.value.nM[A_mask][i]        += m_p*ShapeFunction_pA*MPM_Mesh.Phi.dis.nM[p][i];
        U_n.d_value_dt.nM[A_mask][i]   += m_p*ShapeFunction_pA*MPM_Mesh.Phi.vel.nM[p][i];
        U_n.d2_value_dt2.nM[A_mask][i] += m_p*ShapeFunction_pA*MPM_Mesh.Phi.acc.nM[p][i];
      }
	}

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }


  /*
    Call the LAPACK solver to compute the accelerations and velocities
  */
  int Order = Nnodes_mask*Ndim;
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
    Solve for the velocity
  */
  dgetrs_(&TRANS,&Order,&NRHS,Effective_Mass.nV,&LDA,IPIV,U_n.value.nV,&LDB,&INFO);
  dgetrs_(&TRANS,&Order,&NRHS,Effective_Mass.nV,&LDA,IPIV,U_n.d_value_dt.nV,&LDB,&INFO);
  dgetrs_(&TRANS,&Order,&NRHS,Effective_Mass.nV,&LDA,IPIV,U_n.d2_value_dt2.nV,&LDB,&INFO);
  free(IPIV);

  /*
    Add some usefulll info
  */
  strcpy(U_n.value.Info,"Displacement");
  strcpy(U_n.d_value_dt.Info,"Velocity");
  strcpy(U_n.d2_value_dt2.Info,"Acceleration");

  return U_n;
}

/**************************************************************/

static Nodal_Field initialise_Nodal_Increments(
  Nodal_Field U_n,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Newmark_parameters Params,
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

  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;
  double alpha_4 = Params.alpha_4;
  double alpha_5 = Params.alpha_5;
  double alpha_6 = Params.alpha_6;

  double D_U_value_It;
  Nodal_Field D_U;
  D_U.value        = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  D_U.d_value_dt   = allocZ__MatrixLib__(Nnodes_mask,Ndof);
  D_U.d2_value_dt2 = allocZ__MatrixLib__(Nnodes_mask,Ndof);

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
            Assign the boundary condition :
            First remove the value obtained during the projection and compute the value
            using the evolution of the boundary condition
          */
          U_n.value.nM[Id_BCC_mask][k] = 0.0;
          U_n.d_value_dt.nM[Id_BCC_mask][k] = 0.0;
          U_n.d2_value_dt2.nM[Id_BCC_mask][k] = 0.0;

          for(int t = 0 ; t<TimeStep ; t++)
          {
            D_U_value_It = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[t]*(double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];
            U_n.value.nM[Id_BCC_mask][k] += D_U_value_It;                    
            U_n.d2_value_dt2.nM[Id_BCC_mask][k] += alpha_1*D_U_value_It - alpha_2*U_n.d_value_dt.nM[Id_BCC_mask][k] - (alpha_3 + 1)*U_n.d2_value_dt2.nM[Id_BCC_mask][k];
            U_n.d_value_dt.nM[Id_BCC_mask][k]   += alpha_4*D_U_value_It + (alpha_5-1)*U_n.d_value_dt.nM[Id_BCC_mask][k] + alpha_6*U_n.d2_value_dt2.nM[Id_BCC_mask][k];
          }

          /*
            Initialise increments using newmark and the value of the boundary condition
          */
          D_U.value.nM[Id_BCC_mask][k] = FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*(double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];                    
          D_U.d2_value_dt2.nM[Id_BCC_mask][k] = alpha_1*D_U.value.nM[Id_BCC_mask][k] - alpha_2*U_n.d_value_dt.nM[Id_BCC_mask][k] - (alpha_3 + 1)*U_n.d2_value_dt2.nM[Id_BCC_mask][k];
          D_U.d_value_dt.nM[Id_BCC_mask][k]   = alpha_4*D_U.value.nM[Id_BCC_mask][k] + (alpha_5-1)*U_n.d_value_dt.nM[Id_BCC_mask][k] + alpha_6*U_n.d2_value_dt2.nM[Id_BCC_mask][k];

        }
      }
    }    
  }

  /*
    Add some usefulll info
  */
  strcpy(D_U.value.Info,"Displacement");
  strcpy(D_U.d_value_dt.Info,"Velocity");
  strcpy(D_U.d2_value_dt2.Info,"Acceleration");

  return D_U;
}


/**************************************************************/


static void update_Local_State(
  Nodal_Field D_U,
	Mask ActiveNodes,
	Particle MPM_Mesh,
	Mesh FEM_Mesh,
	double TimeStep)
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
  Matrix gradient_p;
  Matrix D_Displacement_Ap;
  Matrix D_Velocity_Ap;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor DF_p;
  Tensor dFdt_n_p;
  Tensor dFdt_n1_p;
  Tensor dt_DF_p;
  Tensor P_p;
  
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
      D_Displacement_Ap = get_set_field__MeshTools__(D_U.value, Nodes_p, ActiveNodes);
      D_Velocity_Ap = get_set_field__MeshTools__(D_U.d_value_dt, Nodes_p, ActiveNodes);

      /*
	Evaluate the shape function gradient in the coordinates of the particle
      */
      gradient_p = compute_dN__MeshTools__(Nodes_p,MPM_Mesh,FEM_Mesh);
	  
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
      	Update the First Piola-Kirchhoff stress tensor (P)
      */
      MatIndx_p = MPM_Mesh.MatIdx[p];
      MatProp_p = MPM_Mesh.Mat[MatIndx_p];
      P_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2);
      
      if(strcmp(MatProp_p.Type,"Neo-Hookean-Wriggers") == 0)
      {
        J_n1_p = I3__TensorLib__(F_n1_p);
        P_p = compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(P_p,F_n1_p,J_n1_p,MatProp_p);
      }
      else if(strcmp(MatProp_p.Type,"Newtonian-Fluid-Compressible") == 0)
      {
        J_n1_p = I3__TensorLib__(F_n1_p);
        P_p = compute_1PK_Stress_Tensor_Newtonian_Fluid(P_p,F_n1_p,dFdt_n1_p,J_n1_p,MatProp_p);
      }
      else
      {
        fprintf(stderr,"%s : %s %s %s \n",
          "Error in update_Local_State()",
          "The material",MatProp_p.Type,"has not been yet implemnented");
        exit(EXIT_FAILURE);
      }
      
      /*
	       Free memory 
      */
      free__MatrixLib__(D_Displacement_Ap);
      free__MatrixLib__(D_Velocity_Ap);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);
	  
    }
  
}

/**************************************************************/

static Matrix compute_Nodal_Forces(
	Mask ActiveNodes,
	Particle MPM_Mesh,
	Mesh FEM_Mesh,
	int TimeStep)
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  Matrix Forces = allocZ__MatrixLib__(Nnodes_mask,Ndim);

  /*
    Add internal forces contribution
  */
  compute_Nodal_Internal_Forces(Forces, ActiveNodes, MPM_Mesh, FEM_Mesh);

  /*
    Add nominal traction contribution
  */
  compute_Nodal_Nominal_traction_Forces(Forces, ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep);

  /*
    Add body forces contribution
  */
  compute_Nodal_Body_Forces(Forces, ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep);

  return Forces;
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
  double V0_p; /* Volume of the particle in the reference configuration */

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
      Take the value of the deformation gradient at t = n 
    */
    F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    transpose_F_n_p = transpose__TensorLib__(F_n_p);

    /*
      Get the first Piola-Kirchhoff stress tensor
    */
    P_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);
    
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
	      Forces.nM[A_mask][i] += InternalForcesDensity_Ap.n[i]*V0_p;
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
          Forces.nM[A_mask][k] -= ShapeFunction_pA*T.n[k]*A0_p;
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

static void compute_Nodal_Body_Forces(
  Matrix Forces,
  Mask ActiveNodes,
	Particle MPM_Mesh,
	Mesh FEM_Mesh,
  int TimeStep)
{
  /* Define auxilar variables */
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int NumBodyForces = MPM_Mesh.NumberBodyForces;
  int NumParticles_i; /* Number of particles with the */
  int NumNodes_p; /* Number of tributary nodes of p */
  int A_mask; /* Index of the node where we apply the body force */
  int idx_A_mask_k; /* Index of the node where we apply the body force */
  int Ap; /* Tributary node A of particle p */
  int p; /* Particle index */
  
  double m_p; /* Mass of the particle */
  Load * B = MPM_Mesh.B; /* List with the load cases */
  Element Nodes_p; /* Element for each particle */
  Matrix ShapeFunction_p; /* Nodal values of the sahpe function */
  double ShapeFunction_pA; /* Evaluation in the node I for the particle p */
  
  Tensor b = alloc__TensorLib__(1); /* Body forces vector */
  
  for(int i = 0 ; i<NumBodyForces ; i++)
    {
      /* Get the number of particles with the body load i */
      NumParticles_i = B[i].NumNodes;

      for(int j = 0 ; j<NumParticles_i ; j++){

	/* Get the index of the Gauss-Point */
	p = B[i].Nodes[j];
	
	/* Get the value of the mass */
	m_p = MPM_Mesh.Phi.mass.nV[p];

	/* Define tributary nodes of the particle */
	NumNodes_p = MPM_Mesh.NumberNodes[p];
	Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);

	/* Compute shape functions */
	ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);


	/* Fill vector of body forces */
	for(int k = 0 ; k<Ndim ; k++){
	  if(B[i].Dir[k]){
	    if( (TimeStep < 0) || (TimeStep > B[i].Value[k].Num)){
	      printf("%s : %s\n",
		     "Error in compute_Nodal_Body_Forces()",
		     "The time step is out of the curve !!");
	      exit(EXIT_FAILURE);
	    }
	    b.n[k] = B[i].Value[k].Fx[TimeStep];
	  }
	}


	/* Get the node of the mesh for the contribution */
	for(int A = 0 ; A<NumNodes_p ; A++)
	  {
	    
	    /* Pass the value of the nodal shape function to a scalar */
	    ShapeFunction_pA = ShapeFunction_p.nV[A];
	    
	    /* Get the node of the mesh for the contribution */
	    Ap = Nodes_p.Connectivity[A];
	    A_mask = ActiveNodes.Nodes2Mask[Ap];
	    
	    /* Compute body forces */
	    for(int k = 0 ; k<Ndim ; k++)
	      {
          Forces.nM[A_mask][k] -= ShapeFunction_pA*b.n[k]*m_p;
	      } 
	    
	  }
	
	/* Free the matrix with the nodal gradient of the element */
	free__MatrixLib__(ShapeFunction_p);
	free(Nodes_p.Connectivity);
	
      }
      
    }

  /*
    Free auxiliar tensor
  */
  free__TensorLib__(b);

  
}


/**********************************************************************/

static Matrix compute_Nodal_Reactions(
  Mesh FEM_Mesh,
  Matrix Forces,
  Mask ActiveNodes)
/*
  Compute the nodal reactions
*/
{
  /* 1º Define auxilar variables */
  int Ndim = NumberDimensions;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

  Matrix Reactions = allocZ__MatrixLib__(Nnodes_mask,Ndim);
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
  Newmark_parameters Params)
{
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Ndim*Nnodes_mask;
  Matrix Acceleration_n1 = allocZ__MatrixLib__(Nnodes_mask,Ndim);
  Matrix Inertial_Forces = allocZ__MatrixLib__(Nnodes_mask,Ndim);
  Matrix Residual = allocZ__MatrixLib__(Nnodes_mask,Ndim);
  double alpha_1 = Params.alpha_1;
  double alpha_2 = Params.alpha_2;
  double alpha_3 = Params.alpha_3;

  /*
    Compute nodal acceleration (Vectorized)
  */
  for(int idx_B = 0 ; idx_B<Order ; idx_B++)
    {
      Acceleration_n1.nV[idx_B] = alpha_1*D_U.value.nV[idx_B] - alpha_2*U_n.d_value_dt.nV[idx_B] - alpha_3*U_n.d2_value_dt2.nV[idx_B];
    }
  /*
    Compute inertial forces (Vectorized)
  */
  for(int idx_A = 0 ; idx_A<Order ; idx_A++)
    {
      for(int idx_B = 0 ; idx_B<Order ; idx_B++)
      {
        Inertial_Forces.nV[idx_A] += Mass.nM[idx_A][idx_B]*Acceleration_n1.nV[idx_B];
      }
    }

  /*
    Compute (-) residual (Vectorized). The minus symbol is due to
    solver purposes. See compute_D_Displacement
  */
  for(int idx_A = 0 ; idx_A<Order ; idx_A++)
    {
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
static bool check_convergence(
  Matrix Residual,
  double TOL,
  int Iter,
  int MaxIter,
  int Step)
{
  bool convergence;
  int Ndim = NumberDimensions;
  int Nnodes_mask = Residual.N_rows;
  int Total_dof = Ndim*Nnodes_mask;
  double Error = 0;
  double Error_relative = 0;

  if(Iter > MaxIter)
    {
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in check_convergence()",
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
  Tensor Stiffness_density_p;

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
        if(strcmp(MatProp_p.Type,"Neo-Hookean-Wriggers") == 0)
        {      
          Stiffness_density_p = compute_stiffness_density_Neo_Hookean_Wriggers(GRADIENT_pA,GRADIENT_pB,F_n1_p,J_p,MatProp_p);
        }
        else if(strcmp(MatProp_p.Type,"Newtonian-Fluid-Compressible") == 0)
        {
          Stiffness_density_p = compute_stiffness_density_Newtonian_Fluid(GRADIENT_pA,GRADIENT_pB,F_n1_p,dFdt_n1_p,J_p,alpha_4,MatProp_p);
        }
        else
        {
          fprintf(stderr,"%s : %s %s %s \n",
            "Error in assemble_Nodal_Tangent_Stiffness()",
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
		        Tangent_Stiffness.nM[A_mask*Ndim+i][B_mask*Ndim+j] += Stiffness_density_p.N[i][j]*V0_p;
		      }
		    }

	      /*
		Free memory 
	      */
	      free__TensorLib__(GRADIENT_pB);
	      free__TensorLib__(Stiffness_density_p);
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

static void solve_non_reducted_system(
  Nodal_Field D_U,
	Matrix Tangent_Stiffness,
	Matrix Effective_Mass,
	Matrix Residual,
  Newmark_parameters Params)
/*
  This function is deboted to update the vector with the nodal displacement by
  solving :
  
  ((2/Dt^2)*M_AB + K_AB)*delta_B + R_A = 0_A 
  
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
  double alpha_1 = Params.alpha_1;
  Matrix K_Global = allocZ__MatrixLib__(Order,Order);

  /*
    Compute the adition of the mass matrix and the tangent stifness matrix
  */
  for(int idx_AB_ij = 0 ; idx_AB_ij<Order*Order ; idx_AB_ij++)
    {
      K_Global.nV[idx_AB_ij] = alpha_1*Effective_Mass.nV[idx_AB_ij] + Tangent_Stiffness.nV[idx_AB_ij];
    }


  /*
    Compute the LU factorization 
  */
  dgetrf_(&Order,&Order,K_Global.nV,&LDA,IPIV,&INFO);

  /*
    Check error messages in the LAPACK LU descompistion  
  */
  if(INFO)
    {
      fprintf(stderr,"%s : %s %s %s \n",
	      "Error in solve_non_reducted_system",
	      "The function",
	      "dgetrf_",
	      "returned an error message !!!" );
      exit(EXIT_FAILURE);
    }

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order,&NRHS,K_Global.nV,&LDA,IPIV,Residual.nV,&LDB,&INFO);
  free(IPIV);

  
  /*
    Check error messages in the LAPACK solver  
  */
  if(INFO)
    {
      fprintf(stderr,"%s : %s %s %s \n",
	      "Error in solve_non_reducted_system",
	      "The function",
	      "dgetrs_",
	      "returned an error message !!!" );
      exit(EXIT_FAILURE);
    }

  /*
    Update 
  */
  for(int idx_A_i = 0 ; idx_A_i < Order ; idx_A_i++)
    {	
      D_U.value.nV[idx_A_i] -= Residual.nV[idx_A_i];
    }

  /*
    Free auxiliar global matrix
  */
  free__MatrixLib__(K_Global);
  
}


/**************************************************************/

static void solve_reducted_system(
  Nodal_Field D_U,
  Matrix Tangent_Stiffness,
  Matrix Effective_Mass,
  Matrix Residual,
  Mask Free_and_Restricted_Dofs,
  Newmark_parameters Params)
/*
  This function is deboted to update the vector with the nodal displacement by
  solving :
  
  ((2/Dt^2)*M_AB + K_AB)*delta_B + R_A = 0_A 
  
*/
{
  int Nnodes_mask = Residual.N_rows;
  int Ndof = NumberDOF;
  int Order = Nnodes_mask*Ndof;
  int Num_Free_dofs = Free_and_Restricted_Dofs.Nactivenodes;
  int Free_A_ij;
  int idx_A_ij, idx_B_ij;
  int Mask_idx_A_ij, Mask_idx_B_ij;

  double alpha_1 = Params.alpha_1;
  /*
    Guyan reduction : static condensation
    Here, we will work directly with vectorized matrix
  */
  Matrix K_Global_FF = allocZ__MatrixLib__(Num_Free_dofs,Num_Free_dofs);
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
        K_Global_FF.nM[Mask_idx_A_ij][Mask_idx_B_ij] = alpha_1*Effective_Mass.nM[idx_A_ij][idx_B_ij] + Tangent_Stiffness.nM[idx_A_ij][idx_B_ij];
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
  dgetrf_(&Order_FF,&Order_FF,K_Global_FF.nV,&LDA,IPIV,&INFO);

  /*
    Check error messages in the LAPACK LU descompistion  
  */
  if(INFO)
    {
      fprintf(stderr,"%s : %s %s %s \n",
        "Error in solve_reducted_system",
        "The function","dgetrf_",
        "returned an error message !!!" );
      exit(EXIT_FAILURE);
    }

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order_FF,&NRHS,K_Global_FF.nV,&LDA,IPIV,Residual_F.nV,&LDB,&INFO);
  free(IPIV);
  
  /*
    Check error messages in the LAPACK solver  
  */
  if(INFO)
    {
      fprintf(stderr,"%s : %s %s %s \n",
        "Error in solve_reducted_system",
        "The function","dgetrs_",
        "returned an error message !!!" );
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

        D_U.value.nV[idx_A_ij] -= Residual_F.nV[Mask_idx_A_ij];
      } 
      
    }

  /*
    Free auxiliar
  */
  free__MatrixLib__(K_Global_FF);
  free__MatrixLib__(Residual_F);

}

/**************************************************************/

static void update_Newmark_Nodal_Increments(
  Nodal_Field D_U,
  Nodal_Field U_n,
  Newmark_parameters Params)
{

  int Nnodes = U_n.value.N_rows;
  int Ndim = NumberDimensions;
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
    for(int i = 0 ; i<Ndim ; i++)
    {
      D_U.d2_value_dt2.nM[A][i] = alpha_1*D_U.value.nM[A][i] - alpha_2*U_n.d_value_dt.nM[A][i] - (alpha_3 + 1)*U_n.d2_value_dt2.nM[A][i];
      D_U.d_value_dt.nM[A][i]   = alpha_4*D_U.value.nM[A][i] + (alpha_5-1)*U_n.d_value_dt.nM[A][i] + alpha_6*U_n.d2_value_dt2.nM[A][i];
    }
  }
}

/**************************************************************/

static void update_Particles(
  Nodal_Field D_U,
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
  int MatIndx_p;
  Material MatProp_p;
  Matrix D_Displacement_Ap;
  Matrix ShapeFunction_p; /* Value of the shape-function in the particle */
  Matrix gradient_p;
  double ShapeFunction_pI; /* Nodal value for the particle */
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor DF_p;
  Tensor dFdt_n_p;
  Tensor dFdt_n1_p;
  double Delta_J_p;
  double rho_n_p;
  Element Nodes_p; /* Element for each particle */
  double D_U_pI;
  double D_V_pI;
  double D_A_pI;
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
      	Take the values of the deformation gradient from the previous step
      */
      F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
      DF_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.DF.nM[p],2);
      dFdt_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n.nM[p],2);
      dFdt_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n1.nM[p],2);
      /*
	     Update density with the jacobian of the increment deformation gradient
      */
      Delta_J_p = I3__TensorLib__(DF_p);
      rho_n_p = MPM_Mesh.Phi.rho.nV[p];
      MPM_Mesh.Phi.rho.nV[p] = rho_n_p/Delta_J_p;

      /*
	     Replace the deformation gradient at t = n with the converged deformation gradient
      */
      for(int i = 0 ; i<Ndim  ; i++)
	     {
	     for(int j = 0 ; j<Ndim  ; j++)
	       {
	         F_n_p.N[i][j] = F_n1_p.N[i][j];
           dFdt_n_p.N[i][j] = dFdt_n1_p.N[i][j];
	       }
	     }

      /* Compute the deformation energy */
      Vol_0_p = MPM_Mesh.Phi.Vol_0.nV[p];
      MatIndx_p = MPM_Mesh.MatIdx[p];
      MatProp_p = MPM_Mesh.Mat[MatIndx_p];
//      MPM_Mesh.Phi.W.nV[p]= finite_strains_internal_energy__Particles__(F_n_p, MatProp_p,Vol_0_p);
      
      /* Iterate over the nodes of the particle */
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
	    Update velocity position and deformation gradient of the particles
	  */
	  for(int i = 0 ; i<Ndim ; i++)
	    {
	      D_U_pI = ShapeFunction_pI*D_U.value.nM[A_mask][i];
	      D_V_pI = ShapeFunction_pI*D_U.d_value_dt.nM[A_mask][i];
	      D_A_pI = ShapeFunction_pI*D_U.d2_value_dt2.nM[A_mask][i];
	      	      
	      MPM_Mesh.Phi.acc.nM[p][i]  += D_A_pI;
	      MPM_Mesh.Phi.vel.nM[p][i]  += D_V_pI;
        MPM_Mesh.Phi.dis.nM[p][i]  += D_U_pI;
	      MPM_Mesh.Phi.x_GC.nM[p][i] += D_U_pI;
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
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Mask ActiveNodes,
  Nodal_Field D_U,
  Nodal_Field U_n,
  Matrix Forces,
  Matrix Reactions,
  Matrix Residual,
  int TimeStep,
  int ResultsTimeStep)
{

  /*
    vtk results
  */
  if(TimeStep % ResultsTimeStep == 0)
  {
    particle_results_vtk__InOutFun__(MPM_Mesh,TimeStep,ResultsTimeStep);

    nodal_results_vtk__InOutFun__(FEM_Mesh, ActiveNodes, Reactions, TimeStep, ResultsTimeStep);
  }

  // /* 
  //   csv results 
  // */
  // for(int i = 0 ; i<Number_Out_nodal_path_csv ; i++)
  // {

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_Velocity)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(Velocity, FEM_Mesh.Coordinates,"Nodal_path_velocity_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_Acceleration)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(Acceleration, FEM_Mesh.Coordinates,"Nodal_path_acceleration_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_D_Displacement)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(D_Displacement, FEM_Mesh.Coordinates,"Nodal_path_displacement_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_Forces)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(Forces, FEM_Mesh.Coordinates,"Nodal_path_forces_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_Reactions)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(Reactions, FEM_Mesh.Coordinates,"Nodal_path_reactions_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }

  //   if(Out_nodal_path_csv[i].Out_csv_nodes_path_Residual)
  //   {
  //     path_nodes_analysis_csv__InOutFun__(Residual, FEM_Mesh.Coordinates,"Nodal_path_residual_csv", ActiveNodes, Out_nodal_path_csv[i], i, TimeStep, DeltaTimeStep);
  //   }
  // }


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
