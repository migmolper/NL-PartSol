#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>

#endif

/*
  Auxiliar functions 
*/
static Matrix compute_Nodal_Effective_Mass(GaussPoint, Mesh, Mask, double);
static Matrix compute_Nodal_Momentum(GaussPoint, Mesh, Mask);
static Matrix compute_Nodal_Velocity(Matrix, Matrix);
static void   update_Local_State(Matrix,Mask,GaussPoint,Mesh,double);
static Matrix compute_Nodal_Forces(Matrix, Mask, GaussPoint, Mesh);
static void   compute_Nodal_Internal_Forces(Matrix,Matrix,Mask,GaussPoint, Mesh);
/* static void   compute_Nodal_Body_Forces(Matrix,Mask,GaussPoint, Mesh, int); */
static Matrix compute_Nodal_Residual(Matrix, Matrix, Matrix, Matrix, double);
static bool   check_convergence(Matrix,double,int,int);
static Matrix assemble_Nodal_Tangent_Stiffness(Mask, GaussPoint, Mesh);
static void   assemble_Nodal_Tangent_Stiffness_Geometric(Matrix, Mask, GaussPoint, Mesh);
static void   assemble_Nodal_Tangent_Stiffness_Material(Matrix, Mask, GaussPoint, Mesh);
static Tensor compute_stiffness_density(Tensor, Tensor, Material);
static Tensor compute_Nodal_Tangent_Stiffness_Material(Tensor,Tensor,Tensor);
static void   update_D_Displacement(Matrix, Matrix, Matrix, Matrix, double);
static void   imposed_displacements(Matrix, Mask, Mesh, int);
static Matrix compute_D_Velocity(Matrix, Matrix,double);
static void   update_Particles(Matrix, Matrix, GaussPoint, Mesh, Mask, double);

/**************************************************************/

void U_Discrete_Energy_Momentum(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)
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
  Matrix Momentum;
  Matrix Velocity;
  Matrix D_Displacement;
  Matrix D_Velocity;
  Matrix Residual;
  Mask ActiveNodes;
  double TOL = 0.000000000001;
  double epsilon = 1.0;
  double DeltaTimeStep;
  bool Convergence;
  int Iter = 0;
  int MaxIter = 500;

  /*
    Time step is defined at the init of the simulation throught the
    CFL condition. Notice that for this kind of solver, CFL confition is
    not required to be satisfied. The only purpose of it is to use the existing
    software interfase.
   */
  DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);

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
      print_Status("Third step : Compute nodal momentum ... WORKING",TimeStep);
      /*
	Compute the nodal value of the momentum
      */
      Momentum = compute_Nodal_Momentum(MPM_Mesh, FEM_Mesh, ActiveNodes);
      print_Status("DONE !!!",TimeStep);
      
      print_Status("*************************************************",TimeStep);
      print_Status("Four step : Compute nodal velocity ... WORKING",TimeStep);
      /*
	Compute the nodal valocities with the effective mass matrix and the nodal momentum
      */
      Velocity = compute_Nodal_Velocity(Effective_Mass, Momentum);
      /* imposed_velocty(); */
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Five step : Compute equilibrium ... WORKING",TimeStep);
      /*
	Set to zero the increment of velocity and displacement for nodal values
      */
      D_Displacement = allocZ__MatrixLib__(Ndim,Nactivenodes);  
      D_Velocity = allocZ__MatrixLib__(Ndim,Nactivenodes);
      
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
	  update_Local_State(D_Displacement,ActiveNodes,MPM_Mesh,FEM_Mesh,DeltaTimeStep);

	  /*
	    Compute the nodal forces
	  */
	  Forces = compute_Nodal_Forces(D_Displacement,ActiveNodes,MPM_Mesh,FEM_Mesh);

	 	  
	  /*
	    Compute the numerical residual to check the equilibrium
	  */
	  Residual = compute_Nodal_Residual(Velocity,Forces,D_Displacement,
					    Effective_Mass,DeltaTimeStep);
	  
	  /*
	    If the norm of the residual for each nodal value is below a tolerace
	    skip
	  */
	  Convergence = check_convergence(Residual,TOL,Iter,MaxIter);

	  /*
	    If not, solve the linearized equilibrium to compute the next step
	  */
	  if(Convergence == false)
	    {

	      /*
		Assemble the tangent stiffness matrix as a sum of the material tangent
		stiffness and the geometrical tangent stiffness.
	      */
	      Tangent_Stiffness = assemble_Nodal_Tangent_Stiffness(ActiveNodes,
								   MPM_Mesh,FEM_Mesh);
	      
	      /*
		Solve the resulting equation 
	      */
	      update_D_Displacement(D_Displacement,Tangent_Stiffness,
				    Effective_Mass,Residual,DeltaTimeStep);

	      /*
		Impose dirichlet boundary conditions over the increment of
		displacement
	      */
	      imposed_displacements(D_Displacement, ActiveNodes, FEM_Mesh, TimeStep);

	      /* PrintMatrix(Residual,Ndim,Nactivenodes); */

	      /*
		Update the iteration number
	      */
	      Iter++;

	      free__MatrixLib__(Forces);
	      free__MatrixLib__(Residual);
	      free__MatrixLib__(Tangent_Stiffness);
	    }
            
	}
      
      print_iteration(TimeStep,Iter);
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Six step : Compute increment of velocity ... WORKING",TimeStep);
      /*
	Once the equilibrium is reached, obtain the increment of nodal velocity
      */
      D_Velocity = compute_D_Velocity(Velocity,D_Displacement,DeltaTimeStep);
      print_Status("DONE !!!",TimeStep);

      
      print_Status("*************************************************",TimeStep);
      print_Status("Seven step : Update particles lagrangian ... WORKING",TimeStep);
      /*
	Update Lagrangians with D_Displacement
      */
      update_Particles(D_Displacement,D_Velocity,MPM_Mesh,FEM_Mesh,ActiveNodes,DeltaTimeStep);
      print_Status("DONE !!!",TimeStep);
      /*
	Reload the connectivity information for each particle
      */
      local_search__Particles__(MPM_Mesh,FEM_Mesh);
      print_Status("DONE !!!",TimeStep);
      
      /*
	Outputs
      */
      if(TimeStep % ResultsTimeStep == 0)
	{
	  /* WriteVtk_FEM("Mesh",FEM_Mesh,R_I,TimeStep); */
	  WriteVtk_MPM("MPM_VALUES",MPM_Mesh,"ALL",TimeStep,ResultsTimeStep);
	}

      print_Status("*************************************************",TimeStep);
      print_Status("Eight step : Reset nodal values ... WORKING",TimeStep);
      /*
	Free memory.
      */
      free__MatrixLib__(Effective_Mass); 
      free__MatrixLib__(Velocity);
      free__MatrixLib__(D_Velocity);
      free__MatrixLib__(D_Displacement);
      free__MatrixLib__(Forces);
      free__MatrixLib__(Residual);
      free(ActiveNodes.Mask2Nodes);
      free(ActiveNodes.Nodes2Mask);
      
      print_Status("DONE !!!",TimeStep);

    }
  
}

/**************************************************************/

static Matrix compute_Nodal_Effective_Mass(GaussPoint MPM_Mesh,
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
  Matrix Effective_MassMatrix = allocZ__MatrixLib__(Nnodes_mask*Ndof, Nnodes_mask*Ndof);

  /* Define and allocate the lumped mass matrix */
  Matrix Lumped_MassMatrix = allocZ__MatrixLib__(Nnodes_mask*Ndof, 1);

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
	      idx_A_mask_i = A_mask + i*Nnodes_mask;
	      Lumped_MassMatrix.nV[idx_A_mask_i] += m_A_p;	      
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
		  idx_AB_mask_i = B_mask + (i*Order + i)*Nnodes_mask + A_mask*Order;
		  
		  Effective_MassMatrix.nV[idx_AB_mask_i] += m_AB_p;
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
  for(int A = 0 ; A<Nnodes_mask ; A++)
    {
      for(int B = 0 ; B<Nnodes_mask ; B++)
	{
	  for(int i = 0 ; i<Ndof ; i++)
	    {
	      idx_AB_mask_i = B + (i*Ndof*Nnodes_mask + i)*Nnodes_mask + A*Ndof*Nnodes_mask;
	      idx_A_mask_i = A + i*Nnodes_mask;
	      
	      Effective_MassMatrix.nV[idx_AB_mask_i] =
		(1-epsilon)*Effective_MassMatrix.nV[idx_AB_mask_i] +
		(A == B)*epsilon*Lumped_MassMatrix.nV[idx_A_mask_i];
	    }
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

static Matrix compute_Nodal_Momentum(GaussPoint MPM_Mesh,
				     Mesh FEM_Mesh,
				     Mask ActiveNodes)
/*

*/
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
  /* Mass of the particle */
  double m_p;
  /* Element for each particle */
  Element Nodes_p;

 /* Define and allocate the momentum vector */
  Matrix Momentum = allocZ__MatrixLib__(Ndim,Nnodes_mask);
    
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
      
	  /* Nodal momentum */
	  for(int i = 0 ; i<Ndim ; i++)
	    {
	      idx_A_mask_i = A_mask + i*Nnodes_mask;
	      Momentum.nV[idx_A_mask_i] += m_p*ShapeFunction_pA*MPM_Mesh.Phi.vel.nM[p][i];
	    }
	}

      /* Free the value of the shape functions */
      free__MatrixLib__(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }

  /*
     Add some usefulll info
  */
  strcpy(Momentum.Info,"Nodal-Momentum");
 
  return Momentum;
}

/**************************************************************/

static Matrix compute_Nodal_Velocity(Matrix Mass,
				     Matrix Momentum)
/* 
   Call the LAPACK solver to compute the nodal velocity. The operation is linearized and
   all the dof split the velocity array in n components like :
   | M 0 |   |V.x|   | p.x |
   | 0 M | * |V.y| = | p.y |
*/
{

  Matrix Velocity;

  int Ndim = NumberDimensions;
  int Nnodes = Mass.N_rows;
  int Order = Nnodes;
  int LDA   = Nnodes;
  int LDB = Nnodes;
  char  TRANS = 'N'; /* (No transpose) */
  int   INFO= 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;

  /*
    Compute the LU factorization 
  */
  dgetrf_(&Order,&Order,Mass.nV,&LDA,IPIV,&INFO);

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order,&NRHS,Mass.nV,&LDA,IPIV,Momentum.nV,&LDB,&INFO);
  free(IPIV);

  Velocity = Momentum;

  /*
    Add some usefulll info
  */
  strcpy(Velocity.Info,"Nodal-Velocity");

  return Velocity;
}

/**************************************************************/

static void update_Local_State(Matrix D_Displacement,
			       Mask ActiveNodes,
			       GaussPoint MPM_Mesh,
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
  double J_p;  
  Element Nodes_p;
  Material MatProp_p;
  Matrix gradient_p;
  Matrix D_Displacement_Ap;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor f_n1_p;
  Tensor S_p;
  
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
      
      /*
      	Compute the increment of the deformation gradient
      */
      f_n1_p = increment_Deformation_Gradient__Particles__(D_Displacement_Ap,gradient_p);

      /*
	Update the deformation gradient in t = n + 1 with the information
	from t = n and the increment of deformation gradient.
      */  
      update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, f_n1_p);
      
      /*
      	Update the second Piola-Kirchhoff stress tensor (S) with an apropiate
	integration rule.
      */
      MatIndx_p = MPM_Mesh.MatIdx[p];
      MatProp_p = MPM_Mesh.Mat[MatIndx_p];
      S_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2);      
      S_p = average_strain_integration_Stress__Particles__(S_p,F_n1_p,F_n_p,MatProp_p);
      
      /*
	Free memory 
      */
      free__TensorLib__(f_n1_p);
      free__MatrixLib__(D_Displacement_Ap);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);
	  
    }
  
}

/**************************************************************/

static Matrix compute_Nodal_Forces(Matrix D_Displacement,
				   Mask ActiveNodes,
				   GaussPoint MPM_Mesh,
				   Mesh FEM_Mesh)
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  Matrix Forces = allocZ__MatrixLib__(Ndim,Nnodes_mask);

  /*
    Add internal forces contribution
  */
  compute_Nodal_Internal_Forces(Forces, D_Displacement, ActiveNodes, MPM_Mesh, FEM_Mesh);

  /* compute_Nodal_Body_Forces(D_Displacement, ActiveNodes, MPM_Mesh, FEM_Mesh, TimeStep); */
  
  return Forces;
}


/**************************************************************/

static void compute_Nodal_Internal_Forces(Matrix Forces,
					  Matrix D_Displacement,
					  Mask ActiveNodes,
					  GaussPoint MPM_Mesh,
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
  Tensor S_p; /* Second Piola-Kirchhoff Stress tensor */
  Tensor InternalForcesDensity_Ap;

  Element Nodes_p; /* List of nodes for particle */
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pA;
  Tensor GRADIENT_pA;
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor F_n12_p;
  Tensor transpose_F_n_p;
  
  double m_p; /* Mass of the Gauss-Point */
  double rho_p; /* Density of the Gauss-Point */
  double V0_p; /* Volume of the Gauss-Point */
  double J_p; /* Jacobian of the deformation gradient */

  /*
    Loop in the particles 
  */
  for(int p = 0 ; p<Np ; p++)
    {
      
      /*
	Get the value of the density 
      */
      rho_p = MPM_Mesh.Phi.rho.nV[p];

      /*
	Get the value of the mass 
      */
      m_p = MPM_Mesh.Phi.mass.nV[p];

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
      F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
      F_n12_p = Convex_combination__TensorLib__(F_n1_p,F_n_p,0.5);
      transpose_F_n_p = transpose__TensorLib__(F_n_p);

      /*
	Compute the first Piola-Kirchhoff stress tensor
      */
      S_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);
      P_p = matrix_product__TensorLib__(F_n12_p, S_p);

      /*
	Compute the volume of the particle in the reference configuration 
      */
      J_p = I3__TensorLib__(F_n_p);
      V0_p = (1/J_p)*(m_p/rho_p);
    
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
	      idx_A_mask_i = A_mask + i*Nnodes_mask;
	      Forces.nV[idx_A_mask_i] += InternalForcesDensity_Ap.n[i]*V0_p;
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
      free__TensorLib__(F_n12_p);
      free__TensorLib__(transpose_F_n_p);
      free__TensorLib__(P_p);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);
    }
   
}

/**************************************************************/

/* static void compute_Nodal_Body_Forces(Matrix D_Displacement, */
/* 				      Mask ActiveNodes, */
/* 				      GaussPoint MPM_Mesh, */
/* 				      Mesh FEM_Mesh, */
/* 				      int TimeStep) */
/* { */
  
/* } */


/**************************************************************/

static Matrix compute_Nodal_Residual(Matrix Velocity,
				     Matrix Forces,
				     Matrix D_Displacement,
				     Matrix Mass,
				     double Dt)
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = Velocity.N_cols;
  int Order = Ndim*Nnodes_mask;
  Matrix Acceleration = allocZ__MatrixLib__(Ndim,Nnodes_mask);
  Matrix Inertial_Forces = allocZ__MatrixLib__(Ndim,Nnodes_mask);
  Matrix Residual = allocZ__MatrixLib__(Ndim,Nnodes_mask);
  int idx_AB;

  /*
    Compute nodal acceleration (Vectorized)
  */
  for(int idx_B = 0 ; idx_B<Order ; idx_B++)
    {
      Acceleration.nV[idx_B] = (2/DSQR(Dt))*(D_Displacement.nV[idx_B] - Dt*Velocity.nV[idx_B]);
    }
  /*
    Compute inertial forces (Vectorized)
  */
  for(int idx_A = 0 ; idx_A<Order ; idx_A++)
    {
      for(int idx_B = 0 ; idx_B<Order ; idx_B++)
	{
	  idx_AB = idx_A*Order + idx_B;
	  Inertial_Forces.nV[idx_A] += Mass.nV[idx_AB]*Acceleration.nV[idx_B];
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
  free__MatrixLib__(Acceleration);
  free__MatrixLib__(Inertial_Forces);

  
  return Residual;
}

/**************************************************************/
static bool check_convergence(Matrix Residual,
			      double TOL,
			      int Iter,
			      int MaxIter)
{
  bool convergence;
  int Ndim = NumberDimensions;
  int Nnodes_mask = Residual.N_cols;
  double Error_A;

  if(Iter > MaxIter)
    {
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in U_Discrete_Energy_Momentum()",
	      "Convergence not reached in the maximum number of iterations");
      exit(EXIT_FAILURE);
    }
  else
    {
      for(int A = 0 ; A<Nnodes_mask ; A++)
	{
	  Error_A = 0;
	  for(int i = 0 ; i<Ndim ; i++)
	    {
	      Error_A += DSQR(Residual.nM[i][A]);
	    }
	  Error_A = pow(Error_A,0.5);
	  if(Error_A > TOL)
	    {
	      return false;
	    }
	}

      return true; 
    }
}

/**************************************************************/

static Matrix assemble_Nodal_Tangent_Stiffness(Mask ActiveNodes,
					       GaussPoint MPM_Mesh,
					       Mesh FEM_Mesh)

/*
  This function computes the tangent stiffness matrix as a combination
  of the geometrical stiffness matrix and the material stiffness matrix. 

  | K_xx |  K_xy |   | K_geom_xx |    0      |     | K_mat_xx | K_mat_xy |
  ---------------- = -------------------------  +  ----------------------- 
  | K_yx |  K_yy |   |	   0	 | K_geom_yy |	   | K_mat_yx | K_mat_y  |

*/
{

  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ndof = NumberDOF;
  int Order = Ndof*Nnodes_mask;

  Matrix Tangent_Stiffness = allocZ__MatrixLib__(Order, Order);
  
  /*
    Compute terms related to the geometric non-linearities.
  */
  assemble_Nodal_Tangent_Stiffness_Geometric(Tangent_Stiffness,ActiveNodes,MPM_Mesh,FEM_Mesh);
  
  /*
    Compute term related to the material non-linearities
   */
  assemble_Nodal_Tangent_Stiffness_Material(Tangent_Stiffness,ActiveNodes,MPM_Mesh,FEM_Mesh);

  return Tangent_Stiffness;
}


/**************************************************************/

static void assemble_Nodal_Tangent_Stiffness_Geometric(Matrix Tangent_Stiffness,
						       Mask ActiveNodes,
						       GaussPoint MPM_Mesh,
						       Mesh FEM_Mesh)
/*
  Introduce the geometric contribution G_AB to the full stiffness matrix K_AB
*/
{
  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Ndim*Nnodes_mask;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int Bp;
  int B_mask;
  int NumNodes_p;
  int idx_AB_mask_i;

  Tensor S_p; /* Second Piola-Kirchhoff Stress tensor */

  Element Nodes_p; /* List of nodes for particle */
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pA;
  Tensor GRADIENT_pA;
  Tensor gradient_pB;
  Tensor GRADIENT_pB;
  Tensor F_n_p;
  Tensor transpose_F_n_p;
  Tensor GRADIENT_pA_o_GRADIENT_pB;
  
  double m_p; /* Mass of the Gauss-Point */
  double rho_p; /* Density of the Gauss-Point */
  double V0_p; /* Volume of the Gauss-Point */
  double J_p; /* Jacobian of the deformation gradient */
  double Geometric_AB_p; /* Geometric contribution */

  /*
    Loop in the particles 
  */
  for(int p = 0 ; p<Np ; p++)
    {

      /*
      	Get the value of the density 
      */
      rho_p = MPM_Mesh.Phi.rho.nV[p];

      /*
	Get the value of the mass 
      */
      m_p = MPM_Mesh.Phi.mass.nV[p];

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
	Take the values of the deformation gradient ant t = n and  
	the transpose of the deformation gradient.
      */
      F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      transpose_F_n_p = transpose__TensorLib__(F_n_p);

      /*
	Compute the first Piola-Kirchhoff stress tensor
      */
      S_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);

      /*
	Compute the volume of the particle in the reference configuration 
      */
      J_p = I3__TensorLib__(F_n_p);
      V0_p = (1/J_p)*(m_p/rho_p);

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
		Compute the dyadic product of both gradients 
	       */
	      GRADIENT_pA_o_GRADIENT_pB = dyadic_Product__TensorLib__(GRADIENT_pA,GRADIENT_pB);

	      /*
		Get the nodal contribution
	       */
	      Geometric_AB_p = inner_product__TensorLib__(S_p, GRADIENT_pA_o_GRADIENT_pB);

	      /*
		Add the geometric contribution
	      */
	      for(int i = 0 ; i<Ndim ; i++)
		{
		  /*
		    Compute the vectorized index
		  */
		  idx_AB_mask_i = B_mask + (i*Order + i)*Nnodes_mask + A_mask*Order;

		  Tangent_Stiffness.nV[idx_AB_mask_i] += Geometric_AB_p*V0_p;
		}

	      /*
		Free memory 
	      */
	      free__TensorLib__(GRADIENT_pB);
	      free__TensorLib__(GRADIENT_pA_o_GRADIENT_pB);
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

}


/**************************************************************/
static void assemble_Nodal_Tangent_Stiffness_Material(Matrix Tangent_Stiffness,
						      Mask ActiveNodes,
						      GaussPoint MPM_Mesh,
						      Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Ndim*Nnodes_mask;
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
  Tensor F_n12_p;
  Tensor transpose_F_n_p;
  Tensor transpose_F_n12_p;
  Tensor C_AB;
  Tensor Material_AB;

  Material MatProp_p;
  double m_p; /* Mass of the Gauss-Point */
  double rho_p; /* Density of the Gauss-Point */
  double V0_p; /* Volume of the Gauss-Point */
  double J_p; /* Jacobian of the deformation gradient */
  double Geometric_AB_p; /* Geometric contribution */

  /*
    Loop in the particles for the assembling process
  */
  for(int p = 0 ; p<Np ; p++)
    {
      /*
      	Get the value of the density 
      */
      rho_p = MPM_Mesh.Phi.rho.nV[p];

      /*
	Get the value of the mass 
      */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      
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
      F_n12_p = Convex_combination__TensorLib__(F_n1_p,F_n_p,0.5);
      transpose_F_n_p = transpose__TensorLib__(F_n_p);
      transpose_F_n12_p = transpose__TensorLib__(F_n12_p);

      /*
	Compute the volume of the particle in the reference configuration 
      */
      J_p = I3__TensorLib__(F_n_p);
      V0_p = (1/J_p)*(m_p/rho_p);

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
		Get the nodal contribution of the material mass matrix
	      */
	      C_AB = compute_stiffness_density(GRADIENT_pA,GRADIENT_pB,MatProp_p);

	      /*
		Compute the nodal matrix with the contribution to each degree of freedom
	       */
	      Material_AB = compute_Nodal_Tangent_Stiffness_Material(F_n12_p,C_AB,
								     transpose_F_n12_p);
	      
	      /*
		Add the geometric contribution to each dof for the assembling process
	      */
	      for(int i = 0 ; i<Ndim ; i++)
		{
		  for(int j = 0 ; j<Ndim ; j++)
		    {
		      idx_AB_mask_ij = B_mask + (i*Order + j)*Nnodes_mask + A_mask*Order;
		      
		      Tangent_Stiffness.nV[idx_AB_mask_ij] += Material_AB.N[i][j]*V0_p;
		    }
		}

	      /*
		Free memory 
	      */
	      free__TensorLib__(GRADIENT_pB);
	      free__TensorLib__(C_AB);
	      free__TensorLib__(Material_AB);
	    }

	  /*
	    Free memory 
	  */
	  free__TensorLib__(GRADIENT_pA);	  
	}
      

      /* 
	 Free memory 
      */
      free__TensorLib__(F_n12_p);
      free__TensorLib__(transpose_F_n_p);
      free__TensorLib__(transpose_F_n12_p);
      free__MatrixLib__(gradient_p);
      free(Nodes_p.Connectivity);

    }

}

/**************************************************************/

static Tensor compute_stiffness_density(Tensor GRADIENT_pA,
					Tensor GRADIENT_pB,
					Material MatProp_p)
{

  Tensor C_AB;

  if(strcmp(MatProp_p.Type,"Saint-Venant-Kirchhoff") == 0)
    {
      C_AB = compute_stiffness_density_Saint_Venant_Kirchhoff(GRADIENT_pA,
							      GRADIENT_pB,
							      MatProp_p);
    }
  else
    {
      fprintf(stderr,"%s : %s %s %s \n",
	      "Error in compute_stiffness_density()",
	      "The material",MatProp_p.Type,"has not been yet implemnented");
      exit(EXIT_FAILURE);
    }

  return C_AB;
}

/**************************************************************/

static Tensor compute_Nodal_Tangent_Stiffness_Material(Tensor F_n12_p,
						       Tensor C_AB,
						       Tensor Ft_beta_p)
{
  int Ndim = NumberDimensions;
  Tensor F_x_C_x_Ft = alloc__TensorLib__(2);
  Tensor C_x_Ft = alloc__TensorLib__(2);
  
  C_x_Ft = matrix_product__TensorLib__(C_AB, Ft_beta_p);

  F_x_C_x_Ft = matrix_product__TensorLib__(F_n12_p, C_x_Ft);

  free__TensorLib__(C_x_Ft);
  
  return F_x_C_x_Ft;
  
}

/**************************************************************/

static void update_D_Displacement(Matrix D_Displacement,
				  Matrix Tangent_Stiffness,
				  Matrix Effective_Mass,
				  Matrix Residual,
				  double Dt)
/*
  This function is deboted to update the vector with the nodal displacement by
  solving :
  
  ((2/Dt^2)*M_AB + K_AB)*delta_B + R_A = 0_A 
  
*/
{
  int Nnodes_mask = Residual.N_cols;
  int Ndof = Residual.N_rows;
  int Order = Nnodes_mask*Ndof;
  int LDA   = Nnodes_mask*Ndof;
  int LDB = Nnodes_mask*Ndof;
  char  TRANS = 'N'; /* (No transpose) */
  int   INFO = 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;
  
  Matrix Global_Matrix = allocZ__MatrixLib__(Order,Order);

  /*
    Compute the adition of the mass matrix and the tangent stifness matrix
   */
  for(int idx_AB_ij ; idx_AB_ij<Order*Order ; idx_AB_ij++)
    {
      Global_Matrix.nV[idx_AB_ij] =
	(2/DSQR(Dt))*Effective_Mass.nV[idx_AB_ij] + Tangent_Stiffness.nV[idx_AB_ij];
    }  

  /*
    Compute the LU factorization 
  */
  dgetrf_(&Order,&Order,Global_Matrix.nV,&LDA,IPIV,&INFO);

  /*
    Check error messages in the LAPACK LU descompistion  
   */
  if(INFO)
    {
      fprintf(stderr,"%s : %s %s %s \n",
	      "Error in update_D_Displacement",
	      "The function",
	      "dgetrf_",
	      "returned an error message !!!" );
      exit(EXIT_FAILURE);
    }

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order,&NRHS,Global_Matrix.nV,&LDA,IPIV,Residual.nV,&LDB,&INFO);
  free(IPIV);
  
  /*
    Check error messages in the LAPACK solver  
  */
  if(INFO)
    {
      fprintf(stderr,"%s : %s %s %s \n",
	      "Error in update_D_Displacement",
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
      D_Displacement.nV[idx_A_i] -= Residual.nV[idx_A_i];
    }

  /*
    Free auxiliar global matrix
   */
  free__MatrixLib__(Global_Matrix);
  
}


/**************************************************************/

static void imposed_displacements(Matrix D_Displacement,
				  Mask ActiveNodes,
				  Mesh FEM_Mesh,
				  int TimeStep)
/*
  Apply the boundary conditions over the nodes 
*/
{

  /* Define auxilar variables */
  int NumBoundaryConditions = FEM_Mesh.Bounds.NumBounds;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int A_BCC; /* Index of the node where we apply the BCC */
  int A_mask_BCC; /* Index of the node where we apply the BCC */
  
  for(int i_boundary = 0 ; i_boundary<NumBoundaryConditions ; i_boundary++)
    {
      
      /*
	Get the number of nodes of this boundary and 
	the number of dimensions where it is applied
      */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i_boundary].NumNodes;
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i_boundary].Dim;
    
    for(int A = 0 ; A<NumNodesBound ; A++)
      {
	
      /* 
	 Get the index of the node and get the mask node value 
      */
      A_BCC = FEM_Mesh.Bounds.BCC_i[i_boundary].Nodes[A];
      A_mask_BCC = ActiveNodes.Nodes2Mask[A_BCC];

      for(int i_dim = 0 ; i_dim<NumDimBound ; i_dim++)
	{
	  
	if(FEM_Mesh.Bounds.BCC_i[i_boundary].Dir[i_dim] == 1)
	  {
	  /*
	    Check if the curve it is on time 
	  */
	  if( (TimeStep < 0) ||
	      (TimeStep > FEM_Mesh.Bounds.BCC_i[i_boundary].Value[i_dim].Num))
	    {
	      printf("%s : %s \n",
		     "Error in imposed_displacements()",
		     "The time step is out of the curve !!");
	      exit(EXIT_FAILURE);
	    }
	  /* 
	     Assign the boundary condition 
	  */
	  D_Displacement.nM[i_dim][A_mask_BCC] =
	    FEM_Mesh.Bounds.BCC_i[i_boundary].Value[i_dim].Fx[TimeStep]*
	    (double)FEM_Mesh.Bounds.BCC_i[i_boundary].Dir[i_dim];
	  
	  }
	}
      }    
    }
  
}

/**************************************************************/


static Matrix compute_D_Velocity(Matrix Velocity,
				 Matrix D_Displacement,
				 double Dt)
{
  int Nnodes_mask = Velocity.N_cols;
  int Ndim = NumberDimensions;
  
  Matrix D_Velocity = allocZ__MatrixLib__(Ndim,Nnodes_mask);

  /*
    Compute the velocity in the midd-point 
   */
  for(int idx_A_i = 0 ; idx_A_i < Nnodes_mask*Ndim ; idx_A_i++)
    {	
      D_Velocity.nV[idx_A_i] = 2*(D_Displacement.nV[idx_A_i]/Dt - Velocity.nV[idx_A_i]);
    }


  return D_Velocity;
}


/**************************************************************/

static void update_Particles(Matrix D_Displacement,
			     Matrix D_Velocity,
			     GaussPoint MPM_Mesh,
			     Mesh FEM_Mesh,
			     Mask ActiveNodes,
			     double DeltaTimeStep)
{
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Ap;
  int A_mask;
  int idx_A_mask_i;
  int idx_ij;
  Matrix D_Displacement_Ap;
  Matrix ShapeFunction_p; /* Value of the shape-function in the particle */
  Matrix gradient_p;
  double ShapeFunction_pI; /* Nodal value for the particle */
  Tensor F_n_p;
  Tensor F_n1_p;
  Tensor f_n1_p;
  double J_p;
  Element Nodes_p; /* Element for each particle */

  /* iterate over the particles */
  for(int p = 0 ; p<Np ; p++)
    {
      
      /* Define element of the particle */
      Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);
      
      /*
	Get the nodal increment of displacement using the mask
      */
      D_Displacement_Ap = get_set_field__MeshTools__(D_Displacement, Nodes_p, ActiveNodes);
      
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

      /*
      	Compute the increment of the deformation gradient
      */
      f_n1_p = increment_Deformation_Gradient__Particles__(D_Displacement_Ap,gradient_p);

      /*
	Update the deformation gradient in t = n + 1
      */  
      update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, f_n1_p);

      /*
	Update density with the jacobian of the increment deformation gradient
       */
      J_p = I3__TensorLib__(f_n1_p);
      MPM_Mesh.Phi.rho.nV[p] = MPM_Mesh.Phi.rho.nV[p]/J_p;

      /*
	Replace the deformation gradient at t = n with the converged deformation gradient
      */
      for(int i = 0 ; i<Ndim  ; i++)
	{
	  for(int j = 0 ; j<Ndim  ; j++)
	    {
	      F_n_p.N[i][j] = F_n1_p.N[i][j];
	    }
	}      
      
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
	      idx_A_mask_i = A_mask + i*Nnodes_mask;
	      MPM_Mesh.Phi.vel.nM[p][i]  += ShapeFunction_pI*D_Velocity.nV[idx_A_mask_i];
	      MPM_Mesh.Phi.x_GC.nM[p][i] += ShapeFunction_pI*D_Displacement.nV[idx_A_mask_i];
	    } 
	}

      /*
	Free memory
       */
      free(Nodes_p.Connectivity);
      free__MatrixLib__(ShapeFunction_p);
      free__MatrixLib__(gradient_p);
      free__TensorLib__(f_n1_p);
    }  
}

/**************************************************************/
