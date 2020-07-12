#include "nl-partsol.h"
#include "lapacke.h"

/*
  Auxiliar functions 
*/
static Mask   generate_NodalMask(Mesh);
static Matrix get_Nodal_Values_for_Particle(Matrix, Element, Mask);
static Matrix compute_Nodal_Effective_Mass(GaussPoint, Mesh, Mask, double);
static Matrix compute_Nodal_Momentum(GaussPoint, Mesh, Mask);
static Matrix compute_Nodal_Velocity(Matrix, Matrix);
static void   update_Local_State(Matrix,Mask,GaussPoint,Mesh,double);
static Matrix compute_Nodal_Internal_Forces(Matrix,Mask,GaussPoint, Mesh);
static Matrix compute_Nodal_Residual(Matrix, Matrix, Matrix, Matrix, double);
static bool   check_convergence(Matrix,double,int,int);
static Matrix assemble_Nodal_Tangent_Stiffness_2D(Mask, GaussPoint, Mesh);
static void   assemble_Nodal_Tangent_Stiffness_Geometric_2D(Matrix,  Mask, GaussPoint, Mesh);
static void   assemble_Nodal_Tangent_Stiffness_Material_2D(Matrix, Mask, GaussPoint, Mesh);
static Tensor compute_Nodal_Tangent_Stiffness_Material(Tensor,Tensor,Tensor);
static void   update_D_Displacement(Matrix, Matrix, Matrix, Matrix, double);
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
  double TOL = 0.0001;
  double epsilon = 1.0;
  double DeltaTimeStep = 1;
  bool Not_Convergence;
  int Iter = 0;
  int MaxIter = 100;

  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {
      /*
	With the active set of nodes generate a mask to help the algorithm to compute
	the equilibrium only in the active nodes
      */
      ActiveNodes = generate_NodalMask(FEM_Mesh);
      Nactivenodes = ActiveNodes.Nactivenodes;

      /*
	Compute the effective mass matrix as a convex combination of the consistent mass
	matrix and the lumped mass matrix.
      */
      Effective_Mass = compute_Nodal_Effective_Mass(MPM_Mesh,FEM_Mesh,
						    ActiveNodes,epsilon);

      /*
	Compute the nodal value of the momentum
      */
      Momentum = compute_Nodal_Momentum(MPM_Mesh, FEM_Mesh, ActiveNodes);

      /*
	Compute the nodal valocities with the effective mass matrix and the nodal momentum
      */
      Velocity = compute_Nodal_Velocity(Effective_Mass, Momentum);

      /*
	Set to zero the increment of velocity and displacement for nodal values
      */
      D_Displacement = MatAllocZ(Ndim,Nactivenodes);  
      D_Velocity = MatAllocZ(Ndim,Nactivenodes);

      /*
	Set the convergence false by default and start the iterations to compute
	the incement of velocity and displacement in the nodes of the mesh
      */
      Not_Convergence = false;
      while(Not_Convergence)
	{

	  /*
	    Compute the stress-strain state for each particle
	  */
	  update_Local_State(D_Displacement,ActiveNodes,
			     MPM_Mesh,FEM_Mesh,DeltaTimeStep);

	  /*
	    Compute the nodal forces
	  */
	  Forces = compute_Nodal_Internal_Forces(D_Displacement,ActiveNodes,
						 MPM_Mesh,FEM_Mesh);

	  /*
	    Compute the numerical residual to check the equilibrium
	  */
	  Residual = compute_Nodal_Residual(Velocity,Forces,D_Displacement,
					    Effective_Mass,DeltaTimeStep);

	  /*
	    If the norm of the residual for each nodal value is below a tolerace
	    skip
	  */
	  Not_Convergence = check_convergence(Residual,TOL,Iter,MaxIter);

	  /*
	    If not, solve the linearized equilibrium to compute the next step
	  */
	  if(Not_Convergence)
	    {

	      /*
		Assemble the tangent stiffness matrix as a sum of the material tangent
		stiffness and the geometrical tangent stiffness.
	      */
	      Tangent_Stiffness = assemble_Nodal_Tangent_Stiffness_2D(ActiveNodes,
								      MPM_Mesh,FEM_Mesh);

	      /*
		Solve the resulting equation 
	      */
	      update_D_Displacement(D_Displacement,Tangent_Stiffness,
				    Effective_Mass,Residual,DeltaTimeStep);

	      /*
		Update the iteration number
	      */
	      Iter++;

	      FreeMat(Forces);
	      FreeMat(Residual);
	      FreeMat(Tangent_Stiffness);
	    }
            
	}

      /*
	Once the equilibrium is reached, obtain the increment of nodal velocity
      */
      D_Velocity = compute_D_Velocity(Velocity,D_Displacement,DeltaTimeStep);

      /*
	Update Lagrangians with D_Displacement
      */
      update_Particles(D_Displacement,D_Velocity,MPM_Mesh,FEM_Mesh,ActiveNodes,DeltaTimeStep);

      /*
	Reload the connectivity information for each particle
      */
      LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);

      /*
	Outputs
      */
      if(TimeStep % ResultsTimeStep == 0)
	{
	  /* WriteVtk_FEM("Mesh",FEM_Mesh,R_I,TimeStep); */
	  WriteVtk_MPM("MPM_VALUES",MPM_Mesh,"ALL",TimeStep,ResultsTimeStep);
	}
  
      /*
	Free memory.
      */
      FreeMat(Effective_Mass); 
      FreeMat(Velocity);
      FreeMat(D_Velocity);
      FreeMat(D_Displacement);
    }
  
}

/**************************************************************/

static Mask generate_NodalMask(Mesh FEM_Mesh)
{
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Nactivenodes = 0;
  int * Nodes2Mask = (int *)Allocate_ArrayZ(Nnodes,sizeof(int));
  ChainPtr Mask2Nodes = NULL;
  Mask M;

  for(int A = 0 ; A<Nnodes ; A++)
    {
      if(FEM_Mesh.NumParticles[A] > 0)
	{
	  Nodes2Mask[A] = Nactivenodes;
	  push_to_Set(&Mask2Nodes,I);
	  Nactivenodes++;
	}
      else
	{
	  Nodes2Mask[A] = - 1;
	}
    }

  M.Nactivenodes = Nactivenodes;
  M.Mask2Nodes = Set_to_Pointer(Mask2Nodes,Nactivenodes);
  M.Nodes2Mask = Nodes2Mask;  
 
  return M;
}

/**************************************************************/

static Matrix get_Nodal_Values_for_Particle(Matrix Field, Element Nodes_p, Mask ActiveNodes)
/*
  This function performs two operations. First takes the nodal connectivity of the particle, 
  and translate it to the mask numeration. Second, generate a Matrix with the nodal values.
  To help in the future computations. Nodal data is substracted in the shape (nodesxndofs).
 */
  {
    int Nnodes = Nodes_p.NumberNodes;
    int Ndof = Field.N_rows;
    Matrix Field_Ap = MatAllocZ(Nnodes,Ndof);
    int Ap;
    int A_mask;

    if(Ndof > 2)
      {
	for(int A = 0 ; A<Nnodes ; A++)
	  {
	    
	    /* 
	       Get the node in the mass matrix with the mask
	    */
	    Ap = Nodes_p.Connectivity[A];
	    A_mask = ActiveNodes.Nodes2Mask[Ap];
	
	    for(int i = 0 ; i<Ndof ; i++)
	      {
		Field_Ap.nM[A][i] = Field.nM[i][A_mask];
	      }
	  }
      }
    else
      {
	for(int A = 0 ; A<Nnodes ; A++)
	  {

	    /* 
	       Get the node in the mass matrix with the mask
	    */
	    Ap = Nodes_p.Connectivity[A];
	    A_mask = ActiveNodes.Nodes2Mask[Ap];
	    
	    Field_Ap.nV[A] = Field.nV[A_mask];
	  }
      }
    
    return Field_Ap;
  }


/**************************************************************/

static Matrix compute_Nodal_Effective_Mass(GaussPoint MPM_Mesh, Mesh FEM_Mesh,
					   Mask ActiveNodes, double epsilon)
/*
  This function computes the effective mass matrix 
 */
{

  int Nnodes = ActiveNodes.Nactivenodes;
  int Ndof = NumberDOF;
  int Np = MPM_Mesh.NumGP;
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
  Matrix Effective_MassMatrix = MatAllocZ(Nnodes*Ndof, Nnodes*Ndof);

  /* Define and allocate the lumped mass matrix */
  Matrix Lumped_MassMatrix = MatAllocZ(Nnodes*Ndof, 1);

  /*
    Iterate over the particles to get the nodal values 
  */
  for(int p = 0 ; p<Np ; p++)
    {

      /*
	Define tributary nodes of the particle 
      */
      Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /* 
	 Evaluate the shape function in the coordinates of the particle
       */
      ShapeFunction_p = compute_ShapeFunction(Nodes_p, MPM_Mesh, FEM_Mesh);
   
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
	      idx_A_mask_i = A_mask + i*Nnodes;
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
		  idx_AB_mask_i = A_mask + i*Nnodes + (B_mask + i*Nnodes)*Ndof*Nnodes;
		  
		  Effective_MassMatrix.nV[idx_AB_mask_i] += m_AB_p;
		}
	  
	    }
	}

      /* Free the value of the shape functions */
      FreeMat(ShapeFunction_p);
      free(Nodes_p.Connectivity);      
      
    }

  /*
    At this point the effective mass matrix coincides with the consistent mass
    matrix. We can tune it by a convecx combination with the lumped mass matrix
  */
  for(int A = 0 ; A<Nnodes ; A++)
    {
      for(int B = 0 ; B<Nnodes ; B++)
	{
	  for(int i = 0 ; i<Ndof ; i++)
	    {
	      idx_AB_mask_i = A + i*Nnodes + (B + i*Nnodes)*Ndof*Nnodes;
	      idx_A_mask_i = A + i*Nnodes;
	      
	      Effective_MassMatrix.nV[idx_AB_mask_i] =
		(1-epsilon)*Effective_MassMatrix.nV[idx_AB_mask_i] +
		(A == B)*epsilon*Lumped_MassMatrix.nV[idx_A_mask_i];
	    }
	}
    }

  /*
    Free lumped mass matrix.
   */
  FreeMat(Lumped_MassMatrix);

  /* 
     Add some usefulll info 
  */
  strcpy(Effective_MassMatrix.Info,"Effective-Mass-Matrix");

  return Effective_MassMatrix; 
}

/**************************************************************/

static Matrix compute_Nodal_Momentum(GaussPoint MPM_Mesh,Mesh FEM_Mesh, Mask ActiveNodes)
/*

*/
{
  int Ndim = NumberDimensions;
  int Nnodes = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;

  /* Evaluation of the particle in the node */
  double ShapeFunction_pA;
  /* Mass of the particle */
  double m_p;
  /* Element for each particle */
  Element Nodes_p;

 /* Define and allocate the momentum vector */
  Matrix Momentum = MatAllocZ(Ndim,Nnodes);
    
  /* Iterate over the particles to get the nodal values */
  for(int p = 0 ; p<Np ; p++)
    {

      /* Define element of the particle */
      Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /* Evaluate the shape function in the coordinates of the particle */
      ShapeFunction_p = compute_ShapeFunction(Nodes_p, MPM_Mesh, FEM_Mesh);

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
	      Momentum.nM[i][A_mask] += m_p*ShapeFunction_pA*MPM_Mesh.Phi.vel.nM[p][i];
	    }
	}

      /* Free the value of the shape functions */
      FreeMat(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }

  /*
     Add some usefulll info
  */
  strcpy(Momentum.Info,"Nodal-Momentum");
 
  return Momentum;
}

/**************************************************************/

static Matrix compute_Nodal_Velocity(Matrix Mass, Matrix Momentum)
/* 
   Call the LAPACK solver to compute the nodal velocity
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
  LAPACK_dgetrf(&Order,&Order,Mass.nV,&LDA,IPIV,&INFO);

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order,&NRHS,Mass.nV,&LDA,IPIV,Momentum.nV,&LDB,&INFO);

  Velocity = Momentum;

  /*
    Add some usefulll info
  */
  strcpy(Velocity.Info,"Nodal-Velocity");

  return Velocity;
}

/**************************************************************/

static void update_Local_State(Matrix D_Displacement, Mask ActiveNodes,
			       GaussPoint MPM_Mesh, Mesh FEM_Mesh,
			       double TimeStep)
{

  /*
    Auxiliar variables
   */
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes = ActiveNodes.Nactivenodes;
  int MatIndx_p;
  int Nnodes_p;
  double J_p;  
  Element Nodes_p;
  Material MatProp_p;
  Matrix gradient_p;
  Matrix D_Displacement_Ap;
  Tensor F_n_p;
  Tensor F_n12_p;
  Tensor F_n1_p;
  Tensor S_p;
  
  /* Loop in the material point set */
  for(int p = 0 ; p<Np ; p++)
    {
      /*
	Define tributary nodes of the particle 
      */
      Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);
      
      /*
	Get the nodal increment of displacement using the mask
       */
      D_Displacement_Ap = get_Nodal_Values_for_Particle(D_Displacement, Nodes_p, ActiveNodes);

      /*
      	 Evaluate the shape function gradient in the coordinates of the particle
       */
      gradient_p = compute_ShapeFunction_gradient(Nodes_p,MPM_Mesh,FEM_Mesh);
	  
      /*
      	Take the values of the deformation gradient from the previous step
      */
      F_n_p  = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);

      /*
      	Compute the value of the deformation gradient at t = n+1
      */
      compute_Strain_Deformation_Gradient_n1(F_n1_p,F_n_p,D_Displacement_Ap,gradient_p);
      
      /*
      	Update the second Piola-Kirchhoff stress tensor (S) with an apropiate
	intgration rule.
      */
      MatIndx_p = MPM_Mesh.MatIdx[p];
      MatProp_p = MPM_Mesh.Mat[MatIndx_p];
      S_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2);      
      S_p = Itegration_Stress_Average_Strain(S_p,F_n1_p,F_n_p,MatProp_p);
      
      /* Free the gradient */
      FreeMat(D_Displacement_Ap);
      FreeMat(gradient_p);
	  
    }
  
}

/**************************************************************/

static Matrix compute_Nodal_Internal_Forces(Matrix D_Displacement,
					    Mask ActiveNodes,
					    GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Nactivenodes = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int Nn;

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

  Matrix Forces = MatAllocZ(Ndim,Nactivenodes);

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
      Nn = MPM_Mesh.NumberNodes[p];
      Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], Nn);

      /*
	Compute gradient of the shape function in each node 
      */
      gradient_p = compute_ShapeFunction_gradient(Nodes_p, MPM_Mesh, FEM_Mesh);
    	  
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
    
      for(int A = 0 ; A<Nn ; A++)
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
	      Forces.nM[i][A_mask] += InternalForcesDensity_Ap.n[i]*V0_p;
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
      FreeMat(gradient_p);
      free(Nodes_p.Connectivity);
    }

  return Forces;
    
}

/**************************************************************/

static Matrix compute_Nodal_Residual(Matrix Velocity, Matrix Forces,
				     Matrix D_Displacement, Matrix Mass,
				     double Dt)
{
  int Ndim = NumberDimensions;
  int Nnodes = Mass.N_rows;
  Matrix Acceleration = MatAllocZ(Ndim,Nnodes);
  Matrix Inertial_Forces = MatAllocZ(Ndim,Nnodes);
  Matrix Residual = MatAllocZ(Ndim,Nnodes);
  int idx_AB;

  /*
    Compute nodal acceleration (Vectorized)
  */
  for(int idx_B = 0 ; idx_B<Ndim*Nnodes ; idx_B++)
    {
      Acceleration.nV[idx_B] = (2/Dt*Dt)*(D_Displacement.nV[idx_B] - Dt*Velocity.nV[idx_B]);
    }

  /*
    Compute inertial forces (Vectorized)
  */
  for(int idx_A = 0 ; idx_A<Ndim*Nnodes ; idx_A++)
    {
      for(int idx_B = 0 ; idx_B<Ndim*Nnodes ; idx_B++)
	{
	  idx_AB = idx_A + idx_B*Ndim*Nnodes;
	  Inertial_Forces.nV[idx_A] += Mass.nV[idx_AB]*Acceleration.nV[idx_B];
	}
    }

  /*
    Compute (-) residual (Vectorized). The minus symbol is due to
    solver purposes. See compute_D_Displacement
  */
  for(int idx_A = 0 ; idx_A<Ndim*Nnodes ; idx_A++)
    {
      Residual.nV[idx_A] = - Inertial_Forces.nV[idx_A] - Forces.nV[idx_A];
    }

  /*
    Free Memory 
   */
  FreeMat(Acceleration);
  FreeMat(Inertial_Forces);

  
  return Residual;
}

/**************************************************************/
static bool check_convergence(Matrix Residual,double TOL,int Iter,int MaxIter)
{
  bool convergence;
  int Ndim = NumberDimensions;
  int Nnodes = Residual.N_cols;
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
      for(int A = 0 ; A<Nnodes ; A++)
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

static Matrix assemble_Nodal_Tangent_Stiffness_2D(Mask ActiveNodes,
						  GaussPoint MPM_Mesh,
						  Mesh FEM_Mesh)
{

  int Nnodes = ActiveNodes.Nactivenodes;
  int Ndof = NumberDOF;

  Matrix Tangent_Stiffness = MatAllocZ(Ndof*Nnodes, Ndof*Nnodes);
  
  /*
    Compute term related to the geometric non-linearities
  */  
  assemble_Nodal_Tangent_Stiffness_Geometric_2D(Tangent_Stiffness,ActiveNodes,
  						MPM_Mesh,FEM_Mesh);
  
  /*
    Compute term related to the material non-linearities
   */
  assemble_Nodal_Tangent_Stiffness_Material_2D(Tangent_Stiffness,ActiveNodes,
  					       MPM_Mesh,FEM_Mesh);

  return Tangent_Stiffness;
}


/**************************************************************/

static void assemble_Nodal_Tangent_Stiffness_Geometric_2D(Matrix Tangent_Stiffness,
							  Mask ActiveNodes,
							  GaussPoint MPM_Mesh,
							  Mesh FEM_Mesh)
/*
  Introduce the geometric contribution G_AB to the full stiffness matrix K_AB
*/
{
  int Ndim = NumberDimensions;
  int Nnodes = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int Bp;
  int B_mask;
  int Nn;
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
      Nn = MPM_Mesh.NumberNodes[p];
      Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], Nn);

      /*
	Compute gradient of the shape function in each node 
      */
      gradient_p = compute_ShapeFunction_gradient(Nodes_p, MPM_Mesh, FEM_Mesh);
    	  
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

      for(int A = 0 ; A<Nn ; A++)
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

	  
	  for(int B = 0 ; B<Nn ; B++)
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
		Add the geometric contribution of the term .xx and .yy
	      */
	      for(int i = 0 ; i<Ndim ; i++)
		{
		  /*
		    Compute the vectorized index
		  */
		  idx_AB_mask_i = A_mask + i*Nnodes + (B_mask + i*Nnodes)*Ndim*Nnodes;

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
      FreeMat(gradient_p);
      free(Nodes_p.Connectivity);
    }

}


/**************************************************************/
static void assemble_Nodal_Tangent_Stiffness_Material_2D(Matrix Tangent_Stiffness,
							 Mask ActiveNodes,
							 GaussPoint MPM_Mesh,
							 Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Nnodes = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int Bp;
  int B_mask;
  int Nn;
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
      Nn = MPM_Mesh.NumberNodes[p];
      Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], Nn);

      /*
	Compute gradient of the shape function in each node 
      */
      gradient_p = compute_ShapeFunction_gradient(Nodes_p, MPM_Mesh, FEM_Mesh);
      
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

      for(int A = 0 ; A<Nn ; A++)
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

	  
	  for(int B = 0 ; B<Nn ; B++)
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
	      C_AB = compute_stiffness_density_Saint_Venant_Kirchhoff(GRADIENT_pA,
								      GRADIENT_pB,
								      MatProp_p);

	      /*
		Compute the nodal matrix with the contribution to each degree of freedom
	       */
	      Material_AB =
		compute_Nodal_Tangent_Stiffness_Material(F_n12_p,C_AB,transpose_F_n12_p);
	      
	      /*
		Add the geometric contribution to each dof for the assembling process
	      */
	      for(int i = 0 ; i<Ndim ; i++)
		{
		  for(int j = 0 ; j<Ndim ; j++)
		    {
		      idx_AB_mask_ij = A_mask + i*Nnodes + (B_mask + j*Nnodes)*Ndim*Nnodes;
		      
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
      FreeMat(gradient_p);
      free(Nodes_p.Connectivity);

    }

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
{
  int Nnodes = Residual.N_cols;
  int Ndof = Residual.N_rows;
  int Order = Nnodes*Ndof;
  int LDA   = Nnodes*Ndof;
  int LDB = Nnodes*Ndof;
  char  TRANS = 'N'; /* (No transpose) */
  int   INFO= 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;
  
  Matrix Global_Matrix = MatAllocZ(Nnodes*Ndof,Nnodes*Ndof);


  for(int idx_AB_ij ; idx_AB_ij<Nnodes*Ndof*Nnodes*Ndof ; idx_AB_ij++)
    {
      Global_Matrix.nV[idx_AB_ij] =
	(2/DSQR(Dt))*Effective_Mass.nV[idx_AB_ij] +
	Tangent_Stiffness.nV[idx_AB_ij];
    }
  

  /*
    Compute the LU factorization 
  */
  LAPACK_dgetrf(&Order,&Order,Global_Matrix.nV,&LDA,IPIV,&INFO);

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order,&NRHS,Global_Matrix.nV,&LDA,IPIV,Residual.nV,&LDB,&INFO);

  /*
    Update 
   */
    for(int idx_A_i = 0 ; idx_A_i < Nnodes*Ndof ; idx_A_i++)
      {
	
	D_Displacement.nV[idx_A_i] += Residual.nV[idx_A_i];
      }
    
}


/**************************************************************/


static Matrix compute_D_Velocity(Matrix Velocity, Matrix D_Displacement, double Dt)
{
  int Nnodes = Velocity.N_cols;
  int Ndim = NumberDimensions;
  
  Matrix D_Velocity = MatAllocZ(Ndim,Nnodes);

  /*
    Compute the velocity in the midd-point 
   */
  for(int idx_A_i = 0 ; idx_A_i < Nnodes*Ndim ; idx_A_i++)
    {	
      D_Velocity.nV[idx_A_i] = 2*(D_Displacement.nV[idx_A_i]/Dt - Velocity.nV[idx_A_i]);
    }


  return D_Velocity;
}


/**************************************************************/

static void update_Particles(Matrix D_Displacement, Matrix D_Velocity,
			     GaussPoint MPM_Mesh,Mesh FEM_Mesh,
			     Mask ActiveNodes, double DeltaTimeStep)
{
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nnodes = ActiveNodes.Nactivenodes;
  int Ap;
  int A_mask;
  int idx_A_mask_i;
  Matrix ShapeFunction_p; /* Value of the shape-function in the particle */
  double ShapeFunction_pI; /* Nodal value for the particle */
  Element Nodes_p; /* Element for each particle */

  /* iterate over the particles */
  for(int p = 0 ; p<Np ; p++)
    {
      
      /* Define element of the particle */
      Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);
      
      /* Evaluate the shape function in the coordinates of the particle */
      ShapeFunction_p = compute_ShapeFunction(Nodes_p, MPM_Mesh, FEM_Mesh);
      
      /* Iterate over the nodes of the particle */
      for(int A = 0; A<Nodes_p.NumberNodes; A++)
	{

	  /*
	    Get the node in the nodal momentum with the mask
	  */
	  Ap = Nodes_p.Connectivity[A];
	  A_mask = ActiveNodes.Nodes2Mask[Ap];
	  
	  /* Evaluate the GP function in the node */
	  ShapeFunction_pI = ShapeFunction_p.nV[A_mask];
	  /* If this node has a null Value of the SHF continue */
	  if(fabs(ShapeFunction_pI) <= TOL_zero){
	    continue;
	  }

	  /*
	    Update velocity and position of the particles
	   */
	  for(int i = 0 ; i<Ndim ; i++)
	    {
	      idx_A_mask_i = A_mask + i*Nnodes;
	      MPM_Mesh.Phi.vel.nM[p][i] += ShapeFunction_pI*D_Velocity.nV[idx_A_mask_i];
	      MPM_Mesh.Phi.x_GC.nM[p][i] += ShapeFunction_pI*D_Displacement.nV[idx_A_mask_i];
	    } 
	}
      
      free(Nodes_p.Connectivity);
      FreeMat(ShapeFunction_p);
    }  
}

/**************************************************************/
