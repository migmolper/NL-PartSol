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
static Matrix compute_Nodal_Internal_Forces(Matrix, Matrix,Mask,GaussPoint, Mesh);
static Matrix compute_Nodal_Residual(Matrix, Matrix, Matrix, Matrix, double);
static bool check_convergence(Matrix,double,int,int);
static void assemble_Nodal_Tangent_Stiffness_2D(Matrix, Matrix, Matrix, Matrix,
					       Mask, GaussPoint, Mesh);
static void assemble_Nodal_Tangent_Stiffness_Geometric_2D(Matrix, Matrix,
							 Mask, GaussPoint, Mesh);
static void assemble_Nodal_Tangent_Stiffness_Material_2D(Matrix, Matrix, Matrix, Matrix,
							Mask, GaussPoint, Mesh);
static Tensor compute_Nodal_Tangent_Stiffness_Material(Tensor,Tensor,Tensor);

/**************************************************************/

void U_Discrete_Energy_Momentum(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)
{

  /*!
    Integer variables 
  */
  int Ndim = NumberDimensions;
  int Nactivenodes;
  /*!
    Auxiliar variables for the solver
  */
  Matrix Effective_Mass;
  Matrix Stiffness_Uxx;
  Matrix Stiffness_Uyy;
  Matrix Stiffness_Uxy;
  Matrix Stiffness_Uyx;
  Matrix Forces;
  Matrix Momentum;
  Matrix Velocity;
  Matrix DeltaU;
  Matrix Residual;
  Mask ActiveNodes;
  double NormResidual;
  double TOL = 0.0001;
  double epsilon = 0.9;
  double DeltaTimeStep = 1;
  bool Not_Convergence;
  int Iter = 0;
  int MaxIter = 100;
  
  ActiveNodes = generate_NodalMask(FEM_Mesh);

  Nactivenodes = ActiveNodes.Nactivenodes;
  
  Effective_Mass = compute_Nodal_Effective_Mass(MPM_Mesh,FEM_Mesh,ActiveNodes,epsilon);

  Momentum = compute_Nodal_Momentum(MPM_Mesh, FEM_Mesh, ActiveNodes);

  Velocity = compute_Nodal_Velocity(Effective_Mass, Momentum);

  DeltaU = MatAllocZ(Ndim,Nactivenodes);

  Not_Convergence = false;

  /* while(Not_Convergence) */
  /*   { */
      
  update_Local_State(DeltaU,ActiveNodes,MPM_Mesh,FEM_Mesh,DeltaTimeStep);

  Forces = MatAllocZ(Ndim,Nactivenodes);
      
  Forces = compute_Nodal_Internal_Forces(Forces,DeltaU,ActiveNodes,MPM_Mesh,FEM_Mesh);
      
  Residual = compute_Nodal_Residual(Velocity,Forces,DeltaU,Effective_Mass,DeltaTimeStep);

  Not_Convergence = check_convergence(Residual,TOL,Iter,MaxIter);

  if(Not_Convergence)
    {

      Stiffness_Uxx = MatAllocZ(Nactivenodes, Nactivenodes);
      Stiffness_Uxy = MatAllocZ(Nactivenodes, Nactivenodes);
      Stiffness_Uyy = MatAllocZ(Nactivenodes, Nactivenodes);
      Stiffness_Uyx = MatAllocZ(Nactivenodes, Nactivenodes);
  
      assemble_Nodal_Tangent_Stiffness_2D(Stiffness_Uxx,Stiffness_Uxy,
					  Stiffness_Uyy,Stiffness_Uyx,
					  ActiveNodes,MPM_Mesh,FEM_Mesh);
	  
      /* update_DeltaU(DeltaU,Stiffness,Residual,TimeStep); */
	  
      Iter++;

      FreeMat(Forces);
      FreeMat(Residual);
      FreeMat(Stiffness_Uxx);
      FreeMat(Stiffness_Uxy);
      FreeMat(Stiffness_Uyx);
      FreeMat(Stiffness_Uyy);
    }
  else
    {
      FreeMat(Residual);
    }
            
  /*   } */
  
  /*!
    Free memory.
  */
  FreeMat(Effective_Mass); 
  FreeMat(Velocity);
  FreeMat(DeltaU);
  
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
  This function computes the consistent mass matrix 
 */
{

  int Nnodes = ActiveNodes.Nactivenodes;
  int Np = MPM_Mesh.NumGP;
  int Ap, Bp, A_mask, B_mask;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;  

  /* Evaluation of the particle in the node */
  double ShapeFunction_pA, ShapeFunction_pB;
  /* Mass of the particle */
  double m_p;
  /* Element for each particle */
  Element Nodes_p;

  /* Define and allocate the consistent mass matrix */
  Matrix Consistent_MassMatrix = MatAllocZ(Nnodes, Nnodes);

  /* Define and allocate the lumped mass matrix */
  Matrix Lumped_MassMatrix = MatAllocZ(Nnodes, 1);

  /* Define the effective mass matrix */
  Matrix Effective_MassMatrix;
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
	    Get the lumped mass matrix
	  */
	  Lumped_MassMatrix.nV[A_mask] += m_p*ShapeFunction_pA;
	  	  
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
		 Get the consistent mass matrix 
	      */
	      Consistent_MassMatrix.nM[A_mask][B_mask] +=
		m_p*ShapeFunction_pA*ShapeFunction_pB;
	  
	    }
	}

      /* Free the value of the shape functions */
      FreeMat(ShapeFunction_p);
      free(Nodes_p.Connectivity);      
      
    }

  /*
    Compute the effective mass matrix as : 
    a CONVEX combination of the consistent and lumped mass matrix
   */
  Effective_MassMatrix = Consistent_MassMatrix;
  for(int A = 0 ; A<Nnodes ; A++)
    {
      for(int B = 0 ; B<Nnodes ; B++)
	{
	  if(A != B){
	    Effective_MassMatrix.nM[A][B] = (1-epsilon)*Effective_MassMatrix.nM[A][B];
	  }
	  else{
	    Effective_MassMatrix.nM[A][B] =
	      (1-epsilon)*Effective_MassMatrix.nM[A][B] + epsilon*Lumped_MassMatrix.nV[A];
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
    
  /* Iterate over the GP to get the nodal values */
  for(int p = 0 ; p<Np ; p++)
    {

      /* Define element of the GP */
      Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /* Evaluate the shape function in the coordinates of the GP */
      ShapeFunction_p = compute_ShapeFunction(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Get the mass of the GP */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      /* Get the nodal mass and mommentum */
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
  int NRHS = Ndim;

  /* Compute the LU factorization */
  LAPACK_dgetrf(&Order,&Order,Mass.nV,&LDA,IPIV,&INFO);

  dgetrs_(&TRANS,&Order,&NRHS,Mass.nV,&LDA,IPIV,Momentum.nV,&LDB,&INFO);

  Velocity = Momentum;

  /*
    Add some usefulll info
  */
  strcpy(Velocity.Info,"Nodal-Velocity");

  return Velocity;
}

/**************************************************************/

static void update_Local_State(Matrix DeltaU, Mask ActiveNodes,
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
  Matrix DeltaU_Ap;
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
      DeltaU_Ap = get_Nodal_Values_for_Particle(DeltaU, Nodes_p, ActiveNodes);

      PrintMatrix(DeltaU_Ap, DeltaU_Ap.N_rows, DeltaU_Ap.N_cols);

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
      compute_Strain_Deformation_Gradient_n1(F_n1_p,F_n_p,DeltaU_Ap,gradient_p);
      
      /*
      	Update the second Piola-Kirchhoff stress tensor (S) with an apropiate
	intgration rule.
      */
      MatIndx_p = MPM_Mesh.MatIdx[p];
      MatProp_p = MPM_Mesh.Mat[MatIndx_p];
      S_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2);      
      S_p = Itegration_Stress_Average_Strain(S_p,F_n1_p,F_n_p,MatProp_p);
      
      /* Free the gradient */
      FreeMat(DeltaU_Ap);
      FreeMat(gradient_p);
	  
    }
  
}

/**************************************************************/

static Matrix compute_Nodal_Internal_Forces(Matrix Forces, Matrix DeltaU,
					    Mask ActiveNodes,
					    GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
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
				     Matrix DeltaU, Matrix Mass,
				     double Dt)
{
  int Ndim = NumberDimensions;
  int Nnodes = Mass.N_rows;
  Matrix Acceleration = MatAllocZ(Ndim,Nnodes);
  Matrix Inertial_Forces = MatAllocZ(Ndim,Nnodes);
  Matrix Residual = MatAllocZ(Ndim,Nnodes);
  int idx_A;
  int idx_B;
  int idx_AB;

  /*
    Compute nodal acceleration (Vectorized)
  */
  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int B = 0 ; B<Nnodes ; B++)
	{
	  idx_B = B + i*Nnodes;
	  Acceleration.nV[idx_B] = (2/Dt*Dt)*(DeltaU.nV[idx_B] - Dt*Velocity.nV[idx_B]);
	}
    }

  /*
    Compute inertial forces (Vectorized)
  */
  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int A = 0 ; A < Nnodes ; A++)
	{
	  idx_A = A + i*Nnodes;	  
	  for(int B = 0 ; B < Nnodes ; B++)
	    {
	      idx_B = B + i*Nnodes;
	      idx_AB = A + B*Nnodes;
	      Inertial_Forces.nV[idx_A] += Mass.nV[idx_AB]*Acceleration.nV[idx_B];
	    }
	}
    }

  /*
    Compute residual (Vectorized)
  */
  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int A = 0 ; A < Nnodes ; A++)
	{
	  idx_A = A + i*Nnodes;
	  Residual.nV[idx_A] = Inertial_Forces.nV[idx_A] + Forces.nV[idx_A];
	}
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

static void assemble_Nodal_Tangent_Stiffness_2D(Matrix Stiffness_Uxx,
						Matrix Stiffness_Uxy,
						Matrix Stiffness_Uyy,
						Matrix Stiffness_Uyx,
						Mask ActiveNodes,
						GaussPoint MPM_Mesh,
						Mesh FEM_Mesh)
{

  int Nnodes = ActiveNodes.Nactivenodes;

  /*
    Compute term related to the geometric non-linearities
   */  
  assemble_Nodal_Tangent_Stiffness_Geometric_2D(Stiffness_Uxx,Stiffness_Uyy,
  						ActiveNodes,
  						MPM_Mesh,FEM_Mesh);
  
  /*
    Compute term related to the material non-linearities
   */
  assemble_Nodal_Tangent_Stiffness_Material_2D(Stiffness_Uxx,
  					       Stiffness_Uxy,
  					       Stiffness_Uyy,
  					       Stiffness_Uyx,
  					       ActiveNodes,
  					       MPM_Mesh,FEM_Mesh);

}


/**************************************************************/

static void assemble_Nodal_Tangent_Stiffness_Geometric_2D(Matrix Stiffness_Uxx,
							  Matrix Stiffness_Uyy,
							  Mask ActiveNodes,
							  GaussPoint MPM_Mesh,
							  Mesh FEM_Mesh)
/*
  Introduce the geometric contribution G_AB to the full stiffness matrix K_AB
*/
{
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int Bp;
  int B_mask;
  int Nn;

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
		Add the geometric contribution
	       */
	      Stiffness_Uxx.nM[A_mask][B_mask] += Geometric_AB_p*V0_p;
	      Stiffness_Uyy.nM[A_mask][B_mask] += Geometric_AB_p*V0_p;

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
static void assemble_Nodal_Tangent_Stiffness_Material_2D(Matrix Stiffness_Uxx,
							 Matrix Stiffness_Uxy,
							 Matrix Stiffness_Uyy,
							 Matrix Stiffness_Uyx,
							 Mask ActiveNodes,
							 GaussPoint MPM_Mesh,
							 Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Ap;
  int A_mask;
  int Bp;
  int B_mask;
  int Nn;
  int MatIndx_p;

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
		Add the geometric contribution
	       */
	      Stiffness_Uxx.nM[A_mask][B_mask] += Material_AB.N[0][0]*V0_p;
	      Stiffness_Uxy.nM[A_mask][B_mask] += Material_AB.N[0][1]*V0_p;
	      Stiffness_Uyy.nM[A_mask][B_mask] += Material_AB.N[1][1]*V0_p;
	      Stiffness_Uyx.nM[A_mask][B_mask] += Material_AB.N[1][0]*V0_p;

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
