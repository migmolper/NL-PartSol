#include "nl-partsol.h"

/*
  Auxiliar functions 
*/
static Mask generate_NodalMask(Mesh);
static Matrix get_Nodal_Values_for_Particle(Matrix, Element, Mask);
static Matrix compute_Nodal_Effective_Mass(GaussPoint, Mesh, Mask, double);
static Matrix compute_Nodal_Momentum(GaussPoint, Mesh, Mask);
static void update_Local_State(Matrix,Mask,GaussPoint,Mesh,double);
static Matrix compute_Nodal_Internal_Forces(Matrix, Matrix,Mask,GaussPoint, Mesh);
static Matrix compute_Nodal_Velocity(Matrix, Matrix);

/**************************************************************/

void U_Discrete_Energy_Momentum(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)
{

  /*!
    Integer variables 
  */
  int Ndim = NumberDimensions;

  /*!
    Auxiliar variables for the solver
  */
  Matrix Effective_Mass;
  Matrix Stiffness;
  Matrix Forces;
  Matrix Momentum;
  Matrix Velocity;
  Matrix DeltaU;
  Matrix Residual;
  Mask ActiveNodes;
  double NormResidual;
  double TOL;
  double epsilon = 0.9;
  bool Not_Convergence;
  int Iter;
  int MaxIter;
  
  ActiveNodes = generate_NodalMask(FEM_Mesh);
  
  Effective_Mass = compute_Nodal_Effective_Mass(MPM_Mesh,FEM_Mesh,ActiveNodes,epsilon);

  Momentum = compute_Nodal_Momentum(MPM_Mesh, FEM_Mesh, ActiveNodes);

  Velocity = compute_Nodal_Velocity(Effective_Mass, Momentum);

  DeltaU = MatAllocZ(2,ActiveNodes.Nactivenodes);

  Not_Convergence = false;

  /* while(Not_Convergence) */
  /*   { */
      
  update_Local_State(DeltaU,ActiveNodes,MPM_Mesh,FEM_Mesh,DeltaTimeStep);

  Forces = MatAllocZ(2,ActiveNodes.Nactivenodes);
      
  Forces = compute_Nodal_Internal_Forces(Forces,DeltaU,ActiveNodes,MPM_Mesh,FEM_Mesh);
      
  /* Residual = compute_Nodal_Residual(Velocity,Forces,DeltaU,Effective_Mass,TimeStep); */

  /*     Not_Convergence = check_convergence(Residual,TOL,Iter,MaxIter); */

  /*     if(Not_Convergence) */
  /* 	{ */
  /* 	  Stiffness = compute_Nodal_Tangent_Stiffness(TimeStep); */
	  
  /* 	  update_DeltaU(DeltaU,Stiffness,Residual,TimeStep); */
	  
  /* 	  Iter++; */

  FreeMat(Forces);
  /* FreeMat(Residual); */
  /* 	  FreeMat(Stiffness); */
  /* 	} */
  /*     else */
  /* 	{ */
  /* 	  FreeMat(Residual);	   */
  /* 	} */
            
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

  for(int I = 0 ; I<Nnodes ; I++)
    {
      if(FEM_Mesh.NumParticles[I] > 0)
	{
	  Nodes2Mask[I] = Nactivenodes;
	  push_to_Set(&Mask2Nodes,I);
	  Nactivenodes++;
	}
      else
	{
	  Nodes2Mask[I] = - 1;
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
    Matrix Field_Ip = MatAllocZ(Nnodes,Ndof);
    int Ip;
    int I_mask;

    if(Ndof > 2)
      {
	for(int I = 0 ; I<Nnodes ; I++)
	  {
	    
	    /* 
	       Get the node in the mass matrix with the mask
	    */
	    Ip = Nodes_p.Connectivity[I];
	    I_mask = ActiveNodes.Nodes2Mask[Ip];
	
	    for(int i = 0 ; i<Ndof ; i++)
	      {
		Field_Ip.nM[I][i] = Field.nM[i][I_mask];
	      }
	  }
      }
    else
      {
	for(int I = 0 ; I<Nnodes ; I++)
	  {

	    /* 
	       Get the node in the mass matrix with the mask
	    */
	    Ip = Nodes_p.Connectivity[I];
	    I_mask = ActiveNodes.Nodes2Mask[Ip];
	    
	    Field_Ip.nV[I] = Field.nV[I_mask];
	  }
      }
    
    return Field_Ip;
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
  int Ip, Jp, I_mask, J_mask;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;  

  /* Evaluation of the particle in the node */
  double ShapeFunction_pI, ShapeFunction_pJ;
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
      

      for(int I = 0 ; I<Nodes_p.NumberNodes ; I++)
	{

	  /* 
	     Get the node in the mass matrix with the mask
	  */
	  Ip = Nodes_p.Connectivity[I];
	  I_mask = ActiveNodes.Nodes2Mask[Ip];

	  /* Get the value of the shape function */
	  ShapeFunction_pI = ShapeFunction_p.nV[I];

	  /*
	    Get the lumped mass matrix
	  */
	  Lumped_MassMatrix.nV[I_mask] += m_p*ShapeFunction_pI;
	  	  
	  for(int J = 0 ; J<Nodes_p.NumberNodes ; J++)
	    {
	      
	      /* 
		 Get the node in the mass matrix with the mask
	      */
	      Jp = Nodes_p.Connectivity[J];
	      J_mask = ActiveNodes.Nodes2Mask[Jp];

	      /* Get the value of the shape function */
	      ShapeFunction_pJ = ShapeFunction_p.nV[J];
	  
	      /* 
		 Get the consistent mass matrix 
	      */
	      Consistent_MassMatrix.nM[I_mask][J_mask] +=
		m_p*ShapeFunction_pI*ShapeFunction_pJ;
	  
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
  for(int I = 0 ; I<Nnodes ; I++)
    {
      for(int J = 0 ; J<Nnodes ; J++)
	{
	  if(I != J){
	    Effective_MassMatrix.nM[I][J] = (1-epsilon)*Effective_MassMatrix.nM[I][J];
	  }
	  else{
	    Effective_MassMatrix.nM[I][J] =
	      (1-epsilon)*Effective_MassMatrix.nM[I][J] + epsilon*Lumped_MassMatrix.nV[I];
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
  int Ip;
  int I_mask;

  /* Value of the shape-function */
  Matrix ShapeFunction_p;

  /* Evaluation of the particle in the node */
  double ShapeFunction_pI;
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
      for(int I = 0 ; I<Nodes_p.NumberNodes ; I++)
	{
      
	  /*
	     Get the node in the nodal momentum with the mask
	  */
	  Ip = Nodes_p.Connectivity[I];
	  I_mask = ActiveNodes.Nodes2Mask[Ip];
      
	  /* Evaluate the GP function in the node */
	  ShapeFunction_pI = ShapeFunction_p.nV[I];
      
	  /* Nodal momentum */
	  for(int i = 0 ; i<Ndim ; i++)
	    {
	      Momentum.nM[i][I_mask] += m_p*ShapeFunction_pI*MPM_Mesh.Phi.vel.nM[p][i];
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
  Matrix DeltaU_Ip;
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
      DeltaU_Ip = get_Nodal_Values_for_Particle(DeltaU, Nodes_p, ActiveNodes);

      PrintMatrix(DeltaU_Ip, DeltaU_Ip.N_rows, DeltaU_Ip.N_cols);

      /*
      	 Evaluate the shape function gradient in the coordinates of the particle
       */
      gradient_p = compute_ShapeFunction_gradient(Nodes_p,MPM_Mesh,FEM_Mesh);
	  
      /*
      	Take the values of the deformation gradient from the previous step
      */
      F_n_p  = memory_to_Tensor(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p = memory_to_Tensor(MPM_Mesh.Phi.F_n1.nM[p],2);

      /*
      	Compute the value of the deformation gradient at t = n+1
      */
      compute_Strain_Deformation_Gradient_n1(F_n1_p,F_n_p,DeltaU_Ip,gradient_p);
      
      /*
      	Update the second Piola-Kirchhoff stress tensor (S) with an apropiate
	intgration rule.
      */
      MatIndx_p = MPM_Mesh.MatIdx[p];
      MatProp_p = MPM_Mesh.Mat[MatIndx_p];
      S_p = memory_to_Tensor(MPM_Mesh.Phi.Stress.nM[p],2);      
      S_p = Itegration_Stress_Average_Strain(S_p,F_n1_p,F_n_p,MatProp_p);
      
      /* Free the gradient */
      FreeMat(DeltaU_Ip);
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
  int Ip;
  int I_mask;
  int Nn;

  Tensor P_p; /* First Piola-Kirchhoff Stress tensor */
  Tensor S_p; /* Second Piola-Kirchhoff Stress tensor */
  Tensor InternalForcesDensity_Ip;

  Element Nodes_p; /* List of nodes for particle */
  Matrix gradient_p; /* Shape functions gradients */
  Tensor gradient_pI;
  Tensor GRADIENT_pI;
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
      F_n_p  = memory_to_Tensor(MPM_Mesh.Phi.F_n.nM[p],2);
      F_n1_p = memory_to_Tensor(MPM_Mesh.Phi.F_n1.nM[p],2);
      F_n12_p = compute_midpoint_Tensor(F_n1_p,F_n_p,0.5);
      transpose_F_n_p = get_Transpose_Of(F_n_p);

      /*
	Compute the first Piola-Kirchhoff stress tensor
      */
      S_p = memory_to_Tensor(MPM_Mesh.Phi.Stress.nM[p], 2);
      P_p = get_matrixProduct_Of(F_n12_p, S_p);

      /*
	Compute the volume of the particle in the reference configuration 
      */
      J_p = get_I3_Of(F_n_p);
      V0_p = (1/J_p)*(m_p/rho_p);
    
      for(int I = 0 ; I<Nn ; I++)
	{
      
	  /*
	    Compute the gradient in the reference configuration 
	  */
	  gradient_pI = memory_to_Tensor(gradient_p.nM[I], 1);
	  GRADIENT_pI = get_firstOrderContraction_Of(transpose_F_n_p,gradient_pI);
      
	  /*
	    Compute the nodal forces of the particle 
	  */
	  InternalForcesDensity_Ip = get_firstOrderContraction_Of(P_p, GRADIENT_pI);
      
	  /*
	    Get the node of the mesh for the contribution 
	  */
	  Ip = Nodes_p.Connectivity[I];
	  I_mask = ActiveNodes.Nodes2Mask[Ip];
      
	  /*
	    Asign the nodal forces contribution to the node 
	  */
	  for(int i = 0 ; i<Ndim ; i++)
	    {
	      Forces.nM[i][I_mask] -= InternalForcesDensity_Ip.n[i]*V0_p;
	    }

	  /*
	    Free memory 
	  */
	  free_Tensor(InternalForcesDensity_Ip);
	  free_Tensor(GRADIENT_pI);
	}
        
      /* 
	 Free memory 
      */
      free_Tensor(F_n12_p);
      free_Tensor(transpose_F_n_p);
      free_Tensor(P_p);
      FreeMat(gradient_p);
      free(Nodes_p.Connectivity);
    }

  return Forces;
    
}
/**************************************************************/

static Matrix compute_Nodal_Velocity(Matrix Mass, Matrix Momentum)
/* 
   Call the LAPACK solver to compute the nodal velocity
*/
{

#include "lapacke.h"

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
