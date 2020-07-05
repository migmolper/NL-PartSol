#include "nl-partsol.h"

/*
  Auxiliar functions 
*/
static Mask generate_NodalMask(Mesh);
static Matrix compute_Nodal_Effective_Mass(GaussPoint, Mesh, Mask, double);
static Matrix compute_Nodal_Momentum(GaussPoint, Mesh, Mask);
static Matrix compute_Nodal_Velocity(Matrix, Matrix);
static void update_Local_State(Matrix,Mask,GaussPoint,Mesh,double);

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
  bool Convergence;
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
      
  /*     update_Local_State(DeltaU,ActiveNodes,FEM_Mesh,MPM_Mesh,TimeStep); */
      
  /*     Forces = compute_Nodal_Forces(); */
      
  /*     Residual = compute_Nodal_Residual(Velocity,Forces,DeltaU,Effective_Mass,TimeStep); */

  /*     Not_Convergence = check_convergence(Residual,TOL,Iter,MaxIter); */

  /*     if(Not_Convergence) */
  /* 	{ */
  /* 	  Stiffness = compute_Nodal_Tangent_Stiffness(TimeStep); */
	  
  /* 	  update_DeltaU(DeltaU,Stiffness,Residual,TimeStep); */
	  
  /* 	  Iter++; */

  /* 	  FreeMat(Forces); */
  /* 	  FreeMat(Stiffness); */
  /* 	  FreeMat(Residual); */
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
  FreeMat(Momentum);
  FreeMat(DeltaU);
  FreeMat(Forces);
  
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
  Tensor F_n_p;
  Tensor F_n12_p;
  Tensor F_n1_p;
  Tensor C_n_p;
  Tensor C_n1_p;
  Tensor PK2_p;
  
  /* Loop in the material point set */
  for(int p = 0 ; p<Np ; p++)
    {
      /*
	Define tributary nodes of the particle 
      */
      Nodes_p = get_particle_Set(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

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
      update_DeformationGradient(F_n1_p,F_n_p,DeltaU,gradient_p);

      /*
	Compute the right Cauchy-Green tensor
      */
      C_n_p  = compute_RightCauchyGreen(F_n_p);
      C_n1_p = compute_RightCauchyGreen(F_n1_p);

      /*
	Update the second Piola-Kirchhoff stress tensor (S)
      */
      PK2_p = memory_to_Tensor(MPM_Mesh.Phi.PK2.nM[p],2);
      MatIndx_p = MPM_Mesh.MatIdx[p];
      MatProp_p = MPM_Mesh.Mat[MatIndx_p];
      update_PK2_StressTensor_Avg(PK2_p,C_n_p,C_n1_p,MatProp_p); 
      
      /* Free the gradient */
      FreeMat(gradient_p);
      free_Tensor(C_n_p);
      free_Tensor(C_n1_p);
	  
    }
  
}
