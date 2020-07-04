#include "nl-partsol.h"

/*
  Auxiliar functions 
*/
static Matrix compute_Effective_MassMatrix(GaussPoint, Mesh, Mask, double);
static Mask generate_NodalMask(Mesh FEM_Mesh);

/**************************************************************/

void U_Discrete_Energy_Momentum(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)
{

  /*!
    Integer variables 
  */
  int Ndim = NumberDimensions;

  double epsilon = 1;

  /*!
    Auxiliar variable for the mass and momentum 
  */
  Matrix Effective_MassMatrix;
  strcpy(Effective_MassMatrix.Info,"Effective-Mass-Matrix");

  Mask ActiveNodes = generate_NodalMask(FEM_Mesh);
  
  Effective_MassMatrix = compute_Effective_MassMatrix(MPM_Mesh,FEM_Mesh,
						      ActiveNodes, epsilon);

  PrintMatrix(Effective_MassMatrix, ActiveNodes.Nactivenodes, ActiveNodes.Nactivenodes);
  
  /*!
    Free memory.
   */
  FreeMat(Effective_MassMatrix);
  
}

/**************************************************************/

static Matrix compute_Effective_MassMatrix(GaussPoint MPM_Mesh, Mesh FEM_Mesh,
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
  

  return Effective_MassMatrix; 
}

/**************************************************************/

static Mask generate_NodalMask(Mesh FEM_Mesh)
{
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Nactivenodes = 0;
  int * Nodes2Mask = (int *)Allocate_ArrayZ(Nnodes,sizeof(int));;
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
	  Nodes2Mask[I] = -1;
	}
    }

  M.Nactivenodes = Nactivenodes;
  M.Mask2Nodes = Set_to_Pointer(Mask2Nodes,Nactivenodes);
  M.Nodes2Mask = Nodes2Mask;  
 
  return M;
}

/**************************************************************/

