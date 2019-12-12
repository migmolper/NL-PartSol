#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*******************************************************/

void UpdateGaussPointStrain(GaussPoint MPM_Mesh,
			    Mesh FEM_Mesh,
			    Matrix Mesh_Vel)
/*
  Calcule the particle stress increment :

  \Delta\Epsilon_{ij,p}^{k-1/2} = 
  \frac{\Delta t}{2} \cdot
  \sum_{I=0}^{Nn}(N^{k}_{Ip,j} \cdot
  v_{iI}^{k-1/2} + N^{k}_{Ip,i} \cdot 
  v_{jI}^{k-1/2})
*/
{
  /* Gauss-Point properties */
  Matrix X_GP; /* Element coordinates of the Gauss-Point */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  
  /* Element for each Gauss-Point */
  int GP_NumNodes; /* Number of nodes */
  int * GP_Connect; /* Connectivity of the element */
  Matrix GP_ElemCoord; /* Coordinates of the nodes */
  int GP_I; /* Get the node for the GP */
  
  /* Mesh variables */
  Matrix Elem_Vel; /* Array with the nodal velocities */
  Matrix N_GP; /* Matrix with the nodal shape functions */
  Matrix dNdx_GP; /* Matrix with the nodal derivatives */
  Matrix B; /* B marix to get the deformation */
  Matrix Increment_Strain_GP; /* Vectoriced Strain tensor */
  double Incr_TraceStrain; /* Increment of the trace of the Stress tensor */
 
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
   
    /* 1º Define element of the GP */
    GP_NumNodes = MPM_Mesh.NumberNodes[i];
    GP_Connect = ChainToArray(MPM_Mesh.ListNodes[i],GP_NumNodes);
    
    /* 2º Get the element gradient */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
      /* Fill the poligon */
      GP_ElemCoord = MatAllocZ(GP_NumNodes,NumberDimensions);
      for(int k = 0; k<GP_NumNodes ; k++){
    	/* Get the node for the GP */
    	GP_I = GP_Connect[k];
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  GP_ElemCoord.nM[k][l] =
    	    FEM_Mesh.Coordinates.nM[GP_I][l];
    	}
      }
      /* Get the element coordinates of the GP */
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      /* Evaluate the shape function gradient */
      dNdx_GP = Get_dNdX_Q4(X_GP,GP_ElemCoord);
      /* Free memory */
      FreeMat(GP_ElemCoord);
    }
    else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      /* Generate a matrix with the distances to the nodes */
      Delta_Xip = MatAlloc(GP_NumNodes,2);
      for(int k = 0 ; k<GP_NumNodes ; k++){
    	/* Get the node for the GP */
    	GP_I = GP_Connect[k];
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  Delta_Xip.nM[k][l] =
    	    MPM_Mesh.Phi.x_GC.nM[i][l]-
    	    FEM_Mesh.Coordinates.nM[GP_I][l];
    	}
      }
      /* Get the GP voxel */
      lp.nV = MPM_Mesh.lp.nM[i];
      /* Evaluate the shape function gradient */
      dNdx_GP = dGIMP_2D(Delta_Xip,lp,FEM_Mesh.DeltaX);
      /* Free memory */
      FreeMat(Delta_Xip);
    }
    else if(strcmp(MPM_Mesh.ShapeFunctionGP,"LME") == 0){
      /* Get the distance of the GP to the nodes */
      Delta_Xip = MatAlloc(GP_NumNodes,2);
      for(int k = 0 ; k<GP_NumNodes ; k++){
    	GP_I = GP_Connect[k];
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  Delta_Xip.nM[k][l] =
    	    MPM_Mesh.Phi.x_GC.nM[i][l]-
    	    FEM_Mesh.Coordinates.nM[GP_I][l];
    	}
      }
      /* Evaluate the shape function and it gradient */
      N_GP = LME_pa(Delta_Xip, MPM_Mesh.lambda,
		    FEM_Mesh.DeltaX, MPM_Mesh.Gamma);
      dNdx_GP = LME_dpa(Delta_Xip, N_GP);
      /* Free memory */
      FreeMat(Delta_Xip);
      FreeMat(N_GP);
    }
	    
    /* Calcule the B matrix */
    B = Get_B_GP(dNdx_GP);
    /* Free shape-function derivatives */
    FreeMat(dNdx_GP);

    /* 3º Allocate and array with the velocities of the element */
    Elem_Vel = MatAllocZ(GP_NumNodes*NumberDimensions,1);

    /* 4º Get the nodal velocities in the element */
    for(int k = 0; k<GP_NumNodes; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<NumberDimensions ; l++){
	Elem_Vel.nV[k*NumberDimensions + l] =
	  Mesh_Vel.nM[l][GP_I];
      }
    }
    /* Free data */
    free(GP_Connect);
   
    /* 5º Multiply B by the velocity array and by the time step to get
       the increment stress tensor */
    Increment_Strain_GP = Scalar_prod(B,Elem_Vel);    
    for(int j = 0 ; j<MPM_Mesh.Phi.Strain.N_cols ; j++){
      Increment_Strain_GP.nV[j] *= DeltaTimeStep;
    }

    /* Free memory */
    FreeMat(Elem_Vel);
    FreeMat(B);

    /* 6º Set to zero the trace of the stress tensor */
    Incr_TraceStrain = 0;

    /* 7º Udate the Gauss-Point strain tensor */
    for(int j = 0 ; j<MPM_Mesh.Phi.Strain.N_cols ; j++){
      MPM_Mesh.Phi.Strain.nM[i][j] += Increment_Strain_GP.nV[j];
      /* 7º Get the trace of the stress tensor */
      if(j<NumberDimensions){
	Incr_TraceStrain += Increment_Strain_GP.nV[j];
      }
    }
    FreeMat(Increment_Strain_GP);

    /* 8º Update the density of the GP */
    MPM_Mesh.Mat.rho.nV[i] =
      UpdateGaussPointDensity(MPM_Mesh.Mat.rho.nV[i],
			      Incr_TraceStrain);    
  }
}

/*******************************************************/
