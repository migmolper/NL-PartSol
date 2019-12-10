#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*********************************************************************/

Matrix GetNodalMassMomentum(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{

  /* 0º Variable declaration */

  /* Output */
  Matrix Nodal_FIELDS;
  
  /* Gauss-Point definition */
  Matrix X_GP; /* Local coordinates */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  Matrix N_GP; /* Value of the shape-function */
  double N_GP_I; /* Evaluation of the GP in the node */
  double GP_mass; /* Mass of the GP */

  /* Element for each Gauss-Point */
  int GP_NumNodes; /* Number of nodes */
  int * GP_Connect; /* Connectivity of the element */
  int GP_I;

  /* 1º Allocate the output list of fields */
  Nodal_FIELDS = MatAllocZ(NumberDimensions + 1,FEM_Mesh.NumNodesMesh);
  
  /* 2º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 3º Define element of the GP */
    GP_NumNodes = MPM_Mesh.NumberNodes[i];
    GP_Connect = ChainToArray(MPM_Mesh.ListNodes[i],GP_NumNodes);

    /* 4º Evaluate the shape function in the coordinates of the GP */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      N_GP = Q4(X_GP);
    }
    else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      /* Fill the poligon */
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
      /* Get the lenght voxel */
      lp.nV = MPM_Mesh.Phi.lp.nM[i];
      /* Evaluate the shape function */
      N_GP = GIMP_2D(Delta_Xip,lp,FEM_Mesh.DeltaX);
      /* Free memory */
      FreeMat(Delta_Xip);
    }
   
    /* 5º Get the mass of the GP */
    GP_mass = MPM_Mesh.Phi.mass.nV[i];

    /* 6º Get the nodal mass and mommentum */
    for(int k = 0 ; k<GP_NumNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      /* Evaluate the GP function in the node */
      N_GP_I = N_GP.nV[k];
      /* If this node has a null Value of the SHF continue */
      if(N_GP_I == 0){
	continue;
      }
      /* Nodal mass */
      Nodal_FIELDS.nM[0][GP_I] += GP_mass*N_GP_I;
      /* Nodal momentum */
      for(int l = 0 ; l<NumberDimensions ; l++){
	Nodal_FIELDS.nM[l+1][GP_I] +=
	  GP_mass*MPM_Mesh.Phi.vel.nM[i][l]*N_GP_I;
      }   
    }

    /* 7º Free the value of the shape functions */
    FreeMat(N_GP);
    free(GP_Connect);
  }
 
  return Nodal_FIELDS;
  
}

/*******************************************************/

Matrix GetNodalVelocity(Mesh FEM_Mesh,
			Matrix Nodal_MOMENTUM,
			Matrix Nodal_MASS){
  /*
    Get the nodal velocity using : 
    v_{i,I}^{k-1/2} = \frac{p_{i,I}^{k-1/2}}{m_I^{k}}
    Initialize nodal velocities 
  */
  Matrix Vel_Mesh;
  Vel_Mesh = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
  strcpy(Vel_Mesh.Info,"VELOCITY");
  
  /* 1º Get nodal values of the velocity */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<NumberDimensions ; j++){
      if(Nodal_MASS.nV[i] > 0){
	Vel_Mesh.nM[j][i] = (double)Nodal_MOMENTUM.nM[j][i]/Nodal_MASS.nV[i];
      }
    }    
  }
  
  return Vel_Mesh;
}

/*******************************************************/
