#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*******************************************************/

void UpdateVelocityAndPositionGP(GaussPoint MPM_Mesh,
				 Mesh FEM_Mesh,
				 Matrix Nodal_MASS,
				 Matrix Nodal_MOMENTUM,
				 Matrix Nodal_TOT_FORCES){

  /* 0º Variable declaration */

  /* Gauss-Point definition */
  Matrix X_GP; /* Local coordinates of the Gauss-Point */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  Matrix N_GP; /* Value of the shape-function in the GP */

  /* Element for each Gauss-Point */
  int GP_NumNodes; /* Number of nodes */
  int * GP_Connect; /* Connectivity of the element */
  int GP_I;

  /* Mesh properties */
  double mass_I; /* Value of the nodal mass */
  double N_I_GP; /* Nodal value for the GP */

  /* 1º iterate over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Define element of the GP */
    GP_NumNodes = MPM_Mesh.NumberNodes[i];
    GP_Connect = ChainToArray(MPM_Mesh.ListNodes[i],GP_NumNodes);

    /* 3º Evaluate MPM-Q4 shape function  */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      N_GP = Q4(X_GP);
    }
    /* 4º Evaluate GIMP shape function */
    else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      /* Get the distance of the GP to the nodes */
      Delta_Xip = MatAlloc(GP_NumNodes,2);
      for(int k = 0 ; k<GP_NumNodes ; k++){
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  Delta_Xip.nM[k][l] =
    	    MPM_Mesh.Phi.x_GC.nM[i][l]-
    	    FEM_Mesh.Coordinates.nM[GP_Connect[k]][l];
    	}
      }
      /* Get the GP voxel */
      lp.nV = MPM_Mesh.lp.nM[i];
      /* Evaluate shape function */
      N_GP = GIMP_2D(Delta_Xip,lp,FEM_Mesh.DeltaX);
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
      /* Free memory */
      FreeMat(Delta_Xip);
    }
    
    /* 5º Iterate over the nodes of the element */
    for(int j = 0; j<GP_NumNodes; j++){
      /* Node of the GP */
      GP_I = GP_Connect[j];
      /* Evaluate the GP function in the node */
      N_I_GP = N_GP.nV[j];
      /* If this node has a null Value of the SHF continue */
      if(N_I_GP == 0){
	continue;
      }
      /* Get the nodal mass */
      mass_I = Nodal_MASS.nV[GP_I];
      /* Update GP cuantities with nodal values */
      for(int k = 0 ; k<NumberDimensions ; k++){
	/* Update the GP velocities */
	MPM_Mesh.Phi.vel.nM[i][k] +=
	  DeltaTimeStep*N_I_GP*
	  Nodal_TOT_FORCES.nM[k][GP_I]/mass_I;	
	/* Update the GP position */
	MPM_Mesh.Phi.x_GC.nM[i][k] +=
	  DeltaTimeStep*N_I_GP*
	  Nodal_MOMENTUM.nM[k][GP_I]/mass_I;
      }     
    }
    
    /* 6º Free memory */
    free(GP_Connect);
    FreeMat(N_GP);
    
  }  
}


/*******************************************************/
