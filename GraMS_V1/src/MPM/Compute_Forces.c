#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*******************************************************/

Matrix GetNodalForces(GaussPoint MPM_Mesh, Mesh FEM_Mesh, int TimeStep)
{
 
  /* Properties of each Gauss-Point */
  Matrix X_GP = /* Coordinate for each Gauss-Point */
    MatAssign(NumberDimensions,1,NAN,NULL,NULL); 
  Matrix lp; /* Just for GIMP -> Particle voxel */
  double Beta; /* Tunning parameter for LME */
  Matrix lambda_GP = /* Just for LME -> Lagrange multipliers */
    MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  Matrix Delta_Xip; /* Just for GIMP/LME -> Distance from GP to the nodes */

  /* Mass of the GP */
  double GP_mass; 
  /* Gauss-Point volumen */
  double Vol_GP; 
  /* Damage parameter */
  double ji_GP; 

  /* Mesh properties evaluated in Gauss-Point coords */
  Matrix N_GP; /* Matrix with the nodal shape functions */
  double N_GP_I; /* Evaluation of the GP in the node */
  Matrix dNdx_GP; /* Matrix with the nodal derivatives */
  Matrix B, B_T;

  /* Element for each Gauss-Point */
  int GP_NumNodes; /* Number of nodes */
  int * GP_Connect; /* Connectivity of the element */
  Matrix GP_ElemCoord; /* Coordinates of the nodes */
  int GP_I; /* Node of the GP */

  /* Stress tensor of a Gauss-Point and its divergence */
  Matrix Stress_GP =
    MatAssign(MPM_Mesh.Phi.Stress.N_cols,1,NAN,NULL,NULL);
  Matrix D_Stress_GP;

  /* Total forces */
  Matrix Nodal_TOT_FORCES =
    MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
  strcpy(Nodal_TOT_FORCES.Info,"Nodal_TOT_FORCES");

  /* 1º Fill matrix with the body forces for TimeStep */
  Matrix Body_Forces_t =
    Eval_Body_Forces(MPM_Mesh.B,MPM_Mesh.NumGP,TimeStep);
  
  /* 2º Fill matrix with the contact forces for TimeStep */
  Matrix Contact_Forces_t =
    Eval_Contact_Forces(MPM_Mesh.F,MPM_Mesh.NumGP,TimeStep);

  /* 3º Iterate over all the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 4º Define element of the GP */
    GP_NumNodes = MPM_Mesh.NumberNodes[i];
    GP_Connect = ChainToArray(MPM_Mesh.ListNodes[i],GP_NumNodes);

    /* 5º Evaluate the shape function in the GP */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
      /* Get the coordinates of the element */
      GP_ElemCoord = MatAllocZ(GP_NumNodes,NumberDimensions);
      for(int k = 0; k<GP_NumNodes; k++){
    	GP_I = GP_Connect[k];
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  GP_ElemCoord.nM[k][l] =
    	    FEM_Mesh.Coordinates.nM[GP_I][l];
    	}
      }
      /* Get the natural coordinates of the GP */
      X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
      /* Evaluate the shape function and it gradient */
      N_GP = Q4(X_GP);
      dNdx_GP = Get_dNdX_Q4(X_GP,GP_ElemCoord);
      /* Free memory */
      FreeMat(GP_ElemCoord);
    }
    else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      /* Get the distance of the GP to the nodes */
      Delta_Xip = MatAlloc(GP_NumNodes,NumberDimensions);
      for(int k = 0 ; k<GP_NumNodes ; k++){
    	GP_I = GP_Connect[k];
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  Delta_Xip.nM[k][l] =
    	    MPM_Mesh.Phi.x_GC.nM[i][l]-
    	    FEM_Mesh.Coordinates.nM[GP_I][l];
    	}
      }
      /* Get the voxel of the GP */
      lp.nV = MPM_Mesh.lp.nM[i];
      /* Evaluate the shape function and it gradient */
      N_GP = GIMP_2D(Delta_Xip,lp,FEM_Mesh.DeltaX);
      dNdx_GP = dGIMP_2D(Delta_Xip,lp,FEM_Mesh.DeltaX);
      /* Free memory */
      FreeMat(Delta_Xip);
    }
    else if(strcmp(MPM_Mesh.ShapeFunctionGP,"LME") == 0){
      /* Get the distance of the GP to the nodes */
      Delta_Xip = MatAlloc(GP_NumNodes,NumberDimensions);
      for(int k = 0 ; k<GP_NumNodes ; k++){
    	GP_I = GP_Connect[k];
    	for(int l = 0 ; l<NumberDimensions ; l++){
    	  Delta_Xip.nM[k][l] =
    	    MPM_Mesh.Phi.x_GC.nM[i][l]-
    	    FEM_Mesh.Coordinates.nM[GP_I][l];
    	}
      }
      /* Asign lambda to GP */
      lambda_GP.nV = MPM_Mesh.lambda.nM[i];
      /* Evaluate the shape function and it gradient */
      Beta = MPM_Mesh.Gamma/(FEM_Mesh.DeltaX*FEM_Mesh.DeltaX);
      N_GP = LME_p(Delta_Xip, lambda_GP,Beta);
      dNdx_GP = LME_dp(Delta_Xip, N_GP);
      /* Free memory */
      FreeMat(Delta_Xip);
    }
       
    /* 6º Get the B_T matrix for the derivates */
    B = Get_B_GP(dNdx_GP), FreeMat(dNdx_GP);
    B_T = Transpose_Mat(B), FreeMat(B);
    
    /* 7º Asign to an auxiliar variable the value of the stress tensor */
    Stress_GP.nV = MPM_Mesh.Phi.Stress.nM[i];

    /* 8º Get the divergence stress tensor evaluates in the Gauss-Point 
     and free the B_T matrix */
    D_Stress_GP = Scalar_prod(B_T,Stress_GP), FreeMat(B_T);
    
    /* 9º Calcule the volumen of the Gauss-Point */
    GP_mass = MPM_Mesh.Phi.mass.nV[i];
    Vol_GP = GP_mass/MPM_Mesh.Phi.rho.nV[i];

    /* 10º Damage parameter for the Gauss-point (fracture) */
    ji_GP = MPM_Mesh.Phi.ji.nV[i];

    /* 11º Acumulate this forces to the total array with the internal forces */  
    for(int k = 0; k<GP_NumNodes; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      /* Evaluate the GP function in the node */
      N_GP_I = N_GP.nV[k];
      /* If this node has a null Value of the SHF continue */
      if(N_GP_I == 0) continue;
      /* Loop in the dimensions */
      for(int l = 0; l<NumberDimensions; l++){
	/* 10aº Add the internal forces with 
	 damage variable option */
	Nodal_TOT_FORCES.nM[l][GP_I] -= (1-ji_GP)* 
	  D_Stress_GP.nV[k*NumberDimensions+l]*Vol_GP;
	/* 10bº Add the body forces */
	Nodal_TOT_FORCES.nM[l][GP_I] +=
	  N_GP_I*Body_Forces_t.nM[l][i]*GP_mass;	
	/* 10cº Add the contact forces */
	Nodal_TOT_FORCES.nM[l][GP_I] +=
	  N_GP_I*Contact_Forces_t.nM[l][i]*Vol_GP;	
      }      
    }
    
    /* 12 º Free memory */
    free(GP_Connect), FreeMat(D_Stress_GP), FreeMat(N_GP);

  }

  /* 13º Free memory */
  FreeMat(Contact_Forces_t), FreeMat(Body_Forces_t);
  
  return Nodal_TOT_FORCES;
  
}

/*******************************************************/
