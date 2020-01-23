#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*******************************************************/

Matrix GetNodalForces(GaussPoint MPM_Mesh, Mesh FEM_Mesh, int TimeStep)
{
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
  Matrix N_dNdx_GP; /* Operator Matrix */
  Matrix B, B_T;

  /* Element for each Gauss-Point */
  int GP_NumNodes; /* Number of nodes */
  int * GP_Connect; /* Connectivity of the element */
  /* Matrix GP_ElemCoord; /\* Coordinates of the nodes *\/ */
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
    Eval_Body_Forces(MPM_Mesh.B,MPM_Mesh.NumberBodyForces,
		     MPM_Mesh.NumGP,TimeStep);
  
  /* 2º Fill matrix with the contact forces for TimeStep */
  Matrix Contact_Forces_t =
    Eval_Contact_Forces(MPM_Mesh.F,MPM_Mesh.NumNeumannBC,
			MPM_Mesh.NumGP,TimeStep);

  /* 3º Iterate over all the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 4º Define element of the GP */
    GP_NumNodes = MPM_Mesh.NumberNodes[i];
    GP_Connect = ChainToArray(MPM_Mesh.ListNodes[i],GP_NumNodes);

    /* 5º Evaluate the shape function and its gradient in the GP */
    N_dNdx_GP = Get_Operator("N_dNdx",i,GP_Connect,GP_NumNodes,
			     MPM_Mesh,FEM_Mesh);
    /* Asign values to the pointer structures */
    N_GP = MatAssign(1,N_dNdx_GP.N_cols,NAN,N_dNdx_GP.nM[0],NULL);
    dNdx_GP = MatAssign(2,N_dNdx_GP.N_cols,NAN,NULL,
			(double **)malloc(2*sizeof(double *)));
    dNdx_GP.nM[0] = N_dNdx_GP.nM[1];
    dNdx_GP.nM[1] = N_dNdx_GP.nM[2];
    /* Free the original table container */
    free(N_dNdx_GP.nM);
           
    /* 6º Get the B_T matrix for the derivates */
    B = Get_B_GP(dNdx_GP);
    FreeMat(dNdx_GP);
    B_T = Transpose_Mat(B);
    FreeMat(B);
    
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
      if(fabs(N_GP_I) <= TOL_zero) continue;
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
