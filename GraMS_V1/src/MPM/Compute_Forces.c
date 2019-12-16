#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*******************************************************/

Matrix GetNodalForces(GaussPoint MPM_Mesh, Mesh FEM_Mesh, int TimeStep)
{
  
  /* 0º Auxiliar variable declaration */

  /* Properties of each Gauss-Point */
  Matrix X_GP; /* Coordinate for each Gauss-Point */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  double GP_mass; /* Mass of the GP */
  double Vol_GP; /* Gauss-Point volumen */
  double ji_GP; /* Damage parameter */

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

  /* Stress tensor of a Gauss-Point */
  Matrix StressTensor_GP; 
  StressTensor_GP.N_rows = MPM_Mesh.Phi.Stress.N_cols;
  StressTensor_GP.N_cols = 1;
  StressTensor_GP.nM = NULL;

  /* Internal forces for each node in a element (by a GP)*/  
  Matrix Div_Stress_Tensor;

  /* Total forces */
  Matrix Nodal_TOT_FORCES; 
  Nodal_TOT_FORCES = MatAllocZ(NumberDimensions,FEM_Mesh.NumNodesMesh);
  strcpy(Nodal_TOT_FORCES.Info,"Nodal_TOT_FORCES");

  /* Contact forces */
  Matrix Contact_Forces_t;
  Contact_Forces_t = MatAllocZ(NumberDimensions,MPM_Mesh.NumGP);

  /* Body Forces */
  Matrix Body_Forces_t;
  Body_Forces_t = MatAllocZ(NumberDimensions,MPM_Mesh.NumGP);
  int GP_Force;
  
  /* 1º Fill matrix with the body forces for TimeStep */ 
  if(MPM_Mesh.B.NumLoads>0){
    for(int i = 0 ; i<MPM_Mesh.B.NumLoads; i++){
      for(int j = 0 ; j<MPM_Mesh.B.Load_i[i].NumNodes ; j++){
	GP_Force = MPM_Mesh.B.Load_i[i].Nodes[j];
	for(int k = 0 ; k<MPM_Mesh.B.Load_i[i].Dim ; k++){
	  if( (MPM_Mesh.B.Load_i[i].Dir[k] == 1) ||
	      (MPM_Mesh.B.Load_i[i].Dir[k] == -1)){
	    if( (TimeStep < 0) ||
		(TimeStep > MPM_Mesh.B.Load_i[i].Value[k].Num)){
	      printf("%s : %s\n","Error in GetNodalForces()",
		      "The time step is out of the curve !!");
	      exit(0);
	    }
	    Body_Forces_t.nM[k][GP_Force] +=
	      MPM_Mesh.B.Load_i[i].Value[k].Fx[TimeStep]*
	      (double)MPM_Mesh.B.Load_i[i].Dir[k];
	  }
	}
      }
    }
  }

  
  /* 2º Fill matrix with the contact forces for TimeStep */
  if(MPM_Mesh.F.NumLoads>0){
    for(int i = 0 ; i<MPM_Mesh.F.NumLoads; i++){
      for(int j = 0 ; j<MPM_Mesh.F.Load_i[i].NumNodes ; j++){
	GP_Force = MPM_Mesh.F.Load_i[i].Nodes[j];
	for(int k = 0 ; k<MPM_Mesh.F.Load_i[i].Dim ; k++){
	  if( (MPM_Mesh.F.Load_i[i].Dir[k] == 1) ||
	      (MPM_Mesh.F.Load_i[i].Dir[k] == -1)){
	    if( (TimeStep < 0) ||
		(TimeStep > MPM_Mesh.F.Load_i[i].Value[k].Num)){
	      printf("%s : %s \n","Error in GetNodalForces()",
		     "The time step is out of the curve !!");
	      exit(0);
	    }
	    Contact_Forces_t.nM[k][GP_Force] +=
	      MPM_Mesh.F.Load_i[i].Value[k].Fx[TimeStep]*
	      (double)MPM_Mesh.F.Load_i[i].Dir[k];
	  }
	}
      }
    }
  }

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
      Delta_Xip = MatAlloc(GP_NumNodes,2);
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
    }
       
    /* 6º Get the B_T matrix for the derivates */
    B = Get_B_GP(dNdx_GP);    
    FreeMat(dNdx_GP);
    B_T = Transpose_Mat(B);
    FreeMat(B);
    
    /* 7º Asign to an auxiliar variable the value of the stress tensor */
    StressTensor_GP.nV = MPM_Mesh.Phi.Stress.nM[i];

    /* 8º Get the divergence stress tensor evaluates in the Gauss-Point 
     and free the B_T matrix */
    Div_Stress_Tensor = Scalar_prod(B_T,StressTensor_GP);    
    FreeMat(B_T);
    
    /* 9º Calcule the volumen of the Gauss-Point */
    GP_mass = MPM_Mesh.Phi.mass.nV[i];
    Vol_GP = GP_mass/MPM_Mesh.Phi.rho.nV[i];

    /* Damage parameter for the Gauss-point (fracture) */
    ji_GP = MPM_Mesh.Phi.ji.nV[i];

    /* 10º Acumulate this forces to the total array with the internal forces */  
    for(int k = 0; k<GP_NumNodes; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      /* Evaluate the GP function in the node */
      N_GP_I = N_GP.nV[k];
      /* If this node has a null Value of the SHF continue */
      if(N_GP_I == 0){
	continue;
      }
      for(int l = 0; l<NumberDimensions; l++){
	/* 10aº Add the internal forces with 
	 damage variable option */
	Nodal_TOT_FORCES.nM[l][GP_I] -=
	  (1 - ji_GP)* 
	  Div_Stress_Tensor.nV[k*NumberDimensions + l]*
	  Vol_GP;
	/* 10bº Add the body forces */
	Nodal_TOT_FORCES.nM[l][GP_I] +=
	  N_GP_I*Body_Forces_t.nM[l][i]*
	  GP_mass;	
	/* 10cº Add the contact forces */
	Nodal_TOT_FORCES.nM[l][GP_I] +=
	  N_GP_I*Contact_Forces_t.nM[l][i]*
	  Vol_GP;	
      }      
    }
	
    /* 10dº Free memory */
    free(GP_Connect);
    FreeMat(Div_Stress_Tensor);
    FreeMat(N_GP);

  }

  /* 11º Free memory */
  FreeMat(Contact_Forces_t);
  FreeMat(Body_Forces_t);  
  
  return Nodal_TOT_FORCES;
  
}

/*******************************************************/
