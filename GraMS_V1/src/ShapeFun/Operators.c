#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*********************************************************************/

Matrix Get_N_GP(GaussPoint MPM_Mesh,Mesh FEM_Mesh,
		int * GP_Connect, int GP_NumNodes, int i_GP)
{

  /* Gauss-Point definition */
  Matrix X_GP; /* Local coordinates of the Gauss-Point */
  X_GP.N_rows = NumberDimensions;
  X_GP.N_cols = 1;
  Matrix lp; /* Just for GIMP -> Particle voxel */
  double Beta; /* Tunning parameter for LME */
  Matrix lambda_GP = /* Just for LME -> Lagrange multipliers */
    MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  Matrix Delta_Xip; /* Just for GIMP/LME -> Distance from GP to the nodes */
  Matrix N_GP; /* Value of the shape-function in the GP */
  int GP_I;
    
  /* 3º Evaluate MPM-Q4 shape function  */
  if(strcmp(MPM_Mesh.ShapeFunctionGP,"MPMQ4") == 0){
    X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i_GP];
    N_GP = Q4(X_GP);
  }
  /* 4º Evaluate GIMP shape function */
  else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
    /* Get the distance of the GP to the nodes */
    Delta_Xip = MatAlloc(GP_NumNodes,NumberDimensions);
    for(int k = 0 ; k<GP_NumNodes ; k++){
      for(int l = 0 ; l<NumberDimensions ; l++){
	Delta_Xip.nM[k][l] =
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_Connect[k]][l];
      }
    }
    /* Get the GP voxel */
    lp.nV = MPM_Mesh.lp.nM[i_GP];
    /* Evaluate shape function */
    N_GP = GIMP_2D(Delta_Xip,lp,FEM_Mesh.DeltaX);
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
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }
    /* Asign lambda to GP */
    lambda_GP.nV = MPM_Mesh.lambda.nM[i_GP];
    /* Evaluate the shape function and it gradient */
    Beta = MPM_Mesh.Gamma/(FEM_Mesh.DeltaX*FEM_Mesh.DeltaX);
    N_GP = LME_p(Delta_Xip, lambda_GP,Beta);
    /* Free memory */
    FreeMat(Delta_Xip);
  }

  return N_GP;
}

/*********************************************************************/

Matrix Get_dN_GP(GaussPoint MPM_Mesh,Mesh FEM_Mesh,
		 int * GP_Connect, int GP_NumNodes, int i_GP)
{
  
  Matrix GP_ElemCoord; /* Coordinates of the nodes */
  int GP_I; /* Get the node for the GP */
  
  /* Gauss-Point properties */
  Matrix X_GP = /* Element coordinates of the Gauss-Point */
    MatAssign(NumberDimensions,1,NAN,NULL,NULL); 
  Matrix lp; /* Just for GIMP -> Particle voxel */
  double Beta; /* Tunning parameter for LME */
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  Matrix lambda_GP = /* Just for LME/LME -> Lagrange multipliers */
    MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  Matrix N_GP; /* Matrix with the nodal shape functions */
  Matrix dNdx_GP; /* Matrix with the nodal derivatives */
  
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
    X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i_GP];
    /* Evaluate the shape function gradient */
    dNdx_GP = Get_dNdX_Q4(X_GP,GP_ElemCoord);
    /* Free memory */
    FreeMat(GP_ElemCoord);
  }
  else if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
    /* Generate a matrix with the distances to the nodes */
    Delta_Xip = MatAlloc(GP_NumNodes,NumberDimensions);
    for(int k = 0 ; k<GP_NumNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<NumberDimensions ; l++){
	Delta_Xip.nM[k][l] =
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }
    /* Get the GP voxel */
    lp.nV = MPM_Mesh.lp.nM[i_GP];
    /* Evaluate the shape function gradient */
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
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }      
    /* Asign lambda to GP */
    lambda_GP.nV = MPM_Mesh.lambda.nM[i_GP];
    /* Evaluate the shape function and it gradient */
    Beta = MPM_Mesh.Gamma/(FEM_Mesh.DeltaX*FEM_Mesh.DeltaX);
    N_GP = LME_p(Delta_Xip, lambda_GP,Beta);
    dNdx_GP = LME_dp(Delta_Xip, N_GP);
    /* Free memory */
    FreeMat(Delta_Xip);
    FreeMat(N_GP);
  }

  return dNdx_GP;

}

/*********************************************************************/

Matrix Get_B_GP(Matrix dNdX_GP)
/*
   Get the B matrix (Usual in the classical formulation of 
   the finite element method )
   Inputs:
   - Matrix X_NC_GP : Element coordinates
   - Matrix Element : Coordinates of the element 
   (NumNodesElem x NumberDimensions)

   Outputs : Matrix B
*/
{

  /* 0º Define variables */
  /* Declaration of the output matrix (NdimVecStrain x Nnodes*Ndim) */
  Matrix B_GP;

  /* 1º Select the case to solve */
  switch(NumberDimensions){
    
  case 1:  /* 1D stress tensor */
    
    puts("Error in Get_dNdi_matrix() : 1D cases not implemented yet");
    exit(0);
    
  case 2: /* 2D stress tensor */
    
    /* 2º Allocate the output */
    B_GP = MatAlloc(3,2*dNdX_GP.N_cols);
    
    /* 4º Fill the array with the nodal partial derivation 
       of the reference element */    
    for(int i = 0 ; i<dNdX_GP.N_cols ; i++){
      B_GP.nM[0][2*i] = dNdX_GP.nM[0][i];
      B_GP.nM[1][2*i] = 0;
      B_GP.nM[2][2*i] = dNdX_GP.nM[1][i];
      
      B_GP.nM[0][2*i + 1] = 0;
      B_GP.nM[1][2*i + 1] = dNdX_GP.nM[1][i];
      B_GP.nM[2][2*i + 1] = dNdX_GP.nM[0][i];      
    }

    break;
    
  case 3: /* 3D stress tensor */
    puts("Error in Get_dNdi_matrix() : 3D cases not implemented yet");
    exit(0);
    
  default :
    puts("Error in Get_dNdi_matrix() : Wrong case select");
    exit(0);
  }
  
  return B_GP;
}

/*********************************************************************/

