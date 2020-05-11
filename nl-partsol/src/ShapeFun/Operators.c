#include "nl-partsol.h"

/*********************************************************************/

Matrix compute_ShapeFunction(Element GP_Element,GaussPoint MPM_Mesh,Mesh FEM_Mesh) 
{ 
  int i_GP = GP_Element.i_GP;
  int GP_NumNodes = GP_Element.NumberNodes;
  int * GP_Connect = GP_Element.Connectivity;
  
  Matrix GP_ElemCoord; /* Coordinates of the nodes */
  int GP_I; /* Get the node for the GP */

  int Ndim = NumberDimensions;
  
  /* Gauss-Point properties */
  Matrix X_GP = /* Element coordinates of the Gauss-Point */
    MatAssign(Ndim,1,NAN,NULL,NULL); 
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Beta_GP =  /* Tunning parameter for LME */
    MatAssign(Ndim,1,NAN,NULL,NULL);
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  Matrix lambda_GP = /* Just for LME/LME -> Lagrange multipliers */
    MatAssign(Ndim,1,NAN,NULL,NULL);
  
  Matrix ShapeFunction_p; /* Matrix with the nodal shape functions */
  
  if(strcmp(ShapeFunctionGP,"MPMQ4") == 0){

    /* Fill the poligon */
    GP_ElemCoord = MatAllocZ(GP_NumNodes,Ndim);
    for(int k = 0; k<GP_NumNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	GP_ElemCoord.nM[k][l] =
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }
    
    /* Get the element coordinates of the GP */
    X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i_GP];

    /* Evaluate the shape function */
    ShapeFunction_p = Q4_N(X_GP);
    
    /* Free memory */
    FreeMat(GP_ElemCoord);
  }
  else if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
    /* Generate a matrix with the distances to the nodes */
    Delta_Xip = MatAlloc(GP_NumNodes,Ndim);
    for(int k = 0 ; k<GP_NumNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	Delta_Xip.nM[k][l] =
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }
    
    /* Get the GP voxel */
    lp.nV = MPM_Mesh.lp.nM[i_GP];

    /* Evaluate the shape function */
    ShapeFunction_p = uGIMP_N(Delta_Xip,lp,FEM_Mesh.DeltaX);

    /* Free memory */
    FreeMat(Delta_Xip);
  }
  else if(strcmp(ShapeFunctionGP,"LME") == 0){
    /* Get the distance of the GP to the nodes */
    Delta_Xip = MatAlloc(GP_NumNodes,Ndim);
    for(int k = 0 ; k<GP_NumNodes ; k++){
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	Delta_Xip.nM[k][l] =
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }      
    /* Asign lambda to GP */
    lambda_GP.nV = MPM_Mesh.lambda.nM[i_GP];
    /* Evaluate the shape function and it gradient */
    Beta_GP.nV = MPM_Mesh.Beta.nM[i_GP];
    Beta_GP = LME_Beta(Beta_GP, Delta_Xip, gamma_LME);
    
    /* Evaluate the shape function */
    ShapeFunction_p = LME_p(Delta_Xip, lambda_GP,Beta_GP);
    
    /* Free memory */
    FreeMat(Delta_Xip);
  }
  else{
    printf("%s : %s %s %s \n",
	   "Error in Get_Operator()",
	   "The shape-function ",
	   ShapeFunctionGP,
	   "is not implemented");      
    exit(0);
  }

  return ShapeFunction_p;
}


/*********************************************************************/

Matrix compute_ShapeFunction_Gradient(Element GP_Element,GaussPoint MPM_Mesh,
				      Mesh FEM_Mesh) 
{ 
  int i_GP = GP_Element.i_GP;
  int GP_NumNodes = GP_Element.NumberNodes;
  int * GP_Connect = GP_Element.Connectivity;

  int Ndim = NumberDimensions;
  
  Matrix GP_ElemCoord; /* Coordinates of the nodes */
  int GP_I; /* Get the node for the GP */
  
  /* Gauss-Point properties */
  Matrix X_GP = /* Element coordinates of the Gauss-Point */
    MatAssign(Ndim,1,NAN,NULL,NULL); 
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Beta_GP =  /* Tunning parameter for LME */
    MatAssign(Ndim,1,NAN,NULL,NULL);
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  Matrix lambda_GP = /* Just for LME/LME -> Lagrange multipliers */
    MatAssign(Ndim,1,NAN,NULL,NULL);
  
  Matrix ShapeFunction_p; /* Matrix with the nodal shape functions */
  Matrix Gradient_p; /* Matrix with the nodal derivatives */
  
  if(strcmp(ShapeFunctionGP,"MPMQ4") == 0){
    /* Fill the poligon */
    GP_ElemCoord = MatAllocZ(GP_NumNodes,Ndim);
    
    for(int k = 0; k<GP_NumNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	GP_ElemCoord.nM[k][l] =
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }
    
    /* Get the element coordinates of the GP */
    X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i_GP];

    /* Evaluate the shape function gradient */
    Gradient_p = Q4_dN(X_GP,GP_ElemCoord);
    
    /* Free memory */
    FreeMat(GP_ElemCoord);
  }
  else if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
    /* Generate a matrix with the distances to the nodes */
    Delta_Xip = MatAlloc(GP_NumNodes,Ndim);
    for(int k = 0 ; k<GP_NumNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	Delta_Xip.nM[k][l] =
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }
    
    /* Get the GP voxel */
    lp.nV = MPM_Mesh.lp.nM[i_GP];
    
    /* Evaluate the shape function gradient */
    Gradient_p = uGIMP_dN(Delta_Xip,lp,FEM_Mesh.DeltaX);

    /* Free memory */
    FreeMat(Delta_Xip);
  }
  else if(strcmp(ShapeFunctionGP,"LME") == 0){
    /* Get the distance of the GP to the nodes */
    Delta_Xip = MatAlloc(GP_NumNodes,Ndim);
    for(int k = 0 ; k<GP_NumNodes ; k++){
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	Delta_Xip.nM[k][l] =
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }      
    /* Asign lambda to GP */
    lambda_GP.nV = MPM_Mesh.lambda.nM[i_GP];
    /* Evaluate the shape function and it gradient */
    Beta_GP.nV = MPM_Mesh.Beta.nM[i_GP];
    Beta_GP = LME_Beta(Beta_GP, Delta_Xip, gamma_LME);
    
    /* Evaluate the shape function gradient */
    ShapeFunction_p = LME_p(Delta_Xip, lambda_GP,Beta_GP);
    Gradient_p = LME_dp(Delta_Xip, ShapeFunction_p);
    FreeMat(ShapeFunction_p);
    
    /* Free memory */
    FreeMat(Delta_Xip);
  }
  else{
    printf("%s : %s %s %s \n",
	   "Error in Get_Operator()",
	   "The shape-function ",
	   ShapeFunctionGP,
	   "is not implemented");      
    exit(0);
  }

  return Gradient_p;
}
