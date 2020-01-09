#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../GRAMS/grams.h"


/***********************************************/
/******* 2D cuadrilateral linear element *******/
/***********************************************/

/* (3)     (2)  */
/*  o-------o   */
/*  |       |   */
/*  |       |   */
/*  o-------o   */
/* (0)     (1)  */

void Q4_Initialize(GaussPoint MPM_Mesh, Mesh FEM_Mesh)

{
  /* Variables for the GP coordinates */  
  Matrix X_GC_GP = MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  Matrix X_EC_GP = MatAssign(NumberDimensions,1,NAN,NULL,NULL);

  /* Variables for the poligon */
  int NumVertex;
  int * Poligon_Connectivity;
  Matrix Poligon_Coordinates;
  ChainPtr ListNodes_I;

  /* 1º Set to zero the active/non-active node, and the GPs in each 
     element */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.ActiveNode[i] = 0;
  }
  
  for(int i = 0 ; i<FEM_Mesh.NumElemMesh ; i++){
    FreeChain(FEM_Mesh.GPsElements[i]);
    FEM_Mesh.GPsElements[i] = NULL;
  }

  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Assign the value to this auxiliar pointer */ 
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){

      /* 3º Connectivity of the Poligon */
      NumVertex = FEM_Mesh.NumNodesElem[j];
      Poligon_Connectivity =
	ChainToArray(FEM_Mesh.Connectivity[j],NumVertex);

      /* 4º Get the coordinates of the element */
      Poligon_Coordinates =
	ElemCoordinates(FEM_Mesh,Poligon_Connectivity,NumVertex);
      
      /* 5º Check out if the GP is in the Element */
      if(InOut_Poligon(X_GC_GP,Poligon_Coordinates) == 1){

	/* 6º Asign to the GP a element in the background mesh, just for 
	   searching porpuses */
	MPM_Mesh.Element_id[i] = j;
	PushNodeTop(&FEM_Mesh.GPsElements[j],i);

	/* 7º If the GP is in the element, get its natural coordinates */
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	Get_X_EC_Q4(X_EC_GP,X_GC_GP,Poligon_Coordinates);

	/* 8º Get list of nodes near to the GP */
	FreeChain(MPM_Mesh.ListNodes[i]);
	MPM_Mesh.ListNodes[i] = NULL;
	MPM_Mesh.ListNodes[i] = CopyChain(FEM_Mesh.Connectivity[j]);
	
	/* 9º Active those nodes that interact with the GP */
	ListNodes_I = MPM_Mesh.ListNodes[i];
	while(ListNodes_I != NULL){
	  FEM_Mesh.ActiveNode[ListNodes_I->I] += 1;
	  ListNodes_I = ListNodes_I->next; 
	}

	/* 10º Free memory and go for the next GP */
	free(Poligon_Connectivity);
	FreeMat(Poligon_Coordinates);
	break;
	
      }
      
      /* 11º Free memory */
      free(Poligon_Connectivity);
      FreeMat(Poligon_Coordinates);
      
    } 

  }
  
}

/* Shape functions */
Matrix Q4(Matrix X_e){

  /* Error check */
  if( (fabs(X_e.nV[0]) > 1 ) || (fabs(X_e.nV[1]) > 1 ) ){
    printf("Error in Q4() : Out of the element bounds !!! \n");
    exit(0);
  }    
  
  /* Definition and allocation */
  Matrix N_ref =  MatAlloc(1,4);

  /* Fill the array */
  N_ref.nV[0] = 0.25*(1-X_e.nV[0])*(1-X_e.nV[1]);
  N_ref.nV[1] = 0.25*(1+X_e.nV[0])*(1-X_e.nV[1]);
  N_ref.nV[2] = 0.25*(1+X_e.nV[0])*(1+X_e.nV[1]);
  N_ref.nV[3] = 0.25*(1-X_e.nV[0])*(1+X_e.nV[1]);
  
  return N_ref;
}

/*********************************************************************/

/* Derivatives of the shape functions */
Matrix dQ4(Matrix X_e){

  /* Error check */
  if( (fabs(X_e.nV[0]) > 1 ) || (fabs(X_e.nV[1]) > 1 ) ){
    printf("Error in dQ4() : Out of the element bounds !!! \n");
    exit(0);
  }    
  
  /* Definition and allocation */
  Matrix dNdX_ref = MatAlloc(2,4);
  
  /* Fill the matrix */

  /* Node 1 */
  /* \frac{\partial N0}{\partial \xi} */
  dNdX_ref.nM[0][0] = -0.25*(1-X_e.nV[1]);
  /* \frac{\partial N0}{\partial \eta} */
  dNdX_ref.nM[1][0] = -0.25*(1-X_e.nV[0]); 

  /* Node 2 */
  /* \frac{\partial N1}{\partial \xi} */
  dNdX_ref.nM[0][1] = +0.25*(1-X_e.nV[1]);
  /* \frac{\partial N1}{\partial \eta} */
  dNdX_ref.nM[1][1] = -0.25*(1+X_e.nV[0]);
  
  /* Node 3 */
  /* \frac{\partial N2}{\partial \xi} */
  dNdX_ref.nM[0][2] = +0.25*(1+X_e.nV[1]);
  /* \frac{\partial N2}{\partial \eta} */
  dNdX_ref.nM[1][2] = +0.25*(1+X_e.nV[0]);
  
  /* Node 4 */
  /* \frac{\partial N3}{\partial \xi} */
  dNdX_ref.nM[0][3] = -0.25*(1+X_e.nV[1]);
  /* \frac{\partial N3}{\partial \eta} */
  dNdX_ref.nM[1][3] = +0.25*(1-X_e.nV[0]); 
  
  return dNdX_ref;
}

/*********************************************************************/

/* Jacobian of the transformation for the four-nodes quadrilateral */
Matrix Get_F_Ref_Q4(Matrix X_NC_GP,Matrix X_GC_Nodes)
/*
  Get the jacobian of the transformation of the reference element :

  J = grad(N_{\alpha}) \cdot x_{\alpha}
  
  Inputs :
  - X_g -> This are the coordinates of the nodes (4x2)
  - dNdX_Ref_GP -> Derivative gradient evaluated in the GP (2x4)

  Output :
  - Jacobian -> Jacobian of the reference element
  evaluated in the GP
*/
{
  /* Variable declaration */
  Matrix dNdX_Ref_GP;
  Matrix X_alpha = MatAlloc(2,1);
  Matrix dNdx_alpha = MatAlloc(1,2);
  Matrix F_Ref_alpha;
  Matrix F_Ref = MatAllocZ(2,2);

  /* 1º Evaluate the derivarive of the shape function in the GP */
  dNdX_Ref_GP = dQ4(X_NC_GP); 

  /* 2º Get the F_Ref doing a loop over the nodes of the element */
  for(int i = 0 ; i<4 ; i++){

    /* 3º Fill arrays for the tensorial product */
    for(int j = 0 ; j<2 ; j++){
      X_alpha.nV[j] = X_GC_Nodes.nM[i][j];
      dNdx_alpha.nV[j] = dNdX_Ref_GP.nM[j][i];
    }

    /* 4º Get the nodal contribution */
    F_Ref_alpha = Tensorial_prod(X_alpha,dNdx_alpha);

    /* 5º Increment the reference deformation gradient */
    F_Ref = Incr_Mat(F_Ref, F_Ref_alpha);

    /* 6º Free data of the nodal contribution */
    FreeMat(F_Ref_alpha);
    
  }
  
  /* 7º Free memory */
  FreeMat(dNdX_Ref_GP);
  FreeMat(X_alpha);
  FreeMat(dNdx_alpha);

  /* 8º Output */
  return F_Ref;
}

/*********************************************************************/

/* Element gradient in the real element */
Matrix Get_dNdX_Q4(Matrix X_EC_GP,Matrix Element)
/*
  - Matrix X_EC_GP : Element coordinates of the gauss point
  - Matrix Element : Coordinates of the element (NumNodesElem x NumberDimensions)
*/
{
  
  /* 0º Definition and allocation */
  Matrix dNdX_Ref_GP; /* Derivative of the shape function evaluated in the GP (Ndim x Nnodes) */
  Matrix F_GP; /* Jacobian of the transformation evaluated in the GP (Ndim x Ndim) */
  Matrix F_GP_m1; /* Inverse of the Jacobian */
  Matrix F_GP_Tm1; /* Transpose of the inverse Jacobian */
  Matrix dNdx_GP; /* Derivatives of the shape function evaluates in the GP (Ndim x Ndim) */

  /* 1º Evaluate the gradient of the shape function in the GP */
  dNdX_Ref_GP = dQ4(X_EC_GP);
  
  /* 2º Get the Jacobian of the transformation evaluated in the GP */
  F_GP = Get_F_Ref_Q4(X_EC_GP,Element);
    
  /* 3º Get the inverse of the deformation gradient */
  F_GP_m1 = Get_Inverse(F_GP);
  FreeMat(F_GP);
  
  /* 4º Get the transpose of the inverse of the Jacobian */
  F_GP_Tm1 = Transpose_Mat(F_GP_m1);
  FreeMat(F_GP_m1);
  
  /* 5º Get the gradient of the shape functions in global coordinates */
  dNdx_GP = Scalar_prod(F_GP_Tm1,dNdX_Ref_GP);
  FreeMat(F_GP_Tm1);
  FreeMat(dNdX_Ref_GP);

  /* 6º Return result */
  return dNdx_GP;
}

/*********************************************************************/

/* Global coordinates of the four nodes quadrilateral */
Matrix Get_X_GC_Q4(Matrix X_NC_GP,Matrix X_GC_Nodes)
/*
  This function evaluate the position of the GP in the element,
  and get it global coordiantes    
*/
{
  /* 0º Variable declaration */
  Matrix N_ref;
  Matrix X_GC_GP;

  /* 1º Evaluate the Q4 element in the element coordinates */
  N_ref = Q4(X_NC_GP);

  /* 2º Allocate the output coordinates */
  X_GC_GP = MatAlloc(2,1);

  /* 3º Get the global coordinates for this element coordiantes in this element */
  X_GC_GP.nV[0] =
    N_ref.nV[0]*X_GC_Nodes.nM[0][0] +
    N_ref.nV[1]*X_GC_Nodes.nM[1][0] +
    N_ref.nV[2]*X_GC_Nodes.nM[2][0] +
    N_ref.nV[3]*X_GC_Nodes.nM[3][0];
  
  X_GC_GP.nV[1] =
    N_ref.nV[0]*X_GC_Nodes.nM[0][1] +
    N_ref.nV[1]*X_GC_Nodes.nM[1][1] +
    N_ref.nV[2]*X_GC_Nodes.nM[2][1] +
    N_ref.nV[3]*X_GC_Nodes.nM[3][1];

  /* 4º Free memory */
  FreeMat(N_ref);

  /* 5º Output */
  return X_GC_GP;
 
}

/*********************************************************************/

void Get_X_EC_Q4(Matrix X_EC_GP,
		 Matrix X_GC_GP,
		 Matrix Element_GC_Nod)
/* 
   The function return the natural coordinates of a point 
   inside of the element.
 
   Inputs :
   - Coordinates of the element nodes
   - Initial coordinate of the point to start the search 
   - Derivative function of the element

   Depending of the kind of element, we employ differents types
   of shape functions
*/
{
  
  X_EC_GP = Newton_Rapson(Get_X_GC_Q4,Element_GC_Nod,
			  Get_F_Ref_Q4,Element_GC_Nod,
			  X_GC_GP,X_EC_GP);
}

/*********************************************************************/
