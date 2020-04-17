#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grams.h"


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
  int Ndim = NumberDimensions;
  Matrix X_GC_GP = MatAssign(Ndim,1,NAN,NULL,NULL);
  Matrix X_EC_GP = MatAssign(Ndim,1,NAN,NULL,NULL);

  /* Variables for the poligon */
  int NumVertex;
  int * Poligon_Connectivity;
  Matrix Poligon_Coordinates;
  ChainPtr ListNodes_I;

  /* 1º Set to zero the active/non-active node, and the GPs in each 
     element */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.NumParticles[i] = 0;
  }
  
  for(int i = 0 ; i<FEM_Mesh.NumElemMesh ; i++){
    free_Set(FEM_Mesh.I_particles[i]);
    FEM_Mesh.I_particles[i] = NULL;
  }

  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Assign the value to this auxiliar pointer */ 
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){

      /* 3º Connectivity of the Poligon */
      NumVertex = FEM_Mesh.NumNodesElem[j];
      Poligon_Connectivity =
	Set_to_Pointer(FEM_Mesh.Connectivity[j],NumVertex);

      /* 4º Get the coordinates of the element */
      Poligon_Coordinates = ElemCoordinates(FEM_Mesh,Poligon_Connectivity,
					    NumVertex);
      
      /* 5º Check out if the GP is in the Element */
      if(InOut_Poligon(X_GC_GP,Poligon_Coordinates) == 1){

	/* 6º Asign to the GP a element in the background mesh, just for 
	   searching porpuses */
	MPM_Mesh.I0[i] = j;
	push_to_Set(&FEM_Mesh.I_particles[j],i);

	/* 7º If the GP is in the element, get its natural coordinates */	
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	Q4_X_to_Xi(X_EC_GP,X_GC_GP,Poligon_Coordinates);

	/* 8º Get list of nodes near to the GP */
	free_Set(MPM_Mesh.ListNodes[i]);
	MPM_Mesh.ListNodes[i] = NULL;
	MPM_Mesh.ListNodes[i] = CopyChain(FEM_Mesh.Connectivity[j]);
	/* 9º Active those nodes that interact with the GP */
	ListNodes_I = MPM_Mesh.ListNodes[i];
	while(ListNodes_I != NULL){
	  FEM_Mesh.NumParticles[ListNodes_I->I] += 1;
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
Matrix Q4_N(Matrix X_e){
  
  /* Definition and allocation */
  Matrix N_ref =  MatAlloc(1,4);

  /* Fill the array */
  N_ref.nV[0] = 0.25*DMIN(4,DMAX(0,(1-X_e.nV[0])*(1-X_e.nV[1])));
  N_ref.nV[1] = 0.25*DMIN(4,DMAX(0,(1+X_e.nV[0])*(1-X_e.nV[1])));
  N_ref.nV[2] = 0.25*DMIN(4,DMAX(0,(1+X_e.nV[0])*(1+X_e.nV[1])));
  N_ref.nV[3] = 0.25*DMIN(4,DMAX(0,(1-X_e.nV[0])*(1+X_e.nV[1])));
  
  return N_ref;
}

/*********************************************************************/

/* Derivatives of the shape functions */
Matrix Q4_dN_Ref(Matrix X_e){

  int Ndim = NumberDimensions;
  
  /* Definition and allocation */
  Matrix dNdX_ref = MatAllocZ(4,Ndim);
  
  /* Fill the matrix */

  /* Node 1 */
  dNdX_ref.nM[0][0] = -0.25*DMIN(2,DMAX(0,(1-X_e.nV[1])));
  dNdX_ref.nM[0][1] = -0.25*DMIN(2,DMAX(0,(1-X_e.nV[0]))); 

  /* Node 2 */
  dNdX_ref.nM[1][0] = +0.25*DMIN(2,DMAX(0,(1-X_e.nV[1])));
  dNdX_ref.nM[1][1] = -0.25*DMIN(2,DMAX(0,(1+X_e.nV[0])));
  
  /* Node 3 */
  dNdX_ref.nM[2][0] = +0.25*DMIN(2,DMAX(0,(1+X_e.nV[1])));
  dNdX_ref.nM[2][1] = +0.25*DMIN(2,DMAX(0,(1+X_e.nV[0])));
  
  /* Node 4 */
  dNdX_ref.nM[3][0] = -0.25*DMIN(2,DMAX(0,(1+X_e.nV[1])));
  dNdX_ref.nM[3][1] = +0.25*DMIN(2,DMAX(0,(1-X_e.nV[0]))); 
  
  return dNdX_ref;
}

/*********************************************************************/

/* Jacobian of the transformation for the four-nodes quadrilateral */
Matrix Q4_F_Ref(Matrix X_NC_GP, Matrix X_GC_Nodes)
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
  int Ndim = NumberDimensions;
  /* Variable declaration */
  Matrix dNdX_Ref_GP;
  Tensor X_I;
  Tensor dNdx_I;
  Tensor F_Ref_I;
  Matrix F_Ref = MatAllocZ(Ndim,Ndim);

  /* 1º Evaluate the derivarive of the shape function in the GP */
  dNdX_Ref_GP = Q4_dN_Ref(X_NC_GP); 

  /* 2º Get the F_Ref doing a loop over the nodes of the element */
  for(int I = 0 ; I<4 ; I++){

    /* 3º Fill arrays for the tensorial product */
    X_I = memory_to_Tensor(X_GC_Nodes.nM[I], 1);
    dNdx_I = memory_to_Tensor(dNdX_Ref_GP.nM[I], 1);

    /* 4º Get the nodal contribution */
    F_Ref_I = get_dyadicProduct_Of(X_I,dNdx_I);

    /* 5º Increment the reference deformation gradient */
    for(int i = 0 ; i<Ndim ; i++){
      for(int j = 0 ; j<Ndim ; j++){
	F_Ref.nM[i][j] += F_Ref_I.N[i][j];
      }
    }
    /* 6º Free data of the nodal contribution */
    free_Tensor(F_Ref_I);    
  }
  
  /* 7º Free memory */
  FreeMat(dNdX_Ref_GP);

  /* 8º Output */
  return F_Ref;
}

/*********************************************************************/

/* Element gradient in the real element */
Matrix Q4_dN(Matrix X_EC, Matrix Element)
/*
  - Matrix X_EC_GP : Element coordinates of the gauss point
  - Matrix Element : Coordinates of the element (4 x Ndim)
*/
{
  
  /* Derivatives of the shape function evaluates in the GP (4 x Ndim) */
  Matrix dNdX;
  Matrix dNdX_T;

  /* 1º Evaluate the gradient of the shape function in the GP (4 x Ndim) */
  Matrix dNdX_Ref = Q4_dN_Ref(X_EC);
  
  /* 2º Get the Jacobian of the transformation evaluated in the GP */
  Matrix F = Q4_F_Ref(X_EC,Element);
  Matrix F_m1 = Get_Inverse(F);
  Matrix F_Tm1 = Transpose_Mat(F_m1);
  FreeMat(F);
  FreeMat(F_m1);
  Matrix dNdX_Ref_T = Transpose_Mat(dNdX_Ref);
  FreeMat(dNdX_Ref);
  
  /* 5º Get the gradient of the shape functions in global coordinates */
  dNdX_T = Scalar_prod(F_Tm1, dNdX_Ref_T);

  /* Free memory */
  FreeMat(F_Tm1);  
  FreeMat(dNdX_Ref_T);

  dNdX = Transpose_Mat(dNdX_T);
  FreeMat(dNdX_T);  

  
  /* 6º Return result */
  return dNdX;
}

/*********************************************************************/

/* Global coordinates of the four nodes quadrilateral */
Matrix Q4_Xi_to_X(Matrix Xi, Matrix Element)
/*
  This function evaluate the position of the GP in the element,
  and get it global coordiantes    
*/
{
  int Ndim = NumberDimensions;
  /* 1º Evaluate the Q4 element in the element coordinates */
  Matrix N = Q4_N(Xi);

  /* 2º Allocate the output coordinates */
  Matrix X = MatAllocZ(Ndim,1);

  /* 3º Get the global coordinates for this element coordiantes in this element */
  for(int I = 0 ; I<4 ; I++){
    X.nV[0] += N.nV[I]*Element.nM[I][0];
    X.nV[1] += N.nV[I]*Element.nM[I][1];
  }

  /* 4º Free memory */
  FreeMat(N);

  /* 5º Output */
  return X;
 
}

/*********************************************************************/

void Q4_X_to_Xi(Matrix Xi, Matrix X, Matrix Element)
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
  
  Xi = Newton_Rapson(Q4_Xi_to_X, Element, Q4_F_Ref, Element, X, Xi);
  
}

/*********************************************************************/
