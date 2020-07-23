#include "nl-partsol.h"

/*
  Auxiliar functions
*/
static Matrix F_Ref__Q4__(Matrix,Matrix);
static Matrix Xi_to_X__Q4__(Matrix,Matrix);

/*********************************************************************/

void initialize__Q4__(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{
  /* Variables for the GP coordinates */
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nelem = FEM_Mesh.NumElemMesh;

  /*Definition of auxiliar global and local coordinates */
  Matrix X_p, Xi_p;

  /* Variable with stores the conectivity of the element of the particle */
  ChainPtr Elem_p;

  /* Auxiliar variable to loop in the list of tributary nodes of the particle */
  ChainPtr ListNodes_p;

  /* Matrix with the coordinate of the nodes in the element */
  Matrix CoordElement;
 
  /* Loop over the particles to initialize them */
  for(int p = 0 ; p<Np ; p++){

    /* Asign the number of nodes */
    MPM_Mesh.NumberNodes[p] = 4;
    
    /* Get the global and local coodinates of the particle */ 
    X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
    Xi_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_EC.nM[p]);
    
    /* Check for each element of the mesh */
    for(int i = 0 ; i<Nelem ; i++){

      /* Get connectivity of the element */
      Elem_p = FEM_Mesh.Connectivity[i];
      
      /* 5º Check out if the GP is in the Element */
      if(inout_convex_set__MeshTools__(X_p, Elem_p, FEM_Mesh.Coordinates)){

	/* With the element connectivity get the node close to the particle */
	MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,Elem_p,FEM_Mesh.Coordinates);
	
	/* Asign connectivity */
	MPM_Mesh.ListNodes[p] = copy__SetLib__(Elem_p);
	
	/* Active those nodes that interact with the particle */
	asign_to_nodes__Particles__(p, MPM_Mesh.ListNodes[p], FEM_Mesh);
	
	/* Get the coordinates of the element vertex */
	CoordElement = get_nodes_coordinates__MeshTools__(MPM_Mesh.ListNodes[p],
							  FEM_Mesh.Coordinates);

	/* Compute local coordinates of the particle in this element */
	X_to_Xi__Q4__(Xi_p,X_p,CoordElement);
	
	/* Free coordinates of the element */
	free__MatrixLib__(CoordElement);

	break;
	
      }
      
    }

  }
  
}

/*********************************************************************/

/* Shape functions */
Matrix N__Q4__(Matrix X_e){
  
  /* Definition and allocation */
  Matrix N_ref =  allocZ__MatrixLib__(1,4);

  /* Fill the array */
  N_ref.nV[0] = 0.25*DMIN(4,DMAX(0,(1-X_e.nV[0])*(1-X_e.nV[1])));
  N_ref.nV[1] = 0.25*DMIN(4,DMAX(0,(1+X_e.nV[0])*(1-X_e.nV[1])));
  N_ref.nV[2] = 0.25*DMIN(4,DMAX(0,(1+X_e.nV[0])*(1+X_e.nV[1])));
  N_ref.nV[3] = 0.25*DMIN(4,DMAX(0,(1-X_e.nV[0])*(1+X_e.nV[1])));
  
  return N_ref;
}

/*********************************************************************/

/* Derivatives of the shape functions */
Matrix dN_Ref__Q4__(Matrix X_e){

  int Ndim = NumberDimensions;
  
  /* Definition and allocation */
  Matrix dNdX_ref = allocZ__MatrixLib__(4,Ndim);
  
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
static Matrix F_Ref__Q4__(Matrix X_NC_GP, Matrix X_GC_Nodes)
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
  Matrix F_Ref = allocZ__MatrixLib__(Ndim,Ndim);

  /* 1º Evaluate the derivarive of the shape function in the GP */
  dNdX_Ref_GP = dN_Ref__Q4__(X_NC_GP);

  /* 2º Get the F_Ref doing a loop over the nodes of the element */
  for(int I = 0 ; I<4 ; I++){

    /* 3º Fill arrays for the tensorial product */
    X_I = memory_to_tensor__TensorLib__(X_GC_Nodes.nM[I],1);
    dNdx_I = memory_to_tensor__TensorLib__(dNdX_Ref_GP.nM[I],1);

    /* 4º Get the nodal contribution */
    F_Ref_I = dyadic_Product__TensorLib__(X_I,dNdx_I);

    /* 5º Increment the reference deformation gradient */
    for(int i = 0 ; i<Ndim ; i++){
      for(int j = 0 ; j<Ndim ; j++){
	F_Ref.nM[i][j] += F_Ref_I.N[i][j];
      }
    }
    /* 6º Free data of the nodal contribution */
    free__TensorLib__(F_Ref_I);    
  }
  
  /* 7º Free memory */
  free__MatrixLib__(dNdX_Ref_GP);

  /* 8º Output */
  return F_Ref;
}

/*********************************************************************/

/* Element gradient in the real element */
Matrix dN__Q4__(Matrix X_EC, Matrix Element)
/*
  - Matrix X_EC_GP : Element coordinates of the gauss point
  - Matrix Element : Coordinates of the element (4 x Ndim)
*/
{
  
  /* Derivatives of the shape function evaluates in the GP (4 x Ndim) */
  Matrix dNdX;
  Matrix dNdX_T;
    
  /* 1º Evaluate the gradient of the shape function in the GP (4 x Ndim) */
  Matrix dNdX_Ref = dN_Ref__Q4__(X_EC);

  /* 2º Get the Jacobian of the transformation evaluated in the GP */
  Matrix F = F_Ref__Q4__(X_EC,Element);
  Matrix F_m1 = inverse__MatrixLib__(F);  
  Matrix F_Tm1 = transpose__MatrixLib__(F_m1);
 
  free__MatrixLib__(F);
  free__MatrixLib__(F_m1);
  Matrix dNdX_Ref_T = transpose__MatrixLib__(dNdX_Ref);
  free__MatrixLib__(dNdX_Ref);
  
  /* 5º Get the gradient of the shape functions in global coordinates */
  dNdX_T = matrix_product__MatrixLib__(F_Tm1, dNdX_Ref_T);
  
  /* Free memory */
  free__MatrixLib__(F_Tm1);  
  free__MatrixLib__(dNdX_Ref_T);

  dNdX = transpose__MatrixLib__(dNdX_T);
  free__MatrixLib__(dNdX_T);  

  
  /* 6º Return result */
  return dNdX;
}

/*********************************************************************/

/* Global coordinates of the four nodes quadrilateral */
Matrix Xi_to_X__Q4__(Matrix Xi, Matrix Element)
/*
  This function evaluate the position of the GP in the element,
  and get it global coordiantes    
*/
{
  int Ndim = NumberDimensions;
  /* 1º Evaluate the Q4 element in the element coordinates */
  Matrix N = N__Q4__(Xi);

  /* 2º Allocate the output coordinates */
  Matrix X = allocZ__MatrixLib__(Ndim,1);

  /* 3º Get the global coordinates for this element coordiantes in this element */
  for(int I = 0 ; I<4 ; I++){
    X.nV[0] += N.nV[I]*Element.nM[I][0];
    X.nV[1] += N.nV[I]*Element.nM[I][1];
  }

  /* 4º Free memory */
  free__MatrixLib__(N);

  /* 5º Output */
  return X;
 
}

/*********************************************************************/

void X_to_Xi__Q4__(Matrix Xi, Matrix X, Matrix Element)
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
  Xi = Newton_Rapson(Xi_to_X__Q4__, Element, F_Ref__Q4__, Element, X, Xi);  
}

/*********************************************************************/
