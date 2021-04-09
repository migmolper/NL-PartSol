#include "nl-partsol.h"

/*
  Auxiliar functions
 */
static Matrix F_Ref__T3__(Matrix,Matrix);
static Matrix Xi_to_X__T3__(Matrix,Matrix);

/*********************************************************************/

void initialize__T3__(
  GaussPoint MPM_Mesh, 
  Mesh FEM_Mesh)
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
  for(int p = 0 ; p<Np ; p++)
  {

    /* Asign the number of nodes */
    MPM_Mesh.NumberNodes[p] = 3;
    
    /* Get the global and local coodinates of the particle */ 
    X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
    Xi_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_EC.nM[p]);
    
    /* Check for each element of the mesh */
    for(int i = 0 ; i<Nelem ; i++)
    {

      /* Get connectivity of the element */
      Elem_p = FEM_Mesh.Connectivity[i];
      
      /* 5º Check out if the GP is in the Element */
      if(inout_convex_set__MeshTools__(X_p, Elem_p, FEM_Mesh.Coordinates))
      {

        /* With the element connectivity get the node close to the particle */
        MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,Elem_p,FEM_Mesh.Coordinates);

        /* Asign connectivity */
        MPM_Mesh.ListNodes[p] = copy__SetLib__(Elem_p);

        /* Active those nodes that interact with the particle */
        asign_to_nodes__Particles__(p, MPM_Mesh.ListNodes[p], FEM_Mesh);

        /* Get the coordinates of the element vertex */
        CoordElement = get_nodes_coordinates__MeshTools__(MPM_Mesh.ListNodes[p], FEM_Mesh.Coordinates);

        /* Compute local coordinates of the particle in this element */
        X_to_Xi__T3__(Xi_p,X_p,CoordElement);

        /* Free coordinates of the element */
        free__MatrixLib__(CoordElement);

        break;
  
      }
      
    }

  }
  
}

/*********************************************************************/

Matrix N__T3__(
  Matrix X_e)
{

  /* Definition and allocation */
  Matrix N_ref =  alloc__MatrixLib__(1,3);

  /* Fill the array */
  N_ref.nV[0] = 1 - X_e.nV[0] - X_e.nV[1]; 
  N_ref.nV[1] = X_e.nV[0];
  N_ref.nV[2] = X_e.nV[1];
  
  return N_ref;
}

/*********************************************************************/

Matrix dN_Ref__T3__(
  Matrix X_e)
{
  
  /* Definition and allocation */
  Matrix dNdX_ref = alloc__MatrixLib__(2,3);
  
  /* Fill the matrix */
  /* Node 0 */
  dNdX_ref.nM[0][0] = - 1; /* \frac{\partial N0}{\partial \xi} */
  dNdX_ref.nM[1][0] = - 1; /* \frac{\partial N0}{\partial \eta} */
  /* Node 1 */
  dNdX_ref.nM[0][1] = + 1; /* \frac{\partial N1}{\partial \xi} */
  dNdX_ref.nM[1][1] = + 0; /* \frac{\partial N1}{\partial \eta} */
  /* Node 2 */
  dNdX_ref.nM[0][2] = + 0; /* \frac{\partial N2}{\partial \xi} */
  dNdX_ref.nM[1][2] = + 1; /* \frac{\partial N2}{\partial \eta} */
  
  return dNdX_ref;
}


/*********************************************************************/

/* Deformation gradient of the reference element for the three-nodes triangle */
static Matrix F_Ref__T3__(
  Matrix X_NC_GP,
  Matrix X_GC_Nodes)
/*
  Get the deformation gradient of the transformation of the reference element :

  F = grad(N_{\alpha}) \cdot x_{\alpha}
  
  Inputs :
  - X_NC_GP -> This are the natural coordinates of the GP
  - X_GC_Nodes -> Global coordinates of the element nodes

  Output :
  - F -> Deformation gradient of the reference element
  evaluated in the GP
*/
{
  /* Variable declaration */
  Matrix dNdX_Ref_GP;
  Matrix X_alpha = alloc__MatrixLib__(2,1);
  Matrix dNdx_alpha = alloc__MatrixLib__(1,2);
  Matrix F_Ref_alpha;
  Matrix F_Ref = allocZ__MatrixLib__(2,2);

  /* 1º Evaluate the derivarive of the shape function in the GP */
  dNdX_Ref_GP = dN_Ref__T3__(X_NC_GP); 

  /* 2º Get the F_Ref doing a loop over the nodes of the element */
  for(int i = 0 ; i<3 ; i++)
  {

    /* 3º Fill arrays for the tensorial product */
    for(int j = 0 ; j<2 ; j++)
    {
      X_alpha.nV[j] = X_GC_Nodes.nM[i][j];
      dNdx_alpha.nV[j] = dNdX_Ref_GP.nM[j][i];
    }

    /* 4º Get the nodal contribution */
    F_Ref_alpha = dyadic_product__MatrixLib__(X_alpha,dNdx_alpha);

    /* 5º Increment the reference deformation gradient */
    F_Ref = increment__MatrixLib__(F_Ref, F_Ref_alpha);

    /* 6º Free data of the nodal contribution */
    free__MatrixLib__(F_Ref_alpha);
    
  }
  
  /* 7º Free memory */
  free__MatrixLib__(dNdX_Ref_GP);
  free__MatrixLib__(X_alpha);
  free__MatrixLib__(dNdx_alpha);

  /* 8º Output */
  return F_Ref;
}

/*********************************************************************/

/* Element gradient in the real element */
Matrix dN__T3__(
  Matrix X_EC_GP,
  Matrix Element)
/*
  - Matrix X_EC_GP : Element coordinates of the gauss point
  - Matrix Element : Coordinates of the element (NumNodesElem x NumberDimensions)
*/
{
  
  /* 0º Definition and allocation */
  Matrix dNdX_Ref_GP; /* Derivative of the shape function evaluated in the GP (Ndim x Nnodes) */
  Matrix F_GP; /* Deformation gradient of the transformation evaluated in the GP (Ndim x Ndim) */
  Matrix F_GP_m1; /* Inverse of the deformation gradient */
  Matrix F_GP_Tm1; /* Transpose of the deformation gradient */
  Matrix dNdx_GP; /* Derivatives of the shape function evaluates in the GP (Ndim x Ndim) */

  /* 1º Evaluate the gradient of the shape function in the GP */
  dNdX_Ref_GP = dN_Ref__T3__(X_EC_GP);
  
  /* 2º Get the deformation gradient of the transformation evaluated in the GP */
  F_GP = F_Ref__T3__(X_EC_GP,Element);
    
  /* 3º Get the inverse of the deformation gradient */
  F_GP_m1 = inverse__MatrixLib__(F_GP);
  free__MatrixLib__(F_GP);
  
  /* 4º Get the transpose of the inverse of the deformation gradient */
  F_GP_Tm1 = transpose__MatrixLib__(F_GP_m1);
  free__MatrixLib__(F_GP_m1);
  
  /* 5º Get the gradient of the shape functions in global coordinates */
  dNdx_GP = matrix_product__MatrixLib__(F_GP_Tm1,dNdX_Ref_GP);
  free__MatrixLib__(F_GP_Tm1);
  free__MatrixLib__(dNdX_Ref_GP);

  /* 6º Return result */
  return dNdx_GP;
}

/*********************************************************************/

/* Global coordinates of the four nodes quadrilateral */
static Matrix Xi_to_X__T3__(Matrix X_NC_GP,Matrix X_GC_Nodes)
/*
This function evaluate the position of the GP in the element, and get it global coordiantes    
 */
{
  /* 0º Variable declaration */
  Matrix N_ref;
  Matrix X_GC_GP;

  /* 1º Evaluate the Q4 element in the element coordinates */
  N_ref = N__T3__(X_NC_GP);

  /* 2º Allocate the output coordinates */
  X_GC_GP = alloc__MatrixLib__(2,1);

  /* 3º Get the global coordinates for this element coordiantes in this element */
  X_GC_GP.nV[0] = N_ref.nV[0]*X_GC_Nodes.nM[0][0] + N_ref.nV[1]*X_GC_Nodes.nM[1][0] + N_ref.nV[2]*X_GC_Nodes.nM[2][0];
  
  X_GC_GP.nV[1] = N_ref.nV[0]*X_GC_Nodes.nM[0][1] + N_ref.nV[1]*X_GC_Nodes.nM[1][1] + N_ref.nV[2]*X_GC_Nodes.nM[2][1];

  /* 4º Free memory */
  free__MatrixLib__(N_ref);

  /* 5º Output */
  return X_GC_GP;
 
}

/*********************************************************************/

void X_to_Xi__T3__(
  Matrix X_EC_GP,
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
  
  X_EC_GP = Newton_Rapson(Xi_to_X__T3__,Element_GC_Nod,F_Ref__T3__,Element_GC_Nod,X_GC_GP,X_EC_GP);
}

/*********************************************************************/

void element_to_particles__T3__(
  Matrix X_p,
  Mesh FEM_Mesh,
  int GPxElement)
{

  int Ndim = NumberDimensions;
  int NumElemMesh = FEM_Mesh.NumElemMesh;
  int NumNodesElem = 3;
  Matrix N_GP;
  Matrix Xi_p = allocZ__MatrixLib__(GPxElement,Ndim);
  Matrix Xi_p_j;
  Element Element;
  int Node;

  switch(GPxElement)
  {
    case 1:
    Xi_p.nV[0] = 0.0;
    Xi_p.nV[1] = 0.0;
    break;

    case 3:
    Xi_p.nM[0][0] =  0.16666666666;
    Xi_p.nM[0][1] =  0.16666666666;
    Xi_p.nM[1][0] =  0.66666666666;
    Xi_p.nM[1][1] =  0.16666666666;
    Xi_p.nM[2][0] =  0.16666666666;
    Xi_p.nM[2][1] =  0.66666666666;
    break;

    default :
    fprintf(stderr,"%s : %s \n","Error in element_to_particles__T3__()","Wrong number of particles per element");
    exit(EXIT_FAILURE);
  }

  /* Get the coordinate of the center */
  for(int i = 0 ; i<NumElemMesh ; i++)
  {

    Element = nodal_set__Particles__(i, FEM_Mesh.Connectivity[i],FEM_Mesh.NumNodesElem[i]);

    for(int j = 0 ; j<GPxElement ; j++)
    {
      
      /* Evaluate the shape function in the GP position */
      if(GPxElement == 1)
      {
        Xi_p_j.nV = Xi_p.nM[j]; 
        N_GP = N__T3__(Xi_p_j);
      }
      else
      {
        N_GP = N__T3__(Xi_p);
      }

      for(int k = 0 ; k<NumNodesElem ; k++)
      {
        /* Connectivity of each element */
        Node = Element.Connectivity[k];

        for(int l = 0 ; l<Ndim ; l++)
        {
          X_p.nM[i*GPxElement+j][l] += N_GP.nV[k]*FEM_Mesh.Coordinates.nM[Node][l];
        }

      }

      /* Free value of the shape function in the GP */
      free__MatrixLib__(N_GP);
    }

    free(Element.Connectivity);

  }
  /* Free auxiliar matrix with the coordinates */
  free__MatrixLib__(Xi_p);

}

/*********************************************************************/

double min_DeltaX__T3__(ChainPtr Element_Connectivity, Matrix Coordinates)
{
  /* Auxiliar variables of the function */
  int Ndim = NumberDimensions;
  int NumNodesElem = 3; /* Number of nodes of each element */
  int Node_k;
  Matrix Poligon; /* Element Poligon */
  Matrix X_eval = allocZ__MatrixLib__(1,Ndim); /* Where to evaluate the shape function */
  X_eval.nV[0] = 0.0;
  X_eval.nV[1] = 0.0;
  Matrix dNdx; /* Gradient of the shapefunction for each node */
  double MinElementSize = 10e16;

  /* 
    Fill the triangular element with the coordinates of the nodes
  */
  Poligon = allocZ__MatrixLib__(NumNodesElem,Ndim);

  for(int k = 0; k<NumNodesElem; k++)
  {
    Node_k = Element_Connectivity->I;

    for(int l = 0 ; l<Ndim ; l++)
    {
      Poligon.nM[k][l] = Coordinates.nM[Node_k][l];
    }

    Element_Connectivity = Element_Connectivity->next;

  }

  /*
    Get the gradient of the triangle
  */
  dNdx = dN__T3__(X_eval,Poligon);
      
  /*
    Get the minimum minimum height of the triangle
  */
  for(int j = 0 ; j<NumNodesElem ; j++)
  {
    MinElementSize = DMIN(MinElementSize,1/sqrt(dNdx.nM[0][j]*dNdx.nM[0][j] + dNdx.nM[1][j]*dNdx.nM[1][j]));
  }

  /*
    Free memory
  */
  free__MatrixLib__(Poligon);
  free__MatrixLib__(dNdx);
  free__MatrixLib__(X_eval);

  return MinElementSize;
}

/*********************************************************************/
