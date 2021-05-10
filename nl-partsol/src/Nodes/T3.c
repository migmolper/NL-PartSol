#include "nl-partsol.h"

/*
  Global variables
*/
double Thickness_Plain_Stress;


/*
  Auxiliar functions
 */
static Matrix F_Ref__T3__(Matrix,Matrix);
static Matrix Xi_to_X__T3__(Matrix,Matrix);

/*********************************************************************/

void initialize__T3__(
  Particle MPM_Mesh, 
  Mesh FEM_Mesh)
{
  /* Variables for the GP coordinates */
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nelem = FEM_Mesh.NumElemMesh;

  /*Definition of auxiliar global and local coordinates */
  Matrix X_p, Xi_p;

  /* Variable with stores the conectivity of the element of the particle */
  ChainPtr Elem_p_Connectivity;

  /* Auxiliar variable to loop in the list of tributary nodes of the particle */
  ChainPtr ListNodes_p;

  /* Matrix with the coordinate of the nodes in the element */
  Matrix Elem_p_Coordinates;
 
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
      Elem_p_Connectivity = FEM_Mesh.Connectivity[i];
      Elem_p_Coordinates = get_nodes_coordinates__MeshTools__(Elem_p_Connectivity, FEM_Mesh.Coordinates);
      
      /* 5º Check out if the GP is in the Element */
      if(FEM_Mesh.In_Out_Element(X_p,Elem_p_Coordinates))
      {

        /* With the element connectivity get the node close to the particle */
        MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,Elem_p_Connectivity,FEM_Mesh.Coordinates);

        /* Asign connectivity */
        MPM_Mesh.ListNodes[p] = copy__SetLib__(Elem_p_Connectivity);

        /* Active those nodes that interact with the particle */
        asign_to_nodes__Particles__(p, MPM_Mesh.I0[p], MPM_Mesh.ListNodes[p], FEM_Mesh);

        /* Compute local coordinates of the particle in this element */
        X_to_Xi__T3__(Xi_p,X_p,Elem_p_Coordinates);

        /* Free coordinates of the element */
        free__MatrixLib__(Elem_p_Coordinates);

        break;
  
      }

      /* Free coordinates of the element */
      free__MatrixLib__(Elem_p_Coordinates);
      
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
  Matrix dNdX_ref = alloc__MatrixLib__(3,2);
  
  /* Fill the matrix */
  /* Node 0 */
  dNdX_ref.nM[0][0] = - 1; /* \frac{\partial N0}{\partial \xi} */
  dNdX_ref.nM[0][1] = - 1; /* \frac{\partial N0}{\partial \eta} */
  /* Node 1 */
  dNdX_ref.nM[1][0] = + 1; /* \frac{\partial N1}{\partial \xi} */
  dNdX_ref.nM[1][1] = + 0; /* \frac{\partial N1}{\partial \eta} */
  /* Node 2 */
  dNdX_ref.nM[2][0] = + 0; /* \frac{\partial N2}{\partial \xi} */
  dNdX_ref.nM[2][1] = + 1; /* \frac{\partial N2}{\partial \eta} */
  
  return dNdX_ref;
}


/*********************************************************************/

/* Deformation gradient of the reference element for the three-nodes triangle */
static Matrix F_Ref__T3__(
  Matrix Xi,
  Matrix Element)
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
  int Ndim = NumberDimensions;

  /* Variable declaration */
  Matrix dNdX_Ref_GP;
  Tensor X_I;
  Tensor dNdx_I;
  Tensor F_Ref_I;
  Matrix F_Ref = allocZ__MatrixLib__(Ndim,Ndim);

  /* 1º Evaluate the derivarive of the shape function in the GP */
  dNdX_Ref_GP = dN_Ref__T3__(Xi);

  /* 2º Get the F_Ref doing a loop over the nodes of the element */
  for(int I = 0 ; I<3 ; I++)
  {

    /* 3º Fill arrays for the tensorial product */
    X_I    = memory_to_tensor__TensorLib__(Element.nM[I],1);
    dNdx_I = memory_to_tensor__TensorLib__(dNdX_Ref_GP.nM[I],1);

    /* 4º Get the nodal contribution */
    F_Ref_I = dyadic_Product__TensorLib__(X_I,dNdx_I);

    /* 5º Increment the reference deformation gradient */
    for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
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
Matrix dN__T3__(
  Matrix Xi,
  Matrix Element)
/*
  - Matrix Xi_GP : Element coordinates of the gauss point
  - Matrix Element : Coordinates of the element (8 x Ndim)
*/
{
  
  /* Derivatives of the shape function evaluates in the GP (8 x Ndim) */
  Matrix dNdX;
  Matrix dNdX_T;
    
  /* 1º Evaluate the gradient of the shape function in the GP (8 x Ndim)
  and transpose it */
  Matrix dNdX_Ref   = dN_Ref__T3__(Xi);
  Matrix dNdX_Ref_T = transpose__MatrixLib__(dNdX_Ref);

  /* 2º Get the Jacobian of the transformation evaluated in the GP */
  Matrix F     = F_Ref__T3__(Xi,Element);
  Matrix F_m1  = inverse__MatrixLib__(F);  
  Matrix F_Tm1 = transpose__MatrixLib__(F_m1);
  
  /* 3º Get the gradient of the shape functions in global coordinates */
  dNdX_T = matrix_product__MatrixLib__(F_Tm1, dNdX_Ref_T);
  dNdX   = transpose__MatrixLib__(dNdX_T);

  /* 4º Free memory */
  free__MatrixLib__(F);
  free__MatrixLib__(F_m1);
  free__MatrixLib__(F_Tm1);
  free__MatrixLib__(dNdX_Ref);
  free__MatrixLib__(dNdX_Ref_T);
  free__MatrixLib__(dNdX_T);
  
  /* 5º Return result */
  return dNdX;
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

bool in_out__T3__(
  Matrix X,
  Matrix Element)
{
  bool in_out = false;

  Matrix Xi = allocZ__MatrixLib__(2,1);

  Xi = Newton_Rapson(Xi_to_X__T3__, Element, F_Ref__T3__, Element, X, Xi);

  if((Xi.nV[0] >= 0.0) && 
    (Xi.nV[1]  >= 0.0) && 
    (Xi.nV[1] + Xi.nV[0] -1 <=  0.0))
  {
    in_out = true;    
  }

  free__MatrixLib__(Xi);

  return in_out;
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
  Matrix Xi_p_j = memory_to_matrix__MatrixLib__(1,Ndim,NULL);
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
        N_GP = N__T3__(Xi_p);
      }
      else
      {
        Xi_p_j.nV = Xi_p.nM[j]; 
        N_GP = N__T3__(Xi_p_j);
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
    MinElementSize = DMIN(MinElementSize,1/sqrt(dNdx.nM[j][0]*dNdx.nM[j][0] + dNdx.nM[j][1]*dNdx.nM[j][1]));
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

double volume__T3__(
  Matrix Element)
{
  int Ndim = NumberDimensions;
  double J_i;
  double Vol = 0;

  // Use 4 integration points to compute volume
  double table_w[3] = {1./6.,1./6.,1./6.};
  double table_X[3][2] = 
  {
    {0.600000000000000,0.200000000000000},
    {0.200000000000000,0.600000000000000},
    {0.200000000000000,0.200000000000000}
  };

  Matrix F_i;
  Matrix Xi = allocZ__MatrixLib__(2,1);

  for(int i = 0 ; i<3 ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
       Xi.nV[j] = table_X[i][j];
    }

    // Compute deformation gradient and jacobian of this integration point
    F_i = F_Ref__T3__(Xi,Element);
    J_i = I3__MatrixLib__(F_i);

    // Compute volume contribution
    Vol += fabs(J_i)*table_w[i];

    // Free memory
    free__MatrixLib__(F_i);
  }

  Vol *= Thickness_Plain_Stress;

  // Free memory
  free__MatrixLib__(Xi);

  return Vol;
}

/*********************************************************************/
