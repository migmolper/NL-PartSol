

#include "nl-partsol.h"

/*
  Auxiliar functions
 */
static Matrix F_Ref__T4__(Matrix, Matrix);
static Matrix Xi_to_X__T4__(Matrix, Matrix);

/*********************************************************************/

void initialize__T4__(
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
  for(int p = 0 ; p<Np ; p++){

    /* Asign the number of nodes */
    MPM_Mesh.NumberNodes[p] = 4;
    
    /* Get the global and local coodinates of the particle */ 
    X_p  = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
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
        CoordElement = get_nodes_coordinates__MeshTools__(MPM_Mesh.ListNodes[p],FEM_Mesh.Coordinates);

        /* Compute local coordinates of the particle in this element */
        X_to_Xi__T4__(Xi_p,X_p,CoordElement);

        /* Free coordinates of the element */
        free__MatrixLib__(CoordElement);

        break;
	
      }
      
    }

  }
  
}

/*********************************************************************/

/* Shape functions */
Matrix N__T4__(
  Matrix Xi)
{
  
  /* Definition and allocation */
  Matrix N_ref = alloc__MatrixLib__(1,4);

  /* Fill the array */
  N_ref.nV[0] = DMIN(1,DMAX(0,Xi.nV[0]));
  N_ref.nV[1] = DMIN(1,DMAX(0,Xi.nV[1]));
  N_ref.nV[2] = DMIN(1,DMAX(0,Xi.nV[2]));
  N_ref.nV[3] = DMIN(1,DMAX(0,1.-Xi.nV[0]-Xi.nV[1]-Xi.nV[2]));
    
  return N_ref;
}

/*********************************************************************/

/* Derivatives of the shape functions */
Matrix dN_Ref__T4__(
  Matrix Xi)
{
  int Ndim = NumberDimensions;
  
  /* Definition and allocation */
  Matrix dNdX_ref = allocZ__MatrixLib__(4,Ndim);
  
  /* Fill the matrix */
  dNdX_ref.nM[0][0] =  1.;
  dNdX_ref.nM[1][1] =  1.;
  dNdX_ref.nM[2][2] =  1.;
  dNdX_ref.nM[3][0] = -1.;
  dNdX_ref.nM[3][1] = -1.;
  dNdX_ref.nM[3][2] = -1.;
    
  return dNdX_ref;
}


/*********************************************************************/

/* Deformation gradient of the reference element for the three-nodes triangle */
static Matrix F_Ref__T4__(
  Matrix Xi,
  Matrix Element)
/*
  Get the deformation gradient of the transformation of the reference element :

  F = grad(N_{\alpha}) \cdot x_{\alpha}
  
  Inputs :
  - Xi -> This are the natural coordinates of the GP
  - Element -> Global coordinates of the element nodes

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
  dNdX_Ref_GP = dN_Ref__T4__(Xi);

  /* 2º Get the F_Ref doing a loop over the nodes of the element */
  for(int I = 0 ; I<4 ; I++)
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
Matrix dN__T4__(
  Matrix Xi,
  Matrix Element)
/*
  - Matrix Xi_GP : Element coordinates of the gauss point
  - Matrix Element : Coordinates of the element (NumNodesElem x NumberDimensions)
*/
{
  
  /* Derivatives of the shape function evaluates in the GP (4 x Ndim) */
  Matrix dNdX;
  Matrix dNdX_T;
    
  /* 1º Evaluate the gradient of the shape function in the GP (4 x Ndim)
  and get the transpose */
  Matrix dNdX_Ref   = dN_Ref__T4__(Xi);
  Matrix dNdX_Ref_T = transpose__MatrixLib__(dNdX_Ref);

  /* 2º Get the Jacobian of the transformation evaluated in the GP */
  Matrix F     = F_Ref__T4__(Xi,Element);
  Matrix F_m1  = inverse__MatrixLib__(F);  
  Matrix F_Tm1 = transpose__MatrixLib__(F_m1);
  
  /* 5º Get the gradient of the shape functions in global coordinates */
  dNdX_T = matrix_product__MatrixLib__(F_Tm1, dNdX_Ref_T);
  dNdX   = transpose__MatrixLib__(dNdX_T);
  
  /* Free memory */
  free__MatrixLib__(F);
  free__MatrixLib__(F_m1);
  free__MatrixLib__(F_Tm1);  
  free__MatrixLib__(dNdX_Ref);
  free__MatrixLib__(dNdX_Ref_T);
  free__MatrixLib__(dNdX_T);  

  
  /* 6º Return result */
  return dNdX;
}

/*********************************************************************/

/* Global coordinates of the four nodes quadrilateral */
static Matrix Xi_to_X__T4__(
  Matrix Xi,
  Matrix Element)
/*
  This function evaluate the position of the GP in the element,
  and get it global coordiantes    
*/
{
  int Ndim = NumberDimensions;
  
  /* 1º Evaluate the T4 element in the element coordinates */
  Matrix N = N__T4__(Xi);
  
  /* 2º Allocate the output coordinates */
  Matrix X = allocZ__MatrixLib__(Ndim,1);

  /* 3º Get the global coordinates for this element coordiantes in this element */
  for(int I = 0 ; I<4 ; I++)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      X.nV[i] += N.nV[I]*Element.nM[I][i];
    }
  }

  /* 4º Free memory */
  free__MatrixLib__(N);

  /* 5º Output */
  return X;
 
}

/*********************************************************************/

void X_to_Xi__T4__(
  Matrix Xi,
  Matrix X,
  Matrix Element)
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
  Xi = Newton_Rapson(Xi_to_X__T4__,Element,F_Ref__T4__,Element,X,Xi);
}

/*********************************************************************/


void element_to_particles__T4__(
  Matrix X_p,
  Mesh FEM_Mesh,
  int GPxElement)
/*
  Given a T4 element, transform it in to a set of MPM particles
*/
{
  int Ndim = NumberDimensions;
  int NumElemMesh = FEM_Mesh.NumElemMesh;
  int NumNodesElem = 4;
  Matrix N_GP;
  Matrix Xi_p = allocZ__MatrixLib__(GPxElement,Ndim);
  Matrix Xi_p_j;
  Element Element;
  int Node;

  switch(GPxElement)
    {
    case 1 :
      
      Xi_p.nV[0] = 0.25;
      Xi_p.nV[1] = 0.25;
      Xi_p.nV[2] = 0.25;
      break;
      
    case 4 :      
      Xi_p.nM[0][0] = 0.138196601125010;
      Xi_p.nM[0][1] = 0.138196601125010;
      Xi_p.nM[0][2] = 0.138196601125010;
      Xi_p.nM[1][0] = 0.585410196624968;
      Xi_p.nM[1][1] = 0.138196601125010;
      Xi_p.nM[1][2] = 0.138196601125010;
      Xi_p.nM[2][0] = 0.138196601125010;
      Xi_p.nM[2][1] = 0.585410196624968;
      Xi_p.nM[2][2] = 0.138196601125010;
      Xi_p.nM[3][0] = 0.138196601125010;
      Xi_p.nM[3][1] = 0.138196601125010;
      Xi_p.nM[3][2] = 0.585410196624968;
      break;

    case 10:
      Xi_p.nM[0][0] = 0.108103018168070;
      Xi_p.nM[0][1] = 0.108103018168070;
      Xi_p.nM[0][2] = 0.108103018168070;      
      Xi_p.nM[1][0] = 0.816847572980459;
      Xi_p.nM[1][1] = 0.108103018168070;
      Xi_p.nM[1][2] = 0.108103018168070;
      Xi_p.nM[2][0] = 0.108103018168070;
      Xi_p.nM[2][1] = 0.816847572980459;
      Xi_p.nM[2][2] = 0.108103018168070;
      Xi_p.nM[3][0] = 0.108103018168070;
      Xi_p.nM[3][1] = 0.108103018168070;
      Xi_p.nM[3][2] = 0.816847572980459;
      Xi_p.nM[4][0] = 0.445948490915965;
      Xi_p.nM[4][1] = 0.108103018168070;
      Xi_p.nM[4][2] = 0.108103018168070;    
      Xi_p.nM[5][0] = 0.445948490915965;
      Xi_p.nM[5][1] = 0.445948490915965;
      Xi_p.nM[5][2] = 0.108103018168070;
      Xi_p.nM[6][0] = 0.108103018168070;
      Xi_p.nM[6][1] = 0.445948490915965;
      Xi_p.nM[6][2] = 0.108103018168070;
      Xi_p.nM[7][0] = 0.108103018168070;
      Xi_p.nM[7][1] = 0.108103018168070;
      Xi_p.nM[7][2] = 0.445948490915965;
      Xi_p.nM[8][0] = 0.445948490915965;
      Xi_p.nM[8][1] = 0.108103018168070;
      Xi_p.nM[8][2] = 0.445948490915965;
      Xi_p.nM[9][0] = 0.108103018168070;
      Xi_p.nM[9][1] = 0.445948490915965;
      Xi_p.nM[9][2] = 0.445948490915965;
      break;

    default :
      fprintf(stderr,"%s : %s \n",
	      "Error in element_to_particles__T4__()",
	      "Wrong number of particles per element");
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
        N_GP = N__T4__(Xi_p_j);
      }
      else
      {
        N_GP = N__T4__(Xi_p);
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

double min_DeltaX__T4__(ChainPtr Element_Connectivity, Matrix Coordinates)
{

  /* Auxiliar variables of the function */
  int Ndim = NumberDimensions;
  int Node_k;
  int Node_k1;
  int Node_k2;
  int NumNodesElem = 4;
  int NumSides = 4;
  double sqr_lenght = 0.0;
  double MinElementSize = 10e16;

  /*
    Fill the poligon with vectors
  */
  for(int k = 0; k<NumNodesElem; k++)
  {

    Node_k  = Element_Connectivity->I;
    Node_k1 = Element_Connectivity->next->I;

    sqr_lenght = 0.0;
    
    for(int l = 0 ; l<Ndim ; l++)
    {
      sqr_lenght += DSQR(Coordinates.nM[Node_k1][l] - Coordinates.nM[Node_k][l]);
    }

    MinElementSize = DMIN(MinElementSize,sqrt(sqr_lenght));

    Element_Connectivity = Element_Connectivity->next;

  }

  for(int k = 0; k<2; k++)
  {

    Node_k  = Element_Connectivity->I;
    Node_k2 = Element_Connectivity->next->next->I;

    sqr_lenght = 0.0;
    
    for(int l = 0 ; l<Ndim ; l++)
    {
      sqr_lenght += DSQR(Coordinates.nM[Node_k2][l] - Coordinates.nM[Node_k][l]);
    }

    MinElementSize = DMIN(MinElementSize,sqrt(sqr_lenght));

    Element_Connectivity = Element_Connectivity->next;

  }



  return MinElementSize;
}

/*********************************************************************/
