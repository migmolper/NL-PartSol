/***********************************************/
/******* 3D cuadrilateral linear element *******/
/***********************************************/

/*     (7)      (6)   */
/*       o-------o    */
/*      /.      /|    */
/* (4) / . (5) / |    */
/*    o--o----o..o    */
/*    | / (3) | / (2) */
/*    |/      |/      */
/*    o-------o       */
/*   (0)     (1)      */

#include "nl-partsol.h"

/*
  Auxiliar functions
 */
static Matrix F_Ref__Q8__(Matrix, Matrix);
static Matrix Xi_to_X__Q8__(Matrix, Matrix);

/*********************************************************************/

void initialize__Q8__(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
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
    MPM_Mesh.NumberNodes[p] = 8;
    
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
      CoordElement = get_nodes_coordinates__MeshTools__(MPM_Mesh.ListNodes[p],FEM_Mesh.Coordinates);

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
Matrix N__Q8__(Matrix Xi)
{
  
  /* Definition and allocation */
  Matrix N_ref =  allocZ__MatrixLib__(1,8);
  
  N_ref.nV[0] = 0.125*DMIN(8,DMAX(0,(1.-Xi.nV[0])*(1.-Xi.nV[1])*(1.-Xi.nV[2])));
  N_ref.nV[1] = 0.125*DMIN(8,DMAX(0,(1.+Xi.nV[0])*(1.-Xi.nV[1])*(1.-Xi.nV[2])));
  N_ref.nV[2] = 0.125*DMIN(8,DMAX(0,(1.+Xi.nV[0])*(1.+Xi.nV[1])*(1.-Xi.nV[2])));
  N_ref.nV[3] = 0.125*DMIN(8,DMAX(0,(1.-Xi.nV[0])*(1.+Xi.nV[1])*(1.-Xi.nV[2])));
  N_ref.nV[4] = 0.125*DMIN(8,DMAX(0,(1.-Xi.nV[0])*(1.-Xi.nV[1])*(1.+Xi.nV[2])));
  N_ref.nV[5] = 0.125*DMIN(8,DMAX(0,(1.+Xi.nV[0])*(1.-Xi.nV[1])*(1.+Xi.nV[2])));
  N_ref.nV[6] = 0.125*DMIN(8,DMAX(0,(1.+Xi.nV[0])*(1.+Xi.nV[1])*(1.+Xi.nV[2])));
  N_ref.nV[7] = 0.125*DMIN(8,DMAX(0,(1.-Xi.nV[0])*(1.+Xi.nV[1])*(1.+Xi.nV[2])));
    
  return N_ref;
}

/*********************************************************************/

/* Derivatives of the shape functions */
Matrix dN_Ref__Q8__(Matrix Xi)
{

  int Ndim = NumberDimensions;
  
  /* Definition and allocation */
  Matrix dNdX_ref = allocZ__MatrixLib__(8,Ndim);
  
  /* Fill the matrix */
  dNdX_ref.nM[0][0] = -0.125*DMIN(4,DMAX(0,(1.-Xi.nV[1])*(1.-Xi.nV[2])));
  dNdX_ref.nM[0][1] = -0.125*DMIN(4,DMAX(0,(1.-Xi.nV[0])*(1.-Xi.nV[2])));
  dNdX_ref.nM[0][2] = -0.125*DMIN(4,DMAX(0,(1.-Xi.nV[0])*(1.-Xi.nV[1])));

  dNdX_ref.nM[1][0] =  0.125*DMIN(4,DMAX(0,(1.-Xi.nV[1])*(1.-Xi.nV[2])));
  dNdX_ref.nM[1][1] = -0.125*DMIN(4,DMAX(0,(1.+Xi.nV[0])*(1.-Xi.nV[2])));
  dNdX_ref.nM[1][2] = -0.125*DMIN(4,DMAX(0,(1.+Xi.nV[0])*(1.-Xi.nV[1])));

  dNdX_ref.nM[2][0] =  0.125*DMIN(4,DMAX(0,(1.+Xi.nV[1])*(1.-Xi.nV[2])));
  dNdX_ref.nM[2][1] =  0.125*DMIN(4,DMAX(0,(1.+Xi.nV[0])*(1.-Xi.nV[2])));
  dNdX_ref.nM[2][2] = -0.125*DMIN(4,DMAX(0,(1.+Xi.nV[0])*(1.+Xi.nV[1])));

  dNdX_ref.nM[3][0] = -0.125*DMIN(4,DMAX(0,(1.+Xi.nV[1])*(1.-Xi.nV[2])));
  dNdX_ref.nM[3][1] =  0.125*DMIN(4,DMAX(0,(1.-Xi.nV[0])*(1.-Xi.nV[2])));
  dNdX_ref.nM[3][2] = -0.125*DMIN(4,DMAX(0,(1.-Xi.nV[0])*(1.+Xi.nV[1])));

  dNdX_ref.nM[4][0] = -0.125*DMIN(4,DMAX(0,(1.-Xi.nV[1])*(1.+Xi.nV[2])));
  dNdX_ref.nM[4][1] = -0.125*DMIN(4,DMAX(0,(1.-Xi.nV[0])*(1.+Xi.nV[2])));
  dNdX_ref.nM[4][2] =  0.125*DMIN(4,DMAX(0,(1.-Xi.nV[0])*(1.-Xi.nV[1])));

  dNdX_ref.nM[5][0] =  0.125*DMIN(4,DMAX(0,(1.-Xi.nV[1])*(1.+Xi.nV[2])));
  dNdX_ref.nM[5][1] = -0.125*DMIN(4,DMAX(0,(1.+Xi.nV[0])*(1.+Xi.nV[2])));
  dNdX_ref.nM[5][2] =  0.125*DMIN(4,DMAX(0,(1.+Xi.nV[0])*(1.-Xi.nV[1])));

  dNdX_ref.nM[6][0] =  0.125*DMIN(4,DMAX(0,(1.+Xi.nV[1])*(1.+Xi.nV[2])));
  dNdX_ref.nM[6][1] =  0.125*DMIN(4,DMAX(0,(1.+Xi.nV[0])*(1.+Xi.nV[2])));
  dNdX_ref.nM[6][2] =  0.125*DMIN(4,DMAX(0,(1.+Xi.nV[0])*(1.+Xi.nV[1])));

  dNdX_ref.nM[7][0] = -0.125*DMIN(4,DMAX(0,(1.+Xi.nV[1])*(1.+Xi.nV[2])));
  dNdX_ref.nM[7][1] =  0.125*DMIN(4,DMAX(0,(1.-Xi.nV[0])*(1.+Xi.nV[2])));
  dNdX_ref.nM[7][2] =  0.125*DMIN(4,DMAX(0,(1.-Xi.nV[0])*(1.+Xi.nV[1])));
    
  return dNdX_ref;
}

/*********************************************************************/

/* Jacobian of the transformation for the four-nodes quadrilateral */
static Matrix F_Ref__Q8__(Matrix Xi, Matrix Element)
/*
  Get the jacobian of the transformation of the reference element :

  J = grad(N_{\alpha}) \cdot x_{\alpha}
  
  Inputs :
  - X_g -> This are the coordinates of the nodes (8x3)
  - dNdX_Ref_GP -> Derivative gradient evaluated in the GP (3x8)

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
  dNdX_Ref_GP = dN_Ref__Q8__(Xi);

  /* 2º Get the F_Ref doing a loop over the nodes of the element */
  for(int I = 0 ; I<8 ; I++)
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
Matrix dN__Q8__(Matrix Xi, Matrix Element)
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
  Matrix dNdX_Ref   = dN_Ref__Q8__(Xi);
  Matrix dNdX_Ref_T = transpose__MatrixLib__(dNdX_Ref);

  /* 2º Get the Jacobian of the transformation evaluated in the GP */
  Matrix F     = F_Ref__Q8__(Xi,Element);
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
static Matrix Xi_to_X__Q8__(Matrix Xi, Matrix Element)
/*
  This function evaluate the position of the GP in the element,
  and get it global coordiantes    
*/
{
  int Ndim = NumberDimensions;

  /* 1º Evaluate the Q8 element in the element coordinates */
  Matrix N = N__Q8__(Xi);

  /* 2º Allocate the output coordinates */
  Matrix X = allocZ__MatrixLib__(Ndim,1);

  /* 3º Get the global coordinates for this element coordiantes in this element */
  for(int I = 0 ; I<8 ; I++)
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

void X_to_Xi__Q8__(Matrix Xi, Matrix X, Matrix Element)
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
  
  Xi = Newton_Rapson(Xi_to_X__Q8__, Element, F_Ref__Q8__, Element, X, Xi);
  
}

/*********************************************************************/


Matrix element_to_particles__Q8__(int Num_GP)
/*
  Given a Q8 element, transform it in to a set of MPM particles
*/
{
  int Ndim = NumberDimensions;
  Matrix Xi = allocZ__MatrixLib__(Num_GP,Ndim);
  double a;

  switch(Num_GP)
    {
    case 1 :
      a = 0;
      Xi.nV[0] = a;
      Xi.nV[1] = a;
      Xi.nV[2] = a;
      break;
      
    case 8 :
      a = 0.577350269189626;
      
      Xi.nM[0][0] = - a;
      Xi.nM[0][1] = - a;
      Xi.nM[0][2] = - a;

      Xi.nM[1][0] =   a;
      Xi.nM[1][1] = - a;
      Xi.nM[1][2] = - a;

      Xi.nM[2][0] =   a;
      Xi.nM[2][1] =   a;
      Xi.nM[2][2] = - a;

      Xi.nM[3][0] = - a;
      Xi.nM[3][1] =   a;
      Xi.nM[3][2] = - a;

      Xi.nM[4][0] = - a;
      Xi.nM[4][1] = - a;
      Xi.nM[4][2] =   a;

      Xi.nM[5][0] =   a;
      Xi.nM[5][1] = - a;
      Xi.nM[5][2] =   a;

      Xi.nM[6][0] =   a;
      Xi.nM[6][1] =   a;
      Xi.nM[6][2] =   a;

      Xi.nM[7][0] = - a;
      Xi.nM[7][1] =   a;
      Xi.nM[7][2] =   a;      
      break;

    case 27:
      a = 0.774596669241483;

      Xi.nM[0][0] = - a;
      Xi.nM[0][1] = - a;
      Xi.nM[0][2] = - a;

      Xi.nM[1][0] =   a;
      Xi.nM[1][1] = - a;
      Xi.nM[1][2] = - a;

      Xi.nM[2][0] =   a;
      Xi.nM[2][1] =   a;
      Xi.nM[2][2] = - a;

      Xi.nM[3][0] = - a;
      Xi.nM[3][1] =   a;
      Xi.nM[3][2] = - a;

      Xi.nM[4][0] = - a;
      Xi.nM[4][1] = - a;
      Xi.nM[4][2] =   a;

      Xi.nM[5][0] =   a;
      Xi.nM[5][1] = - a;
      Xi.nM[5][2] =   a;

      Xi.nM[6][0] =   a;
      Xi.nM[6][1] =   a;
      Xi.nM[6][2] =   a;

      Xi.nM[7][0] = - a;
      Xi.nM[7][1] =   a;
      Xi.nM[7][2] =   a;

      Xi.nM[8][0] =   0;
      Xi.nM[8][1] = - a;
      Xi.nM[8][2] = - a;

      Xi.nM[9][0] =   a;
      Xi.nM[9][1] =   0;
      Xi.nM[9][2] = - a;

      Xi.nM[10][0] =   0;
      Xi.nM[10][1] =   a;
      Xi.nM[10][2] = - a;

      Xi.nM[11][0] = - a;
      Xi.nM[11][1] =   0;
      Xi.nM[11][2] = - a;

      Xi.nM[12][0] = - a;
      Xi.nM[12][1] = - a;
      Xi.nM[12][2] =   0;

      Xi.nM[13][0] =   a;
      Xi.nM[13][1] = - a;
      Xi.nM[13][2] =   0;

      Xi.nM[14][0] =   a;
      Xi.nM[14][1] =   a;
      Xi.nM[14][2] =   0;

      Xi.nM[15][0] = - a;
      Xi.nM[15][1] =   a;
      Xi.nM[15][2] =   0;

      Xi.nM[16][0] =   0;
      Xi.nM[16][1] = - a;
      Xi.nM[16][2] =   a;

      Xi.nM[17][0] =   a;
      Xi.nM[17][1] =   0;
      Xi.nM[17][2] =   a;

      Xi.nM[18][0] =   0;
      Xi.nM[18][1] =   a;
      Xi.nM[18][2] =   a;

      Xi.nM[19][0] = - a;
      Xi.nM[19][1] =   0;
      Xi.nM[19][2] = - a;

      Xi.nM[20][0] =   0;
      Xi.nM[20][1] =   0;
      Xi.nM[20][2] = - a;

      Xi.nM[21][0] =   0;
      Xi.nM[21][1] = - a;
      Xi.nM[21][2] =   0;

      Xi.nM[22][0] =   a;
      Xi.nM[22][1] =   0;
      Xi.nM[22][2] =   0;

      Xi.nM[23][0] =   0;
      Xi.nM[23][1] =   a;
      Xi.nM[23][2] =   0;

      Xi.nM[24][0] = - a;
      Xi.nM[24][1] =   0;
      Xi.nM[24][2] =   0;

      Xi.nM[25][0] =   0;
      Xi.nM[25][1] =   0;
      Xi.nM[25][2] =   a;

      Xi.nM[26][0] =   0;
      Xi.nM[26][1] =   0;
      Xi.nM[26][2] =   0;
      break;

    default :
      fprintf(stderr,"%s : %s \n",
	      "Error in element_to_particles__Q8__()",
	      "Wrong number of particles per element");
      exit(EXIT_FAILURE);
    }

  return Xi;
}


/*********************************************************************/
