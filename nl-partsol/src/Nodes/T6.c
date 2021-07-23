#include "nl-partsol.h"

/*
  Global variables
*/
double Thickness_Plain_Stress;


/*
  Auxiliar functions
 */
static Matrix F_Ref__T6__(Matrix,Matrix);
static void   X_to_Xi__T6__(Matrix,Matrix,Matrix);
static Matrix Xi_to_X__T6__(Matrix,Matrix);
static ChainPtr tributary__T6__(int,Matrix,ChainPtr,Matrix);

/*********************************************************************/

void initialize__T6__(
  Particle MPM_Mesh, 
  Mesh FEM_Mesh)
{
  /* Variables for the GP coordinates */
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nelem = FEM_Mesh.NumElemMesh;

  bool Init_p;

  /*Definition of auxiliar global and local coordinates */
  Matrix X_p, Xi_p;

  /* Variable with stores the conectivity of the element of the particle */
  ChainPtr Elem_p_Connectivity;

  /* Auxiliar variable to loop in the list of tributary nodes of the particle */
  ChainPtr ListNodes_p;

  /* Matrix with the coordinate of the nodes in the element */
  Matrix Elem_p_Coordinates;

  /* Matrix with the coordinate of the nodes in the internal element */
  Matrix Internal_Elem_p_Coordinates;
 
  /* Loop over the particles to initialize them */
  for(int p = 0 ; p<Np ; p++)
  {

    /* Supose that the particle was not initilise */
    Init_p = false;

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
      if(in_out__T6__(X_p,Elem_p_Coordinates) == true)
      {

        /* Particle will be initilise */
        Init_p = true;

        // Assign the index of the element
        MPM_Mesh.Element_p[p] = i;

        /* With the element connectivity get the node close to the particle */
        MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,Elem_p_Connectivity,FEM_Mesh.Coordinates);

        /* Asign connectivity */
        MPM_Mesh.ListNodes[p] = tributary__T6__(p,X_p,Elem_p_Connectivity,Elem_p_Coordinates);

        /* Active those nodes that interact with the particle */
        asign_to_nodes__Particles__(p, MPM_Mesh.Element_p[p], MPM_Mesh.I0[p], MPM_Mesh.ListNodes[p], FEM_Mesh);

        Internal_Elem_p_Coordinates = get_nodes_coordinates__MeshTools__(MPM_Mesh.ListNodes[p], FEM_Mesh.Coordinates);

        /* Compute local coordinates of the particle in this element */
        X_to_Xi__T6__(Xi_p,X_p,Internal_Elem_p_Coordinates);

        /* Free coordinates of the element */
        free__MatrixLib__(Elem_p_Coordinates);
        free__MatrixLib__(Internal_Elem_p_Coordinates);

        break;
  
      }

      /* Free coordinates of the element */
      free__MatrixLib__(Elem_p_Coordinates);
      
    }

    if(!Init_p)
    {
      fprintf(stderr,"%s : %s %i\n",
        "Error in initialize__T6__()",
        "The search algorithm was unable to find particle",p);
      exit(EXIT_FAILURE);
    }

  }
  
}

/*********************************************************************/

Matrix N__T6__(
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

Matrix dN_Ref__T6__(
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
static Matrix F_Ref__T6__(
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
  dNdX_Ref_GP = dN_Ref__T6__(Xi);

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
Matrix dN__T6__(
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
  Matrix dNdX_Ref   = dN_Ref__T6__(Xi);
  Matrix dNdX_Ref_T = transpose__MatrixLib__(dNdX_Ref);

  /* 2º Get the Jacobian of the transformation evaluated in the GP */
  Matrix F     = F_Ref__T6__(Xi,Element);
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
static Matrix Xi_to_X__T6__(
  Matrix Xi,
  Matrix Element)
/*
This function evaluate the position of the GP in the element, and get it global coordiantes    
 */
{
  int Ndim = NumberDimensions;
  /* 1º Evaluate the Q4 element in the element coordinates */
  Matrix N = N__T6__(Xi);

  /* 2º Allocate the output coordinates */
  Matrix X = allocZ__MatrixLib__(Ndim,1);

  /* 3º Get the global coordinates for this element coordiantes in this element */
  for(int I = 0 ; I<3 ; I++)
  {
    X.nV[0] += N.nV[I]*Element.nM[I][0];
    X.nV[1] += N.nV[I]*Element.nM[I][1];
  }

  /* 4º Free memory */
  free__MatrixLib__(N);

  /* 5º Output */
  return X;
 
}

/*********************************************************************/

static void X_to_Xi__T6__(
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
  
  X_EC_GP = Newton_Rapson(Xi_to_X__T6__,Element_GC_Nod,F_Ref__T6__,Element_GC_Nod,X_GC_GP,X_EC_GP);
}

/*********************************************************************/

bool in_out__T6__(
  Matrix X,
  Matrix Element)
{
  double min[2] = {Element.nM[0][0],Element.nM[0][1]};
  double max[2] = {Element.nM[0][0],Element.nM[0][1]};

  for(int a = 1 ; a<6 ; a++)
  {
    for(int i = 0 ; i<2; i++)
    {
      min[i] = DMIN(min[i],Element.nM[a][i]);
      max[i] = DMAX(max[i],Element.nM[a][i]);
    }
  }

  Matrix Xi;
  Matrix Macro_Element;

  // Check if it is inside
  if((X.nV[0] <= max[0]) && 
    (X.nV[0] >= min[0]) && 
    (X.nV[1] <= max[1]) && 
    (X.nV[1] >= min[1]))
  {

    Xi = allocZ__MatrixLib__(2,1);
    Macro_Element = allocZ__MatrixLib__(3,2);

    Macro_Element.nM[0][0] = Element.nM[3][0];
    Macro_Element.nM[0][1] = Element.nM[3][1];

    Macro_Element.nM[1][0] = Element.nM[4][0];
    Macro_Element.nM[1][1] = Element.nM[4][1];
 
    Macro_Element.nM[2][0] = Element.nM[5][0];
    Macro_Element.nM[2][1] = Element.nM[5][1];

    X_to_Xi__T6__(Xi, X, Macro_Element);

    if((Xi.nV[0] >= 0.0) && 
      (Xi.nV[1]  >= 0.0) && 
      (Xi.nV[1] + Xi.nV[0] -1 <=  0.0))
    {
      free__MatrixLib__(Xi);
      free__MatrixLib__(Macro_Element);
      return true;    
    }
    else
    {
      free__MatrixLib__(Xi);
      free__MatrixLib__(Macro_Element);
      return false;
    }

  }
  else
  {
    return false;
  }

}


/*********************************************************************/

void element_to_particles__T6__(
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
    fprintf(stderr,"%s : %s \n","Error in element_to_particles__T6__()","Wrong number of particles per element");
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
        N_GP = N__T6__(Xi_p);
      }
      else
      {
        Xi_p_j.nV = Xi_p.nM[j]; 
        N_GP = N__T6__(Xi_p_j);
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

double min_DeltaX__T6__(ChainPtr Element_Connectivity, Matrix Coordinates)
{
  /* Auxiliar variables of the function */
  int Ndim = NumberDimensions;
  int NumNodesElem = 6; /* Number of nodes of each element */
  int NumNodesInternalElements = 3;
  int NumInternalElements = 4;
  int Node_i;
  int Node_k;
  ChainPtr * Internal_Element_Connectivity = (ChainPtr *)malloc(NumInternalElements*sizeof(ChainPtr));
  for(int i = 0; i<NumInternalElements ; i++)
  {
    Internal_Element_Connectivity[i] = NULL;
  }
  Matrix Poligon; /* Element Poligon */
  Matrix X_eval = allocZ__MatrixLib__(1,Ndim); /* Where to evaluate the shape function */
  X_eval.nV[0] = 0.0;
  X_eval.nV[1] = 0.0;
  Matrix dNdx; /* Gradient of the shapefunction for each node */
  double MinElementSize = 10e16;

  /* 
    Fill the triangular element with the coordinates of the nodes
  */
  Poligon = allocZ__MatrixLib__(NumNodesInternalElements,Ndim);

  for(int i = 0; i < NumNodesElem; i++)
  {

    Node_i = Element_Connectivity->I;

    switch(i) 
    {
      case 0 :
      push__SetLib__ (&Internal_Element_Connectivity[0], Node_i);
      push__SetLib__ (&Internal_Element_Connectivity[1], Node_i);
      push__SetLib__ (&Internal_Element_Connectivity[3], Node_i);
      break; 

      case 1 :
      push__SetLib__ (&Internal_Element_Connectivity[1], Node_i);
      push__SetLib__ (&Internal_Element_Connectivity[2], Node_i);
      push__SetLib__ (&Internal_Element_Connectivity[3], Node_i);
      break;

      case 2 :
      push__SetLib__ (&Internal_Element_Connectivity[0], Node_i);
      push__SetLib__ (&Internal_Element_Connectivity[1], Node_i);
      push__SetLib__ (&Internal_Element_Connectivity[2], Node_i);
      break; 

      case 3 :
      push__SetLib__ (&Internal_Element_Connectivity[3], Node_i);
      break;

      case 4 :
      push__SetLib__ (&Internal_Element_Connectivity[2], Node_i);

      break;

      case 5 :
      push__SetLib__ (&Internal_Element_Connectivity[0], Node_i); 

      break;
    }

    Element_Connectivity = Element_Connectivity->next;
  }

  /*
    Loop inside of the element and check each individual internal element
  */
  for(int i = 0 ; i<NumInternalElements ; i++)
  {

    /* 
      Fill the internal triangular element with the coordinates of the nodes
    */
    for(int k = 0; k<NumNodesInternalElements; k++)
    {
      Node_k = Internal_Element_Connectivity[i]->I;

      for(int l = 0 ; l<Ndim ; l++)
      {
        Poligon.nM[k][l] = Coordinates.nM[Node_k][l];
      }

      Internal_Element_Connectivity[i] = Internal_Element_Connectivity[i]->next;
    }

    /*
      Get the gradient of the triangle
    */
    dNdx = dN__T6__(X_eval,Poligon);

    /*
      Get the minimum minimum height of the triangle
    */
    for(int j = 0 ; j<NumNodesInternalElements ; j++)
    {
      MinElementSize = DMIN(MinElementSize,1/sqrt(dNdx.nM[j][0]*dNdx.nM[j][0] + dNdx.nM[j][1]*dNdx.nM[j][1]));
    }

    /*
      Free memory
    */
    free__MatrixLib__(dNdx);

  }

  /*
    Free memory
  */
  free_table__SetLib__(Internal_Element_Connectivity,NumInternalElements);
  free__MatrixLib__(Poligon);
  free__MatrixLib__(X_eval);

  return MinElementSize;
}

/*********************************************************************/

double volume__T6__(
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
    F_i = F_Ref__T6__(Xi,Element);
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


void local_search__T6__(Particle MPM_Mesh, Mesh FEM_Mesh)
{

  // Number of dimensions
  int Ndim = NumberDimensions;
  int Num_Particles_Node_i;
  int Num_Particles_Element_i;
  // Velocity and position of the particle
  Matrix Xi_p;
  Matrix X_p;
  Matrix V_p;
  // Previous closest node to the particle
  int I0_p_old; 
  // New closest node to the particle
  int I0_p_new;
  // List of nodes that interact with the particle
  ChainPtr Connectivity_p;

  ChainPtr Elem_p_Connectivity;
  // Coordinates of the nodes of the Element 
  Matrix Internal_Elem_p_Coordinates;
  Matrix Elem_p_Coordinates;
  // List of nodes close to the node I0_p 
  ChainPtr Locality_I0;

  // Set to zero the active/non-active node, and the GPs in each element */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++)
  {
    Num_Particles_Node_i = FEM_Mesh.Num_Particles_Node[i];

    if(Num_Particles_Node_i != 0)
    {
      FEM_Mesh.Num_Particles_Node[i] = 0;
      free__SetLib__(&FEM_Mesh.List_Particles_Node[i]);
    }

    FEM_Mesh.ActiveNode[i] = false;

  }

  for(int i = 0 ; i<FEM_Mesh.NumElemMesh ; i++)
  {
    Num_Particles_Element_i = FEM_Mesh.Num_Particles_Element[i];

    if(Num_Particles_Element_i != 0)
    {
      free__SetLib__(&FEM_Mesh.List_Particles_Element[i]); 
      FEM_Mesh.Num_Particles_Element[i] = 0;
      FEM_Mesh.Vol_patch_n[i] = 0.0;
      FEM_Mesh.Vol_patch_n1[i] = 0.0;
    }
  }

  // Loop over the particles
  for(int p = 0 ; p<MPM_Mesh.NumGP ; p++)
  {

    // Get the global coordinates and velocity of the particle
    X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
    V_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.vel.nM[p]);
    Xi_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_EC.nM[p]);

    // Check if the particle is static or is in movement
    if(norm__MatrixLib__(V_p,2) > 0)
    {

      // Get the index of the node close to the particle
      I0_p_old = MPM_Mesh.I0[p];

      // Get nodes close to the node I0_p
      Locality_I0 = FEM_Mesh.NodalLocality_0[I0_p_old];

      // Compute node close to the particle
      I0_p_new = get_closest_node__MeshTools__(X_p,Locality_I0,FEM_Mesh.Coordinates);

      // Update the tributary nodes of each particle
      MPM_Mesh.Element_p[p] = search_particle_in_surrounding_elements__Particles__(p,X_p,FEM_Mesh.NodeNeighbour[MPM_Mesh.I0[p]],FEM_Mesh);
  
      if(MPM_Mesh.Element_p[p] == -999)
      {
        fprintf(stderr,"%s : %s %i \n",
        "Error in local_search__T6__ -> search_particle_in_surrounding_elements__Particles__",
        "Not posible to find the particle",p);
        exit(EXIT_FAILURE);
      }

      // Free previous connectivity
      free__SetLib__(&MPM_Mesh.ListNodes[p]);
      MPM_Mesh.ListNodes[p] = NULL;  

      // Get connectivity of the element 
      Elem_p_Connectivity = FEM_Mesh.Connectivity[MPM_Mesh.Element_p[p]];
      Elem_p_Coordinates = get_nodes_coordinates__MeshTools__(Elem_p_Connectivity, FEM_Mesh.Coordinates);

      // Asign new connectivity
      MPM_Mesh.ListNodes[p] = tributary__T6__(p,X_p,Elem_p_Connectivity,Elem_p_Coordinates);

      Internal_Elem_p_Coordinates = get_nodes_coordinates__MeshTools__(MPM_Mesh.ListNodes[p], FEM_Mesh.Coordinates);

      /* Compute local coordinates of the particle in this element */
      X_to_Xi__T6__(Xi_p,X_p,Internal_Elem_p_Coordinates);

      /* Free coordinates of the element */
      free__MatrixLib__(Elem_p_Coordinates);
      free__MatrixLib__(Internal_Elem_p_Coordinates);

      // Activate the nodes near the particle
      Connectivity_p = MPM_Mesh.ListNodes[p];
      while(Connectivity_p != NULL)
      {
        if(FEM_Mesh.ActiveNode[Connectivity_p->I] == false)
        {
          FEM_Mesh.ActiveNode[Connectivity_p->I] = true;
        }

        Connectivity_p = Connectivity_p->next;
      }

      // Active those nodes that interact with the particle
      asign_to_nodes__Particles__(p, MPM_Mesh.Element_p[p], MPM_Mesh.I0[p], MPM_Mesh.ListNodes[p], FEM_Mesh);
      
   }
    else
    {
      // Activate the nodes near the particle
      Connectivity_p = MPM_Mesh.ListNodes[p];
      while(Connectivity_p != NULL)
      {
        if(FEM_Mesh.ActiveNode[Connectivity_p->I] == false)
        {
          FEM_Mesh.ActiveNode[Connectivity_p->I] = true;
        }

        Connectivity_p = Connectivity_p->next;
      }

      // Active those nodes that interact with the particle
      asign_to_nodes__Particles__(p, MPM_Mesh.Element_p[p], MPM_Mesh.I0[p], MPM_Mesh.ListNodes[p], FEM_Mesh);
    }

  }


}

/*********************************************************************/

static ChainPtr tributary__T6__(
  int p,
  Matrix X,
  ChainPtr Elem_Connectivity,
  Matrix Elem_Coordinates)
{

  /* Define output */
  ChainPtr Triburary_Nodes = NULL;

  Matrix Xi = allocZ__MatrixLib__(2,1);
  Matrix Macro_Element = allocZ__MatrixLib__(3,2);

  Macro_Element.nM[0][0] = Elem_Coordinates.nM[5][0];
  Macro_Element.nM[0][1] = Elem_Coordinates.nM[5][1];

  Macro_Element.nM[1][0] = Elem_Coordinates.nM[4][0];
  Macro_Element.nM[1][1] = Elem_Coordinates.nM[4][1];
 
  Macro_Element.nM[2][0] = Elem_Coordinates.nM[3][0];
  Macro_Element.nM[2][1] = Elem_Coordinates.nM[3][1];

  X_to_Xi__T6__(Xi, X, Macro_Element);

  int * List_Nodes = set_to_memory__SetLib__(Elem_Connectivity, 6);

  if((Xi.nV[0]  >= 0.0) && 
      (Xi.nV[1] >= 0.0) && 
      (Xi.nV[1] + Xi.nV[0] - 0.5 <=  0.0))
  {
      push__SetLib__ (&Triburary_Nodes, List_Nodes[5]);
      push__SetLib__ (&Triburary_Nodes, List_Nodes[2]);
      push__SetLib__ (&Triburary_Nodes, List_Nodes[0]);
  }
  else if((Xi.nV[0]  >= 0.5) && 
    (Xi.nV[1] >= 0.0) && 
    (Xi.nV[1] + Xi.nV[0] - 1.0 <=  0.0))
  {
      push__SetLib__ (&Triburary_Nodes, List_Nodes[4]);
      push__SetLib__ (&Triburary_Nodes, List_Nodes[2]);
      push__SetLib__ (&Triburary_Nodes, List_Nodes[1]);
  }
  else if((Xi.nV[0]  >= 0.0) && 
    (Xi.nV[1] >= 0.5) && 
    (Xi.nV[1] + Xi.nV[0] - 1.0 <=  0.0))
  {
      push__SetLib__ (&Triburary_Nodes, List_Nodes[3]);
      push__SetLib__ (&Triburary_Nodes, List_Nodes[1]);
      push__SetLib__ (&Triburary_Nodes, List_Nodes[0]);
  }
  else if((Xi.nV[0]  <= 0.5) && 
    (Xi.nV[1] <= 0.5) && 
    (Xi.nV[1] + Xi.nV[0] - 0.5 >=  0.0))
  {
      push__SetLib__ (&Triburary_Nodes, List_Nodes[2]);
      push__SetLib__ (&Triburary_Nodes, List_Nodes[1]);
      push__SetLib__ (&Triburary_Nodes, List_Nodes[0]);
  }
  else
  {
    printf(" %s : %s \n",
     "Error in tributary__T6__",
     "Particle outside of the element");
    exit(EXIT_FAILURE);
  }

  /*
    Free memory
  */
  free(List_Nodes);
  free__MatrixLib__(Xi);
  free__MatrixLib__(Macro_Element);

  return Triburary_Nodes;

}

/*********************************************************************/
