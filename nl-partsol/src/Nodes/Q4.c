#include "nl-partsol.h"


/* 
  Call global variables
*/
double Thickness_Plain_Stress;

/*
  Auxiliar functions
*/
static Matrix F_Ref__Q4__(Matrix,Matrix);
static Matrix Xi_to_X__Q4__(Matrix,Matrix);
static void   X_to_Xi__Q4__(Matrix,Matrix,Matrix);

/*********************************************************************/

void initialize__Q4__(
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
 
  ChainPtr Locality_I0; // List of nodes close to the node I0_p
 
  /* Loop over the particles to initialize them */
  for(int p = 0 ; p<Np ; p++)
  {

    /* Supose that the particle was not initilise */
    Init_p = false;

    /* Asign the number of nodes */
    MPM_Mesh.NumberNodes[p] = 4;
    
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
      if(in_out__Q4__(X_p,Elem_p_Coordinates) == true)
      {

        /* Particle will be initilise */
        Init_p = true;

        // Assign the index of the element
        MPM_Mesh.Element_p[p] = i;

        /* With the element connectivity get the node close to the particle */
        MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,Elem_p_Connectivity,FEM_Mesh.Coordinates);

        /* Asign connectivity */
        MPM_Mesh.ListNodes[p] = copy__SetLib__(Elem_p_Connectivity);

        /* Active those nodes that interact with the particle */
        asign_to_nodes__Particles__(p, MPM_Mesh.Element_p[p], MPM_Mesh.I0[p], MPM_Mesh.ListNodes[p], FEM_Mesh);

        /* Compute local coordinates of the particle in this element */
        X_to_Xi__Q4__(Xi_p,X_p,Elem_p_Coordinates);

        /* Free coordinates of the element */
        free__MatrixLib__(Elem_p_Coordinates);

        /* 
          Activate the nodes near the particle
        */
        Locality_I0 = FEM_Mesh.NodalLocality_0[MPM_Mesh.I0[p]];

        while(Locality_I0 != NULL)
        {
          if(FEM_Mesh.ActiveNode[Locality_I0->Idx] == false)
          {
            FEM_Mesh.ActiveNode[Locality_I0->Idx] = true;
          }

          Locality_I0 = Locality_I0->next;
        }


        break;
	
      }

      /* Free coordinates of the element */
      free__MatrixLib__(Elem_p_Coordinates);
    }

    if(!Init_p)
    {
      fprintf(stderr,"%s : %s %i\n",
        "Error in initialize__Q4__()",
        "The search algorithm was unable to find particle",p);
      exit(EXIT_FAILURE);
    }

  }
  
}

/*********************************************************************/

/* Shape functions */
Matrix N__Q4__(
  Matrix X_e)
{
  
  /* Definition and allocation */
  Matrix N_ref =  allocZ__MatrixLib__(1,4);

  /* Fill the array */
  N_ref.nV[0] = 0.25*(1-X_e.nV[0])*(1-X_e.nV[1]);
  N_ref.nV[1] = 0.25*(1+X_e.nV[0])*(1-X_e.nV[1]);
  N_ref.nV[2] = 0.25*(1+X_e.nV[0])*(1+X_e.nV[1]);
  N_ref.nV[3] = 0.25*(1-X_e.nV[0])*(1+X_e.nV[1]);
  
  return N_ref;
}

/*********************************************************************/

/* Derivatives of the shape functions */
Matrix dN_Ref__Q4__(
  Matrix X_e)
{

  int Ndim = NumberDimensions;
  
  /* Definition and allocation */
  Matrix dNdX_ref = allocZ__MatrixLib__(4,Ndim);
  
  /* Fill the matrix */

  /* Node 1 */
  dNdX_ref.nM[0][0] = -0.25*(1-X_e.nV[1]);
  dNdX_ref.nM[0][1] = -0.25*(1-X_e.nV[0]); 

  /* Node 2 */
  dNdX_ref.nM[1][0] = +0.25*(1-X_e.nV[1]);
  dNdX_ref.nM[1][1] = -0.25*(1+X_e.nV[0]);
  
  /* Node 3 */
  dNdX_ref.nM[2][0] = +0.25*(1+X_e.nV[1]);
  dNdX_ref.nM[2][1] = +0.25*(1+X_e.nV[0]);
  
  /* Node 4 */
  dNdX_ref.nM[3][0] = -0.25*(1+X_e.nV[1]);
  dNdX_ref.nM[3][1] = +0.25*(1-X_e.nV[0]); 
  
  return dNdX_ref;
}

/*********************************************************************/

/* Jacobian of the transformation for the four-nodes quadrilateral */
static Matrix F_Ref__Q4__(
  Matrix X_NC_GP, 
  Matrix X_GC_Nodes)
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
  for(int I = 0 ; I<4 ; I++)
  {

    /* 3º Fill arrays for the tensorial product */
    X_I = memory_to_tensor__TensorLib__(X_GC_Nodes.nM[I],1);
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
Matrix dN__Q4__(
  Matrix Xi,
  Matrix Element)
/*
  - Matrix X_EC_GP : Element coordinates of the gauss point
  - Matrix Element : Coordinates of the element (4 x Ndim)
*/
{
  
  /* Derivatives of the shape function evaluates in the GP (4 x Ndim) */
  Matrix dNdX;
  Matrix dNdX_T;
    
  /* 1º Evaluate the gradient of the shape function in the GP (4 x Ndim)
  and get the transpose */
  Matrix dNdX_Ref   = dN_Ref__Q4__(Xi);
  Matrix dNdX_Ref_T = transpose__MatrixLib__(dNdX_Ref);

  /* 2º Get the Jacobian of the transformation evaluated in the GP */
  Matrix F     = F_Ref__Q4__(Xi,Element);
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
Matrix Xi_to_X__Q4__(
  Matrix Xi,
  Matrix Element)
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

static void X_to_Xi__Q4__(
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
  Xi = Newton_Rapson(Xi_to_X__Q4__, Element, F_Ref__Q4__, Element, X, Xi);  
}

/*********************************************************************/

bool in_out__Q4__(
  Matrix X,
  Matrix Element)
{

  double min[2] = {Element.nM[0][0],Element.nM[0][1]};
  double max[2] = {Element.nM[0][0],Element.nM[0][1]};

  for(int a = 1 ; a<4 ; a++)
  {
    for(int i = 0 ; i<2; i++)
    {
      min[i] = DMIN(min[i],Element.nM[a][i]);
      max[i] = DMAX(max[i],Element.nM[a][i]);
    }
  }

  Matrix Xi;

  // Check if it is inside
  if((X.nV[0] <= max[0]) && 
    (X.nV[0] >= min[0]) && 
    (X.nV[1] <= max[1]) && 
    (X.nV[1] >= min[1]))
  {

    Xi = allocZ__MatrixLib__(2,1);

    X_to_Xi__Q4__(Xi, X, Element);

    if((fabs(Xi.nV[0]) <= 1.0) && (fabs(Xi.nV[1]) <= 1.0))
    {
      free__MatrixLib__(Xi);
      return true;    
    }
    else
    {
      free__MatrixLib__(Xi);
      return false;
    }

  }
  else
  {
    return false;
  }

}

/*********************************************************************/

void element_to_particles__Q4__(
  Matrix X_p,
  Mesh FEM_Mesh,
  int GPxElement)
{

  int Ndim = NumberDimensions;
  int NumElemMesh = FEM_Mesh.NumElemMesh;
  int NumNodesElem = 4;
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

    case 4:
    Xi_p.nM[0][0] =   1./sqrt(3.0);
    Xi_p.nM[0][1] =   1./sqrt(3.0);
    Xi_p.nM[1][0] =   1./sqrt(3.0);
    Xi_p.nM[1][1] = - 1./sqrt(3.0);
    Xi_p.nM[2][0] = - 1./sqrt(3.0);
    Xi_p.nM[2][1] =   1./sqrt(3.0);
    Xi_p.nM[3][0] = - 1./sqrt(3.0);
    Xi_p.nM[3][1] = - 1./sqrt(3.0);
    break;

    case 9:
    Xi_p.nM[0][0] =   0.0;
    Xi_p.nM[0][1] =   0.0;

    Xi_p.nM[1][0] =   0.6666666666666;
    Xi_p.nM[1][1] =   0.0;

    Xi_p.nM[2][0] =   0.6666666666666;
    Xi_p.nM[2][1] =   0.6666666666666;

    Xi_p.nM[3][0] =   0.0;
    Xi_p.nM[3][1] =   0.6666666666666;

    Xi_p.nM[4][0] = - 0.6666666666666;
    Xi_p.nM[4][1] =   0.6666666666666;

    Xi_p.nM[5][0] = - 0.6666666666666;
    Xi_p.nM[5][1] =   0.0;

    Xi_p.nM[6][0] = - 0.6666666666666;
    Xi_p.nM[6][1] = - 0.6666666666666;

    Xi_p.nM[7][0] =   0.0;
    Xi_p.nM[7][1] = - 0.6666666666666;

    Xi_p.nM[8][0] =   0.6666666666666;
    Xi_p.nM[8][1] = - 0.6666666666666;
    break;

    default :
    fprintf(stderr,"%s : %s \n","Error in element_to_particles__Q4__()","Wrong number of particles per element");
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
        N_GP = N__Q4__(Xi_p);
      }
      else
      {
        Xi_p_j.nV = Xi_p.nM[j]; 
        N_GP = N__Q4__(Xi_p_j);
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

double min_DeltaX__Q4__(ChainPtr Element_Connectivity, Matrix Coordinates)
{

  /* Auxiliar variables of the function */
  int Ndim = NumberDimensions;
  int Node_k;
  int Node_k1;
  int NumNodesElem = 4;
  int NumSides = 4;
  double sqr_lenght = 0.0;
  double MinElementSize = 10e16;

  /*
    Fill the poligon with vectors
  */
  for(int k = 0; k<NumNodesElem; k++)
  {

    Node_k = Element_Connectivity->Idx;
    Node_k1 = Element_Connectivity->next->Idx;

    sqr_lenght = 0.0;
    
    for(int l = 0 ; l<Ndim ; l++)
    {
      sqr_lenght += DSQR(Coordinates.nM[Node_k1][l] - Coordinates.nM[Node_k][l]);
    }

    MinElementSize = DMIN(MinElementSize,sqrt(sqr_lenght));

    Element_Connectivity = Element_Connectivity->next;

  }

  return MinElementSize;
}

/*********************************************************************/

double volume__Q4__(
  Matrix Element)
{
  int Ndim = NumberDimensions;
  double J_i;
  double Vol = 0;

  // Use 4 integration points to compute volume
  double table_w[4] = {1.,1.,1.,1.};
  double table_X[4][2] = 
  {
    {-0.577350269200000,-0.577350269200000},
    {0.577350269200000,-0.577350269200000},
    {-0.577350269200000,0.577350269200000},
    {0.577350269200000,0.577350269200000}
  };

  Matrix F_i;
  Matrix Xi = allocZ__MatrixLib__(2,1);

  for(int i = 0 ; i<4 ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
       Xi.nV[j] = table_X[i][j];
    }

    // Compute deformation gradient and jacobian of this integration point
    F_i = F_Ref__Q4__(Xi,Element);
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


void local_search__Q4__(Particle MPM_Mesh, Mesh FEM_Mesh)
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
  // Coordinates of the nodes of the Element 
  Matrix CoordElement;
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

  for(int i = 0 ; i<FEM_Mesh.Num_Patch_Mesh ; i++)
  {
//    Num_Particles_Element_i = FEM_Mesh.Num_Particles_Element[i];

//    if(Num_Particles_Element_i != 0)
//    {
//      free__SetLib__(&FEM_Mesh.List_Particles_Element[i]); 
//      FEM_Mesh.Num_Particles_Element[i] = 0;
      FEM_Mesh.Vol_Patch_n[i] = 0.0;
      FEM_Mesh.Vol_Patch_n1[i] = 0.0;
//    }
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

      MPM_Mesh.I0[p] = I0_p_new;

      // Update the tributary nodes of each particle
      MPM_Mesh.Element_p[p] = search_particle_in_surrounding_elements__Particles__(p,X_p,FEM_Mesh.NodeNeighbour[I0_p_new],FEM_Mesh);
  
      if(MPM_Mesh.Element_p[p] == -999)
      {
        fprintf(stderr,"%s : %s %i \n",
        "Error in local_search__Q4__ -> search_particle_in_surrounding_elements__Particles__",
        "Not posible to find the particle",p);
        exit(EXIT_FAILURE);
      }

      // Free previous connectivity
      free__SetLib__(&MPM_Mesh.ListNodes[p]);
      MPM_Mesh.ListNodes[p] = NULL;  
        
      // Asign new connectivity
      MPM_Mesh.ListNodes[p] = copy__SetLib__(FEM_Mesh.Connectivity[MPM_Mesh.Element_p[p]]);

      // Get the coordinates of the element vertex
      CoordElement = get_nodes_coordinates__MeshTools__(MPM_Mesh.ListNodes[p],FEM_Mesh.Coordinates);

      // Compute local coordinates of the particle in this element
      X_to_Xi__Q4__(Xi_p,X_p,CoordElement);

      // Free coordinates of the element
      free__MatrixLib__(CoordElement);

      // Activate the nodes near the particle
      Connectivity_p = MPM_Mesh.ListNodes[p];
      while(Connectivity_p != NULL)
      {
        if(FEM_Mesh.ActiveNode[Connectivity_p->Idx] == false)
        {
          FEM_Mesh.ActiveNode[Connectivity_p->Idx] = true;
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
        if(FEM_Mesh.ActiveNode[Connectivity_p->Idx] == false)
        {
          FEM_Mesh.ActiveNode[Connectivity_p->Idx] = true;
        }

        Connectivity_p = Connectivity_p->next;
      }

      // Active those nodes that interact with the particle
      asign_to_nodes__Particles__(p, MPM_Mesh.Element_p[p], MPM_Mesh.I0[p], MPM_Mesh.ListNodes[p], FEM_Mesh);
    }

  }

}

/*********************************************************************/
