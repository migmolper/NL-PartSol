#include "nl-partsol.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))


/*********************************************************************/

void initial_position__Particles__(Matrix X_p, Mesh FEM_Mesh, int GPxElement)
/*
 * 
 */
{

  if(strcmp(FEM_Mesh.TypeElem,"Triangle") == 0)
  {
    element_to_particles__T3__(X_p, FEM_Mesh, GPxElement);
  }
  else if(strcmp(FEM_Mesh.TypeElem,"Quadrilateral") == 0)
  {
    element_to_particles__Q4__(X_p, FEM_Mesh,GPxElement);
  }
  else if(strcmp(FEM_Mesh.TypeElem,"Tetrahedra") == 0)
  {
    element_to_particles__T4__(X_p, FEM_Mesh,GPxElement);
  }
  else if(strcmp(FEM_Mesh.TypeElem,"Hexahedra") == 0)
  {
    element_to_particles__H8__(X_p, FEM_Mesh,GPxElement);
  }
  else
  {
    fprintf(stderr,"%s : %s \n",
      "Error in initial_position__Particles__()",
      "Wrong type of shape function");
    exit(EXIT_FAILURE);
  }
  
}

/*********************************************************************/

int search_particle_in_surrounding_elements__Particles__(
  int p,
  Matrix X_p,
  ChainPtr ListElement,
  Mesh FEM_Mesh)
/*

*/
{
  ChainPtr Ixd = NULL;
  int I_element = -999;
  int Nn; /* Number of nodes of the element */
  ChainPtr Nodes;
  Matrix Element_Coordinates;

  Ixd = ListElement;
  
  while(Ixd != NULL)
  {

    Nn = FEM_Mesh.NumNodesElem[Ixd->I];
    Nodes = FEM_Mesh.Connectivity[Ixd->I];

    Element_Coordinates = get_nodes_coordinates__MeshTools__(Nodes, FEM_Mesh.Coordinates);

    /* Check if the particle is in the element */
    if(FEM_Mesh.In_Out_Element(X_p,Element_Coordinates))
    {
      free__MatrixLib__(Element_Coordinates);
      I_element = Ixd->I;
      break; 
    }
    
    /* Cycle */
    free__MatrixLib__(Element_Coordinates);
    Ixd = Ixd->next;

  }
  
  return I_element;
}

/*********************************************************************/

void get_particle_tributary_nodes(Particle MPM_Mesh, Mesh FEM_Mesh, int p)
{

  int Ndim = NumberDimensions;

  /* Indetify the node close to the particle */
  int I0 = MPM_Mesh.I0[p];

  /* Define auxiliar matrix for local/global coordinates */
  Matrix X_p  = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
  Matrix Xi_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_EC.nM[p]);
  /* Lis of elements near the particle */
  ChainPtr Elements_Near_I0 = FEM_Mesh.NodeNeighbour[I0];
  /* Index of the element (Q4/uGIMP) */
  int IdxElement;
  /* Coordinates of the nodes of the Element */
  Matrix CoordElement;
  
  /* 6º Assign the new connectivity of the GP */
  if(strcmp(ShapeFunctionGP,"FEM") == 0)
  {

    /* Get the index of the element */
    IdxElement = search_particle_in_surrounding_elements__Particles__(p,X_p,Elements_Near_I0,FEM_Mesh);
    if(IdxElement != -999)
    {
      /* Free previous list of tributary nodes to the particle */
      free__SetLib__(&MPM_Mesh.ListNodes[p]);
      MPM_Mesh.ListNodes[p] = NULL;    
      
      /* Asign new connectivity */
      MPM_Mesh.ListNodes[p] = copy__SetLib__(FEM_Mesh.Connectivity[IdxElement]);
      
      /* Get the coordinates of the element vertex */
      CoordElement = get_nodes_coordinates__MeshTools__(MPM_Mesh.ListNodes[p],FEM_Mesh.Coordinates);
      
      /* Compute local coordinates of the particle in this element */
      FEM_Mesh.X_to_Xi(Xi_p,X_p,CoordElement);

      /* Free coordinates of the element */
      free__MatrixLib__(CoordElement);
    }
    else
    {
      /* Get the coordinates of the element vertex */
      CoordElement = get_nodes_coordinates__MeshTools__(MPM_Mesh.ListNodes[p],FEM_Mesh.Coordinates);
      /* Compute local coordinates of the particle in this element */
      FEM_Mesh.X_to_Xi(Xi_p,X_p,CoordElement);
      /* Free coordinates of the element */
      free__MatrixLib__(CoordElement);      
    }
    
  }

  else if(strcmp(ShapeFunctionGP,"uGIMP") == 0)
  {

    /* Auxiliar variables for GIMP */
    Matrix lp = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.lp.nM[p]);
    
    /* Get the index of the element */
    IdxElement = search_particle_in_surrounding_elements__Particles__(p,X_p,Elements_Near_I0,FEM_Mesh);

    if(IdxElement != -999)
    {
      /* Get the coordinates of the element vertex */
      CoordElement = get_nodes_coordinates__MeshTools__(FEM_Mesh.Connectivity[IdxElement],
							FEM_Mesh.Coordinates);
      /* Compute local coordinates of the particle in this element */
      X_to_Xi__Q4__(Xi_p,X_p,CoordElement);
      /* Free coordinates of the element */
      free__MatrixLib__(CoordElement);   
      /* Free previous list of tributary nodes to the particle */
      free__SetLib__(&MPM_Mesh.ListNodes[p]);
      /* Asign new connectivity */
      MPM_Mesh.ListNodes[p] = tributary__GIMP__(Xi_p,IdxElement,lp,FEM_Mesh);  
      /* Calculate number of nodes */
      MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);
    }
    
  }

  else if(strcmp(ShapeFunctionGP,"LME") == 0)
  {
    /* 
      Auxiliar variables for LME
    */
    Matrix Metric_p; // Define a metric tensor
    Matrix Delta_Xip; // Distance from particles to the nodes
    Matrix lambda_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.lambda.nM[p]);
    Tensor F_p; // Particle deformation gradient
    double Beta_p = MPM_Mesh.Beta.nV[p]; // Thermalization parameter

    /*
      Compute the metric tensor
    */
    F_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    Metric_p = metric__LME__(F_p);

    /*
      Free previous list of tributary nodes to the particle
    */
    free__SetLib__(&MPM_Mesh.ListNodes[p]);
    MPM_Mesh.ListNodes[p] = NULL;


    /*
      Calculate the new connectivity with the previous value of beta
    */
    MPM_Mesh.ListNodes[p] = tributary__LME__(p,X_p,Metric_p,Beta_p,I0,FEM_Mesh);

    /*
      Calculate number of nodes
    */
    MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);

    /* 
      Generate nodal distance list
    */
    Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p], X_p, FEM_Mesh.Coordinates);    	      

    /*
      Compute the thermalization parameter for the new set of nodes
      and update it
    */
    Beta_p = beta__LME__(gamma_LME, FEM_Mesh.h_avg[MPM_Mesh.I0[p]]);
    MPM_Mesh.Beta.nV[p] = Beta_p;

    /*
      Compute the lagrange multiplier of the new shape functions
    */
    MPM_Mesh.update_lambda(p, Delta_Xip, lambda_p, Metric_p, Beta_p);
    
    /* Free memory */
    free__MatrixLib__(Metric_p);
    free__MatrixLib__(Delta_Xip);
  }

}

/*********************************************************************/

Element nodal_set__Particles__(int i_GP, ChainPtr ListNodes, int NumNodes)
{

  /* Define new element */
  Element GP_Element;

  /* Fill element */
  GP_Element.i_GP = i_GP;
  GP_Element.NumberNodes = NumNodes;
  GP_Element.Connectivity = set_to_memory__SetLib__(ListNodes,NumNodes);

  return GP_Element;
}

/*********************************************************************/

void asign_to_nodes__Particles__(int p, int E_p, int I0, ChainPtr ListNodes_p, Mesh FEM_Mesh)
{

  /*!
   * Assign particle to an element of the background mesh
   * */
  push__SetLib__(&FEM_Mesh.List_Particles_Element[E_p],p);
  FEM_Mesh.Num_Particles_Element[E_p] += 1;
  
  /*!
   * Assign particle to a node of the background mesh
   * */
  push__SetLib__(&FEM_Mesh.List_Particles_Node[I0],p);
  FEM_Mesh.Num_Particles_Node[I0] += 1;
  
}

/*********************************************************************/

