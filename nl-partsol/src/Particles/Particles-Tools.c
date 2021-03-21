#include "nl-partsol.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

/*********************************************************************/

void initial_position__Particles__(Matrix X_p, Mesh FEM_Mesh, int GPxElement)
/*
 * 
 */
{
  int Ndim = NumberDimensions;
  int NumElemMesh = FEM_Mesh.NumElemMesh;
  Matrix N_GP;
  Matrix Xi_p = allocZ__MatrixLib__(GPxElement,Ndim);
  Matrix Xi_p_j;
  Element Element;
  int Node;

  switch(GPxElement){
  case 1 :
    if(strcmp(FEM_Mesh.TypeElem,"Quadrilateral") == 0){
      /* Centred GP */
      Xi_p.nV[0] = 0.0;
      Xi_p.nV[1] = 0.0;
      /* Evaluate the shape function */
      N_GP = N__Q4__(Xi_p);
      /* Get the coordinate of the center */
      for(int i = 0 ; i<NumElemMesh ; i++){
	Element = nodal_set__Particles__(i,FEM_Mesh.Connectivity[i],
					 FEM_Mesh.NumNodesElem[i]);
	
	for(int k = 0 ; k<4 ; k++){
	  Node = Element.Connectivity[k];
	  for(int l = 0 ; l<2 ; l++){
	    X_p.nM[i][l] += N_GP.nV[k]*FEM_Mesh.Coordinates.nM[Node][l];
	  }
	}
	
	free(Element.Connectivity);
      }
      /* Free auxiliar matrix with the coordinates */
      free__MatrixLib__(Xi_p);
      /* Free value of the shape function in the GP */
      free__MatrixLib__(N_GP);
    }
    else if(strcmp(FEM_Mesh.TypeElem,"Triangle") == 0){
      Xi_p.nV[0] = (double)1/3;
      Xi_p.nV[1] = (double)1/3;
      /* Evaluate the shape function */
      N_GP = N__T3__(Xi_p);
      /* Get the coordinate of the center */
      for(int i = 0 ; i<NumElemMesh ; i++){
	
	Element = nodal_set__Particles__(i, FEM_Mesh.Connectivity[i],
					 FEM_Mesh.NumNodesElem[i]);
	for(int k = 0 ; k<3 ; k++){
	  Node = Element.Connectivity[k];
	  for(int l = 0 ; l<2 ; l++){
	    X_p.nM[i][l] += N_GP.nV[k]*FEM_Mesh.Coordinates.nM[Node][l];
	  }
	}
	
	free(Element.Connectivity);
      }
      /* Free auxiliar matrix with the coordinates */
      free__MatrixLib__(Xi_p);
      /* Free value of the shape function in the GP */
      free__MatrixLib__(N_GP);
    }
    break;

  case 3:
    if(strcmp(FEM_Mesh.TypeElem,"Triangle") == 0){
      /* Centred GP */
      Xi_p.nM[0][0] =  0.16666666666;
      Xi_p.nM[0][1] =  0.16666666666;
      Xi_p.nM[1][0] =  0.66666666666;
      Xi_p.nM[1][1] =  0.16666666666;
      Xi_p.nM[2][0] =  0.16666666666;
      Xi_p.nM[2][1] =  0.66666666666;
      /* Get the coordinate of the center */
      for(int i = 0 ; i<NumElemMesh ; i++){
	
	Element = nodal_set__Particles__(i, FEM_Mesh.Connectivity[i],
					 FEM_Mesh.NumNodesElem[i]);
	for(int j = 0 ; j<GPxElement ; j++){
	  /* Evaluate the shape function in the GP position */
	  Xi_p_j.nV = Xi_p.nM[j]; 
	  N_GP = N__T3__(Xi_p_j);
	  for(int k = 0 ; k<3 ; k++){
	    /* Connectivity of each element */
	    Node = Element.Connectivity[k];
	    for(int l = 0 ; l<Ndim ; l++){
	      X_p.nM[i*GPxElement+j][l] += 
		N_GP.nV[k]*FEM_Mesh.Coordinates.nM[Node][l];
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
    break;
    
  case 4:
    if(strcmp(FEM_Mesh.TypeElem,"Quadrilateral") == 0){
      /* Xi_p.nM[0][0] = (double)1/pow(3,0.5); */
      /* Xi_p.nM[0][1] = (double)1/pow(3,0.5); */
      /* Xi_p.nM[1][0] = (double)1/pow(3,0.5); */
      /* Xi_p.nM[1][1] = (double)-1/pow(3,0.5); */
      /* Xi_p.nM[2][0] = (double)-1/pow(3,0.5); */
      /* Xi_p.nM[2][1] = (double)1/pow(3,0.5); */
      /* Xi_p.nM[3][0] = (double)-1/pow(3,0.5); */
      /* Xi_p.nM[3][1] = (double)-1/pow(3,0.5); */
      Xi_p.nM[0][0] =   0.5;
      Xi_p.nM[0][1] =   0.5;
      Xi_p.nM[1][0] =   0.5;
      Xi_p.nM[1][1] = - 0.5;
      Xi_p.nM[2][0] = - 0.5;
      Xi_p.nM[2][1] =   0.5;
      Xi_p.nM[3][0] = - 0.5;
      Xi_p.nM[3][1] = - 0.5;
      /* Get the coordinate of the center */
      for(int i = 0 ; i<NumElemMesh ; i++){
	
	Element = nodal_set__Particles__(i, FEM_Mesh.Connectivity[i],
					 FEM_Mesh.NumNodesElem[i]);
	for(int j = 0 ; j<GPxElement ; j++){
	  /* Evaluate the shape function in the GP position */
	  Xi_p_j.nV = Xi_p.nM[j]; 
	  N_GP = N__Q4__(Xi_p_j);
	  for(int k = 0 ; k<4 ; k++){
	    /* Connectivity of each element */
	    Node = Element.Connectivity[k];
	    for(int l = 0 ; l<Ndim ; l++){
	      X_p.nM[i*GPxElement+j][l] += 
		N_GP.nV[k]*FEM_Mesh.Coordinates.nM[Node][l];
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
    break;
  default :
    fprintf(stderr,"%s : %s \n",
	    "Error in initial_position__Particles__()",
	    "Wrong number of gauss point per element");
    exit(EXIT_FAILURE);
  }
  
}

/*********************************************************************/

int inout_element__Particles__(int p, Matrix X_p, ChainPtr ListElement, Mesh FEM_Mesh)
/*

*/
{
  ChainPtr Ixd = NULL;
  int I_element = -999;
  int Nn; /* Number of nodes of the element */
  ChainPtr Nodes;

  Ixd = ListElement;
  
  while(Ixd!=NULL){

    Nn = FEM_Mesh.NumNodesElem[Ixd->I];
    Nodes = FEM_Mesh.Connectivity[Ixd->I];

    /* Check if the particle is in the element */
    if(inout_convex_set__MeshTools__(X_p, Nodes, FEM_Mesh.Coordinates))
      {
	      I_element = Ixd->I;
	      break;
      }
    
    /* Cycle */
    Ixd = Ixd->next;

  }

  if(I_element == -999)
    {
      fprintf(stderr,"%s : %s %i \n",
	      "Error in inout_element__Particles__()",
	      "Not posible to find the particle",p);
    }
  
  
  return I_element;
}

/*********************************************************************/

void get_particle_tributary_nodes(GaussPoint MPM_Mesh, Mesh FEM_Mesh, int p){

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
  if(strcmp(ShapeFunctionGP,"MPMQ4") == 0){

    /* Get the index of the element */
    IdxElement = inout_element__Particles__(p,X_p,Elements_Near_I0,FEM_Mesh);
    if(IdxElement != -999){
      /* Free previous list of tributary nodes to the particle */
      free__SetLib__(&MPM_Mesh.ListNodes[p]);
      MPM_Mesh.ListNodes[p] = NULL;    
      /* Asign new connectivity */
      MPM_Mesh.ListNodes[p] = copy__SetLib__(FEM_Mesh.Connectivity[IdxElement]);
      /* Get the coordinates of the element vertex */
      CoordElement = get_nodes_coordinates__MeshTools__(MPM_Mesh.ListNodes[p],
							FEM_Mesh.Coordinates);
      /* Compute local coordinates of the particle in this element */
      X_to_Xi__Q4__(Xi_p,X_p,CoordElement);
      /* Free coordinates of the element */
      free__MatrixLib__(CoordElement);
    }
    else{
      /* Get the coordinates of the element vertex */
      CoordElement = get_nodes_coordinates__MeshTools__(MPM_Mesh.ListNodes[p],
							FEM_Mesh.Coordinates);
      /* Compute local coordinates of the particle in this element */
      X_to_Xi__Q4__(Xi_p,X_p,CoordElement);
      /* Free coordinates of the element */
      free__MatrixLib__(CoordElement);      
    }
    
  }
  else if(strcmp(ShapeFunctionGP,"uGIMP") == 0)
  {

    /* Auxiliar variables for GIMP */
    Matrix lp = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.lp.nM[p]);
    
    /* Get the index of the element */
    IdxElement = inout_element__Particles__(p,X_p,Elements_Near_I0,FEM_Mesh);

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
    Matrix Metric_p = metric__LME__(); // Define a metric tensor
    Matrix Delta_Xip; // Distance from particles to the nodes
    Matrix lambda_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.lambda.nM[p]);
    double Beta_p = MPM_Mesh.Beta.nV[p]; // Thermalization parameter

    /*
      Free previous list of tributary nodes to the particle
    */
    free__SetLib__(&MPM_Mesh.ListNodes[p]);
    MPM_Mesh.ListNodes[p] = NULL;

    /*
      Calculate the new connectivity with the previous value of beta
    */
    MPM_Mesh.ListNodes[p] = tributary__LME__(X_p,Beta_p,I0,FEM_Mesh);

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
    Beta_p = beta__LME__(Delta_Xip, gamma_LME);
    MPM_Mesh.Beta.nV[p] = Beta_p;

    /*
      Compute the lagrange multiplier of the new shape functions
    */
    lambda_p = lambda__LME__(Delta_Xip, lambda_p, Metric_p, Beta_p);
    
    /* Free memory */
    free__MatrixLib__(Metric_p);
    free__MatrixLib__(Delta_Xip);
  }

}


/*********************************************************************/

void local_search__Particles__(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
/*
  Search the closest node to the particle based in its previous position.
*/
{

  /* Number of dimensions */
  int Ndim = NumberDimensions;
  /* Velocity and position of the particle */
  Matrix X_p;
  Matrix V_p;
  /* Previous closest node to the particle */
  int I0_p;
  /* List of nodes close to the node I0_p */
  ChainPtr Locality_I0;

  /* Set to zero the active/non-active node, and the GPs in each element */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.NumParticles[i] = 0;
  }
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    free__SetLib__(&FEM_Mesh.I_particles[i]);
  }

  /* Loop over the particles */
  for(int p = 0 ; p<MPM_Mesh.NumGP ; p++){

    /* Get the global coordinates and velocity of the particle */
    X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
    V_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.vel.nM[p]);

    /* Check if the particle is static or is in movement */
    if(norm__MatrixLib__(V_p,2) > 0){

      /* Get the index of the node close to the particle */
      I0_p = MPM_Mesh.I0[p];

      /* Get nodes close to the node I0_p */
      Locality_I0 = FEM_Mesh.NodalLocality[I0_p];

      /* Update the index of the node close to the particle */
      MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,Locality_I0,FEM_Mesh.Coordinates);

      /* Update the tributary nodes of each particle */
      get_particle_tributary_nodes(MPM_Mesh,FEM_Mesh,p);

      /* Active those nodes that interact with the particle */
      asign_to_nodes__Particles__(p, MPM_Mesh.ListNodes[p], FEM_Mesh);
      
    }
    else{
      /* Active those nodes that interact with the particle */
      asign_to_nodes__Particles__(p, MPM_Mesh.ListNodes[p], FEM_Mesh);
    }

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

void asign_to_nodes__Particles__(int p, ChainPtr ListNodes_p, Mesh FEM_Mesh)
{
  
  /* Auxiliar variable to loop in the list of tributary nodes of the particle */
  ChainPtr Nodes_p = NULL;
  
  Nodes_p = ListNodes_p;
  while(Nodes_p != NULL)
    {
      FEM_Mesh.NumParticles[Nodes_p->I] += 1;
      push__SetLib__(&FEM_Mesh.I_particles[Nodes_p->I],p);
      Nodes_p = Nodes_p->next; 
    }

}

/*********************************************************************/

