#include "nl-partsol.h"

/*
  Auxiliar functions
 */
static Matrix   voxel__GIMP__(Matrix, double);
static double   Sip__GIMP__(double, double, double);
static double   dSip__GIMP__(double, double, double);

/*********************************************************************/


void initialize__GIMP__(GaussPoint MPM_Mesh, Mesh FEM_Mesh){

  int Ndim = NumberDimensions;

  /* Variables for the GP coordinates */  
  Matrix X_p;
  Matrix Xi_p;
  Matrix l_p;
  
  double A_p;
  double rho_p;
  double th_p;
  double m_p;
  int Mat_p;
  bool Init_p;

  /* Variable with stores the conectivity of the element of the particle */
  ChainPtr Elem_p;
  
  /* Variable with stores the coordinates of the element of the particle */
  Matrix Poligon_Coordinates;

  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
   
  /* Get some properties for each particle */ 
  X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[i]);
  Xi_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_EC.nM[i]);
  l_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.lp.nM[i]);
  rho_p = MPM_Mesh.Phi.rho.nV[i];
  m_p = MPM_Mesh.Phi.mass.nV[i];
  th_p = MPM_Mesh.Mat[MPM_Mesh.MatIdx[i]].thickness;
  A_p = m_p/(rho_p*th_p);

  /* Supose that the particle was not initilise */
  Init_p = false;

  for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){

    /* Get the element properties */
    Elem_p = FEM_Mesh.Connectivity[j];
      
    /* 5º Check out if the GP is in the Element */
    if(InOut_Element(X_p, Elem_p, FEM_Mesh.Coordinates)){

      /* Particle will be initilise */
      Init_p = true;
	
      /* With the element connectivity get the node close to the particle */
      MPM_Mesh.I0[i] = get_closest_node_to(X_p,Elem_p,FEM_Mesh.Coordinates);

      /* Get the coordinates of the element */
      Poligon_Coordinates = ElemCoordinates(FEM_Mesh.Connectivity[j],
					    FEM_Mesh.Coordinates);

      /* If the GP is in the element, get its natural coordinates */
      Xi_p.nV = MPM_Mesh.Phi.x_EC.nM[i];
      X_to_Xi__Q4__(Xi_p,X_p,Poligon_Coordinates);

      /* Initialize the voxel lenght */
      l_p = voxel__GIMP__(l_p, A_p);

      /* Get the initial connectivity of the particle */
      MPM_Mesh.ListNodes[i] = tributary__GIMP__(Xi_p,j,l_p,FEM_Mesh);

      /* Measure the size of the connectivity */
      MPM_Mesh.NumberNodes[i] = lenght__SetLib__(MPM_Mesh.ListNodes[i]);
	
      /* Active those nodes that interact with the particle */
      asign_particle_to_nodes(i, MPM_Mesh.ListNodes[i], FEM_Mesh);

      /* Free memory */
      free__MatrixLib__(Poligon_Coordinates);
	
      break;
    }
      
      
  }
  if(!Init_p)
    {
      fprintf(stderr,"%s : %s %i %s \n",
	      "Error in initialize__GIMP__",
	      "Particle",i,"Was not initialise");
      exit(EXIT_FAILURE);
    }

  }
  
}

/*********************************************************************/

static Matrix voxel__GIMP__(Matrix l_p, double Volume_p)
/*!
  Get the initial voxel lenght
  Clarify that the volume measure depends if the problem is 2D (area)
  or 3D (volume).
*/
{

  int Ndim = NumberDimensions; 

  /* Fill Voxel */  
  for(int j =0 ; j<Ndim ; j++)
    {
      l_p.nV[j] = 0.5*pow(Volume_p,(double)1/Ndim);
    }

  return l_p;
  
}

/*********************************************************************/

/* Uniform GIMP shape function */
static double Sip__GIMP__(double L, double lp, double Delta_xp){

  /* Evaluation of the shape function */

  if ((-lp < Delta_xp) && (Delta_xp <= lp)){
    return 1 - 0.5*(DSQR(Delta_xp) + lp*lp)*(double)1/(L*lp);
  }
  else if (((-L-lp) < Delta_xp) && (Delta_xp <= (-L+lp))){
    return (double)(0.25/(L*lp))*DSQR(L+lp+Delta_xp);
  }
  else if (((L-lp) < Delta_xp) && (Delta_xp <= (L+lp))){
    return (double)(0.25/(L*lp))*DSQR(L+lp-Delta_xp);
  }
  else if (((-L+lp) < Delta_xp) && (Delta_xp <= -lp)){
    return 1 + (double)Delta_xp/L;
  }
  else if ((lp < Delta_xp) && (Delta_xp <= (L-lp))){
    return 1 - (double)Delta_xp/L;
  }
  else {
    return (double)0.0;
  }
}

/*********************************************************************/

/* Uniform GIMP derivative shape function */
static double dSip__GIMP__(double L, double lp, double Delta_xp){

  /* Evaluation of the shape function */
  if (((-L-lp) < Delta_xp) && (Delta_xp <= (-L+lp))){
    return (double)(0.5/(L*lp))*(L+lp+Delta_xp);
  }
  else if (((-L+lp) < Delta_xp) && (Delta_xp <= -lp)){
    return (double)1/L;
  }
  else if ((-lp < Delta_xp) && (Delta_xp <= lp)){
    return -(double)Delta_xp/(L*lp);
  }
  else if ((lp < Delta_xp) && (Delta_xp <= (L-lp))){
    return -(double)1/L;
  }
  else if (((L-lp) < Delta_xp) && (Delta_xp <= (L+lp))){
    return -(double)(0.5/(L*lp))*(L+lp-Delta_xp);
  }
  else {
    return (double)0.0;
  }
}

/*********************************************************************/

/* Uniform GIMP shape function 2D */
Matrix N__GIMP__(Matrix Delta_Xp, Matrix lp, double L){

  /* 1º Variable declaration */
  int Ndim = NumberDimensions;
  int Nnodes = Delta_Xp.N_rows;
  Matrix S_Ip = allocZ__MatrixLib__(1,Nnodes);

  /* 2º Fill the shape function array */
  switch(Ndim){
  case 2 :
    for(int i = 0 ; i<Nnodes ; i++){
      /* 3º Shape function in this node */
      S_Ip.nV[i] =
	Sip__GIMP__(L, lp.nV[0], Delta_Xp.nM[i][0])*
	Sip__GIMP__(L, lp.nV[1], Delta_Xp.nM[i][1]);
    }
    break;
  case 3:
    for(int i = 0 ; i<Nnodes ; i++){
      /* 3º Shape function in this node */
      S_Ip.nV[i] =
	Sip__GIMP__(L, lp.nV[0], Delta_Xp.nM[i][0])*
	Sip__GIMP__(L, lp.nV[1], Delta_Xp.nM[i][1])*
	Sip__GIMP__(L, lp.nV[2], Delta_Xp.nM[i][2]);
    }
    break;
  }

  /* 4º Output */
  return S_Ip;
}

/*********************************************************************/

/* Uniform GIMP derivative shape function 2D */
Matrix dN__GIMP__(Matrix Delta_xp, Matrix lp, double L){

  /* 1º Variable declaration */
  int Ndim = NumberDimensions;
  int Nnodes = Delta_xp.N_rows;
  Matrix dS_Ip = allocZ__MatrixLib__(Nnodes,Ndim); 
  
  /* 2º Fill the shape function array */
  switch(Ndim){    
  case 2 :
    for(int i = 0 ; i<Nnodes ; i++){
      /* 3º Gradient of the shape function for each node*/
      dS_Ip.nM[i][0] =
	dSip__GIMP__(L, lp.nV[0], Delta_xp.nM[i][0])*
	Sip__GIMP__(L, lp.nV[1], Delta_xp.nM[i][1]);
      dS_Ip.nM[i][1] =
	Sip__GIMP__(L, lp.nV[0], Delta_xp.nM[i][0])*
	dSip__GIMP__(L, lp.nV[1], Delta_xp.nM[i][1]);
    }
    break;
  case 3 :
    for(int i = 0 ; i<Nnodes ; i++){
      /* 3º Gradient of the shape function for each node*/
      dS_Ip.nM[i][0] =
	dSip__GIMP__(L, lp.nV[0], Delta_xp.nM[i][0])*
	Sip__GIMP__(L, lp.nV[1], Delta_xp.nM[i][1])*
	Sip__GIMP__(L, lp.nV[2], Delta_xp.nM[i][2]);
      dS_Ip.nM[i][1] =
	Sip__GIMP__(L, lp.nV[0], Delta_xp.nM[i][0])*
	dSip__GIMP__(L, lp.nV[1], Delta_xp.nM[i][1])*
	Sip__GIMP__(L, lp.nV[2], Delta_xp.nM[i][2]);
      dS_Ip.nM[i][2] =
	Sip__GIMP__(L, lp.nV[0], Delta_xp.nM[i][0])*
	Sip__GIMP__(L, lp.nV[1], Delta_xp.nM[i][1])*
	dSip__GIMP__(L, lp.nV[2], Delta_xp.nM[i][2]);      
    }
    break;
  }
  /* 4º Output */
  return dS_Ip;
}

/*********************************************************************/

ChainPtr tributary__GIMP__(Matrix Xi_p,int Elem_GP,Matrix lp,Mesh FEM_Mesh)
{

  ChainPtr * Table_Nodes = NULL;
  ChainPtr Triburary_Nodes = NULL;
  ChainPtr ChainElements = NULL;
  int * Tributary_Elements;
  int * NodesElem;
  int Elem_i;
  int Num_Elem;
  double Dist[2];

  /* Get the reference distance measured from the center of the element */
  for(int i = 0 ; i<2; i++){
    Dist[i] = 1 - 2*lp.nV[i]/FEM_Mesh.DeltaX;
  }

  /* Check if I am in the central area */
  if ((fabs(Xi_p.nV[0]) <= Dist[0]) &&
      (fabs(Xi_p.nV[1]) <= Dist[1])){
    Triburary_Nodes = copy__SetLib__(FEM_Mesh.Connectivity[Elem_GP]);
  }    
  /* Check if I am in the 1º Quadrant */
  else if((Xi_p.nV[0]>=0) &&
	  (Xi_p.nV[1]>=0)){
    /* Create an array with the nodes of the element */
    NodesElem = set_to_memory__SetLib__(FEM_Mesh.Connectivity[Elem_GP],4);
    if((fabs(Xi_p.nV[0]) >= Dist[0]) &&
       (fabs(Xi_p.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */       
      ChainElements = intersection__SetLib__(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					     FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = lenght__SetLib__(ChainElements);
      Tributary_Elements = set_to_memory__SetLib__(ChainElements,Num_Elem);
      free__SetLib__(&ChainElements);

       /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(Xi_p.nV[0]) <= Dist[0]) &&
	    (fabs(Xi_p.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = intersection__SetLib__(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					     FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Num_Elem = lenght__SetLib__(ChainElements);
      Tributary_Elements = set_to_memory__SetLib__(ChainElements,Num_Elem);
      free__SetLib__(&ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(Xi_p.nV[0]) >= Dist[0]) &&
	    (fabs(Xi_p.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[2]];
      Num_Elem = lenght__SetLib__(ChainElements);
      Tributary_Elements = set_to_memory__SetLib__(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    /* Free memory */
    free(NodesElem);    
  }  
  /* Check if I am in the 2º Quadrant */
  else if((Xi_p.nV[0]<=0) &&
	  (Xi_p.nV[1]>=0)){
    /* Create an array with the nodes of the element */
    NodesElem = set_to_memory__SetLib__(FEM_Mesh.Connectivity[Elem_GP],4);

    if((fabs(Xi_p.nV[0]) <= Dist[0]) &&
       (fabs(Xi_p.nV[1]) >= Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = intersection__SetLib__(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					     FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Num_Elem = lenght__SetLib__(ChainElements);      
      Tributary_Elements =  set_to_memory__SetLib__(ChainElements,Num_Elem);
      free__SetLib__(&ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(Xi_p.nV[0]) >= Dist[0]) &&
	    (fabs(Xi_p.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = intersection__SetLib__(FEM_Mesh.NodeNeighbour[NodesElem[3]],
					     FEM_Mesh.NodeNeighbour[NodesElem[0]]);
      Num_Elem = lenght__SetLib__(ChainElements);
      Tributary_Elements = set_to_memory__SetLib__(ChainElements,Num_Elem);
      free__SetLib__(&ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(Xi_p.nV[0]) >= Dist[0]) &&
	    (fabs(Xi_p.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[3]];
      Num_Elem = lenght__SetLib__(ChainElements);
      Tributary_Elements =  set_to_memory__SetLib__(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }

    /* Free memory */
    free(NodesElem);
    
  }  
  /* Check if I am in the 3º Quadrant */
  else if((Xi_p.nV[0]<=0) &&
	  (Xi_p.nV[1]<=0)){
    /* Create an array with the nodes of the element */
    NodesElem = set_to_memory__SetLib__(FEM_Mesh.Connectivity[Elem_GP],4);   
    if((fabs(Xi_p.nV[0]) >= Dist[0]) &&
       (fabs(Xi_p.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = intersection__SetLib__(FEM_Mesh.NodeNeighbour[NodesElem[3]],
					     FEM_Mesh.NodeNeighbour[NodesElem[0]]);
      Num_Elem = lenght__SetLib__(ChainElements);
      Tributary_Elements = set_to_memory__SetLib__(ChainElements,Num_Elem);
      free__SetLib__(&ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(Xi_p.nV[0]) <= Dist[0]) &&
	    (fabs(Xi_p.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */
      ChainElements = intersection__SetLib__(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					     FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = lenght__SetLib__(ChainElements);
      Tributary_Elements = set_to_memory__SetLib__(ChainElements,Num_Elem);
      free__SetLib__(&ChainElements);
      
      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(Xi_p.nV[1]) >= Dist[1]) &&
	    (fabs(Xi_p.nV[0]) >= Dist[0])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[0]];
      Num_Elem = lenght__SetLib__(ChainElements);
      Tributary_Elements = set_to_memory__SetLib__(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    /* Free memory */
    free(NodesElem);    
  }  
  /* Check if it I am the 4º Quadrant */
  else if((Xi_p.nV[0]>=0) &&
	  (Xi_p.nV[1]<=0)){
    /* Create an array with the nodes of the element */
    NodesElem = set_to_memory__SetLib__(FEM_Mesh.Connectivity[Elem_GP],4);

    if((fabs(Xi_p.nV[0]) <= Dist[0]) &&
       (fabs(Xi_p.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = intersection__SetLib__(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					     FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = lenght__SetLib__(ChainElements);
      Tributary_Elements =  set_to_memory__SetLib__(ChainElements,Num_Elem);
      free__SetLib__(&ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(Xi_p.nV[0]) >= Dist[0]) &&
	    (fabs(Xi_p.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = intersection__SetLib__(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					     FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = lenght__SetLib__(ChainElements);
      Tributary_Elements = set_to_memory__SetLib__(ChainElements,Num_Elem);
      free__SetLib__(&ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    else if((fabs(Xi_p.nV[0]) >= Dist[0]) &&
	    (fabs(Xi_p.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[1]];
      Num_Elem = lenght__SetLib__(ChainElements);
      Tributary_Elements = set_to_memory__SetLib__(ChainElements,Num_Elem);
      
      /* Iterate in the list and select the union of the sets of nodes */
      Table_Nodes = malloc(Num_Elem*sizeof(ChainPtr));
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Table_Nodes[i] = FEM_Mesh.Connectivity[Elem_i];
      }
      Triburary_Nodes = union__SetLib__(Table_Nodes,Num_Elem);
      /* Free memory */
      free(Tributary_Elements);
      free(Table_Nodes);
      Table_Nodes = NULL;
      
    }
    /* Free memory */
    free(NodesElem);
  }
  else{
    printf("%s : %s \n",
	   "Error in Tributary_Nodes_GIMP",
	   "Unlocated GP in the element");
    printf("%s : (%f;%f) \n",
	   "Natural coordinates of the GP",
	   Xi_p.nV[0],Xi_p.nV[1]);
    exit(EXIT_FAILURE);
  }
 
  return Triburary_Nodes;
}

  
/*   /\* Chain with the tributary elements, this is the list of element near the */
/*      gauss point, including where it is *\/ */

/*   /\* Iterate in the list and select the union of the sets of nodes *\/ */
/*   Table_Elem = malloc(NumNodesElem*sizeof(ChainPtr)); */
/*   for(int i = 0 ; i<NumNodesElem ; i++){ */
/*     Table_Elem[i] = FEM_Mesh.NodeNeighbour[NodesElem[i]]; */
/*   } */
/*   Triburary_Elements = union__SetLib__(Table_Elem,NumNodesElem); */
/*   /\* Free memory *\/ */
/*   free(NodesElem); */
/*   free(Table_Elem); */
/*   Table_Elem = NULL; */
  
/*   /\* List with the tributary nodes *\/ */
/*   Num_Elem = lenght__SetLib__(Triburary_Elements); */
/*   List_Elements = set_to_memory__SetLib__(Triburary_Elements,Num_Elem); */

/*   /\* Free the chain with the tributary elements *\/ */
/*   free__SetLib__(Triburary_Elements); */
  
/*   /\* Fill the chain with the preliminary tributary nodes *\/ */
/*   Table_ElemNodes = malloc(Num_Elem*sizeof(ChainPtr)); */
/*   for(int i = 0 ; i<Num_Elem ; i++){ */
/*     Table_ElemNodes[i] = FEM_Mesh.Connectivity[List_Elements[i]]; */
/*   } */
  
/*   List_Nodes = union__SetLib__(Table_ElemNodes,Num_Elem); */
  
/*   /\* Free the array wit the list of tributary elements *\/ */
/*   free(List_Elements); */
/*   free(Table_ElemNodes); */
/*   Table_ElemNodes = NULL; */
  
/*   /\* Initialize the iterator to iterate over the list of tributary nodes *\/ */
/*   iPtr = List_Nodes; */

/*   /\* Loop over the chain with the tributary nodes *\/ */
/*   while(iPtr != NULL){ */

/*     /\* Assign to a pointer the coordinates of the nodes *\/ */
/*     X_I.nV = FEM_Mesh.Coordinates.nM[iPtr->I]; */

/*     /\* Get a vector from the GP to the node *\/ */
/*     Distance = substraction__MatrixLib__(X_GP,X_I); */

/*     /\* If the node is near the GP push in the chain *\/ */
/*     if(norm__MatrixLib__(Distance,2) <= Ra){ */
/*       PushNodeTop(&Triburary_Nodes,iPtr->I); */
/*     } */

/*     /\* Free memory of the distrance vector *\/ */
/*     free__MatrixLib__(Distance); */

/*     /\* Update pointer index *\/ */
/*     iPtr = iPtr->next; */
/*   } */
/*   /\* Free memory *\/ */
/*   free__SetLib__(&List_Nodes); */
  
/*   return Triburary_Nodes; */
/* } */


