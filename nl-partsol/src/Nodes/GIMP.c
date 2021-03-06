#include "nl-partsol.h"

/*
  Global variables
*/
double Thickness_Plain_Stress;

/*
  Auxiliar functions
 */
static Matrix voxel__GIMP__(Matrix, double);
static double Sip__GIMP__(double, double, double);
static double dSip__GIMP__(double, double, double);

/*********************************************************************/


void initialize__GIMP__(
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
{

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

  /* Variables to store the conectivity and coordinates of the element of the particle */
  ChainPtr Elem_p_Connectivity;
  Matrix Elem_p_Coordinates;

  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++)
  {
   
  /* Get some properties for each particle */ 
  X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[i]);
  Xi_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_EC.nM[i]);
  l_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.lp.nM[i]);
  rho_p = MPM_Mesh.Phi.rho.nV[i];
  m_p = MPM_Mesh.Phi.mass.nV[i];
  th_p = Thickness_Plain_Stress;
  A_p = m_p/(rho_p*th_p);

  /* Supose that the particle was not initilise */
  Init_p = false;

  for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++)
  {

    /* Get the element properties */
    Elem_p_Connectivity = FEM_Mesh.Connectivity[j];
    Elem_p_Coordinates = get_nodes_coordinates__MeshTools__(Elem_p_Connectivity, FEM_Mesh.Coordinates);
      
    /* Check out if the GP is in the Element */
    if(FEM_Mesh.In_Out_Element(X_p,Elem_p_Coordinates))
    {

      /* Particle will be initilise */
      Init_p = true;
	
      /* With the element connectivity get the node close to the particle */
      MPM_Mesh.I0[i] = get_closest_node__MeshTools__(X_p,Elem_p_Connectivity,FEM_Mesh.Coordinates);

      /* If the GP is in the element, get its natural coordinates */
      Xi_p.nV = MPM_Mesh.Phi.x_EC.nM[i];
      X_to_Xi__Q4__(Xi_p,X_p,Elem_p_Coordinates);

      /* Initialize the voxel lenght */
      l_p = voxel__GIMP__(l_p, A_p);

      /* Get the initial connectivity of the particle */
      MPM_Mesh.ListNodes[i] = tributary__GIMP__(Xi_p,j,l_p,FEM_Mesh);

      /* Measure the size of the connectivity */
      MPM_Mesh.NumberNodes[i] = lenght__SetLib__(MPM_Mesh.ListNodes[i]);
	
      /* Active those nodes that interact with the particle */
      asign_to_nodes__Particles__(i, j, MPM_Mesh.I0[i], MPM_Mesh.ListNodes[i], FEM_Mesh);

      /* Free memory */
      free__MatrixLib__(Elem_p_Coordinates);
	
      break;
    }

    /* Free memory */
    free__MatrixLib__(Elem_p_Coordinates);  
      
  }
  if(!Init_p)
  {
    fprintf(stderr,"%s : %s %i %s \n","Error in initialize__GIMP__","Particle",i,"Was not initialise");
    exit(EXIT_FAILURE);
  }

  }
  
}

/*********************************************************************/

static Matrix voxel__GIMP__(
  Matrix l_p,
  double Volume_p)
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
static double Sip__GIMP__(
  double L,
  double lp, double Delta_xp)
{

  /* Evaluation of the shape function */

  if ((-lp < Delta_xp) && (Delta_xp <= lp))
  {
    return 1 - 0.5*(DSQR(Delta_xp) + lp*lp)*(double)1/(L*lp);
  }
  else if (((-L-lp) < Delta_xp) && (Delta_xp <= (-L+lp)))
  {
    return (double)(0.25/(L*lp))*DSQR(L+lp+Delta_xp);
  }
  else if (((L-lp) < Delta_xp) && (Delta_xp <= (L+lp)))
  {
    return (double)(0.25/(L*lp))*DSQR(L+lp-Delta_xp);
  }
  else if (((-L+lp) < Delta_xp) && (Delta_xp <= -lp))
  {
    return 1 + (double)Delta_xp/L;
  }
  else if ((lp < Delta_xp) && (Delta_xp <= (L-lp)))
  {
    return 1 - (double)Delta_xp/L;
  }
  else 
  {
    return (double)0.0;
  }
}

/*********************************************************************/

/* Uniform GIMP derivative shape function */
static double dSip__GIMP__(
  double L, 
  double lp, 
  double Delta_xp)
{

  /* Evaluation of the shape function */
  if (((-L-lp) < Delta_xp) && (Delta_xp <= (-L+lp)))
  {
    return (double)(0.5/(L*lp))*(L+lp+Delta_xp);
  }
  else if (((-L+lp) < Delta_xp) && (Delta_xp <= -lp))
  {
    return (double)1/L;
  }
  else if ((-lp < Delta_xp) && (Delta_xp <= lp))
  {
    return -(double)Delta_xp/(L*lp);
  }
  else if ((lp < Delta_xp) && (Delta_xp <= (L-lp)))
  {
    return -(double)1/L;
  }
  else if (((L-lp) < Delta_xp) && (Delta_xp <= (L+lp)))
  {
    return -(double)(0.5/(L*lp))*(L+lp-Delta_xp);
  }
  else 
  {
    return (double)0.0;
  }
}

/*********************************************************************/

Matrix N__GIMP__(
  Matrix Delta_Xp, 
  Matrix lp, 
  double L)
{

  /* Variable declaration */
  int Ndim = NumberDimensions;
  int Nnodes = Delta_Xp.N_rows;
  Matrix S_Ip = allocZ__MatrixLib__(1,Nnodes);

  for(int i = 0 ; i<Nnodes ; i++)
  {

    // Initialise to one for the node i
    S_Ip.nV[i] = 1.0;

    for(int j = 0 ; j<Ndim ; j++)
    {
      S_Ip.nV[i] *= Sip__GIMP__(L,lp.nV[j],Delta_Xp.nM[i][j]);
    }
  }

  return S_Ip;
}

/*********************************************************************/


Matrix dN__GIMP__(
  Matrix Delta_xp,
  Matrix lp,
  double L)
{

  /* Variable declaration */
  int Ndim = NumberDimensions;
  int Nnodes = Delta_xp.N_rows;
  Matrix dS_Ip = allocZ__MatrixLib__(Nnodes,Ndim); 
  
  /* 2º Fill the shape function array */
  for(int i = 0 ; i<Nnodes ; i++)
  {

    // Initialise to one in each direction for the node i
    for(int j = 0 ; j<Ndim ; j++)
    {
      dS_Ip.nM[i][j] = 1.0;
    }

    for(int j = 0 ; j<Ndim ; j++)
    {
      if(i == j)
      {
        dS_Ip.nM[i][j] *= dSip__GIMP__(L,lp.nV[j],Delta_xp.nM[i][j]);
      }
      else
      {
        dS_Ip.nM[i][j] *= Sip__GIMP__(L,lp.nV[j],Delta_xp.nM[i][j]);
      }
    }
  } 


  return dS_Ip;
}

/*********************************************************************/

/*
  This function is a mess, sorry for that
*/
ChainPtr tributary__GIMP__(
  Matrix Xi_p,
  int Elem_GP,
  Matrix lp,
  Mesh FEM_Mesh)
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
