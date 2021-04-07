#include "nl-partsol.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))


// Define some auxiliar functions
static void get_sourrounding_elements(Mesh);
static void fill_nodal_locality(Mesh);
static ChainPtr node_I_locality(int, Mesh);
static ChainPtr ring_search_nodal_locality(ChainPtr *, ChainPtr, Mesh);
static double mesh_size(Mesh);

/**********************************************************************/

Mesh GramsBox(char * Name_File)
/*    

      GramsBox (Type=GID,File=FEM_Mesh.msh) { 
      left = GramsBoundary {
      BcDirichlet U NULL
      BcDirichlet V NULL	
      }
      right = GramsBoundary {
      BcDirichlet U NULL
      BcDirichlet V NULL	
      }
      top = GramsBoundary {
      BcDirichlet U NULL
      BcDirichlet V NULL
      }
      bottom = GramsBoundary {
      BcDirichlet U NULL
      BcDirichlet V NULL
      }
      }
*/
{
  /* Asign material library to an auxiliar variable */
  Mesh FEM_Mesh;

  /* Simulation file */
  FILE * Sim_dat;
  char Route_Mesh[MAXC] = {0};
  char FileMeshRoute[MAXC];
   
  /* Parse line of GramsBox */
  int Num_GramsBox;
  char Line_GramsBox[MAXC] = {0};
  char * Parse_GramsBox[MAXW] = {NULL};
  int Num_GramsBox_Prop;
  char * Parse_Mesh_id[MAXW] = {NULL};

  /* Parse lines of GramsBoundary */
  int Num_GramsBoundary;
  char Line_GramsBoundary[MAXC] = {0};
  char * Parse_GramsBoundary[MAXW] = {NULL};
  int NumBounds;

  /* Auxiliar variable for status */
  char * STATUS_LINE;
    
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsBox()",
	    "Incorrect lecture of",
	    Name_File);
    exit(EXIT_FAILURE);
  }

  /* Generate route */
  generate_route(Route_Mesh,Name_File);
  printf("%s : %s \n","Read mesh from", Route_Mesh);

  /* Read GramsBox line  */
  while( fgets(Line_GramsBox, sizeof(Line_GramsBox), Sim_dat) != NULL ){
    
    /* Read the line with the space as separators */
    Num_GramsBox = parse (Parse_GramsBox, Line_GramsBox," \n\t");
    
    if (Num_GramsBox < 0){
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsBox ()",
	      "Parser failed");
      exit(EXIT_FAILURE);
    }

    /* Find GramsBox line */
    if ((Num_GramsBox > 0) &&
    	(strcmp(Parse_GramsBox[0],"GramsBox") == 0 ) &&
	((strcmp(Parse_GramsBox[2],"{") == 0))){

      /* Read the type of the mesh and the name of the file */
      Num_GramsBox_Prop = parse (Parse_Mesh_id, Parse_GramsBox[1],"(=,)");
      if( (Num_GramsBox_Prop != 4) ||
	  (strcmp(Parse_Mesh_id[0],"Type") != 0) ||
	  (strcmp(Parse_Mesh_id[2],"File") != 0) ){
	fprintf(stderr,"%s : %s \n",
		"Error in GramsBox()",
		"Use this format -> (Type=str,File=str) !!!");
	exit(EXIT_FAILURE);
      }

      /* Read GID-type mesh */
      if(strcmp(Parse_Mesh_id[1],"GID") == 0)
      {
	/* Read file with the mesh */
	sprintf(FileMeshRoute,"%s%s",Route_Mesh,Parse_Mesh_id[3]);

	puts("*************************************************");
	printf(" \t %s : \n \t %s \n",
	       "* Read GID mesh in",FileMeshRoute);
	FEM_Mesh = ReadGidMesh__MeshTools__(FileMeshRoute);
      }
      else
      {
        fprintf(stderr,"%s : %s %s \n","Error in GramsBox(Type=*, )","Unrecognized mesh Type",Parse_Mesh_id[1]);
      }

      /* By default the number of BCCs is 0 */
      NumBounds = 0;

      /* Check format GramsBox (Type=,File=) { } */
      if((Num_GramsBox == 4) &&
	 (strcmp(Parse_GramsBox[3],"}") == 0)){
	puts("*************************************************");
	printf(" \t %s : \n \t %i \n",
	       "* Number of boundaries",0);
	break;
      }

      /* Initial line */
      STATUS_LINE =
	fgets(Line_GramsBoundary,sizeof(Line_GramsBoundary),Sim_dat);
      if(STATUS_LINE == NULL){
	fprintf(stderr,"%s : %s \n",
		"Error in GramsBox()",
		"Unspected EOF !!!");
	exit(EXIT_FAILURE);
      }
      Num_GramsBoundary = parse(Parse_GramsBoundary,Line_GramsBoundary," =\t\n");

      /* Check format GramsBox (Type=,File=) { 
	 }
      */
      if((Num_GramsBoundary > 0) &&
	 (strcmp(Parse_GramsBoundary[0],"}") == 0)){
	puts("*************************************************");
	printf(" \t %s : \n \t %i \n",
	       "* Number of boundaries",0);
	break;
      }

      /* Loop to count the number of boundaries */
      while(STATUS_LINE != NULL){
	if((Num_GramsBoundary > 0) &&
	   (strcmp(Parse_GramsBoundary[0],"GramsBoundary") == 0)){
	  NumBounds++;
	}
	/* Continue reading to the end*/
	STATUS_LINE = fgets(Line_GramsBoundary,
			    sizeof(Line_GramsBoundary),Sim_dat);
	Num_GramsBoundary = parse(Parse_GramsBoundary,
			       Line_GramsBoundary," =\t\n");	
      }
      
    }
  }

  /* Close file */
  fclose(Sim_dat);

  /**************************************************/
  /********** Set the minimum mesh size *************/
  /**************************************************/
  FEM_Mesh.DeltaX = mesh_size(FEM_Mesh);

  puts("*************************************************");
  printf(" \t %s : \n \t %f \n","* Mesh size",FEM_Mesh.DeltaX);

  /**************************************************/	
  /*** Generate nodal connectivity of the mesh : ****/
  /******* list of elements near to a node ***********/
  /**************************************************/
  FEM_Mesh.NumNeighbour =  (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh,sizeof(int)); 
  FEM_Mesh.NodeNeighbour = (ChainPtr *)malloc(FEM_Mesh.NumNodesMesh*sizeof(ChainPtr));
  get_sourrounding_elements(FEM_Mesh);
  
  /*
    Initialize nodal locality of each node :

    1º Define a pointer with the number of nodes close to each node

    2ª Define a table of sets with the list of nodes close to a node

    3ª Fill the table with the nodal locality and the pointer with the
      number of nodes close to each node
  */

  FEM_Mesh.SizeNodalLocality = (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh,sizeof(int));
  FEM_Mesh.NodalLocality = (ChainPtr *)malloc(FEM_Mesh.NumNodesMesh*sizeof(ChainPtr));
  fill_nodal_locality(FEM_Mesh);
  
  /**************************************************/	
  /** Initialize particle connectivity of each node */
  /**************************************************/
  FEM_Mesh.NumParticles = (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh,sizeof(int));
  FEM_Mesh.I_particles = (ChainPtr *)malloc(FEM_Mesh.NumNodesMesh*sizeof(ChainPtr));


  /**************************************************/  
  /******** Allocate array with the boundaries ******/
  /**************************************************/
  if(NumBounds > 0){
    puts("*************************************************");
    printf(" \t %s (%i) : \n","* Boundary conditions",NumBounds);
    FEM_Mesh.Bounds.NumBounds = NumBounds;
    
    if(strcmp(Formulation,"-u") == 0)
    {
      FEM_Mesh.Bounds = GramsBoundary(Name_File,NumBounds);  
    }
    else if(strcmp(Formulation,"-upw") == 0)
    {
      FEM_Mesh.Bounds = Read_upw_Dirichlet_Boundary_Conditions__InOutFun__(Name_File,NumBounds);
    }
    
  }


  return FEM_Mesh;
}

/*********************************************************************/

static double mesh_size(Mesh FEM_Mesh)
/*
  Function to get the minimum mesh size.
*/
{

  /* Auxiliar variables of the function */
  int NumElemMesh = FEM_Mesh.NumElemMesh;
  int NumNodesElem; /* Number of nodes of each element */
  int * Connectivity; /* Connectivity of the element */
  Matrix Poligon; /* Element Poligon */
  Matrix X_eval = allocZ__MatrixLib__(1,2); /* Where to evaluate the shape function */
  X_eval.nV[0] = 0;
  X_eval.nV[1] = 0;
  Matrix dNdx; /* Gradient of the shapefunction for each node */
  double MinElementSize_aux;
  double MinElementSize = 10e16;

  /* 1º Loop over the elements in the mesh */
  for(int i = 0 ; i<NumElemMesh ; i++)
  {

    /* 2º Connectivity of the Poligon */
    NumNodesElem = FEM_Mesh.NumNodesElem[i];
    Connectivity = set_to_memory__SetLib__(FEM_Mesh.Connectivity[i],NumNodesElem);
    
    /* 4º Get the gradient of the element for each node */
    if((NumNodesElem == 3) && (NumberDimensions == 2))
    { 
      /* Fill the triangular element with the coordinates of the nodes */
      Poligon = allocZ__MatrixLib__(3,2);
      for(int k = 0; k<3; k++)
      {
        for(int l = 0 ; l<2 ; l++)
        {
          Poligon.nM[k][l] = FEM_Mesh.Coordinates.nM[Connectivity[k]][l];
        }
      }

      /* Get the gradient of the triangle */
      dNdx = dN__T3__(X_eval,Poligon);
      free__MatrixLib__(Poligon);
      
      /* Get the minimum minimum height of the triangle */
      for(int j = 0 ; j<3 ; j++)
      {
        MinElementSize_aux = 1/pow(dNdx.nM[0][j]*dNdx.nM[0][j] + dNdx.nM[1][j]*dNdx.nM[1][j],0.5);
        MinElementSize = DMIN(MinElementSize,MinElementSize_aux);
      }
      /* Free memory */
      free__MatrixLib__(dNdx);
      
    }
    else if((NumNodesElem == 4) && (NumberDimensions == 2))
    { 
      /* Fill the quadrilateral element with the coordinates of the nodes */
      Poligon = allocZ__MatrixLib__(4,2);

      /* Fill the poligon with vectors */
      for(int k = 0; k<3; k++)
      {
        for(int l = 0 ; l<2 ; l++)
        {
          Poligon.nM[k][l] = FEM_Mesh.Coordinates.nM[Connectivity[k+1]][l] - FEM_Mesh.Coordinates.nM[Connectivity[k]][l];
        }
      }

      for(int l = 0 ; l<2 ; l++)
      {
        Poligon.nM[3][l] = FEM_Mesh.Coordinates.nM[Connectivity[0]][l] - FEM_Mesh.Coordinates.nM[Connectivity[3]][l];
      }
      
      /* Get the minimum minimum height of the triangle */
      for(int k = 0 ; k<4 ; k++)
      {
        MinElementSize_aux = pow(Poligon.nM[k][0]*Poligon.nM[k][0] + Poligon.nM[k][1]*Poligon.nM[k][1] , 0.5);
        MinElementSize = DMIN(MinElementSize,MinElementSize_aux);
      }

      /* Free memory */
      free__MatrixLib__(Poligon);

    }
    else
    {
      printf("%s : %s %i %s \n",
       "Error in mesh_size","Element with ",NumNodesElem,"nodes is not implemented !!!" );
      exit(EXIT_FAILURE);
    }

    /* Free memory */
    free(Connectivity);
    
  }

  /* Free memory */
  free__MatrixLib__(X_eval);

  return MinElementSize;

}

/*********************************************************************/

static void get_sourrounding_elements(Mesh FEM_Mesh)
{

  /* Variable declaration */
  int * Element_Connectivity;
  int NumNodesElem;
  
  /* 1º Start the search of neighbour for each node */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++)
  {
    /* 2º Loop over all the elements in the mesh */
    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++)
    {
      NumNodesElem = FEM_Mesh.NumNodesElem[j];
      Element_Connectivity = set_to_memory__SetLib__(FEM_Mesh.Connectivity[j],NumNodesElem);

      /* 3º Loop over the all the node in an element */
      for(int k = 0 ; k<NumNodesElem ; k++)
      {
        /* 4º If my node belong to the element */
        if(Element_Connectivity[k] == i)
        {

          /* 5º Introduce the element in the chain */
          push__SetLib__(&FEM_Mesh.NodeNeighbour[i], j);
          
          /* 6º Update the counter */
          FEM_Mesh.NumNeighbour[i] += 1;    
        }

      }

      /* Free memory */
      free(Element_Connectivity);
    }
    
  }
  
}

/*********************************************************************/

static void fill_nodal_locality(Mesh FEM_Mesh)
{
  /*
    Auxiliar variables for the nodal neighborhood reconstruction
  */
  int Num_nodal_rings = 3; // Number of search rings
  int k_nodal_ring; // Current search ring
  ChainPtr Search_Set; // Auxiliar set for recursive search

  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++)
  {

    if (Num_nodal_rings == 1)
    {
      FEM_Mesh.NodalLocality[i] = node_I_locality(i, FEM_Mesh);
    }
    else
    {

      k_nodal_ring = 0;
      Search_Set = NULL;
      push__SetLib__(&Search_Set, i);

      while(k_nodal_ring < Num_nodal_rings)
      {

        Search_Set = ring_search_nodal_locality(&FEM_Mesh.NodalLocality[i], Search_Set, FEM_Mesh);

        k_nodal_ring++;
      }

    }

    FEM_Mesh.SizeNodalLocality[i] = lenght__SetLib__(FEM_Mesh.NodalLocality[i]);

  }

}

/*********************************************************************/

static ChainPtr node_I_locality(int I, Mesh FEM_Mesh)
{

  /* Define output */
  ChainPtr Nodes = NULL;

  /* Number of elements sourronding the node */
  int NumNeighbour = FEM_Mesh.NumNeighbour[I];
  
  /* Index of the elements sourronding the node */
  int * NodeNeighbour = set_to_memory__SetLib__(FEM_Mesh.NodeNeighbour[I],NumNeighbour);
  
  /* Table with the nodes of each element */
  ChainPtr * Table_ElemNodes = malloc(NumNeighbour*sizeof(ChainPtr));

  /* Fill each position of the table with a list of nodes in the element */
  for(int i = 0 ; i<NumNeighbour ; i++)
  {
    Table_ElemNodes[i] = FEM_Mesh.Connectivity[NodeNeighbour[i]];
  }
  
  /* Free table with elements */
  free(NodeNeighbour);
  
  /* Get the union of this nodes */
  Nodes = union__SetLib__(Table_ElemNodes,NumNeighbour);
  
  /* Free table with the nodes of each elements */
  free(Table_ElemNodes);

  /* Return nodes close to the node I */
  return Nodes;
}

/*********************************************************************/

static ChainPtr ring_search_nodal_locality(ChainPtr * Set_k, ChainPtr Search_Set, Mesh FEM_Mesh)
{

  /*
    Variables
  */
  ChainPtr aux_Set = NULL;
  ChainPtr new_Search_Set = NULL;
  ChainPtr i_Search_Set = Search_Set; 
  
  /*
    Loop in the search set
  */
  while (i_Search_Set != NULL)
  { 

    /*
      For each node in the set, get the closest node
    */
    aux_Set = node_I_locality(i_Search_Set->I, FEM_Mesh);

    /*
      Loop in the closest set of nodes of the original set
    */
    while (aux_Set != NULL)
    {
      /*
        Update the set with the new node and the search set
      */
      if (!inout__SetLib__(*Set_k, aux_Set->I))
      {
        push__SetLib__(Set_k, aux_Set->I);
        push__SetLib__(&new_Search_Set, aux_Set->I);
      }

      /*
        Update second iterator
      */
      aux_Set = aux_Set->next; 

    }

    /*  
      Free auxiliar set
    */
    free__SetLib__(&aux_Set);

    /*
      Update first iterator
    */
    i_Search_Set = i_Search_Set->next; 

  }

  /*
    Free old search list
  */
  free__SetLib__(&Search_Set);

  return new_Search_Set;
}

/*********************************************************************/
