#include "nl-partsol.h"

Mesh ReadGidMesh(char * MeshName)
{
  /* Output element mesh */
  Mesh GID_Mesh;

  /* Pointer to the file */
  FILE * MeshFile;

  /* Variable to save the line of a text */
  char line[MAXC] = {0};
  char line_coords[MAXC] = {0};
  char line_connectivity[MAXC] = {0};
  char * Aux_line;
  
  /* Char pointer to store the words in a line */
  char * words[MAXW] = {NULL};
  char * read_coords[MAXW] = {NULL};
  char * read_connectivity[MAXW] = {NULL};

  /* Count number of lines */
  int Num_line = 0;
  
  /* Number of words in a text line */
  int nwords;
  int ncoords;
  int nconnectivity;

  /* Variables to fill the Connectivity */
  int NumNodesElem;
  int * ConnectivityElem;
   
  /* Open the mesh file and store in the FILE * variable */
  MeshFile = fopen(MeshName,"r");  
  /* If the file is corrupted, put a wrong message */
  if( MeshFile == NULL )
    {
      perror(MeshName);
    }

  /* Initialize some mesh parameters */
  GID_Mesh.NumNodesMesh = 0; /* Set the number of nodes in the mesh to zero */
  GID_Mesh.NumElemMesh = 0; /* Set the number of elements in the mesh to zero */
    
  /* Read the first line with the header */
  Aux_line = fgets(line, sizeof line, MeshFile);
  if(Aux_line != NULL)
    {
      nwords = parse (words, line, " \n\t");
      Num_line++;
    }
  else
    {
      printf("Error in ReadGidMesh() : during line %i !!! \n",Num_line);
      exit(EXIT_FAILURE);
    }
  
  /* Element properties of the mesh */
  if (strcmp(words[0],"MESH") == 0)
    {
      /* Check if the first word is the keyword MESH */
      if( nwords == 7)
	{
	  /* ElemType */
	  strcpy(GID_Mesh.TypeElem,words[4]);      
	  /* Number of nodes per element */
	  NumNodesElem = atoi(words[6]);
	  ConnectivityElem = (int *)Allocate_Array(NumNodesElem,sizeof(int));
	  /* Shape functions and its derivatives:
	     Here we only pass by reference the function, the output Matrix
	     is allocated insede of N_ref() and dNdX_ref(), we also set the
	     number of dimensions */
	  if((strcmp(GID_Mesh.TypeElem,"Linear") == 0) &&
	     (NumNodesElem == 2))
	    {
	      GID_Mesh.Dimension = 1;
	      GID_Mesh.N_ref = L2;
	      GID_Mesh.dNdX_ref = dL2;
	    }
	  if((strcmp(GID_Mesh.TypeElem,"Quadrilateral") == 0) &&
	     (NumNodesElem == 4))
	    {
	      GID_Mesh.Dimension = 2;
	      GID_Mesh.N_ref = Q4_N;
	      GID_Mesh.dNdX_ref = Q4_dN_Ref;
	    }
	  if( (strcmp(GID_Mesh.TypeElem,"Triangle") == 0) &&
	      (NumNodesElem == 3) ){
	    GID_Mesh.Dimension = 2;
	    GID_Mesh.N_ref = T3;
	    GID_Mesh.dNdX_ref = dT3;
	  }	
	}
      else
	{
	  perror("Error in ReadGidMesh()");
	  printf("Wrong number of mesh properties : %i ! \n",nwords);
	  exit(EXIT_FAILURE);
	}
    }
  
  /* Check the format for the nodes */
  Aux_line = fgets(line, sizeof line, MeshFile);
  if(Aux_line != NULL)
    {
      nwords = parse (words, line, " \n\t");
      if(strcmp(words[0],"Coordinates") != 0)
	{
	  puts("Error in ReadGidMesh() : Format error in the input mesh !!!");
	  exit(EXIT_FAILURE);
	}
      Num_line++;
    }
  else
    {
      printf("Error in ReadGidMesh() : during line %i !!! \n",Num_line);
      exit(EXIT_FAILURE);
    }
  
  /* Count number of nodes */
  Aux_line = fgets(line, sizeof line, MeshFile);
  if(Aux_line != NULL)
    {
      nwords = parse (words, line, " \n\t");
      Num_line++;
    }
  else
    {
      printf("Error in ReadGidMesh() : during line %i !!! \n",Num_line);
      exit(EXIT_FAILURE);
    }

  for( ; strcmp(words[1],"Coordinates") != 0;
       Aux_line = fgets(line, sizeof line, MeshFile),
	 parse (words, line, " \n\t"))
    {
      if(Aux_line != NULL)
	{
	  GID_Mesh.NumNodesMesh++;
	  Num_line++;
	}
      else
	{
	  printf("Error in ReadGidMesh() : during line %i !!! \n",Num_line);
	  exit(EXIT_FAILURE);
	}
    }

  /* White space */
  fgets(line, sizeof line, MeshFile);

  /* Check the format for the elements */
  fgets(line, sizeof line, MeshFile);
  nwords = parse (words, line, " \n\t");
  if(strcmp(words[0],"Elements") != 0){
    puts("Error in ReadGidMesh() : Format error in the input mesh !!!");
    exit(EXIT_FAILURE);
  }

  /* Count number of elements */
  fgets(line, sizeof line, MeshFile);
  nwords = parse (words, line, " \n\t");
  for( ; strcmp(words[1],"Elements") != 0;
       fgets(line, sizeof line, MeshFile),
	 parse (words, line, " \n\t"))
    {
      GID_Mesh.NumElemMesh++;
    }
  
  /* At the end close the mesh file */
  fclose(MeshFile);

  /***************************************************************************/
  /****************************** Allocate data ******************************/
  /***************************************************************************/
  
  /* Allocate the mesh data */
  GID_Mesh.Coordinates = alloc__MatrixLib__(GID_Mesh.NumNodesMesh,GID_Mesh.Dimension);
  GID_Mesh.Connectivity = (ChainPtr *)
    malloc(GID_Mesh.NumElemMesh*sizeof(ChainPtr));
  if(GID_Mesh.Connectivity == NULL)
    {
      puts("Error in Chain declaration");
      exit(EXIT_FAILURE);
    }
  for(int i = 0 ; i<GID_Mesh.NumElemMesh ; i++)
    {
      GID_Mesh.Connectivity[i] = NULL;  
    }

  GID_Mesh.NodeNeighbour = (ChainPtr *)
    malloc(GID_Mesh.NumNodesMesh*sizeof(ChainPtr));
  if(GID_Mesh.NodeNeighbour == NULL)
    {
      puts("Error in Chain declaration");
      exit(EXIT_FAILURE);
    }
  for(int i = 0 ; i<GID_Mesh.NumNodesMesh ; i++)
    {
      GID_Mesh.NodeNeighbour[i] = NULL;  
    }
  
  GID_Mesh.NumNodesElem =  (int *)
    Allocate_ArrayZ(GID_Mesh.NumElemMesh,sizeof(int));

  GID_Mesh.NumParticles = (int *)
    Allocate_ArrayZ(GID_Mesh.NumNodesMesh,sizeof(int));

  GID_Mesh.NumNeighbour =  (int *)
    Allocate_ArrayZ(GID_Mesh.NumNodesMesh,sizeof(int)); 

  /***************************************************************************/
  /***************************************************************************/
  /***************************************************************************/

  
  /* Open the file again and read the mesh */
  MeshFile = fopen(MeshName,"r");
  
  /* Header */
  fgets(line, sizeof line, MeshFile);

  /* Fill the coordinates data */
  fgets(line, sizeof line, MeshFile);
  /* Split the line in word using the space character as a delimiter */
  nwords = parse(words, line, " \n\t");
  if(strcmp(words[0],"Coordinates") == 0)
    {
      for(int i = 0 ; i<GID_Mesh.NumNodesMesh ; i++)
	{
	  fgets(line_coords, sizeof line_coords, MeshFile);
	  ncoords = parse(read_coords, line_coords," \n\t");
	  if(ncoords == 4)
	    {
	      for(int j = 0 ; j<GID_Mesh.Dimension ; j++)
		{
		  GID_Mesh.Coordinates.nM[i][j] = atof(read_coords[j+1]);
		}
	    }
	  else
	    {
	      printf("Check the node : %i \n", atoi(read_coords[0]) - 1);
	      puts("Error in ReadGidMesh() : Format error in the input mesh !!!");
	      exit(EXIT_FAILURE);
	    }
	}
    }
  
  /* End Coordinates */
  fgets(line, sizeof line, MeshFile);
  /* White line */
  fgets(line, sizeof line, MeshFile);
  
  /* Fill the Connectivity data */
  fgets(line, sizeof line, MeshFile);
  /* Split the line in word using the space character as a delimiter */
  nwords = parse(words, line, " \n\t");
  if(strcmp(words[0],"Elements") == 0)
    {
      for(int i = 0 ; i<GID_Mesh.NumElemMesh ; i++)
	{
	  GID_Mesh.NumNodesElem[i] = NumNodesElem;
	  fgets(line_connectivity, sizeof line_connectivity, MeshFile);
	  nconnectivity = parse(read_connectivity, line_connectivity," \n\t");
	  if(nconnectivity == 1 + NumNodesElem)
	    {
	      for(int j = 0 ; j<NumNodesElem ; j++)
		{
		  ConnectivityElem[j] = atoi(read_connectivity[j+1]) - 1;
		}
	      GID_Mesh.Connectivity[i] = Pointer_to_Set(ConnectivityElem,NumNodesElem);
	    }
	  else
	    {
	      printf("Check the element : %i \n", atoi(read_connectivity[0]) - 1);
	      puts("Error in ReadGidMesh() : Format error in the input mesh !!!");
	      exit(EXIT_FAILURE);
	    }
	}
      free(ConnectivityElem);
    }

  /* End Elements */
  fgets(line, sizeof line, MeshFile);

  /* Close mesh file */
  fclose(MeshFile);
  
  /***************************************************************************/
  /***************************************************************************/
  /***************************************************************************/
 
  
  return GID_Mesh;
}
