
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

Mesh ReadGidMesh(char * MeshName)
/*

 */  
{
  /* Output element mesh */
  Mesh GID_Mesh;

  /* Pointer to the file */
  FILE * MeshFile;

  /* Variable to save the line of a text */
  char line[MAXC] = {0};
  char line_coords[MAXC] = {0};
  char line_connectivity[MAXC] = {0};
        
  /* Char pointer to store the words in a line */
  char * words[MAXW] = {NULL};
  char * read_coords[MAXW] = {NULL};
  char * read_connectivity[MAXW] = {NULL};
  
  /* Number of words in a text line */
  int nwords;
  int ncoords;
  int nconnectivity;
  int Count_Connectivity = 0;
  int Count_Coordinates = 0;
  /* Some auxiliar variables */
  int Element_i, Nodes_i;
  int Repeat_Nod;
  int NumNodesBound;
  int * NodesBound_aux;
  int aux_RIGHT = 0; 
  int aux_TOP = 0;
  int aux_LEFT = 0;
  int aux_BOTTOM = 0;
  /* Set to zero X and Y min and max values of the mesh */
  double MAX_X = 0;
  double MAX_Y = 0;
  double MIN_X = 0;
  double MIN_Y = 0;   
  /* Parser dictionay stuff */
  ParserDictionary Dict = InitParserDictionary();
  char * delims = Dict.sep[6];

  /* Screen message */
  printf("************************************************* \n");
  printf("Begin of set mesh properties !!! \n");

  printf(" * Begin of read mesh file : %s \n",MeshName);
  printf("\t -> Reading mesh ...\n");
   
  /* 
     Read line and store in a char array (max 80 characters) and 
     split the line in word using the space character as a delimiter.
     In the first line we have some properties of the mesh
  */

  /* Open the mesh file and store in the FILE * variable */
  MeshFile = fopen(MeshName,"r");  
  /* If the file is corrupted, put a wrong message */
  if( MeshFile == NULL ){
    perror(MeshName);
  }

  /* Initialize some mesh parameters */
  GID_Mesh.NumNodesMesh = 0; /* Set the number of nodes in the mesh to zero */
  GID_Mesh.NumElemMesh = 0; /* Set the number of elements in the mesh to zero */
  GID_Mesh.NumTOP = 0; /* Set to zero the number of nodes in the boundary of the mesh */
  GID_Mesh.NumBOTTOM = 0;
  GID_Mesh.NumRIGHT = 0;
  GID_Mesh.NumLEFT= 0;

  /* Read the file line by line */
  while( fgets(line, sizeof line, MeshFile) != NULL ){

    /* Parse line */
    nwords = parse (words, line, " \n\t");
  
    /* Element properties of the mesh */
    if ( strcmp(words[0],"MESH") == 0 ){
      /* Check if the first word is the keyword MESH */
      if( nwords == 7){
	/* ElemType */
	strcpy(GID_Mesh.TypeElem,words[4]);      
	/* Number of nodes per element */
	GID_Mesh.NumNodesElem = atoi(words[6]);
	/* Shape functions and its derivatives:
	   Here we only pass by reference the function, the output Matrix
	   is allocated insede of N_ref() and dNdX_ref(), we also set the
	   number of dimensions */
	if( (strcmp(GID_Mesh.TypeElem,"Linear") == 0) &&
	    (GID_Mesh.NumNodesElem == 2) ){
	  GID_Mesh.Dimension = 1;
	  GID_Mesh.N_ref = L2;
	  GID_Mesh.dNdX_ref = dL2;
	}
	if( (strcmp(GID_Mesh.TypeElem,"Quadrilateral") == 0) &&
	    (GID_Mesh.NumNodesElem == 4) ){
	  GID_Mesh.Dimension = 2;
	  GID_Mesh.N_ref = Q4;
	  GID_Mesh.dNdX_ref = dQ4;
	}
	if( (strcmp(GID_Mesh.TypeElem,"Triangle") == 0) &&
	    (GID_Mesh.NumNodesElem == 3) ){
	  GID_Mesh.Dimension = 2;
	  GID_Mesh.N_ref = T3;
	  GID_Mesh.dNdX_ref = dT3;
	}	
      }
      else{
	perror("Error in ReadGidMesh()");
	printf("Wrong number of mesh properties : %i ! \n",nwords);
	exit(0);
      }
    }

    /* /\* Initit to count the coordinates *\/ */
    /* if(strcmp(words[0],"Coordinates") == 0){ */
    /*   Count_Coordinates = 1; */
    /* } */
    /* /\* Set the number of nodes coordinates of the mesh *\/ */
    /* while(Count_Coordinates){ */
    /*   /\* Read and parse the first line *\/ */
    /*   fgets(line_coords, sizeof line_coords, MeshFile); */
    /*   ncoords = parse (read_coords, line_coords, " \t\n"); */
    /*   printf("%s \n",read_coords[0]); */
    /*   /\* Add number of nodes to the mesh *\/ */
    /*   if(ncoords == 4){ */
    /* 	GID_Mesh.NumNodesMesh++; */
    /*   } */
    /*   /\* Out of the loop *\/ */
      
    /*   if(strcmp(read_coords[0],"End") == 0){ */
    /* 	Count_Coordinates = 0; */
    /* 	printf("paso \n"); */
    /* 	break; */
    /*   } */
    /* } */

    for( strcmp(words[0],"Coordinates");
	 strcmp(words[1],"Coordinates") == 0;
	 fgets(line, sizeof line, MeshFile)){
      printf("%s \n",words[0]);
      nwords = parse(words, line, delims);       
     }

    /* /\* Initit to count the connectivity *\/ */
    /* if(strcmp(words[0],"Elements") == 0){ */
    /*   Count_Connectivity = 1; */
    /* } */
    /* /\* Read the connectivity *\/ */
    /* while(Count_Connectivity){ */
    /*   /\* Read and parse the first line *\/ */
    /*   fgets(line_connectivity, sizeof line_connectivity, MeshFile); */
    /*   nconnectivity = parse(read_connectivity, line_connectivity," \n\t"); */
    /*   /\* Check the format *\/ */
    /*   if(nconnectivity == 1+GID_Mesh.NumNodesElem){ */
    /* 	/\* Add number of elements to the mesh *\/ */
    /* 	GID_Mesh.NumElemMesh++; */
    /*   } */
    /*   /\* Stop counting the connectivity *\/ */
    /*   if(strcmp(read_connectivity[0],"End") == 0){ */
    /* 	Count_Connectivity = 0; */
    /*   } */
    /* } */
    
  }
  
  /* At the end close the mesh file */
  fclose(MeshFile);

  printf("%i %i \n",GID_Mesh.NumNodesMesh,GID_Mesh.NumElemMesh);

  
  exit(0);
 
  /* Allocate the mesh data */
  GID_Mesh.Coordinates = MatAlloc(GID_Mesh.NumNodesMesh,
  				  GID_Mesh.Dimension);
  GID_Mesh.Connectivity = (int **)
    Allocate_Matrix(GID_Mesh.NumElemMesh,
  		    GID_Mesh.NumNodesElem,sizeof(int));
  GID_Mesh.ActiveNode = (int *)
    Allocate_ArrayZ(GID_Mesh.NumNodesMesh,sizeof(int));
  
  /* Open the file again and read the mesh */
  MeshFile = fopen(MeshName,"r");
  
  /* Read line and store in a char array (max 80 characters)*/
  fgets(line, sizeof line, MeshFile);
  /* Split the line in word using the space character as a delimiter */
  nwords = parse(words, line, delims);

  /* Read line and store in a char array (max 80 characters)*/
  fgets(line, sizeof line, MeshFile);
  /* Split the line in word using the space character as a delimiter */
  nwords = parse(words, line, delims);
  if(strcmp(words[0],"Coordinates") == 0){
    for(int i = 0 ; i<GID_Mesh.NumNodesMesh ; i++){
      fgets(line_coords, sizeof line_coords, MeshFile);
      ncoords = parse(read_coords, line_coords," \n\t");
      if(ncoords == 4){
  	Nodes_i = atoi(read_coords[0]) - 1;
  	for(int j = 0 ; j<GID_Mesh.Dimension ; j++){
  	  GID_Mesh.Coordinates.nM[Nodes_i][j] = atof(read_coords[j+1]);
  	}
      }
      else{
  	printf("Check the node : %i \n", atoi(read_coords[0]) - 1);
  	puts("Error in ReadGidMesh() : Format error in the input mesh !!!");
  	exit(0);
      }
    }
  }
  
  /* Read line and store in a char array (max 80 characters)*/
  fgets(line, sizeof line, MeshFile);
  /* Split the line in word using the space character as a delimiter */
  nwords = parse(words, line, delims);

    
  /* Read line and store in a char array (max 80 characters)*/
  fgets(line, sizeof line, MeshFile);
  /* Split the line in word using the space character as a delimiter */
  nwords = parse(words, line, delims);

  fgets(line, sizeof line, MeshFile);
  nwords = parse(words, line, delims);

  for(strcmp(words[0],"Elements") == 0 ;
      strcmp(words[0],"End") != 0 ;
      fgets(line, sizeof line, MeshFile),
  	nwords = parse(words, line, delims) ){
    if(nwords>2){
      Element_i = atoi(words[0]) - 1;
      for(Nodes_i = 0 ; Nodes_i<GID_Mesh.NumNodesElem ; Nodes_i++){
  	GID_Mesh.Connectivity[Element_i][Nodes_i] = atoi(words[Nodes_i+1]) - 1;
  	printf("%i ",atoi(words[Nodes_i+1]) - 1);
      }
      printf("\n");
    }
        
  }
  
  /* Find the nodes in the boundary */
  switch(NumberDimensions){
    
  case 1: /******************** 1D mesh ********************/
    /* In a 1D mesh we only have two nodes in the boundary */
    GID_Mesh.NumLEFT = 1;
    GID_Mesh.NumRIGHT = 1;
    /* Allocate the size of the array with the nodes */
    GID_Mesh.LEFT = (int *)Allocate_ArrayZ(1,sizeof(int));
    GID_Mesh.RIGHT = (int *)Allocate_ArrayZ(1,sizeof(int));
    /* Set the nodes of the boundaries */
    GID_Mesh.LEFT[0] = 0;
    GID_Mesh.RIGHT[0] = 1-GID_Mesh.NumNodesMesh;
    break; /******************** 2D mesh ********************/
    
  case 2: /******************** 2D mesh ********************/

    if(GID_Mesh.NumNodesElem == 4){ /* Quadrilateral elements */
      /* 0º Allocate an array of zeros to assign a 1 to those nodes in the boundary */
      NodesBound_aux = (int *)Allocate_ArrayZ(GID_Mesh.NumNodesMesh,sizeof(int));
      /* 1º Set to zero the number of nodes in the boundary */
      NumNodesBound = 0;
      /* 2º Iterate over the nodes to fin the nodes in the boundary */
      for(int i = 0 ; i<GID_Mesh.NumNodesMesh ; i++){

  	/* 3º Get the max values of the boundary */
  	MAX_X = MAXVAL(MAX_X,GID_Mesh.Coordinates.nM[i][0]);
  	MAX_Y = MAXVAL(MAX_Y,GID_Mesh.Coordinates.nM[i][1]);
  	MIN_X = MINVAL(MIN_X,GID_Mesh.Coordinates.nM[i][0]);
  	MIN_Y = MINVAL(MIN_Y,GID_Mesh.Coordinates.nM[i][1]);
	
  	/* 4º Set the counter to zero */
  	Repeat_Nod = 0;
  	/* 5º Loop over the connectivity mesh */
  	for(int j = 0 ; j<GID_Mesh.NumElemMesh ; j++){
  	  for(int k = 0 ; k<GID_Mesh.NumNodesElem ; k++){
  	    if(GID_Mesh.Connectivity[j][k] == i){
  	      Repeat_Nod++;
  	    }
  	  }
  	}
  	/* 6º Add this element to the boundary */
  	if (Repeat_Nod < 4){
  	  NodesBound_aux[i] = 1;
  	  NumNodesBound++;
  	}
      }
      
      /* 7º Fill the array GID_Mesh.NodesBound */
      for(int i = 0 ; i<GID_Mesh.NumNodesMesh ; i++){
  	if(NodesBound_aux[i] == 1){

  	  if(GID_Mesh.Coordinates.nM[i][0] == MAX_X){
  	    GID_Mesh.NumRIGHT ++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][1] == MAX_Y){
  	    GID_Mesh.NumTOP ++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][0] == MIN_X){
  	    GID_Mesh.NumLEFT ++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][1] == MIN_Y){
  	    GID_Mesh.NumBOTTOM ++;
  	  }
	    
  	}
      }
      
      /* Allocate the arrays with the boundary nodes */
      GID_Mesh.TOP = (int *)Allocate_ArrayZ(GID_Mesh.NumTOP,sizeof(int));
      GID_Mesh.BOTTOM = (int *)Allocate_ArrayZ(GID_Mesh.NumBOTTOM,sizeof(int));
      GID_Mesh.LEFT = (int *)Allocate_ArrayZ(GID_Mesh.NumLEFT,sizeof(int));
      GID_Mesh.RIGHT = (int *)Allocate_ArrayZ(GID_Mesh.NumRIGHT,sizeof(int));

      /* Fill the arrays  */
      for(int i = 0 ; i<GID_Mesh.NumNodesMesh ; i++){
  	if(NodesBound_aux[i] == 1){

  	  if(GID_Mesh.Coordinates.nM[i][0] == MAX_X){
  	    GID_Mesh.RIGHT[aux_RIGHT] = i;
  	    aux_RIGHT++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][1] == MAX_Y){
  	    GID_Mesh.TOP[aux_TOP] = i;
  	    aux_TOP++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][0] == MIN_X){
  	    GID_Mesh.LEFT[aux_LEFT] = i;
  	    aux_LEFT++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][1] == MIN_Y){
  	    GID_Mesh.BOTTOM[aux_BOTTOM] = i;
  	    aux_BOTTOM++;
  	  }
	    
  	}
      }

      
    } /* Quadrilateral elements */
    if(GID_Mesh.NumNodesElem == 3){ /* Triangular elements */
      puts("Error in ReadGidMesh() : Boundary nodes localization for T3 not implemented yet !");
    } /* Triangular elements */
    
    
    free(NodesBound_aux);

    break; /******************** 2D mesh ********************/
    
  case 3: /******************** 3D mesh ********************/
    printf("************************************************* \n");
    puts("Error in ReadGidMesh() : 3D cases not implemented yet !");
    printf("************************************************* \n");
    exit(0);
    break; /******************** 2D mesh ********************/
    
  default :
    printf("************************************************* \n");
    puts("Error in ReadGidMesh() : Wrong number of dimensions !");
    printf("************************************************* \n");
    exit(0);
  }

  /* Close mesh file */
  fclose(MeshFile);
  printf(" * End of read mesh file \n");

  return GID_Mesh;
}
