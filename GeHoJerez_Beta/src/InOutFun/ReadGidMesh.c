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
  /* Some auxiliar variables */
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

  /* Screen message */
  printf("************************************************* \n");
  printf("Begin of set mesh properties !!! \n");

  printf(" * Begin of read mesh file : %s \n",MeshName);
  printf("\t -> Reading mesh ...\n");
   
  /***************************************************************************/ 
  /*     Read line and store in a char array (max 80 characters) and         */ 
  /*   split the line in word using the space character as a delimiter.      */
  /*     In the first line we have some properties of the mesh               */
  /***************************************************************************/

  /* Open the mesh file and store in the FILE * variable */
  MeshFile = fopen(MeshName,"r");  
  /* If the file is corrupted, put a wrong message */
  if( MeshFile == NULL ){
    perror(MeshName);
  }

  /* Initialize some mesh parameters */
  GID_Mesh.NumNodesMesh = 0; /* Set the number of nodes in the mesh to zero */
  GID_Mesh.NumElemMesh = 0; /* Set the number of elements in the mesh to zero */

  /* Set to zero the number of nodes in the boundary of the mesh */
  GID_Mesh.TOP.NumNodes = 0;
  GID_Mesh.BOTTOM.NumNodes = 0;
  GID_Mesh.RIGHT.NumNodes = 0;
  GID_Mesh.LEFT.NumNodes = 0;
  
    
  /* Read the first line with the header */
  Aux_line = fgets(line, sizeof line, MeshFile);
  if(Aux_line != NULL){
    nwords = parse (words, line, " \n\t");
    Num_line++;
  }
  else{
    printf("Error in ReadGidMesh() : during line %i !!! \n",Num_line);
    exit(0);
  }  
  
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
  
  /* Check the format for the nodes */
  Aux_line = fgets(line, sizeof line, MeshFile);
  if(Aux_line != NULL){
    nwords = parse (words, line, " \n\t");
    if(strcmp(words[0],"Coordinates") != 0){
      puts("Error in ReadGidMesh() : Format error in the input mesh !!!");
      exit(0);
    }
    Num_line++;
  }
  else{
    printf("Error in ReadGidMesh() : during line %i !!! \n",Num_line);
    exit(0);
  }
  
  /* Count number of nodes */
  Aux_line = fgets(line, sizeof line, MeshFile);
  if(Aux_line != NULL){
    nwords = parse (words, line, " \n\t");
    Num_line++;
  }
  else{
    printf("Error in ReadGidMesh() : during line %i !!! \n",Num_line);
    exit(0);
  }

  for( ; strcmp(words[1],"Coordinates") != 0;
       Aux_line = fgets(line, sizeof line, MeshFile),
	 parse (words, line, " \n\t")){
    if(Aux_line != NULL){
      GID_Mesh.NumNodesMesh++;
      Num_line++;
    }
    else{
      printf("Error in ReadGidMesh() : during line %i !!! \n",Num_line);
      exit(0);
    }
  }

  /* White space */
  fgets(line, sizeof line, MeshFile);

  /* Check the format for the elements */
  fgets(line, sizeof line, MeshFile);
  nwords = parse (words, line, " \n\t");
  if(strcmp(words[0],"Elements") != 0){
    puts("Error in ReadGidMesh() : Format error in the input mesh !!!");
    exit(0);
  }

  /* Count number of elements */
  fgets(line, sizeof line, MeshFile);
  nwords = parse (words, line, " \n\t");
  for( ; strcmp(words[1],"Elements") != 0;
       fgets(line, sizeof line, MeshFile),
	 parse (words, line, " \n\t")){
    GID_Mesh.NumElemMesh++;
  }
  
  /* At the end close the mesh file */
  fclose(MeshFile);

  /***************************************************************************/
  /****************************** Allocate data ******************************/
  /***************************************************************************/
  
  /* Allocate the mesh data */
  GID_Mesh.Coordinates = MatAlloc(GID_Mesh.NumNodesMesh,
  				  GID_Mesh.Dimension);
  GID_Mesh.Connectivity = (int **)
    Allocate_Matrix(GID_Mesh.NumElemMesh,
  		    GID_Mesh.NumNodesElem,sizeof(int));
  GID_Mesh.ActiveNode = (int *)
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
  if(strcmp(words[0],"Coordinates") == 0){
    for(int i = 0 ; i<GID_Mesh.NumNodesMesh ; i++){
      fgets(line_coords, sizeof line_coords, MeshFile);
      ncoords = parse(read_coords, line_coords," \n\t");
      if(ncoords == 4){
  	for(int j = 0 ; j<GID_Mesh.Dimension ; j++){
  	  GID_Mesh.Coordinates.nM[i][j] = atof(read_coords[j+1]);
  	}
      }
      else{
  	printf("Check the node : %i \n", atoi(read_coords[0]) - 1);
  	puts("Error in ReadGidMesh() : Format error in the input mesh !!!");
  	exit(0);
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
  if(strcmp(words[0],"Elements") == 0){
    for(int i = 0 ; i<GID_Mesh.NumElemMesh ; i++){
      fgets(line_connectivity, sizeof line_connectivity, MeshFile);
      nconnectivity = parse(read_connectivity, line_connectivity," \n\t");
      if(nconnectivity == 1 + GID_Mesh.NumNodesElem){
  	for(int j = 0 ; j<GID_Mesh.NumNodesElem ; j++){
  	  GID_Mesh.Connectivity[i][j] = atoi(read_connectivity[j+1]) - 1;
  	}
      }
      else{
  	printf("Check the element : %i \n", atoi(read_connectivity[0]) - 1);
  	puts("Error in ReadGidMesh() : Format error in the input mesh !!!");
  	exit(0);
      }
    }
  }

  /* End Elements */
  fgets(line, sizeof line, MeshFile);

  /* Close mesh file */
  fclose(MeshFile);
  printf(" * End of read mesh file \n");

  /***************************************************************************/
  /***************************************************************************/
  /***************************************************************************/

  
  /* Find the nodes in the boundary */
  switch(NumberDimensions){
    
  case 1: /******************** 1D mesh ********************/
    /* In a 1D mesh we only have two nodes in the boundary */
    GID_Mesh.LEFT.NumNodes = 1;
    GID_Mesh.RIGHT.NumNodes = 1;
    /* Allocate the size of the array with the nodes */
    GID_Mesh.LEFT.Nodes = (int *)Allocate_ArrayZ(1,sizeof(int));
    GID_Mesh.RIGHT.Nodes = (int *)Allocate_ArrayZ(1,sizeof(int));
    /* Set the nodes of the boundaries */
    GID_Mesh.LEFT.Nodes[0] = 0;
    GID_Mesh.RIGHT.Nodes[0] = 1-GID_Mesh.NumNodesMesh;
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
  	    GID_Mesh.RIGHT.NumNodes ++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][1] == MAX_Y){
  	    GID_Mesh.TOP.NumNodes ++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][0] == MIN_X){
  	    GID_Mesh.LEFT.NumNodes ++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][1] == MIN_Y){
  	    GID_Mesh.BOTTOM.NumNodes ++;
  	  }
	    
  	}
      }
      
      /* Allocate the arrays with the boundary nodes */
      GID_Mesh.TOP.Nodes =
	(int *)Allocate_ArrayZ(GID_Mesh.TOP.NumNodes,sizeof(int));
      GID_Mesh.BOTTOM.Nodes =
	(int *)Allocate_ArrayZ(GID_Mesh.BOTTOM.NumNodes,sizeof(int));
      GID_Mesh.LEFT.Nodes =
	(int *)Allocate_ArrayZ(GID_Mesh.LEFT.NumNodes,sizeof(int));
      GID_Mesh.RIGHT.Nodes =
	(int *)Allocate_ArrayZ(GID_Mesh.RIGHT.NumNodes,sizeof(int));

      /* Fill the arrays  */
      for(int i = 0 ; i<GID_Mesh.NumNodesMesh ; i++){
  	if(NodesBound_aux[i] == 1){

  	  if(GID_Mesh.Coordinates.nM[i][0] == MAX_X){
  	    GID_Mesh.RIGHT.Nodes[aux_RIGHT] = i;
  	    aux_RIGHT++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][1] == MAX_Y){
  	    GID_Mesh.TOP.Nodes[aux_TOP] = i;
  	    aux_TOP++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][0] == MIN_X){
  	    GID_Mesh.LEFT.Nodes[aux_LEFT] = i;
  	    aux_LEFT++;
  	  }
  	  if(GID_Mesh.Coordinates.nM[i][1] == MIN_Y){
  	    GID_Mesh.BOTTOM.Nodes[aux_BOTTOM] = i;
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

  
  return GID_Mesh;
}
