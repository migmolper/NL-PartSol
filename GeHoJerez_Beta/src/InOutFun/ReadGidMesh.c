
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

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
  /* Char pointer to store the words in a line */
  char * words[MAXW] = {NULL};
  /* Number of words in a text line */
  int nwords;
  /* Some auxiliar variables */
  int aux, Element_i, Nodes_i;
  int Repeat_Nod;
  int * NodesBound_aux;
  /* Parser dictionay stuff */
  ParserDictionary Dict = InitParserDictionary();
  char * delims = Dict.sep[6];

  /* Screen message */
  printf("************************************************* \n");
  printf("Begin of set mesh properties !!! \n");

  printf(" * Begin of read mesh file : %s \n",MeshName);
  printf("\t -> Reading mesh ...\n");

  /* Open the mesh file and store in the FILE * variable */
  MeshFile = fopen(MeshName,"r");
  
  /* If the file is corrupted, put a wrong message */
  if( MeshFile == NULL ){
    perror(MeshName);
  }

  /* Set the number of nodes in the mesh to zero */
  GID_Mesh.NumNodesMesh = 0;
  /* Set the number of elements in the mesh to zero */
  GID_Mesh.NumElemMesh = 0;
  /* Set to zero the number of nodes in the boundary of the mesh */
  GID_Mesh.NumNodesBound = 0;
  
  /* 
     Read line and store in a char array (max 80 characters) and 
     split the line in word using the space character as a delimiter.
     In the first line we have some properties of the mesh
  */
  fgets(line, sizeof line, MeshFile);
  nwords = parse (words, line, " \n");

  /* Element properties of the mesh */
  if ( strcmp(words[0],"MESH") == 0 ){
    /* Check if the first word is the keyword MESH */
    if( nwords == 7){
      /* ElemType */
      strcpy(GID_Mesh.TypeElem,words[4]);      
      /* Number of nodes per element */
      GID_Mesh.NumNodesElem = atoi(words[6]);      
    }
    else{
      perror("Error in ReadGidMesh()");
      printf("Wrong number of mesh properties : %i ! \n",nwords);
      exit(0);
    }
  } /* end if */

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

  fgets(line, sizeof line, MeshFile);
  nwords = parse(words, line, delims);
  /* Set the number of nodes coordinates of the mesh */
  for(strcmp(words[0],"Coordinates") == 0 ;
      strcmp(words[0],"End") != 0 ;
      fgets(line, sizeof line, MeshFile),
	nwords = parse(words, line,delims) ){
    ++GID_Mesh.NumNodesMesh;
  }
  
  /* Correct the number of nodes in the mesh : The prevoius algorithim counts 
     also the header ("Coordinates") and the footer ("End Coordiantes") */
  GID_Mesh.NumNodesMesh -= 1;

  /* Continuing reading the mesh */
  fgets(line, sizeof line, MeshFile);
  nwords = parse(words, line, delims);
  fgets(line, sizeof line, MeshFile);
  nwords = parse(words, line, delims);

  /* Set the number of elements of the mesh */
  for(strcmp(words[0],"Elements") == 0 ;
      strcmp(words[0],"End") != 0 ;
      fgets(line, sizeof line, MeshFile),
	nwords = parse (words, line, delims) ){
    ++GID_Mesh.NumElemMesh;
  }
  /* Correct the number of elements in the mesh : The prevoius algorithim counts 
     also the header ("Elements") and the footer ("End Elements") */
  GID_Mesh.NumElemMesh -= 1;
      
  /* At the end close the mesh file */
  fclose(MeshFile);
 
  /* Allocate the mesh data */
  GID_Mesh.Coordinates = MatAlloc(GID_Mesh.NumNodesMesh,GID_Mesh.Dimension);
  /* GID_Mesh.Coordinates =  (double **) Allocate_Matrix(GID_Mesh.NumNodesMesh,3,sizeof(double)); */  
  GID_Mesh.Connectivity = (int **)
    Allocate_Matrix(GID_Mesh.NumElemMesh,
		    GID_Mesh.NumNodesElem,sizeof(int));
  GID_Mesh.ActiveElem = (int *)
    Allocate_ArrayZ(GID_Mesh.NumElemMesh,sizeof(int));
  
  /* Open the file again and read the mesh */
  MeshFile = fopen(MeshName,"r");
  
  /* Read line and store in a char array (max 80 characters)*/
  fgets(line, sizeof line, MeshFile);
  /* Split the line in word using the space character as a delimiter */
  nwords = parse(words, line, delims);

  for(strcmp(words[0],"Coordinates") == 0 ;
      strcmp(words[0],"End") != 0 ;
      fgets(line, sizeof line, MeshFile),
	nwords = parse(words, line, delims) ){
    if(nwords == 4){
      Nodes_i = atoi(words[0]) - 1;
      for(int i = 0 ; i<GID_Mesh.Dimension ; i++){
	GID_Mesh.Coordinates.nM[Nodes_i][i] = atof(words[1+i]);
      }
    }
  }

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
      }
    }
        
  }
  
  /* Find the nodes in the boundary */
  switch(NumberDimensions){
  case 1:/* 1D mesh */
    /* In a 1D mesh we only have two nodes in the boundary */
    GID_Mesh.NumNodesBound = 2;
    /* Allocate the size of the array with the nodes */
    GID_Mesh.NodesBound = (int *)Allocate_Array(2,sizeof(int));
    /* Set the boundary nodes */
    GID_Mesh.NodesBound[0] = 0;
    GID_Mesh.NodesBound[1] = GID_Mesh.NumNodesMesh-1;
    break;
    
  case 2: /* 2D mesh (quadrangular/triangular) */
    NodesBound_aux = (int *)Allocate_Array(GID_Mesh.NumNodesMesh,sizeof(int));
    GID_Mesh.NumNodesBound = 0;
    
    for(int i = 0 ; i<GID_Mesh.NumNodesMesh ; i++){

      /* Set the counter to zero */
      Repeat_Nod = 0;
      
      /* Loop over the connectivity mesh */
      for(int j = 0 ; j<GID_Mesh.NumElemMesh ; j++){
	for(int k = 0 ; k<GID_Mesh.NumNodesElem ; k++){
	  if(GID_Mesh.Connectivity[j][k] == i){
	    Repeat_Nod++;
	  }
	}
      }
      /* Add this elemente to the boundary */
      if (Repeat_Nod < 4){
	NodesBound_aux[i] = 1;
	GID_Mesh.NumNodesBound++;
      }
      
    }
    
    /* Allocate the array with in the index of the nodal boundaries */
    GID_Mesh.NodesBound = (int *)Allocate_Array(GID_Mesh.NumNodesBound,sizeof(int));

    /* Fill the array GID_Mesh.NodesBound */
    aux = 0;
    for(int i = 0 ; i<GID_Mesh.NumNodesMesh ; i++){
      if(NodesBound_aux[i] == 1){
	GID_Mesh.NodesBound[aux] = i;
	aux++;
      }
    }
    
    free(NodesBound_aux);
    
    break;
  case 3:
    break;
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
