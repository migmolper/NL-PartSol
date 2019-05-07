
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

Element ReadGidMesh(char * MeshName)
/*

 */  
{
  /* Output element mesh */
  Element Elem;
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
  /* Parser dictionay stuff */
  ParserDictionary Dict = InitParserDictionary();
  char * delims = Dict.sep[6];

  /* Screen message */
  printf("Begin of read mesh file : %s \n",MeshName);

  /* Open the mesh file and store in the FILE * variable */
  MeshFile = fopen(MeshName,"r");
  
  /* If the file is corrupted, put a wrong message */
  if( MeshFile == NULL ){
    perror(MeshName);
  }

  /* Set the number of nodes in the mesh to zero */
  Elem.NumNodesMesh = 0;
  /* Set the number of elements in the mesh to zero */
  Elem.NumElemMesh = 0;
  
  /* 
     Read line and store in a char array (max 80 characters) and 
     split the line in word using the space character as a delimiter.
     In the first line we have some properties of the mesh
  */
  fgets(line, sizeof line, MeshFile);
  nwords = parse (words, line, delims);

  /* Element properties of the mesh */
  if ( strcmp(words[0],"MESH") == 0 ){
    if( nwords == 7){
      /* ElemType */
      strcpy(Elem.TypeElem,words[4]);
      /* Number of nodes per element */
      Elem.NumNodesElem = atoi(words[6]); 
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
  if( (strcmp(Elem.TypeElem,"Linear") == 0) &&
      (Elem.NumNodesElem == 2) ){
    Elem.Dimension = 1;
    Elem.N_ref = L2;
    Elem.dNdX_ref = dL2;
  }
  if( (strcmp(Elem.TypeElem,"Linear") == 0) &&
      (Elem.NumNodesElem == 4) ){
    Elem.Dimension = 2;
    Elem.N_ref = Q4;
    Elem.dNdX_ref = dQ4;
  }

  fgets(line, sizeof line, MeshFile);
  nwords = parse(words, line, delims);
  /* Set the number of nodes coordinates of the mesh */
  for(strcmp(words[0],"Coordinates") == 0 ;
      strcmp(words[0],"End") != 0 ;
      fgets(line, sizeof line, MeshFile),
	nwords = parse(words, line,delims) ){
    ++Elem.NumNodesMesh;
  }
  
  /* Correct the number of nodes in the mesh : The prevoius algorithim counts 
     also the header ("Coordinates") and the footer ("End Coordiantes") */
  Elem.NumNodesMesh -= 1;

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
    ++Elem.NumElemMesh;
  }
  /* Correct the number of elements in the mesh : The prevoius algorithim counts 
     also the header ("Elements") and the footer ("End Elements") */
  Elem.NumElemMesh -= 1;
      
  /* At the end close the mesh file */
  fclose(MeshFile);
 
  /* Allocate the mesh data */
  Elem.Coordinates = (double **)
    Allocate_Matrix(Elem.NumNodesMesh,
		    3,sizeof(double));
  Elem.Connectivity = (int **)
    Allocate_Matrix(Elem.NumElemMesh,
		    Elem.NumNodesElem,sizeof(int));
  Elem.ActiveElem = (int *)
    Allocate_ArrayZ(Elem.NumElemMesh,sizeof(int));
  
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
      Elem.Coordinates[Nodes_i][0] = atof(words[1]); 
      Elem.Coordinates[Nodes_i][1] = atof(words[2]);
      Elem.Coordinates[Nodes_i][2] = atof(words[3]);
    }
  }

  /* Read line and store in a char array (max 80 characters)*/
  fgets(line, sizeof line, MeshFile);
  /* Split the line in word using the space character as a delimiter */
  nwords = parse(words, line, delims);

  fgets(line, sizeof line, MeshFile);
  nwords = parse(words, line, delims);

  printf("%s \n",words[0]);
  for(strcmp(words[0],"Elements") == 0 ;
      strcmp(words[0],"End") != 0 ;
      fgets(line, sizeof line, MeshFile),
	nwords = parse(words, line, delims) ){
    if(nwords>2){
      Element_i = atoi(words[0]) - 1;
      for(Nodes_i = 0 ; Nodes_i<Elem.NumNodesElem ; Nodes_i++){
	Elem.Connectivity[Element_i][Nodes_i] = atoi(words[Nodes_i+1]);
      }
    }
        
  }
  
  fclose(MeshFile);

  puts("End of read mesh");

  return Elem;
  
}
