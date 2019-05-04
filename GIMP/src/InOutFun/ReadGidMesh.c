
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

void ReadGidMesh(char * MeshName)
/*

 */  
{
  
  FILE * MeshFile;
  int nwords;
  int aux;
  int Element_i,Nodes_i;  
  char line[MAXC] = {0};
  char * words[MAXW] = {NULL};
  ParserDictionary Dict = InitParserDictionary();
  char * delims = Dict.sep[6];

  printf("Begin of read mesh file : %s \n",MeshName);
  
  MeshFile = fopen(MeshName,"r");
  
  /* If the file is corrupted, put a wrong message */
  if( MeshFile == NULL ){
    perror(MeshName);
  }
  
  MeshProp.Nnodes = 0;
  MeshProp.Nelem = 0;
  
  /* 
     Read line and store in a char array (max 80 characters) and 
     split the line in word using the space character as a delimiter.
     In the first line we have some properties of the mesh
  */

  fgets(line, sizeof line, MeshFile);
  nwords = parse (words, line, delims);

  /* Mesh properties */
  if ( strcmp(words[0],"MESH") == 0 ){
    if( nwords == 7){
      ElemProp.Dimension = atoi(words[2]); /* Number of dimensions */
      strcpy(ElemProp.Type,words[4]); /* ElemType */
      ElemProp.Nnodes = atoi(words[6]); /* Number of nodes per element */
    }
    else{
      perror("Error in ReadGidMesh()");
      printf("Wrong number of mesh properties : %i ! \n",nwords);
      exit(0);
    }
  } /* end if */
  
  /* Set the number of nodes coordinates of the mesh */
  for(strcmp(words[0],"Coordinates") == 0 ;
      strcmp(words[0],"End") != 0 ;
      fgets(line, sizeof line, MeshFile),
	nwords = parse(words, line,delims) ){
    ++MeshProp.Nnodes;
  }
  /* Correct the number of nodes in the mesh : The prevoius algorithim counts 
     also the header ("Coordinates") and the footer ("End Coordiantes") */
  MeshProp.Nnodes -= 2;

  /* Continuing reading the mesh */
  fgets(line, sizeof line, MeshFile);
  nwords = parse(words, line, delims);
  
  /* Set the number of elements of the mesh */
  for(strcmp(words[0],"Elements") == 0 ;
      strcmp(words[0],"End") != 0 ;
      fgets(line, sizeof line, MeshFile),
	nwords = parse (words, line, delims) ){
    ++MeshProp.Nelem;
  }
  /* Correct the number of elements in the mesh : The prevoius algorithim counts 
     also the header ("Elements") and the footer ("End Elements") */
  MeshProp.Nelem -= 2;
      
  /* At the end close the mesh file */
  fclose(MeshFile);
 
  /* Allocate the mesh data */
  MeshProp.Coordinates = (double **)
    Allocate_Matrix(MeshProp.Nnodes,3,sizeof(double));
  MeshProp.Connectivity =
    (int **)Allocate_Matrix(MeshProp.Nelem,ElemProp.Nnodes,sizeof(int));
  
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
      MeshProp.Coordinates[Nodes_i][0] = atof(words[1]); 
      MeshProp.Coordinates[Nodes_i][1] = atof(words[2]);
      MeshProp.Coordinates[Nodes_i][2] = atof(words[3]);
    }
  }

  /* Read line and store in a char array (max 80 characters)*/
  fgets(line, sizeof line, MeshFile);
  /* Split the line in word using the space character as a delimiter */
  nwords = parse(words, line, delims);
  
  for(strcmp(words[0],"Elements") == 0 ;
      strcmp(words[0],"End") != 0 ;
      fgets(line, sizeof line, MeshFile),
	nwords = parse(words, line, delims) ){
    if(nwords == (ElemProp.Nnodes + 1)){
      Element_i = atoi(words[0]) - 1;
      for(Nodes_i = 0 ; Nodes_i<ElemProp.Nnodes ; Nodes_i++){
	MeshProp.Connectivity[Element_i][Nodes_i] = atoi(words[Nodes_i+1]);
      }
    }
  }
  
  fclose(MeshFile);

  puts("End of read mesh");
  
}
