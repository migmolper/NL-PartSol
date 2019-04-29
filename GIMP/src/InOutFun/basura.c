
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"


void ReadGidMesh(char * MeshName,
		 ElementProperties ElemProp,
		 MeshProperties MeshProp,
		 ParserDictionary Dict)

/*

*/
  
{
  
  FILE * MeshFile;
  char line[80];
  int nwords;
  int aux;
  char * words[20];
  
  MeshFile  = fopen(MeshName,"r");
  
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
  nwords = getwords(line, words, Dict.ascii_sep[6], 20);
    
  
  /* Mesh properties */
  if ( strcmp(words[0],"MESH") == 0 ){
    if( nwords == 7){
      ElemProp.Dimension = atoi(words[2]);
      //	  ElemProp.Type = words[4];
      ElementPrElemPropoperties.Nnodes = atoi(words[6]);
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
	nwords = getwords(line, words, Dict.ascii_sep[6], 20)
      ){
    ++MeshProp.Nnodes;
  }
  /* Correct the number of nodes in the mesh : The prevoius algorithim counts 
     also the header ("Coordinates") and the footer ("End Coordiantes") */
  MeshProperties.Nnodes -= 2;
  
  /* Continuing reading the mesh */      
  fgets(line, sizeof line, MeshFile);
  nwords = getwords(line, words, Dict.ascii_sep[6], 20);
  
  /* Set the number of elements of the mesh */
  for(strcmp(words[0],"Elements") == 0 ;
      strcmp(words[0],"End") != 0 ;
      fgets(line, sizeof line, MeshFile),
	nwords = getwords(line, words, Dict.ascii_sep[6], 20) ){
    ++MeshProp.Nelem;
  }
  /* Correct the number of elements in the mesh : The prevoius algorithim counts 
     also the header ("Elements") and the footer ("End Elements") */
  MeshProperties.Nelem -= 2;
  
    
  /* At the end close the mesh file */
  fclose(MeshFile);
  
  }
/* If the file is corrupted, put a wrong message */
 else{
   perror("Error in ReadGidMesh() : Wrong mesh file !");
 }


/* Open the file again and read the mesh */
MeshFile = fopen(MeshName,"r");

if( MeshFile != NULL ){
  
  char line[80];
  int nwords;
  char * words[20];
  
  while ( feof(MeshFile) == 0 ){

    /* Read line and store in a char array (max 80 characters)*/
    fgets(line, sizeof line, MeshFile);
    /* Split the line in word using the space character as a delimiter */
    nwords = getwords(line, words,Dict.ascii_sep[6], 20);


    for(strcmp(words[0],"Coordinates") == 0 ;
	strcmp(words[0],"End") != 0 ;
	fgets(line, sizeof line, MeshFile),
	  nwords = getwords(line, words, Dict.ascii_sep[6], 20)
	){
      if(nwords == 4){	  
	printf("Nodo [%i] : [%f ; %f ; %f] \n",
	       atoi(words[0]),
	       atof(words[1]),
	       atof(words[2]),
	       atof(words[3]));	  
      }
    }

    /* Read line and store in a char array (max 80 characters)*/
    fgets(line, sizeof line, MeshFile);
    /* Split the line in word using the space character as a delimiter */
    nwords = getwords(line, words,ascii_sep[6], 20);
    
    for(strcmp(words[0],"Elements") == 0 ;
	strcmp(words[0],"End") != 0 ;
	fgets(line, sizeof line, MeshFile),
	  nwords = getwords(line, words,Dict.ascii_sep[6], 20) ){
      if(nwords == 10){
	printf("Nodo [%i] : ",atoi(words[0]));
	for(int i = 1 ; i<=ElementProperties.Nnodes ; i++){
	  printf(" %i ",atoi(words[i]));
	}
      }
      printf("\n");
    }

      
  }
  
  fclose(MeshFile);

  puts("End of read mesh");
  
 }
