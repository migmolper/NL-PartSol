#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))


/**********************************************************************/

Load * GramsBodyForces(char * Name_File,GaussPoint GP_Mesh)
/*
  GramsBodyForces (Nodes=ListNodes.txt) {
  b.x Load_x.txt
  b.y Load_x.txt
  }
*/
{
  /* Define new load case for the body forces */
  Load * B = (Load *)Allocate_Array(GP_Mesh.NumberBodyForces,
				    sizeof(Load));
  
  /* Simulation file */
  FILE * Sim_dat;

  /* Number of words */
  int Num_words_line;
  int Num_words_route;

  /* Parse lines of GramsBodyForces */
  char Line_GramsBodyForces[MAXC] = {0};
  char * Parse_GramsBodyForces[MAXW] = {NULL};
  char * Parse_Nodes[MAXW] = {NULL};

  /* Boundaries iterator */
  int IndexLoad = 0;

  /* Parse file name with the list of nodes */
  char * Name_File_Copy = malloc(strlen(Name_File)); 
  char * Name_Parse[MAXW] = {NULL};
  char Route_Nodes[MAXC] = {0};
  char FileNodesRoute[MAXC];
  ChainPtr Chain_Nodes = NULL;

  /* Parse GramsBodyForces properties */
  char Line_Properties[MAXC] = {0};
  char * Parse_Properties[MAXW] = {NULL};
  char FileLoadRoute[MAXC];

  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsBodyForces()",
	    "Incorrect lecture of",
	    Name_File);
    exit(0);
  }

  /* Generate route */
  strcpy(Name_File_Copy, Name_File);
  Num_words_route = parse(Name_Parse,Name_File_Copy,"(/)");
  strcat(Route_Nodes,"./");
  for(int i = 0 ; i<Num_words_route-1 ; i++){
    strcat(Route_Nodes, Name_Parse[i]);
    strcat(Route_Nodes,"/");
  }

  /* Read GramsBodyForces line  */
  while( fgets(Line_GramsBodyForces, sizeof(Line_GramsBodyForces), Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    Num_words_line = parse (Parse_GramsBodyForces, Line_GramsBodyForces," \n\t");
    if (Num_words_line < 0){
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsBodyForces ()",
	      "Parser failed");
      exit(0);
    }

    /* Find GramsBodyForces line */
    if ((Num_words_line > 0) &&
    	(strcmp(Parse_GramsBodyForces[0],"GramsBodyForces") == 0 ) &&
	((strcmp(Parse_GramsBodyForces[2],"{") == 0))){

      /* File with the nodes */
      Num_words_line = parse (Parse_Nodes, Parse_GramsBodyForces[1],"(=)");
      if( (Num_words_line != 2) ||
	  (strcmp(Parse_Nodes[0],"Nodes") != 0)){
	fprintf(stderr,"%s : %s \n",
		"Error in GramsBodyForces()",
		"Use this format -> (Nodes=str) !!!");
	exit(0);
      }
      
      /* Read file with the nodes */
      sprintf(FileNodesRoute,"%s%s",Route_Nodes,Parse_Nodes[1]);

      /* Get an array with the nodes */
      Chain_Nodes = File2Chain(FileNodesRoute);
      B[IndexLoad].NumNodes = get_Lenght_Set(Chain_Nodes);
      B[IndexLoad].Nodes = Set_to_Pointer(Chain_Nodes,B[IndexLoad].NumNodes);
      free_Set(Chain_Nodes);

      /* Number of dimensions of the BCC */
      B[IndexLoad].Dim = NumberDOF;
      /* Direction of the BCC */
      B[IndexLoad].Dir =
	(int *)Allocate_ArrayZ(NumberDOF,sizeof(int));
      /* Curve for each dimension */
      B[IndexLoad].Value =
	(Curve *)Allocate_Array(NumberDOF,sizeof(Curve));
      /* Information of the BCC */
      if(strcmp(Formulation,"-V") == 0){
  	/* Name of the BCC */
  	strcpy(B[IndexLoad].Info,"Velocity");
      }

      /* Read properties line  */
      while( fgets(Line_Properties, sizeof(Line_Properties), Sim_dat) != NULL ){

	/* Read the line with the space as separators */
	Num_words_line = parse (Parse_Properties, Line_Properties," \n\t");
	if (Num_words_line < 0){
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsBodyForces ()",
		  "Parser failed");
	  exit(0);
	}
	if((Num_words_line > 0) &&
	   (strcmp(Parse_Properties[0],"}") != 0 )){
	    if(strcmp(Formulation,"-V") == 0){ 
	      if(strcmp(Parse_Properties[0],"b.x") == 0){
		/* Fill the direction of the BCC */
		B[IndexLoad].Dir[0] = 1;
		/* Read the curve to impose the boundary condition */
		if(strcmp(Parse_Properties[1],"NULL") != 0){
		  sprintf(FileLoadRoute,"%s%s",Route_Nodes,Parse_Properties[1]);
		  B[IndexLoad].Value[0] =
		    ReadCurve(FileLoadRoute);
		  puts("*************************************************");
		  printf(" \t %s (%i) : \n \t %s %s \n",
			 "* Body BC",
			 B[IndexLoad].NumNodes,
			 Parse_Properties[1],FileLoadRoute);
		}
		else{
		  B[IndexLoad].Dir[0] = 0;
		}
	      }
	      else if(strcmp(Parse_Properties[0],"b.y") == 0){
		/* Fill the direction of the BCC */
		B[IndexLoad].Dir[1] = 1;
		/* Read the curve to impose the boundary condition */
		if(strcmp(Parse_Properties[1],"NULL") != 0){
		  sprintf(FileLoadRoute,"%s%s",Route_Nodes,Parse_Properties[1]);
		  B[IndexLoad].Value[1] = ReadCurve(FileLoadRoute);
		  puts("*************************************************");
		  printf(" \t %s (%i) : \n \t %s %s \n",
			 "* Body BC",
			 B[IndexLoad].NumNodes,
			 Parse_Properties[1],FileLoadRoute);
		}
		else{
		  B[IndexLoad].Dir[1] = 0;
		}
	      }
	      else if(strcmp(Parse_Properties[0],"b.z") == 0){
		/* Fill the direction of the BCC */
		B[IndexLoad].Dir[2] = 1;
		/* Read the curve to impose the boundary condition */
		if(strcmp(Parse_Properties[1],"NULL") != 0){
		  sprintf(FileLoadRoute,"%s%s",Route_Nodes,Parse_Properties[1]);
		  B[IndexLoad].Value[2] = ReadCurve(FileLoadRoute);
		  puts("*************************************************");
		  printf(" \t %s (%i) : \n \t %s %s \n",
			 "* Body BC",
			 B[IndexLoad].NumNodes,
			 Parse_Properties[1],FileLoadRoute);
		}
		else{
		  B[IndexLoad].Dir[2] = 0;
		}
	      }
	      else{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsBodyForces()", "Use -> b.x, b.y, b.z");
		exit(0);
	      }
	    }
	}
	else{
	  break;
	}
	
      }
      
      /* Increment the index */
      IndexLoad++;
    }

  }
  return B;
}

/**********************************************************************/

