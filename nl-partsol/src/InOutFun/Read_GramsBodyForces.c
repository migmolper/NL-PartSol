#include "nl-partsol.h"

/*
  Call global variables
*/
int NumberDOF;
char * Formulation;

/**********************************************************************/

Load * GramsBodyForces(char * Name_File, int NumBodyForces, int GPxElement)
{
  /* Define new load case for the body forces */
  Load * B = (Load *)Allocate_Array(NumBodyForces, sizeof(Load));
  
  /* Simulation file */
  FILE * Sim_dat;

  /* Number of words */
  int Num_words_line;

  /* Parse lines of GramsBodyForces */
  char Line_GramsBodyForces[MAXC] = {0};
  char * Parse_GramsBodyForces[MAXW] = {NULL};
  char * Parse_Nodes[MAXW] = {NULL};

  /* Boundaries iterator */
  int IndexLoad = 0;

  /* Parse file name with the list of nodes */
  char Route_Nodes[MAXC] = {0};
  char FileNodesRoute[MAXC];

  /* Array */
  int Num_Nodes;
  ChainPtr Chain_Nodes = NULL;
  int * Array_Nodes;
  
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
    exit(EXIT_FAILURE);
  }

  /* Generate route */
  generate_route(Route_Nodes,Name_File);

  /* Read GramsBodyForces line  */
  while(fgets(Line_GramsBodyForces,sizeof(Line_GramsBodyForces),Sim_dat) != NULL){

    /* Read the line with the space as separators */
    Num_words_line = parse (Parse_GramsBodyForces, Line_GramsBodyForces," \n\t");
    if (Num_words_line < 0){
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsBodyForces ()",
	      "Parser failed");
      exit(EXIT_FAILURE);
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
	exit(EXIT_FAILURE);
      }
      
      /* Read file with the nodes */
      sprintf(FileNodesRoute,"%s%s",Route_Nodes,Parse_Nodes[1]);

      /* Get an array with the nodes */
      Chain_Nodes = File2Chain(FileNodesRoute);
      Num_Nodes = get_Lenght_Set(Chain_Nodes);
      Array_Nodes = Set_to_Pointer(Chain_Nodes,Num_Nodes);
      free_Set(&Chain_Nodes);

      /* Fill GPs */
      B[IndexLoad].NumNodes = Num_Nodes*GPxElement; 
      B[IndexLoad].Nodes =
	(int *)Allocate_ArrayZ(B[IndexLoad].NumNodes,sizeof(int));
      for(int i = 0 ; i<Num_Nodes ; i++){
	for(int j = 0 ; j<GPxElement ; j++){
	  B[IndexLoad].Nodes[i*GPxElement+j] = Array_Nodes[i]*GPxElement+j;
	}
      }

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
	  exit(EXIT_FAILURE);
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
		exit(EXIT_FAILURE);
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

