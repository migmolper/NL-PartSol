#include "nl-partsol.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))


void GramsInitials(char * Name_File, GaussPoint GP_Mesh, int GPxElement)
/*
  GramsInitials (Nodes=ListInit.txt) {
  Value=[5.0,0.0,0.0]
  }
*/
{

  /* Number of dimensions */
  int Ndim = NumberDimensions;
  
  /* Simulation file */
  FILE * Sim_dat;
  
  /* Parser num chars */
  int Num_words_parse;

  /* Special variables GramsInitials */
  char Line_GramsInitials[MAXC] = {0}; 
  char * Parse_GramsInitials[MAXW] = {NULL};

  /* Parse file name with the list of nodes */
  char * Parse_Init_Nodes[MAXW] = {NULL};
  char * Name_File_Copy = malloc(strlen(Name_File)); 
  char Route_Nodes[MAXC] = {0};
  char FileNodesRoute[MAXC];

  /* Array */
  int Num_Nodes;
  ChainPtr Chain_Nodes = NULL;
  int * Array_Nodes;

  /* Parse priperties of the initial condition function */
  char Line_Init_Prop[MAXC] = {0}; 
  char * Parse_Init_Prop[MAXW] = {NULL};

  /* */
  char * IC_value[MAXW] = {NULL};

  /* Auxiliar variable for status */
  char * STATUS_LINE;

  /* Initial message */  
  puts("*************************************************");
  printf(" \t %s : \n\t %s \n",
	 "* Read initial conditions in ",
	 Name_File);
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsInitials()",
	    "Incorrect lecture of",
	    Name_File);
    exit(EXIT_FAILURE);
  }

  /* Generate route */
  generate_route(Route_Nodes,Name_File);

  /* Read the file line by line */
  while(fgets(Line_GramsInitials,sizeof(Line_GramsInitials),Sim_dat) != NULL){

    /* Read the line with the space as separators */
    Num_words_parse = parse(Parse_GramsInitials,Line_GramsInitials," \n\t");
    if (Num_words_parse < 0){
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsInitials()",
	      "Parser failed");
      exit(EXIT_FAILURE);
    }

    if ((Num_words_parse > 0) &&
	(strcmp(Parse_GramsInitials[0],"GramsInitials") == 0)){

      /* Read temporal integrator scheme */
      Num_words_parse = parse(Parse_Init_Nodes,Parse_GramsInitials[1],"(=)");
      if( (Num_words_parse != 2) ||
	  (strcmp(Parse_Init_Nodes[0],"Nodes") != 0)){
	fprintf(stderr,"%s : %s \n",
		"Error in GramsInitials()",
		"Use this format -> (Nodes=str) !!!");
	exit(EXIT_FAILURE);
      }

      /* Read file with the nodes */
      sprintf(FileNodesRoute,"%s%s",Route_Nodes,Parse_Init_Nodes[1]);

      /* Get an array with the nodes */
      Chain_Nodes = File2Chain(FileNodesRoute);
      Num_Nodes = get_Lenght_Set(Chain_Nodes);
      Array_Nodes = Set_to_Pointer(Chain_Nodes,Num_Nodes);
      free_Set(Chain_Nodes);
	
      /* Look for the curly brace { */
      if(strcmp(Parse_GramsInitials[2],"{") == 0){
	/* Initial line */
	STATUS_LINE = fgets(Line_Init_Prop,
			    sizeof(Line_Init_Prop),
			    Sim_dat);
	if(STATUS_LINE == NULL){
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsInitials()",
		  "Unspected EOF !!!");
	  exit(EXIT_FAILURE);	
	}
	Num_words_parse = parse(Parse_Init_Prop,Line_Init_Prop," =\t\n");
	if(strcmp(Parse_Init_Prop[0],"}") == 0){
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsInitials()",
		  "No initial value defined !!!");
	  exit(EXIT_FAILURE);
	}
	while(STATUS_LINE != NULL){
	  
	  if(Num_words_parse != 2){
	    fprintf(stderr,"%s : %s \n",
		    "Error in GramsInitials()",
		    "Use this format -> Value = [,] !!!");
	    exit(EXIT_FAILURE);
	  }

	  /* Fill the initial conditions */
 	  if(strcmp(Parse_Init_Prop[0],"Value") == 0){

	    printf("\t -> %s : %s \n",
		   "For nodes in",FileNodesRoute);
	    printf("\t    %s : %s = %s \n",
		   "Initial conditions",
		   Formulation,Parse_Init_Prop[1]);
	    
	    Num_words_parse = parse(IC_value,Parse_Init_Prop[1],"[,]");

	    /* Check the number of dimensions */
	    if(Num_words_parse != Ndim){
	      fprintf(stderr,"%s : %s %i-D %s\n",
		      "Error in GramsInitials()",
		      "This is a",Ndim,"problem");
	      exit(EXIT_FAILURE);
	    }

	    /* Fill it IC for velocity formulation */
	    if(strcmp(Formulation,"-V") == 0){
	      for(int i = 0 ; i<Num_Nodes ; i++){
		for(int j = 0 ; j<GPxElement ; j++){
		  for(int k = 0 ; k<Ndim ; k++){
		    GP_Mesh.Phi.vel.nM[Array_Nodes[i]*GPxElement+j][k] =
		    atof(IC_value[k]);
		  }
		}
	      }
	    }
	    
	  }
	  else{
	    fprintf(stderr,"%s : %s %s \n",
		    "Error in GramsInitials()",
		    "Undefined",Parse_Init_Prop[0]);
	    exit(EXIT_FAILURE);
	  }
	  /* Read next line and check */
	  STATUS_LINE = fgets(Line_Init_Prop,
			      sizeof(Line_Init_Prop),
			      Sim_dat);
	  Num_words_parse = parse(Parse_Init_Prop,Line_Init_Prop," =\t\n");
	  if(strcmp(Parse_Init_Prop[0],"}") == 0){
	    break;
	  }
	}
	if(STATUS_LINE == NULL){
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsInitials()",
		  "you forget to put a } !!!");
	  exit(EXIT_FAILURE);	  
	}
      }
      else{
	fprintf(stderr,"%s : %s \n",
		"Error in GramsInitials()",
		"Use this format -> GramsInitials (Type=string) { !!!");
	exit(EXIT_FAILURE);
      }
    }
  }

  /* Free array nodes */
  free(Array_Nodes);

  /* Free data */
  free(Name_File_Copy);

  /* Close .dat file */
  fclose(Sim_dat);
}
