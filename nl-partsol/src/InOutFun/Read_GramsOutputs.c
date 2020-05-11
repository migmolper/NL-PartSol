#include "nl-partsol.h"

/**********************************************************************/

void GramsOutputs(char * Name_File)
/*
  Example : 
  GramsOutputs (i=100) {
  DIR=test/Sulsky_MPM	
  }
*/
{
  /* Simulation file */
  FILE * Sim_dat;

  /* Temporal integator */
  int Aux_Out_id;
  char * Parse_Out_id[MAXW] = {NULL};

  /* Temporal integrator properties */
  char Line_Out_Prop[MAXC] = {0};
  char * Parse_Out_Prop[MAXW] = {NULL};

  /* Parse file for the route */
  char Route_Outs[MAXC] = {0};

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  /* Auxiliar variable for status */
  char * STATUS_LINE;

  /* Initial message */  
  puts("*************************************************");
  printf(" \t %s : \n\t %s \n",
	 "* Read Outputs properties ",
	 Name_File);
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	   "Error in GramsOutputs()",
	   "Incorrect lecture of",
	   Name_File);
    exit(0);
  }

  /* Generate route */
  generate_route(Route_Outs,Name_File);
    
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");
    if (nkwords < 0){
      fprintf(stderr,"%s : %s \n",
	     "Error in GramsOutputs()",
	     "Parser failed");
      exit(0);
    }

    if ((nkwords > 0) &&
	(strcmp(kwords[0],"GramsOutputs") == 0 )){

      /* Read temporal integrator scheme */
      Aux_Out_id = parse (Parse_Out_id, kwords[1],"(=)");
      if( (Aux_Out_id != 2) ||
	  (strcmp(Parse_Out_id[0],"i") != 0)){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsOutputs()",
	       "Use this format -> (i=int) !!!");
	exit(0);
      }
      ResultsTimeStep = atoi(Parse_Out_id[1]);
      if(ResultsTimeStep > NumTimeStep){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsOutputs()",
	       "The result time step should be less than the total time steps !!!");
	exit(0);	
      }
      printf("\t -> %s : %i \n","Output values each",ResultsTimeStep);

      /* Look for the curly brace { */
      if(strcmp(kwords[2],"{") == 0){
	/* Initial line */
	STATUS_LINE = fgets(Line_Out_Prop,
			    sizeof(Line_Out_Prop),
			    Sim_dat);
	if(STATUS_LINE == NULL){
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsOutputs()",
		  "Unspected EOF !!!");
	  exit(0);	
	}
	Aux_Out_id = parse(Parse_Out_Prop,Line_Out_Prop," =\t\n");
	if(strcmp(Parse_Out_Prop[0],"}") == 0){
	  /* Check output dir */
	  if(OutputDir == NULL){
	    fprintf(stderr,"%s : %s \n",
		    "Error in GramsOutputs()",
		    "Non output dir defined !!!");
	    exit(0);
	  }
	  break;
	}
	while(STATUS_LINE != NULL){
	  
	  if(Aux_Out_id != 2){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsOutputs()",
		   "Use this format -> Propertie = value !!!");
	    exit(0);
	  }

 	  if(strcmp(Parse_Out_Prop[0],"DIR") == 0){
	    sprintf(OutputDir,"%s%s",Route_Outs,Parse_Out_Prop[1]);
	    printf("\t -> %s : %s \n","Output directory",OutputDir);
	  }
	  else{
	    fprintf(stderr,"%s : %s %s \n",
		   "Error in GramsOutputs()",
		   "Undefined",Parse_Out_Prop[0]);
	    exit(0);
	  }
	  /* Read next line and check */
	  STATUS_LINE = fgets(Line_Out_Prop,
			      sizeof(Line_Out_Prop),
			      Sim_dat);
	  Aux_Out_id = parse(Parse_Out_Prop,Line_Out_Prop," =\t\n");
	  if(strcmp(Parse_Out_Prop[0],"}") == 0){
	    break;
	  }
	}
	if(STATUS_LINE == NULL){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsOutputs()",
	       "you forget to put a } !!!");
	exit(0);	  
	}

	if(strcmp(Parse_Out_Prop[0],"}") == 0){
	  /* Check output dir */
	  if(OutputDir == NULL){
	    fprintf(stderr,"%s : %s \n",
		    "Error in GramsOutputs()",
		    "Non output dir defined !!!");
	    exit(0);
	  }
	  break;
	}
      }
      else{
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsOutputs()",
	       "Use this format -> GramsOutputs (Type=string) { !!!");
	exit(0);
      }
    }
  }

  /* Close .dat file */
  /* Final message */
  fclose(Sim_dat);

}


