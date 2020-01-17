#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "../GRAMS/grams.h"

/**********************************************************************/

void GramsTime(char * Name_File)
/*
Example : 
GramsTime(Scheme=FE){
	CEL=70.71
	CFL=0.6
	N=4000
}
*/
{

  /* Simulation file */
  FILE * Sim_dat;

  /* Temporal integator */
  int Aux_Temp_id;
  char * Parse_Temp_id[MAXW] = {NULL};

  /* Temporal integrator properties */
  char Line_Temp_Prop[MAXC] = {0};
  char * Parse_Temp_Prop[MAXW] = {NULL};

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  /* Auxiliar variable for status */
  char * STATUS_LINE;

  /* Initial message */  
  puts("*************************************************");
  printf(" \t %s : \n\t %s \n",
	 "* Read time integration properties in ",
	 Name_File);
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	   "Error in GramsTime()",
	   "Incorrect lecture of",
	   Name_File);
    exit(0);
  }
  
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");
    if (nkwords < 0){
      fprintf(stderr,"%s : %s \n",
	     "Error in GramsTime()",
	     "Parser failed");
      exit(0);
    }

    if (strcmp(kwords[0],"GramsTime") == 0 ){

      /* Read temporal integrator scheme */
      Aux_Temp_id = parse (Parse_Temp_id, kwords[1],"(=)");
      if( (Aux_Temp_id != 2) ||
	  (strcmp(Parse_Temp_id[0],"Type") != 0)){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsTime()",
	       "Use this format -> (Type=string) !!!");
	exit(0);
      }
      TimeIntegration = Parse_Temp_id[1];
      printf("\t -> %s : %s \n","Time integrator",TimeIntegration);

      /* Set to default all it properties */
      CFL=0.8;
      CEL=NAN;
      NumTimeStep=0;

      /* Look for the curly brace { */
      if(strcmp(kwords[2],"{") == 0){
	/* Initial line */
	STATUS_LINE = fgets(Line_Temp_Prop,
			    sizeof(Line_Temp_Prop),
			    Sim_dat);
	if(STATUS_LINE == NULL){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsTime()",
		   "Unspected EOF !!!");
	    exit(0);	
	}
	Aux_Temp_id = parse(Parse_Temp_Prop,Line_Temp_Prop," =\t\n");
	if(strcmp(Parse_Temp_Prop[0],"}") == 0){

	  /* Check celerity */
	  if(CEL != CEL){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsTime()",
		   "Celerity parameter required !!!");
	    exit(0);
	  }
	  /* Check number of time steps */
	  
	  if(NumTimeStep == 0){
	    fprintf(stderr,"%s : %s \n",
	   "Error in GramsTime()",
		   "Number of time steps <= 0 !!!");
	    exit(0);
	  }
	  break;
	}
	while(STATUS_LINE != NULL){

	  if(Aux_Temp_id != 2){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsTime()",
		   "Use this format -> Propertie = value !!!");
	    exit(0);
	  }

 	  if(strcmp(Parse_Temp_Prop[0],"CEL") == 0){
	    CEL = atof(Parse_Temp_Prop[1]);
	    printf("\t -> %s : %f \n","Celerity",CEL);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"CFL") == 0){
	    CFL = atof(Parse_Temp_Prop[1]);
	    printf("\t -> %s : %f \n","CFL condition",CFL);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"N") == 0){
	    NumTimeStep = atoi(Parse_Temp_Prop[1]);
	    printf("\t -> %s : %i \n","Number of time-steps",NumTimeStep);
	  }  
	  else{
	    fprintf(stderr,"%s : %s %s \n",
		   "Error in GramsTime()",
		   "Undefined",Parse_Temp_Prop[0]);
	    exit(0);
	  }
	  /* Read next line and check */
	  STATUS_LINE = fgets(Line_Temp_Prop,
			      sizeof(Line_Temp_Prop),
			      Sim_dat);
	  Aux_Temp_id = parse(Parse_Temp_Prop,Line_Temp_Prop," =\t\n");
	  if(strcmp(Parse_Temp_Prop[0],"}") == 0){
	    break;
	  }
	}
	if(STATUS_LINE == NULL){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsTime()",
	       "you forget to put a } !!!");
	exit(0);	  
	}

	if(strcmp(Parse_Temp_Prop[0],"}") == 0){

	  /* Check celerity */
	  if(CEL != CEL){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsTime()",
		   "Celerity parameter required !!!");
	    exit(0);
	  }
	  /* Check number of time steps */
	  
	  if(NumTimeStep == 0){
	    fprintf(stderr,"%s : %s \n",
	   "Error in GramsTime()",
		   "Number of time steps <= 0 !!!");
	    exit(0);
	  }
	  break;
	}
      }
      else{
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsTime()",
	       "Use this format -> GramsTime (Type=string) { !!!");
	exit(0);
      }
    }
  }

  /* Close .dat file */
  /* Final message */
  fclose(Sim_dat);

}
