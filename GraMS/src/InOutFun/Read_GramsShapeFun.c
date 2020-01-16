#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "../GRAMS/grams.h"

/**********************************************************************/

void GramsShapeFun(char * Name_File)
/*
  Example : 
  GramsShapeFun (Type=MPMQ4) {
  }
  GramsShapeFun (Type=uGIMP) {
  }
  GramsShapeFun (Type=LME) {
  gamma=2.3
  TOL_lambda=10e-6
  }
*/
{
  /* Simulation file */
  FILE * Sim_dat;

  /* Temporal integator */
  int Aux_Shf_id;
  char * Parse_Shf_id[MAXW] = {NULL};

  /* Temporal integrator properties */
  char Line_Shf_Prop[MAXC] = {0};
  char * Parse_Shf_Prop[MAXW] = {NULL};

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  /* Auxiliar variable for status */
  char * STATUS_LINE;

  /* Initial message */  
  puts("*************************************************");
  printf(" \t %s : \n\t %s \n",
	 "* Read shape function properties ",
	 Name_File);
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	   "Error in GramsShapeFun()",
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
	     "Error in GramsShapeFun()",
	     "Parser failed");
      exit(0);
    }

    if (strcmp(kwords[0],"GramsShapeFun") == 0 ){

      /* Read temporal integrator scheme */
      Aux_Shf_id = parse (Parse_Shf_id, kwords[1],"(=)");
      if( (Aux_Shf_id != 2) ||
	  (strcmp(Parse_Shf_id[0],"Type") != 0)){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsShapeFun()",
	       "Use this format -> (Type=string) !!!");
	exit(0);
      }
      ShapeFunctionGP = Parse_Shf_id[1];
      printf("\t -> %s : %s \n","Type of shape function",ShapeFunctionGP);

      /* Set to default all it properties */
      gamma_LME=0;
      TOL_lambda=10e-6;

      /* Look for the curly brace { */
      if(strcmp(kwords[2],"{") == 0){
	/* Initial line */
	STATUS_LINE = fgets(Line_Shf_Prop,
			    sizeof(Line_Shf_Prop),
			    Sim_dat);
	if(STATUS_LINE == NULL){
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsShapeFun()",
		  "Unspected EOF !!!");
	  exit(0);	
	}
	Aux_Shf_id = parse(Parse_Shf_Prop,Line_Shf_Prop," =\t\n");
	if(strcmp(Parse_Shf_Prop[0],"}") == 0){
	  /* Check gamma */
	  if((strcmp(ShapeFunctionGP,"LME") == 0) &&
	     (gamma_LME == 0)){
	    fprintf(stderr,"%s : %s \n",
		    "Error in GramsShapeFun()",
		    "gamma parameter required for LME !!!");
	    exit(0);
	  }
	  /* Check number of time steps */
	  break;
	}
	while(STATUS_LINE != NULL){
	  
	  if(Aux_Shf_id != 2){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsShapeFun()",
		   "Use this format -> Propertie = value !!!");
	    exit(0);
	  }

 	  if(strcmp(Parse_Shf_Prop[0],"gamma") == 0){
	    gamma_LME = atof(Parse_Shf_Prop[1]);
	  }
	  else if(strcmp(Parse_Shf_Prop[0],"TOL_lambda") == 0){
	    TOL_lambda = atof(Parse_Shf_Prop[1]);
	  }
	  else{
	    fprintf(stderr,"%s : %s %s \n",
		   "Error in GramsShapeFun()",
		   "Undefined",Parse_Shf_Prop[0]);
	    exit(0);
	  }
	  /* Read next line and check */
	  STATUS_LINE = fgets(Line_Shf_Prop,
			      sizeof(Line_Shf_Prop),
			      Sim_dat);
	  Aux_Shf_id = parse(Parse_Shf_Prop,Line_Shf_Prop," =\t\n");
	  if(strcmp(Parse_Shf_Prop[0],"}") == 0){
	    break;
	  }
	}
	if(STATUS_LINE == NULL){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsShapeFun()",
	       "you forget to put a } !!!");
	exit(0);	  
	}

	if(strcmp(Parse_Shf_Prop[0],"}") == 0){
	  /* Check gamma */
	  if((strcmp(ShapeFunctionGP,"LME") == 0) &&
	     (gamma_LME == 0)){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsShapeFun()",
		   "gamma parameter required for LME !!!");
	    exit(0);
	  }
	  /* Check number of time steps */
	  break;
	}
      }
      else{
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsShapeFun()",
	       "Use this format -> GramsShapeFun (Type=string) { !!!");
	exit(0);
      }
    }
  }

  /* Close .dat file */
  /* Final message */
  fclose(Sim_dat);

}

