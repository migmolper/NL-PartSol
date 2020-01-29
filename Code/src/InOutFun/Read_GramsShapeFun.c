#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "grams.h"

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
  int Num_GramsShapefun_Type;
  char * Parse_GramsShapefun_Type[MAXW] = {NULL};

  /* Temporal integrator properties */
  int Num_ShapeFun_Prop;
  char Line_Shf_Prop[MAXC] = {0};
  char * Parse_Shf_Prop[MAXW] = {NULL};

  /* Special variables for GramsShapeFun */
  int Num_GramsShapeFun;
  char Line_GramsShapeFun[MAXC] = {0}; 
  char * Parse_GramsShapeFun[MAXW] = {NULL};

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
  while( fgets(Line_GramsShapeFun, sizeof Line_GramsShapeFun, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    Num_GramsShapeFun = parse (Parse_GramsShapeFun, Line_GramsShapeFun," \n\t");
    if (Num_GramsShapeFun < 0){
      fprintf(stderr,"%s : %s \n",
	     "Error in GramsShapeFun()",
	     "Parser failed");
      exit(0);
    }

    if ((Num_GramsShapeFun > 0) &&
	(strcmp(Parse_GramsShapeFun[0],"GramsShapeFun") == 0)){

      /* Read Type of shape function */
      Num_GramsShapefun_Type =
	parse(Parse_GramsShapefun_Type, Parse_GramsShapeFun[1],"(=)");
      if( (Num_GramsShapefun_Type != 2) ||
	  (strcmp(Parse_GramsShapefun_Type[0],"Type") != 0)){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsShapeFun()",
	       "Use this format -> (Type=string) !!!");
	exit(0);
      }
      
      ShapeFunctionGP = Parse_GramsShapefun_Type[1];
      printf("\t -> %s : %s \n","Type of shape function",ShapeFunctionGP);

      /* Set to default all it properties */
      gamma_LME=0;
      TOL_lambda=10e-6;

      /* Look for the curly brace { */
      if((Num_GramsShapeFun>=3) &&
	 (strcmp(Parse_GramsShapeFun[2],"{") == 0)){

	/* In case GramsShapeFun (Type=*) {} */
	if((Num_GramsShapeFun == 4) &&
	   (strcmp(Parse_GramsShapeFun[3],"}") == 0)){
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
	
	/* Initial line */
	STATUS_LINE = fgets(Line_Shf_Prop,sizeof(Line_Shf_Prop),Sim_dat);
	if(STATUS_LINE == NULL){
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsShapeFun()",
		  "Unspected EOF !!!");
	  exit(0);	
	}
	Num_ShapeFun_Prop = parse(Parse_Shf_Prop,Line_Shf_Prop," =\t\n");

	/* In case GramsShapeFun (Type=*) {
	   } */
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

	/* Loop */
	while(STATUS_LINE != NULL){
	  
	  if(Num_ShapeFun_Prop != 2){
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
	  Num_ShapeFun_Prop = parse(Parse_Shf_Prop,Line_Shf_Prop," =\t\n");
	  if((Num_ShapeFun_Prop>0) &&
	     (strcmp(Parse_Shf_Prop[0],"}") == 0)){
	    break;
	  }
	}
	if(STATUS_LINE == NULL){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsShapeFun()",
	       "you forget to put a } !!!");
	exit(0);	  
	}

	if((Num_ShapeFun_Prop>0) &&
	   (strcmp(Parse_Shf_Prop[0],"}") == 0)){
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

