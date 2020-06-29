#include "nl-partsol.h"

/*
  Call global functions
*/
double CFL; /* Courant number (0-1) */
double DeltaTimeStep;
double SpectralRadius;
int NumTimeStep;

/**********************************************************************/

void GramsTime(char * Name_File)
/*
Example : 
GramsTime(Scheme=FE){
	CFL=0.6
	N=4000
}
*/
{
  /* Check the number of GramsTime calls */
  bool Is_GramsTime = false;
  int Counter_GramsTime = 0;

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
    exit(EXIT_FAILURE);
  }
  
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");
    if (nkwords < 0){
      fprintf(stderr,"%s : %s \n",
	     "Error in GramsTime()",
	     "Parser failed");
      exit(EXIT_FAILURE);
    }

    if ((nkwords > 0) &&
	(strcmp(kwords[0],"GramsTime") == 0)){

      /* Set to true the boolean flag */
      Is_GramsTime = true;
      Counter_GramsTime++;

      /* Read temporal integrator scheme */
      Aux_Temp_id = parse (Parse_Temp_id, kwords[1],"(=)");
      if( (Aux_Temp_id != 2) ||
	  (strcmp(Parse_Temp_id[0],"Type") != 0)){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsTime()",
	       "Use this format -> (Type=string) !!!");
	exit(EXIT_FAILURE);
      }
      TimeIntegrationScheme = Parse_Temp_id[1];
      printf("\t -> %s : %s \n","Time integrator",TimeIntegrationScheme);

      /* Set to default all it properties */
      CFL=0.8;
      NumTimeStep=0;
      SpectralRadius=0.6;

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
	    exit(EXIT_FAILURE);
	}
	Aux_Temp_id = parse(Parse_Temp_Prop,Line_Temp_Prop," =\t\n");
	if(strcmp(Parse_Temp_Prop[0],"}") == 0){

	  /* Check number of time steps */	  
	  if(NumTimeStep == 0){
	    fprintf(stderr,"%s : %s \n",
	   "Error in GramsTime()",
		   "Number of time steps <= 0 !!!");
	    exit(EXIT_FAILURE);
	  }
	  break;
	}
	while(STATUS_LINE != NULL){

	  if(Aux_Temp_id != 2){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsTime()",
		   "Use this format -> Propertie = value !!!");
	    exit(EXIT_FAILURE);
	  }

	  if(strcmp(Parse_Temp_Prop[0],"CFL") == 0){
	    CFL = atof(Parse_Temp_Prop[1]);
	    printf("\t -> %s : %f \n","CFL condition",CFL);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"N") == 0){
	    NumTimeStep = atoi(Parse_Temp_Prop[1]);
	    printf("\t -> %s : %i \n","Number of time-steps",NumTimeStep);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"rb") == 0){
	    SpectralRadius = atof(Parse_Temp_Prop[1]);
	    printf("\t -> %s : %f \n","Spectral radio",SpectralRadius);
	  } 
	  else{
	    fprintf(stderr,"%s : %s %s \n",
		   "Error in GramsTime()",
		   "Undefined",Parse_Temp_Prop[0]);
	    exit(EXIT_FAILURE);
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
	exit(EXIT_FAILURE);	  
	}

	if(strcmp(Parse_Temp_Prop[0],"}") == 0){

	  /* Check number of time steps */	  
	  if(NumTimeStep == 0){
	    fprintf(stderr,"%s : %s \n",
	   "Error in GramsTime()",
		   "Number of time steps <= 0 !!!");
	    exit(EXIT_FAILURE);
	  }
	  break;
	}
      }
      else{
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsTime()",
	       "Use this format -> GramsTime (Type=string) { !!!");
	exit(EXIT_FAILURE);
      }
    }
  }

  if(Is_GramsTime == false){
    fprintf(stderr,"%s : %s \n",
	    "Error in GramsSolid2D()",
	    "GramsTime no defined");
    exit(EXIT_FAILURE);
  }
  if(Counter_GramsTime != 1){
    fprintf(stderr,"%s : %s \n",
	    "Error in GramsSolid2D()",
	    "More than one call to GramsTime");
    exit(EXIT_FAILURE);
  }

  /* Close .dat file */
  /* Final message */
  fclose(Sim_dat);

}
