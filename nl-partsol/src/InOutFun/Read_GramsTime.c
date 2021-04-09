#include "nl-partsol.h"


/*
  Call global variables
*/
int NumTimeStep;
double SpectralRadius;
double CFL;
double epsilon_Mass_Matrix; 
double beta_Newmark_beta;   
double gamma_Newmark_beta;
double TOL_Newmark_beta;

#ifdef _WIN32
static char * delimiters_1 = " \r\n\t";
static char * delimiters_2 = " =\t\r\n"; 
#else
static char * delimiters_1 = " \n\t";
static char * delimiters_2 = " =\t\n"; 
#endif
static char * delimiters_3 = "(=)";

/**********************************************************************/

void Solver_selector__InOutFun__(char * Name_File)
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
  printf(" \t %s : %s \n","* Read time integration properties in",Name_File);
  
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
    nkwords = parse (kwords, line,delimiters_1);
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
      Aux_Temp_id = parse (Parse_Temp_id, kwords[1],delimiters_3);
      if( (Aux_Temp_id != 2) ||
	  (strcmp(Parse_Temp_id[0],"Type") != 0)){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsTime()",
	       "Use this format -> (Type=string) !!!");
	exit(EXIT_FAILURE);
      }
      TimeIntegrationScheme = Parse_Temp_id[1];
      printf("\t \t -> %s : %s \n","Time integrator",TimeIntegrationScheme);
     
      /* Set to default all it properties */
      CFL=0.8;
      NumTimeStep=0;
      SpectralRadius=0.6;
      epsilon_Mass_Matrix = 0.0;
      beta_Newmark_beta = 0.25;
      gamma_Newmark_beta = 0.5;
      TOL_Newmark_beta = 0.000000001;

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
	Aux_Temp_id = parse(Parse_Temp_Prop,Line_Temp_Prop,delimiters_2);
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
	    printf("\t \t -> %s : %f \n","CFL condition",CFL);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"N") == 0){
	    NumTimeStep = atoi(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %i \n","Number of time-steps",NumTimeStep);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"rb") == 0){
	    SpectralRadius = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Spectral radio",SpectralRadius);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Epsilon") == 0){
	    epsilon_Mass_Matrix = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Epsilon",epsilon_Mass_Matrix);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Beta-Newmark") == 0){
	    beta_Newmark_beta = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Beta-Newmark",beta_Newmark_beta);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Gamma-Newmark") == 0){
	    gamma_Newmark_beta = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Gamma-Newmark",gamma_Newmark_beta);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"TOL-Newmark") == 0){
	    TOL_Newmark_beta = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Tolerance Newmark",TOL_Newmark_beta);
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
	  Aux_Temp_id = parse(Parse_Temp_Prop,Line_Temp_Prop,delimiters_2);
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
