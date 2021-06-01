#include "nl-partsol.h"



#ifdef _WIN32
static char * delimiters_1 = " \r\n\t";
static char * delimiters_2 = " =\t\r\n"; 
#else
static char * delimiters_1 = " \n\t";
static char * delimiters_2 = " =\t\n"; 
#endif
static char * delimiters_3 = "(=)";

/**********************************************************************/

Time_Int_Params Solver_selector__InOutFun__(char * Name_File, double DeltaX)
/*
Example : 
GramsTime(Scheme=FE){
	CFL=0.6
	N=4000
}
*/
{

  Time_Int_Params Parameters;


  bool Is_NumTimeStep = false;
  bool Is_FinalTime = false;
  bool Is_Cel = false;

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
      Parameters.CFL = 0.8;
      Parameters.Cel = 0.0;
      Parameters.InitialTimeStep = 0;
      Parameters.NumTimeStep = 0;
      Parameters.FinalTime = 0.0;
      Parameters.epsilon_Mass_Matrix = 0.0;

      Parameters.TOL_Conserving_Energy_Momentum = 1E-10;

      Parameters.rb_Generalized_alpha = 0.6;
      Parameters.TOL_Generalized_alpha = 1E-10;

      Parameters.beta_Newmark_beta = 0.25;
      Parameters.gamma_Newmark_beta = 0.5;
      Parameters.TOL_Newmark_beta = 1E-10;

      Parameters.MaxIter = 10;

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
	  break;
	}
	while(STATUS_LINE != NULL)
	{

	  if(Aux_Temp_id != 2){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsTime()",
		   "Use this format -> Propertie = value !!!");
	    exit(EXIT_FAILURE);
	  }

	  if(strcmp(Parse_Temp_Prop[0],"CFL") == 0)
	  {
	    Parameters.CFL = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","CFL condition",Parameters.CFL);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Cel") == 0)
	  {
	  	Is_Cel = true;
	  	Parameters.Cel = atof(Parse_Temp_Prop[1]);
	  	printf("\t \t -> %s : %f \n","Celerity",Parameters.Cel);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Tend") == 0)
	  {
	  	Is_FinalTime = true;
	  	Parameters.FinalTime = atof(Parse_Temp_Prop[1]);
	  	printf("\t \t -> %s : %e \n","Final time step",Parameters.FinalTime);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"i0") == 0)
	  {
	  	Parameters.InitialTimeStep = atoi(Parse_Temp_Prop[1]);
	  	printf("\t \t -> %s : %i \n","Initial time",Parameters.InitialTimeStep);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"N") == 0)
	  {
	  	Is_NumTimeStep = true;
	    Parameters.NumTimeStep = atoi(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %i \n","Number of time-steps",Parameters.NumTimeStep);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Epsilon") == 0)
	  {
	    Parameters.epsilon_Mass_Matrix = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Epsilon",Parameters.epsilon_Mass_Matrix);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"rb-Generalized-alpha") == 0)
	  {
	    Parameters.rb_Generalized_alpha = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Spectral radio Generalized-alpha",Parameters.rb_Generalized_alpha);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"TOL-Generalized-alpha") == 0)
	  {
	  	Parameters.TOL_Generalized_alpha = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Tolerance Generalized-alpha",Parameters.TOL_Generalized_alpha);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Beta-Newmark-beta") == 0)
	  {
	    Parameters.beta_Newmark_beta = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Beta Newmark-beta",Parameters.beta_Newmark_beta);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Gamma-Newmark-beta") == 0)
	  {
	    Parameters.gamma_Newmark_beta = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Gamma Newmark-beta",Parameters.gamma_Newmark_beta);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"TOL-Newmark-beta") == 0)
	  {
	    Parameters.TOL_Newmark_beta = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Tolerance Newmark-beta",Parameters.TOL_Newmark_beta);
	  } 
	  else if(strcmp(Parse_Temp_Prop[0],"Max-Iter") == 0)
	  {
	  	Parameters.MaxIter = atoi(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %i \n","Max number of interations",Parameters.MaxIter);
	  }
	  else
	  {
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
	if(STATUS_LINE == NULL)
	{
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsTime()",
	       "you forget to put a } !!!");
	exit(EXIT_FAILURE);	  
	}

	if(strcmp(Parse_Temp_Prop[0],"}") == 0)
	{
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

  if(Is_GramsTime == false)
  {
    fprintf(stderr,"%s : %s \n",
	    "Error in GramsSolid()",
	    "GramsTime no defined");
    exit(EXIT_FAILURE);
  }

  if(Counter_GramsTime != 1)
  {
    fprintf(stderr,"%s : %s \n",
	    "Error in GramsSolid2D()",
	    "More than one call to GramsTime");
    exit(EXIT_FAILURE);
  }

  if((Is_NumTimeStep == false) && (Is_FinalTime == true) && (Is_Cel == true))
  {
	Parameters.NumTimeStep = (int)(Parameters.FinalTime*Parameters.Cel*DeltaX)/(Parameters.CFL);
	printf("\t \t -> %s : %i \n","Number of time-steps",Parameters.NumTimeStep);
  }
  else if((Is_NumTimeStep == true) && (Is_FinalTime == false) && (Is_Cel == true))
  {
	Parameters.FinalTime = (int)(Parameters.NumTimeStep*Parameters.CFL*DeltaX/Parameters.Cel);
	printf("\t \t -> %s : %e \n","Final time step",Parameters.FinalTime);
  }
  else
  {
	exit(EXIT_FAILURE);
  }



  /* Close .dat file */
  /* Final message */
  fclose(Sim_dat);

  return Parameters;
}
