#include "nl-partsol.h"


/*
  Auxiliar functions and variables
*/
#ifdef _WIN32
static char * delimiters_1 = " \r\n\t";
static char * delimiters_2 = " =\t\r\n"; 
#else
static char * delimiters_1 = " \n\t";
static char * delimiters_2 = " =\t\n"; 
#endif
static char * delimiters_3 = "(=)";


typedef struct
{
  bool Is_NumTimeStep;
  bool Is_CFL;
  bool Is_Cel;
  bool Is_Solver;
  bool Is_Epsilon_Mass_Matrix;  
  int Counter_Solver;

  /*!
  * Check Newmark parameters
  */
  bool Is_beta_Newmark_beta;   
  bool Is_gamma_Newmark_beta;
  bool Is_TOL_Newmark_beta;

  /*!
   * Check alpha parameters 
   */
  bool Is_rb_Generalized_alpha;
  bool Is_TOL_Generalized_alpha;
  
  /*!
  * Check Conserving Energy-Momentum parameters 
  */
  bool Is_TOL_Conserving_Energy_Momentum;

} Check_Params;

static Time_Int_Params Initialise_Parameters();
static Check_Params Initialise_Check_Params();
static void check_Solver(Time_Int_Params,Check_Params);

/**********************************************************************/

Time_Int_Params Solver_selector__InOutFun__(char * Name_File)
{

  /* Set to default some parameters */
  Time_Int_Params Parameters = Initialise_Parameters();
  Check_Params ChkParam = Initialise_Check_Params();

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
  if(Sim_dat==NULL)
  {
    fprintf(stderr,"%s : \n\t %s %s","Error in NLPS-Solver()","Incorrect lecture of",Name_File);
    exit(EXIT_FAILURE);
  }
  
  /* Read the file line by line */
  while(fgets(line, sizeof line, Sim_dat) != NULL)
  {

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line,delimiters_1);
    if(nkwords < 0)
	{
      fprintf(stderr,"%s : %s \n","Error in NLPS-Solver()","Parser failed");
      exit(EXIT_FAILURE);
    }

    if((nkwords > 0) 
	&& (strcmp(kwords[0],"NLPS-Solver") == 0))
	{

      /* Set to true the boolean flag */
      ChkParam.Is_Solver = true;
      ChkParam.Counter_Solver++;

      /* Read temporal integrator scheme */
      Aux_Temp_id = parse (Parse_Temp_id, kwords[1],delimiters_3);
      if((Aux_Temp_id != 2) 
	  || (strcmp(Parse_Temp_id[0],"Type") != 0))
	  {
		fprintf(stderr,"%s : %s \n","Error in NLPS-Solver()","Use this format -> (Type=string) !!!");
		exit(EXIT_FAILURE);
      }

	  strcpy(Parameters.TimeIntegrationScheme,Parse_Temp_id[1]);
    
      /* Look for the curly brace { */
      if(strcmp(kwords[2],"{") == 0)
	  {
		  /* Initial line */
		  STATUS_LINE = fgets(Line_Temp_Prop,sizeof(Line_Temp_Prop),Sim_dat);
		  if(STATUS_LINE == NULL)
		  {
			fprintf(stderr,"%s : %s \n","Error in NLPS-Solver()","Unspected EOF !!!");
			exit(EXIT_FAILURE);
		  }
		  
		  Aux_Temp_id = parse(Parse_Temp_Prop,Line_Temp_Prop,delimiters_2);
		  if(strcmp(Parse_Temp_Prop[0],"}") == 0)
		  {
			  break;
		  }
		  
		  while(STATUS_LINE != NULL)
		  {
			if(Aux_Temp_id != 2)
			{
				fprintf(stderr,"%s : %s \n","Error in NLPS-Solver()","Use this format -> Propertie = value !!!");
				exit(EXIT_FAILURE);
			}

	  if(strcmp(Parse_Temp_Prop[0],"CFL") == 0)
	  {
		  ChkParam.Is_CFL = true;
	    Parameters.CFL = atof(Parse_Temp_Prop[1]);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Cel") == 0)
	  {
	  	ChkParam.Is_Cel = true;
	  	Parameters.Cel = atof(Parse_Temp_Prop[1]);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"i0") == 0)
	  {
	  	Parameters.InitialTimeStep = atoi(Parse_Temp_Prop[1]);
	  	printf("\t \t -> %s : %i \n","Initial time",Parameters.InitialTimeStep);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"N") == 0)
	  {
	  	ChkParam.Is_NumTimeStep = true;
	    Parameters.NumTimeStep = atoi(Parse_Temp_Prop[1]);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Epsilon") == 0)
	  {
	    ChkParam.Is_Epsilon_Mass_Matrix = true;
	    Parameters.epsilon_Mass_Matrix = atof(Parse_Temp_Prop[1]);
	    printf("\t \t -> %s : %f \n","Epsilon",Parameters.epsilon_Mass_Matrix);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"rb-Generalized-alpha") == 0)
	  {
		ChkParam.Is_rb_Generalized_alpha = true;
	    Parameters.rb_Generalized_alpha = atof(Parse_Temp_Prop[1]);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"TOL-Generalized-alpha") == 0)
	  {
		ChkParam.Is_TOL_Generalized_alpha = true;
	  	Parameters.TOL_Generalized_alpha = atof(Parse_Temp_Prop[1]);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Beta-Newmark-beta") == 0)
	  {
		ChkParam.Is_beta_Newmark_beta = true;
	    Parameters.beta_Newmark_beta = atof(Parse_Temp_Prop[1]);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"Gamma-Newmark-beta") == 0)
	  {
		ChkParam.Is_gamma_Newmark_beta = true;
	    Parameters.gamma_Newmark_beta = atof(Parse_Temp_Prop[1]);
	  }
	  else if(strcmp(Parse_Temp_Prop[0],"TOL-Newmark-beta") == 0)
	  {
		ChkParam.Is_TOL_Newmark_beta = true;
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
	    fprintf(stderr,"%s : %s %s \n","Error in NLPS-Solver()","Undefined",Parse_Temp_Prop[0]);
	    exit(EXIT_FAILURE);
	  }
	  /* Read next line and check */
	  STATUS_LINE = fgets(Line_Temp_Prop,sizeof(Line_Temp_Prop),Sim_dat);
	  Aux_Temp_id = parse(Parse_Temp_Prop,Line_Temp_Prop,delimiters_2);
	  if(strcmp(Parse_Temp_Prop[0],"}") == 0)
	  {
	    break;
	  }
	}
	if(STATUS_LINE == NULL)
	{
		fprintf(stderr,"%s : %s \n","Error in NLPS-Solver()","you forget to put a } !!!");
		exit(EXIT_FAILURE);	  
	}

	if(strcmp(Parse_Temp_Prop[0],"}") == 0)
	{
	  break;
	}
      }
      else
	  {
		fprintf(stderr,"%s : %s \n","Error in NLPS-Solver()","Use this format -> NLPS-Solver (Type=string) { !!!");
		exit(EXIT_FAILURE);
      }
    }
  }

  /* Close .dat file */
  fclose(Sim_dat);

  check_Solver(Parameters,ChkParam);

  return Parameters;
}

/**********************************************************************/

static Check_Params Initialise_Check_Params()
{
  Check_Params ChkParam;

  ChkParam.Is_NumTimeStep = false;
  ChkParam.Is_Cel = false;
  ChkParam.Is_CFL = false;
  ChkParam.Is_Solver = false;
  ChkParam.Is_Epsilon_Mass_Matrix = false;
  ChkParam.Counter_Solver = 0;

  ChkParam.Is_beta_Newmark_beta = false;   
  ChkParam.Is_gamma_Newmark_beta = false;
  ChkParam.Is_TOL_Newmark_beta = false;

  ChkParam.Is_rb_Generalized_alpha = false;
  ChkParam.Is_TOL_Generalized_alpha = false;

  ChkParam.Is_TOL_Conserving_Energy_Momentum = false;

  return ChkParam;
}

/**********************************************************************/

static Time_Int_Params Initialise_Parameters()
{
  Time_Int_Params Parameters;
  
  Parameters.CFL = 0.8;
  Parameters.Cel = 0.0;
  Parameters.InitialTimeStep = 0;
  Parameters.NumTimeStep = 0;
  Parameters.FinalTime = 0.0;
  Parameters.epsilon_Mass_Matrix = 1.0;
  Parameters.TOL_Conserving_Energy_Momentum = 1E-10;
  Parameters.rb_Generalized_alpha = 0.6;
  Parameters.TOL_Generalized_alpha = 1E-10;
  Parameters.beta_Newmark_beta = 0.25;
  Parameters.gamma_Newmark_beta = 0.5;
  Parameters.TOL_Newmark_beta = 1E-10;
  Parameters.MaxIter = 10;

  return Parameters;
}



/**********************************************************************/

static void check_Solver(
	Time_Int_Params Parameters, 
	Check_Params ChkParam)
{
	
	if(ChkParam.Is_Solver == false)
	{
		fprintf(stderr,"%s : %s \n","Error in NLPS-Solver()","NLPS-Solver no defined");
		exit(EXIT_FAILURE);
	}
	
	if(ChkParam.Counter_Solver != 1)
	{
		fprintf(stderr,"%s : %s \n","Error in NLPS-Solver()","More than one call to NLPS-Solver");
		exit(EXIT_FAILURE);
	}


  if(ChkParam.Is_NumTimeStep 
  && ChkParam.Is_Cel
  && ChkParam.Is_CFL)
  {
    printf("\t \t -> %s : %i \n","Number of TimeStep",Parameters.NumTimeStep);
    printf("\t \t -> %s : %f \n","Celerity",Parameters.Cel);
    printf("\t \t -> %s : %f \n","Courant condition",Parameters.CFL);

    if(strcmp(Parameters.TimeIntegrationScheme,"Newmark-beta-Finite-Strains") == 0)
    {
		if(ChkParam.Is_beta_Newmark_beta  
		&& ChkParam.Is_gamma_Newmark_beta
		&& ChkParam.Is_TOL_Newmark_beta
		&& ChkParam.Is_Epsilon_Mass_Matrix)
		{
			printf("\t \t -> %s : %s \n","Time integrator",Parameters.TimeIntegrationScheme);
			printf("\t \t \t -> %s : %f \n","beta",Parameters.beta_Newmark_beta);
			printf("\t \t \t -> %s : %f \n","gamma",Parameters.gamma_Newmark_beta);
			printf("\t \t \t -> %s : %e \n","epsilon",Parameters.epsilon_Mass_Matrix);
			printf("\t \t \t -> %s : %e \n","TOL",Parameters.TOL_Newmark_beta);
		}
		else
		{
			fprintf(stderr,"%s : %s \n",
			"Error in NLPS-Solver()","Some parameter for the Newmark-beta-Finite-Strains");
			fputs(ChkParam.Is_beta_Newmark_beta ? "beta : true \n" : "beta : false \n", stdout);
			fputs(ChkParam.Is_gamma_Newmark_beta ? "gamma : true \n" : "gamma : false \n", stdout);
			fputs(ChkParam.Is_Epsilon_Mass_Matrix ? "epsilon : true \n" : "epsilon : false \n", stdout);
			fputs(ChkParam.Is_TOL_Newmark_beta ? "TOL : true \n" : "TOL : false \n", stdout);
			exit(EXIT_FAILURE);
		}
	}
	else if(strcmp(Parameters.TimeIntegrationScheme,"Generalized-alpha") == 0)
	{
		if(ChkParam.Is_rb_Generalized_alpha 
		&& ChkParam.Is_TOL_Generalized_alpha)
		{
			printf("\t \t -> %s : %s \n","Time integrator",Parameters.TimeIntegrationScheme);
			printf("\t \t \t -> %s : %f \n","radius",Parameters.rb_Generalized_alpha);
			printf("\t \t \t -> %s : %e \n","TOL",Parameters.TOL_Generalized_alpha);
		}
		else
		{
			fprintf(stderr,"%s : %s \n",
			"Error in NLPS-Solver()","Some parameter for the Generalized-alpha");
			fputs(ChkParam.Is_rb_Generalized_alpha ? "radius : true \n" : "radius : false \n", stdout);
			fputs(ChkParam.Is_TOL_Generalized_alpha ? "TOL : true \n" : "TOL : false \n", stdout);
			exit(EXIT_FAILURE);
		}

	}
	else if(strcmp(Parameters.TimeIntegrationScheme,"Discrete-Energy-Momentum") == 0)
	{
		if(ChkParam.Is_rb_Generalized_alpha 
		&& ChkParam.Is_TOL_Generalized_alpha
		&& ChkParam.Is_Epsilon_Mass_Matrix)
		{
			printf("\t \t -> %s : %s \n","Time integrator",Parameters.TimeIntegrationScheme);
			printf("\t \t \t -> %s : %f \n","radius",Parameters.rb_Generalized_alpha);
			printf("\t \t \t -> %s : %e \n","epsilon",Parameters.epsilon_Mass_Matrix);
			printf("\t \t \t -> %s : %e \n","TOL",Parameters.TOL_Conserving_Energy_Momentum);
		}
		else
		{
			fprintf(stderr,"%s : %s \n",
			"Error in NLPS-Solver()","Some parameter for the Discrete-Energy-Momentum");
			fputs(ChkParam.Is_rb_Generalized_alpha ? "radius : true \n" : "radius : false \n", stdout);
			fputs(ChkParam.Is_Epsilon_Mass_Matrix ? "epsilon : true \n" : "epsilon : false \n", stdout);
			fputs(ChkParam.Is_TOL_Conserving_Energy_Momentum ? "TOL : true \n" : "TOL : false \n", stdout);
			exit(EXIT_FAILURE);
		}

	}
	else
	{
		printf("\t \t -> %s : %s \n","Time integrator",Parameters.TimeIntegrationScheme);
	}
  }
  else
  {
    fprintf(stderr,"%s : %s \n",
      "Error in NLPS-Solver()","Some parameter is missed");
    fputs(ChkParam.Is_NumTimeStep ? "Number of TimeStep : true \n" : "Number of TimeStep : false \n", stdout);
    fputs(ChkParam.Is_Cel ? "Celerity : true \n" : "Celerity : false \n", stdout);
    fputs(ChkParam.Is_CFL ? "Courant condition : true \n" : "Courant condition : false \n", stdout);
    exit(EXIT_FAILURE);
  }
}

/**********************************************************************/