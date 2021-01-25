#include "nl-partsol.h"


/*
  Call global variables
*/
double TOL_Radial_Returning;
int Max_Iterations_Radial_Returning;


/*
  Local structures
*/
typedef struct
{
  int   Idx;
  char  Model [100];

} Param_Index_and_Model;

typedef struct
{

  bool Is_Soil_Idx;
  bool Is_Water_Idx;
  bool Is_Permeability;
  bool Is_phi_s_0;
  bool Is_phi_f_0;

} Check_Mixture;

/*
  Auxiliar functions and variables
*/

#ifdef _WIN32
static char * delimiters_1 = " ,()\r\n\t";
static char * delimiters_2 = " =\t\r\n"; 
#else
static char * delimiters_1 = " ,()\n\t";
static char * delimiters_2 = " =\t\n"; 
#endif
static char * delimiters_3 = "=";

static char Error_message[MAXW];

static Param_Index_and_Model Read_Index_and_Model(char *,char *);
static Mixture  Define_Mixture(FILE *, Param_Index_and_Model);
static bool Activate_Options(char *, char *);
static Check_Mixture Initialise_Check_Mixture();
static void check_Borja_2004_Soil_water_Mixture(Mixture, Check_Mixture);
static void standard_error();
static void standard_output(char *);
static FILE * Open_and_Check_simulation_file(char *);

/**********************************************************************/

Mixture * Read_Soil_Water_Mixtures__InOutFun__(char * SimulationFile, int Number_Soil_Water_Mixtures)
/*

Define-Mixture (idx=0, Model=Borja-2004)
{
  Soil-Idx=0
  Fluid-Idx=1
  Permeability-Skeleton=0.1
  Reference-Volume-Fraction-Soil=0.58
  Reference-Volume-Fraction-Fluid=0.42
}

*/
{
	/* Simulation file */
	FILE * Sim_dat;

  /* Variables for reading purposes */
  char line[MAXC] = {0};
  char * kwords[MAXW] = {NULL};
  int nkwords;

	/* 
    Index for the mixtures 
  */
	int idx = 0;
  
  /* 
    Auxiliar parameter 
  */
  Param_Index_and_Model Index_and_Model;

	/* 
    Allocate table with the mixture 
  */
  Mixture * List_Mixtures = (Mixture *)malloc(Number_Soil_Water_Mixtures*sizeof(Mixture));
  if(List_Mixtures == NULL)
  {
  	sprintf(Error_message,"%s","Memory error for table of mixtures");
		standard_error();
  }

  /* 
    Open and check file
  */
  Sim_dat = Open_and_Check_simulation_file(SimulationFile);

  while(fgets(line, sizeof line, Sim_dat) != NULL)
  {

    	nkwords = parse (kwords, line, delimiters_1);
    	if (nkwords < 0)
    	{
        sprintf(Error_message,"%s","Parser failed");
    		standard_error();
    	}

    	/* 
        Read Initial-nodal-values 
      */
    	if ((nkwords > 0) && (strcmp(kwords[0],"Define-Mixture") == 0 ))
    	{

    		/* 
          Read index and model 
        */
    		Index_and_Model = Read_Index_and_Model(kwords[1],kwords[2]);

    		/* 
          Read Mixture properties and asign to the Mixture 
        */
    		List_Mixtures[idx] = Define_Mixture(Sim_dat,Index_and_Model);

    		idx++;
    	}

	}

  /* Close  file */
  fclose(Sim_dat);

	return List_Mixtures;
}

/**********************************************************************/

static Param_Index_and_Model Read_Index_and_Model(char * String_Index, char * String_Model)
{

  int Num_Idx;
  char * Idx_pars[MAXW] = {NULL};

  int Num_Model;
  char * Model_pars[MAXW] = {NULL};

  /* Define outputs */
  Param_Index_and_Model Parameters;

  Num_Idx   = parse (Idx_pars,String_Index,delimiters_3);
  Num_Model = parse (Model_pars,String_Model,delimiters_3);

  /* Fill index of the Mixture */
  if(Num_Idx == 2)
  {
    Parameters.Idx = atoi(Idx_pars[1]);
  }
  else if(Num_Idx == 1)
  {
    Parameters.Idx = atoi(Idx_pars[0]);
  }    

  /* Fill model information */
  if(Num_Model == 2)
  {
    strcpy(Parameters.Model,Model_pars[1]);
  }
  else if(Num_Model == 1)
  {
    strcpy(Parameters.Model,Model_pars[0]);
  }

  /* Return outputs */
  return Parameters;
}

/***************************************************************************/

static Mixture Define_Mixture(FILE * Simulation_file, 
                                Param_Index_and_Model Index_and_Model)
{

  int Ndim = NumberDimensions;

  /* Define outputs */
  Mixture New_Mixture;

  /* Variables for reading purposes */
  char Parameter_line[MAXC] = {0};
  char * Parameter_pars[MAXW] = {NULL};
  int Parser_status;

  /* Check variables for sintax */
  bool Is_Open = false;
  bool Is_Close = false;
  Check_Mixture ChkMix = Initialise_Check_Mixture();

  /*
    Auxiliar scalar to store the permeability tensor (isotropic)
  */
  double Permeability_Isotropic;

  while(fgets(Parameter_line, sizeof(Parameter_line), Simulation_file) != NULL)
  {
    /* Parse line */    
    Parser_status = parse(Parameter_pars,Parameter_line,delimiters_2);

    if((strcmp(Parameter_pars[0],"{") == 0) && (Parser_status == 1))
    {
      Is_Open = true;
    }
    else if(strcmp(Parameter_pars[0],"Soil-Idx") == 0)
    {
      ChkMix.Is_Soil_Idx = true;
      New_Mixture.Soil_Idx = atoi(Parameter_pars[1]);
    }
    else if(strcmp(Parameter_pars[0],"Fluid-Idx") == 0)
    {
      ChkMix.Is_Water_Idx = true;
      New_Mixture.Water_Idx = atoi(Parameter_pars[1]);
    }
    else if(strcmp(Parameter_pars[0],"Permeability-Skeleton") == 0)
    {
      ChkMix.Is_Permeability = true;
      Permeability_Isotropic = atof(Parameter_pars[1]);
      New_Mixture.Permeability = alloc__TensorLib__(2);
      for(int i = 0 ; i<Ndim ; i++)
      {
        New_Mixture.Permeability.N[i][i] = Permeability_Isotropic;
      }
    }
    else if(strcmp(Parameter_pars[0],"Reference-Volume-Fraction-Soil") == 0)
    {
      ChkMix.Is_phi_s_0 = true;
      New_Mixture.phi_s_0 = atof(Parameter_pars[1]);
    } 
    else if(strcmp(Parameter_pars[0],"Reference-Volume-Fraction-Fluid") == 0)
    {
      ChkMix.Is_phi_f_0 = true;
      New_Mixture.phi_f_0 = atof(Parameter_pars[1]);
    }
    else if((strcmp(Parameter_pars[0],"}") == 0) && (Parser_status == 1))
    {
        Is_Close = true;
        break;
    }
    else if(Parser_status > 0)
    {
      sprintf(Error_message,"%s %s","Undefined",Parameter_pars[0]);
      standard_error(); 
    }
  
  }

  strcpy(New_Mixture.Type,Index_and_Model.Model);

  if(strcmp(New_Mixture.Type,"Borja-2004") == 0)
  {
    check_Borja_2004_Soil_water_Mixture(New_Mixture,ChkMix);
  }
  else
  {
    sprintf(Error_message,"%s","Unrecognized kind of Mixture");
    standard_error();
  }

  /* Return outputs */
  return New_Mixture;
}

/***************************************************************************/

static bool Activate_Options(char * Option, char * status_text)
{
  bool status;

  if(strcmp(status_text,"true") == 0)
  {
    printf("\t -> %s : True \n", Option);
    return true;
  }
  else if(strcmp(status_text,"false") == 0)
  {
    printf("\t -> %s : False \n", Option);
    return false;
  }
  else
  {
    sprintf(Error_message,"The status was %s. Please, use : true/false",status_text);
    standard_error(); 
  }

  return status;
}

/***************************************************************************/

static Check_Mixture Initialise_Check_Mixture()
{
  Check_Mixture ChkMix;

  ChkMix.Is_Soil_Idx = false;
  ChkMix.Is_Water_Idx = false;
  ChkMix.Is_Permeability = false;
  ChkMix.Is_phi_s_0 = false;
  ChkMix.Is_phi_f_0 = false;

  return ChkMix;
}

/**********************************************************************/

static void check_Borja_2004_Soil_water_Mixture(Mixture Mix_particle, Check_Mixture ChkMix)
{

    if(ChkMix.Is_Soil_Idx && 
      ChkMix.Is_Water_Idx && 
      ChkMix.Is_Permeability && 
      ChkMix.Is_phi_s_0 && 
      ChkMix.Is_phi_f_0)
    {
      printf("\t \t -> %s : %i \n","Material model for the Soil phase",Mix_particle.Soil_Idx);
      printf("\t \t -> %s : %i \n","Material model for the fluid phase",Mix_particle.Water_Idx);
      printf("\t \t -> %s : %f \n","Volume fraction (Soil phase)",Mix_particle.phi_s_0);
      printf("\t \t -> %s : %f \n","Volume fraction (fluid phase)",Mix_particle.phi_f_0);
    }
    else
    {
      fprintf(stderr,"%s : %s %s %s\n",
        "Error in Define-Mixture()","Some parameter is missed for",Mix_particle.Type,"Mixture");
      fputs(ChkMix.Is_Soil_Idx  ? "Material model (Soil phase) : true \n" : "Material model (Soil phase) : false \n", stdout);
      fputs(ChkMix.Is_Water_Idx  ? "Material model (fluid phase) : true \n" : "Material model (fluid phase) : false \n", stdout);
      fputs(ChkMix.Is_Permeability  ? "Permeabily Soil skeleton : true \n" : "Permeabily Soil skeleton : false \n", stdout);
      fputs(ChkMix.Is_phi_s_0  ? "Volume fraction (Soil phase) : true \n" : "Volume fraction (Soil phase) : false \n", stdout);
      fputs(ChkMix.Is_phi_f_0  ? "Volume fraction (fluid phase) : true \n" : "Volume fraction (fluid phase) : false \n", stdout);
      exit(EXIT_FAILURE);
    }
}

/**********************************************************************/

static void standard_error()
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in Define-Mixture()",Error_message);
    exit(EXIT_FAILURE);
}

/***************************************************************************/

static void standard_output(char * Status_message)
{
  fprintf(stdout,"%s \n",Status_message);
}

/**********************************************************************/

static FILE * Open_and_Check_simulation_file(char * Name_File)
{
  FILE * Simulation_file = fopen(Name_File,"r");  
  
  if (Simulation_file==NULL)
  {
    sprintf(Error_message,"%s %s","Incorrect lecture of",Name_File);
    standard_error(); 
  }  

  return Simulation_file;
}

/***************************************************************************/