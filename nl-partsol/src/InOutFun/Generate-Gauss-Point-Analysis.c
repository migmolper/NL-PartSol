#include "nl-partsol.h"


/*
  Call global variables
*/
Event * Out_Gauss_Point_evolution_csv;
int Number_Out_Gauss_Point_evolution_csv;

/*
  Auxiliar functions 
*/
static void  standard_error(char *);
static FILE * Open_and_Check_simulation_file(char *);

/**************************************************************/

GaussPoint Generate_Gauss_Point_Analysis__InOutFun__(char * SimulationFile)
{
	GaussPoint PointAnalysis;

  PointAnalysis.NumberMaterials = 1;
  PointAnalysis.Mat = Read_Materials(SimulationFile, 1);

	return PointAnalysis;
}

/**************************************************************/

static void standard_error(char * Error_message)
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in Generate_Gauss_Point_Analysis()",Error_message);
    exit(EXIT_FAILURE);
}

/**************************************************************/

static FILE * Open_and_Check_simulation_file(char * Name_File)
{
  FILE * Simulation_file = fopen(Name_File,"r");  
  char Error_message[MAXW];
  
  if (Simulation_file==NULL)
  {
  	sprintf(Error_message,"%s %s","Incorrect lecture of",Name_File);
	  standard_error(Error_message); 
  }  

  return Simulation_file;
}

/**************************************************************/