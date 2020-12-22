#include "nl-partsol.h"
#include <sys/stat.h>

/*
  Call global variables
*/
int NumTimeStep;
Event * Out_Gauss_Point_evolution_csv;
int Number_Out_Gauss_Point_evolution_csv;


/*
  Local structures
*/

typedef struct
{

  bool Is_Deformation_Gradient;
  bool Is_Right_Cauchy_Green;
  bool Is_Small_Strain;

} Check_Strain;

typedef struct {

  Matrix Value;

  Check_Strain ChkStr;

} Strain_curve;

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

static Strain_curve * Read_Strains_curves(Fields, char *, int);
static Check_Strain Read_Strain_Measure(char *);
static Matrix Read_Curve_File_Name(FILE *);
static void Fill_Strains(Fields, Strain_curve);
static bool Check_Path(char *);
static void standard_error();
static FILE * Open_and_Check_simulation_file(char *);

/***************************************************************************/

GaussPoint Generate_Gauss_Point_Analysis__InOutFun__(char * SimulationFile)
{
  int NumStrainCurves = 1;
  int NumberMaterials = 1;

	GaussPoint PointAnalysis;

  PointAnalysis.NumberMaterials = NumberMaterials;
  PointAnalysis.Mat = Read_Materials__InOutFun__(SimulationFile, NumberMaterials);

  PointAnalysis.Phi = allocate_Fields(NumTimeStep);

  Strain_curve * Strain_Case = Read_Strains_curves(PointAnalysis.Phi,SimulationFile,NumStrainCurves);

  for(int i = 0 ; i<NumStrainCurves ; i++)
  {
    Fill_Strains(PointAnalysis.Phi, Strain_Case[i]);
    free__MatrixLib__(Strain_Case[i].Value);
  }

  exit(1);

	return PointAnalysis;
}

/***************************************************************************/

static Strain_curve * Read_Strains_curves(Fields Phi, char * SimulationFile, int NumStrainCurves)
/*
 Define-Strain-Curve (Var=Deformation-Gradient)
  {
    File=Strain.txt
  }
*/
{
  /* Simulation file */
  FILE * Sim_dat = Open_and_Check_simulation_file(SimulationFile);

  /* Variables for reading purposes */
  char line[MAXC] = {0};
  char * kwords[MAXW] = {NULL};
  int nkwords;

  /* Index for the diferent strain curves */
  int idx = 0;

  Strain_curve * Strain_Case = (Strain_curve *)malloc(NumStrainCurves*sizeof(Strain_curve));

  while(fgets(line, sizeof line, Sim_dat) != NULL)
  {

    /* Read the line with the delimiter_1 */
      nkwords = parse (kwords, line, delimiters_1);
      if (nkwords < 0)
      {
        sprintf(Error_message,"%s","Parser failed");
        standard_error();
      }

      /* Read Define-Strain-Curve */
      if ((nkwords > 0) && (strcmp(kwords[0],"Define-Strain-Curve") == 0 ))
      {
        printf("\t -> Read Strain curve\n");

        Strain_Case[idx].ChkStr = Read_Strain_Measure(kwords[1]);

        Strain_Case[idx].Value = Read_Curve_File_Name(Sim_dat);        

        idx++;
      }

  }

  return Strain_Case;
}

/***************************************************************************/

static Check_Strain Read_Strain_Measure(char * Strain_Measure_char)
{
  Check_Strain ChkStr;
  char * kwords[MAXW] = {NULL};
  int nkwords;

  ChkStr.Is_Deformation_Gradient = false;
  ChkStr.Is_Right_Cauchy_Green = false;
  ChkStr.Is_Small_Strain = false;

  nkwords = parse(kwords, Strain_Measure_char, delimiters_3);

  if(strcmp(kwords[1],"Deformation-Gradient") == 0)
  {
    ChkStr.Is_Deformation_Gradient = true;
    printf("\t \t -> Kind of analysis : %s \n","Deformation-Gradient");
  }
  else if(strcmp(kwords[1],"Right-Cauchy-Green") == 0)
  {
    ChkStr.Is_Right_Cauchy_Green = true;
    printf("\t \t -> Kind of analysis : %s \n","Right-Cauchy-Green");
  }
  else if(strcmp(kwords[1],"Small-Strain") == 0)
  { 
    ChkStr.Is_Small_Strain = true;
    printf("\t \t -> Kind of analysis : %s \n","Small-Strain");
  }
  else
  {
    sprintf(Error_message,"%s","Strain meassure not recognised");
    standard_error();
  }


  return ChkStr;
}

/***************************************************************************/

static Matrix Read_Curve_File_Name(FILE * Simulation_file)
{
  Matrix Value;

  /* Variables for reading purposes */
  char Parameter_line[MAXC] = {0};
  char * Parameter_pars[MAXW] = {NULL};
  int Parser_status;

  char File_Name[MAXW];

  /* Check variables for sintax */
  bool Is_Open = false;
  bool Is_Close = false;
  bool Is_File = false;

  while(fgets(Parameter_line, sizeof(Parameter_line), Simulation_file) != NULL)
  {

    /* Parse line */    
    Parser_status = parse(Parameter_pars,Parameter_line,delimiters_2);

    if((strcmp(Parameter_pars[0],"{") == 0) && (Parser_status == 1))
    {
      Is_Open = true;
    }
    else if(strcmp(Parameter_pars[0],"File") == 0)
    {
      strcpy(File_Name,Parameter_pars[1]);
      Is_File = Check_Path(File_Name); 
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

  
  if(!Is_Open || !Is_Close)
  {
    sprintf(Error_message,"%s","Non balansed statement {}");
    standard_error(); 
  }
  else if(!Is_File)
  {
    sprintf(Error_message,"%s","None file provided");
    standard_error(); 
  }

  Value = Read_Delimited_File__InOutLib__(File_Name);

  return Value;
}

/***************************************************************************/

static void Fill_Strains(Fields Phi, Strain_curve Strain_i)
{

    if(Strain_i.ChkStr.Is_Deformation_Gradient)
    {
      Phi.F_n = increment__MatrixLib__(Phi.F_n, Strain_i.Value);
    }

    else if(Strain_i.ChkStr.Is_Right_Cauchy_Green)
    {
      Phi.Strain = increment__MatrixLib__(Phi.Strain, Strain_i.Value);
    }

    else if(Strain_i.ChkStr.Is_Small_Strain)
    {
      Phi.Strain = increment__MatrixLib__(Phi.Strain, Strain_i.Value);
    }

}

/***************************************************************************/

static bool Check_Path(char * PATH_Name)
{
  struct stat info;
  stat(PATH_Name,&info);
  bool status_check;

  if(S_ISREG(info.st_mode))
  {
    printf("\t \t -> %s : %s \n","Path file",PATH_Name);
    status_check = true;
  }
  else
  {
    sprintf(Error_message,"%s %s %s","Path file",PATH_Name,"does not exists");
    standard_error();
  } 

  return status_check;
}

/***************************************************************************/

static void standard_error()
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in Generate_Gauss_Point_Analysis()",Error_message);
    exit(EXIT_FAILURE);
}

/***************************************************************************/

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

