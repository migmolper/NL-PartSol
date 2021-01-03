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

typedef struct
{

  Matrix Value;
  Check_Strain ChkStr;

} Strain_curve;

typedef struct 
{

  int i_start;
  int i_step;
  int i_end;

} Intervals;

typedef struct 
{

  char File_Name[MAXC];
  char Directory_Name[MAXC];
  bool Is_Current_directory;
  bool Is_Defined_directory;

  bool Out_Stress;
  bool Out_Strain;
  bool Out_Deformation_gradient;
  bool Out_Plastic_Deformation_gradient;
  bool Out_EPS;
  bool Out_Cohesion;

} Parameters;

/*
  Auxiliar functions and variables
*/
#ifdef _WIN32
static char * delimiters_1 = " ()\r\n\t";
static char * delimiters_2 = " =\t\r\n"; 
#else
static char * delimiters_1 = " ()\n\t";
static char * delimiters_2 = " =\t\n"; 
#endif
static char * delimiters_3 = "=";
static char * delimiters_4 = ";";

static char Error_message[MAXW];

static Strain_curve * Read_Strains_curves(Fields, char *, int);
static Check_Strain Read_Strain_Measure(char *);
static Matrix Read_Curve_File_Name(FILE *);
static void Fill_Strains(Fields, Strain_curve);
static int Number_Output_directives(char * );
static void Read_Output_directives(char *);
static Intervals Read_CSV_Intervals(char *);
static Parameters Read_CSV_Parameters(FILE *);
static bool Is_Output_Activate(char *, char *);
static Event Fill_CSV_Output_Directives(Intervals, Parameters);
static bool Check_Output_directory(char *);
static bool Check_File(char *);
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

  Read_Output_directives(SimulationFile);
  

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

  /* Close  file */
  fclose(Sim_dat);

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
      Is_File = Check_File(File_Name); 
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

static int Number_Output_directives(char * Name_File)
{
  int NOD = 0;

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  /* Simulation file */
  FILE * Sim_dat = Open_and_Check_simulation_file(Name_File);
  
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL )
    {
      nkwords = parse (kwords, line, delimiters_1);

      if ((nkwords > 0) && (strcmp(kwords[0],"Particle-variables-evolution-to-csv") == 0 ))
      {
        NOD++;
      }
    }

  fclose(Sim_dat);

  return NOD;
}

/***************************************************************************/

static void Read_Output_directives(char * SimulationFile)
/*
  Particle-variables-evolution-to-csv (i_start=0;i_step=1;i_end=4000)
  {
    File=Strain-Stress.csv
    Directory=Foo

  } 
*/
{
  /* Simulation file */
  FILE * Sim_dat = Open_and_Check_simulation_file(SimulationFile);

  int Num_NOGPE = Number_Output_directives(SimulationFile);
  Number_Out_Gauss_Point_evolution_csv = Num_NOGPE;
  Out_Gauss_Point_evolution_csv = (Event *)malloc(Num_NOGPE*sizeof(Event));

  Intervals CSV_Intervals;
  Parameters CSV_Parameters;

  int idx = 0;

    /* Variables for reading purposes */
  char line[MAXC] = {0};
  char * kwords[MAXW] = {NULL};
  int nkwords;

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
      if ((nkwords > 0) && (strcmp(kwords[0],"Particle-variables-evolution-to-csv") == 0 ))
      {
        printf("\t -> Read Output directive\n");

        CSV_Intervals = Read_CSV_Intervals(kwords[1]);

        CSV_Parameters = Read_CSV_Parameters(Sim_dat);  

        Out_Gauss_Point_evolution_csv[idx] = Fill_CSV_Output_Directives(CSV_Intervals,CSV_Parameters);

        idx++;           

      }

  }

  /* Close  file */
  fclose(Sim_dat);
}

/***************************************************************************/

static Intervals Read_CSV_Intervals(char * Interval_message)
{
  Intervals CSV_Intervals;
  
  int Interval_status_1;
  char * Aux_Parse_1[MAXW] = {NULL};
  int Interval_status_2;
  char * Aux_Parse_2[MAXW] = {NULL};

  /* Set parameters to default */
  CSV_Intervals.i_start = 0;
  CSV_Intervals.i_step = 1;
  CSV_Intervals.i_end = NumTimeStep;

  Interval_status_1 = parse (Aux_Parse_1, Interval_message,delimiters_4);

  /* Check format */
  if(Interval_status_1 > 3)
    {
      sprintf(Error_message,"You have exceded the maximum number of parameters");
      standard_error();
    }
  

  for(int i = 0 ; i<Interval_status_1 ; i++)
  {
    Interval_status_2 = parse (Aux_Parse_2, Aux_Parse_1[i],delimiters_3); 

    if((strcmp(Aux_Parse_2[0],"i_start") == 0) && (Interval_status_2 == 2))
      {
        CSV_Intervals.i_start  = atoi(Aux_Parse_2[1]);
      }
   
      else if((strcmp(Aux_Parse_2[0],"i_step") == 0) && (Interval_status_2 == 2))
      {
        CSV_Intervals.i_step  = atoi(Aux_Parse_2[1]);
      }
  
      else if((strcmp(Aux_Parse_2[0],"i_end") == 0) && (Interval_status_2 == 2))
      {
        CSV_Intervals.i_end  = atoi(Aux_Parse_2[1]);
      }
      else
      {
        sprintf(Error_message,"The statement %s is not recognised",Aux_Parse_2[0]);
        standard_error();
      }
  }
   
    /* Check interval output */
    if(CSV_Intervals.i_end < CSV_Intervals.i_step)
      {
        sprintf(Error_message,"The result interval step should be less than final time");
        standard_error();
      }
    if(CSV_Intervals.i_end < CSV_Intervals.i_start)
      {
        sprintf(Error_message,"The initial result time step should be less than the final result time");
        standard_error();
      }       
    if(CSV_Intervals.i_end > NumTimeStep)
      {
        sprintf(Error_message,"The final result time step should be less than the final time");
        standard_error();
      }

  /* Print some info */
  printf("\t \t -> %s : %i \n","i_start",CSV_Intervals.i_start);
  printf("\t \t -> %s : %i \n","i_step",CSV_Intervals.i_step);
  printf("\t \t -> %s : %i \n","i_end",CSV_Intervals.i_end);

  return CSV_Intervals;
}

/***************************************************************************/

static Parameters Read_CSV_Parameters(FILE * Simulation_file)
{
  Parameters Output_csv;

  /* Set outputs to default */
  Output_csv.Is_Current_directory = true;
  Output_csv.Is_Defined_directory = false;

  /* Variables for reading purposes */
  char Line_Out_Prop[MAXC] = {0};
  char * Parse_Out_Prop[MAXW] = {NULL};
  int Aux_Out_id;

  /* Check variables for sintax */
  bool Is_Open = false;
  bool Is_Close = false;
  bool Is_File = false;

  while(fgets(Line_Out_Prop, sizeof(Line_Out_Prop), Simulation_file) != NULL)
    {
      /* Parse line */    
    Aux_Out_id = parse(Parse_Out_Prop,Line_Out_Prop,delimiters_2);

    if((strcmp(Parse_Out_Prop[0],"{") == 0) && (Aux_Out_id == 1))
    {
      Is_Open = true;
    }
    else if((strcmp(Parse_Out_Prop[0],"Directory") == 0) && (Aux_Out_id == 2)) 
    {
        sprintf(Output_csv.Directory_Name,"./%s/",Parse_Out_Prop[1]);
        printf("\t \t -> %s : %s \n","Directory",Output_csv.Directory_Name);
        Output_csv.Is_Current_directory = false;
        Output_csv.Is_Defined_directory = Check_Output_directory(Output_csv.Directory_Name);
    }
    else if((strcmp(Parse_Out_Prop[0],"File") == 0) && (Aux_Out_id == 2))
    {
        sprintf(Output_csv.File_Name,"%s",Parse_Out_Prop[1]);
        printf("\t \t -> %s : %s \n","File",Parse_Out_Prop[1]);
        Is_File = true;
    }    
    else if(strcmp(Parse_Out_Prop[0],"Out-stress") == 0)
    {
        Output_csv.Out_Stress = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
    }
    else if(strcmp(Parse_Out_Prop[0],"Out-strain") == 0)
    {
        Output_csv.Out_Strain = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
    }
    else if(strcmp(Parse_Out_Prop[0],"Out-deformation-gradient") == 0)
    {
        Output_csv.Out_Deformation_gradient = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
    }
    else if(strcmp(Parse_Out_Prop[0],"Out-plastic-deformation-gradient") == 0)
    {
        Output_csv.Out_Plastic_Deformation_gradient = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
    }
    else if(strcmp(Parse_Out_Prop[0],"Out-EPS") == 0)
    {
        Output_csv.Out_EPS = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
    }
    else if(strcmp(Parse_Out_Prop[0],"Out-Cohesion") == 0)
    {
        Output_csv.Out_Cohesion = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
    }
    else if((strcmp(Parse_Out_Prop[0],"}") == 0) && (Aux_Out_id == 1))
    {
        Is_Close = true;
        break;
    }   
    else if(Aux_Out_id > 0)
      {
        sprintf(Error_message,"%s %s","Undefined",Parse_Out_Prop[0]);
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

  /* Check syntax */

  return Output_csv;
}

/***************************************************************************/

static bool Is_Output_Activate(char * output_field, char * status_text)
{
  bool status;

  if(strcmp(status_text,"true") == 0)
  {
    printf("\t -> %s : true \n", output_field);
    return true;
  }
  else if(strcmp(status_text,"false") == 0)
  {
    printf("\t -> %s : False \n", output_field);
    return false;
  }
  else
  {
    sprintf(Error_message,"The input was %s. Please, use : true/false",status_text);
    standard_error(); 
  }

  return status;
}

/***************************************************************************/

static Event Fill_CSV_Output_Directives(Intervals CSV_Intervals, Parameters CSV_Parameters)
{
  Event CSV_Event;

  /* Set outputs to default */
  CSV_Event.Out_csv_Gauss_Point_evolution_Stress = false;
  CSV_Event.Out_csv_Gauss_Point_evolution_Strain = false;
  CSV_Event.Out_csv_Gauss_Point_evolution_Deformation_gradient = false;
  CSV_Event.Out_csv_Gauss_Point_evolution_Plastic_Deformation_gradient = false;
  CSV_Event.Out_csv_Gauss_Point_evolution_EPS = false;
  CSV_Event.Out_csv_Gauss_Point_evolution_Cohesion = false;


  /* Write name of the ouput file directory */
  if(CSV_Parameters.Is_Current_directory)
  {
    sprintf(CSV_Event.Directory,"./%s.csv",CSV_Parameters.File_Name);
  }
  else if(CSV_Parameters.Is_Defined_directory)
  {
    sprintf(CSV_Event.Directory,"./%s/%s.csv",CSV_Parameters.Directory_Name,CSV_Parameters.File_Name);
  }

  /* Read outputs intervals */
  CSV_Event.i_start = CSV_Intervals.i_start;
  CSV_Event.i_step  = CSV_Intervals.i_step;
  CSV_Event.i_end   = CSV_Intervals.i_end;
  
  /* Select ouput variable */ 
  if(CSV_Parameters.Out_Stress)
  {
    CSV_Event.Out_csv_Gauss_Point_evolution_Stress = true;
  }
  else if(CSV_Parameters.Out_Strain)
  {
    CSV_Event.Out_csv_Gauss_Point_evolution_Strain = true;
  }
  else if(CSV_Parameters.Out_Deformation_gradient)
  {
    CSV_Event.Out_csv_Gauss_Point_evolution_Deformation_gradient = true;
  }
  else if(CSV_Parameters.Out_Plastic_Deformation_gradient)
  {
    CSV_Event.Out_csv_Gauss_Point_evolution_Plastic_Deformation_gradient = true;
  }
  else if(CSV_Parameters.Out_EPS)
  {
    CSV_Event.Out_csv_Gauss_Point_evolution_EPS = true;
  }
  else if(CSV_Parameters.Out_Cohesion)
  {
    CSV_Event.Out_csv_Gauss_Point_evolution_Cohesion = true;
  }


  return CSV_Event;
}


/***************************************************************************/

static bool Check_Output_directory(char * Output_directory)
{
  struct stat info;
  stat(Output_directory,&info);
  char Error_message[MAXW];
  bool status_check;

  if(S_ISDIR(info.st_mode))
  {
    printf("\t -> %s : %s \n","Output directory",Output_directory);
    status_check = true;
  }
  else
  {
    sprintf(Error_message,"%s : %s %s \n","Output directory",Output_directory,"does not exists");
    standard_error(); 
  } 

  return status_check;
}

/***************************************************************************/

static bool Check_File(char * Path_File)
{
  struct stat info;
  stat(Path_File,&info);
  bool status_check;

  if(S_ISREG(info.st_mode))
  {
    printf("\t \t -> %s : %s \n","File",Path_File);
    status_check = true;
  }
  else
  {
    sprintf(Error_message,"%s %s %s","File",Path_File,"does not exists");
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

