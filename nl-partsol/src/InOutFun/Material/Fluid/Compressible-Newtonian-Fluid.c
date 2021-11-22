#include "nl-partsol.h"

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

typedef struct
{
    bool Is_rho; 
    bool Is_Compressibility;
    bool Is_ReferencePressure;
    bool Is_Viscosity;
    bool Is_n_Macdonald_model;
    bool Is_Locking_Control_Fbar; // For incompressible limit
    bool Is_alpha_Fbar;

} Check_Material;

static bool Activate_Options(char *);
static void standard_error();
static Check_Material Initialise_Check_Material();
static void check_Compressible_Newtonian_Fluid_Material(Material,Check_Material,int);

/***************************************************************************/

Material Define_Compressible_Newtonian_Fluid(
  FILE * Simulation_file,
  char * Material_Model,
  int Material_Idx)
{

  int Ndim = NumberDimensions;
  /* Define outputs */
  Material New_Material;

  /* Variables for reading purposes */
  char Parameter_line[MAXC] = {0};
  char * Parameter_pars[MAXW] = {NULL};
  int Parser_status;

  /* Check variables for sintax */
  bool Is_Open = false;
  bool Is_Close = false;
  Check_Material ChkMat = Initialise_Check_Material();

  /* Default parameters */
  New_Material.Locking_Control_Fbar = false;
  New_Material.alpha_Fbar = 0.0;

  while(fgets(Parameter_line, sizeof(Parameter_line), Simulation_file) != NULL)
  {
    /* Parse line */    
    Parser_status = parse(Parameter_pars,Parameter_line,delimiters_2);

    if((strcmp(Parameter_pars[0],"{") == 0) && (Parser_status == 1))
    {
      Is_Open = true;
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"rho") == 0)
    {
      ChkMat.Is_rho = true;
      New_Material.rho = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Compressibility") == 0)
    {
      ChkMat.Is_Compressibility = true;
      New_Material.Compressibility = atof(Parameter_pars[1]);
    }
	  /**************************************************/
	  else if(strcmp(Parameter_pars[0],"Reference-Pressure") == 0)
	  {
	  	ChkMat.Is_ReferencePressure = true;
	  	New_Material.ReferencePressure = atof(Parameter_pars[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parameter_pars[0],"Viscosity") == 0)
	  {
	  	ChkMat.Is_Viscosity = true;
	  	New_Material.Viscosity = atof(Parameter_pars[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parameter_pars[0],"Macdonald-parameter") == 0)
	  {
	  	ChkMat.Is_n_Macdonald_model = true;
	  	New_Material.n_Macdonald_model = atof(Parameter_pars[1]);
	  }
    else if((strcmp(Parameter_pars[0],"}") == 0) && (Parser_status == 1))
    {
        Is_Close = true;
        break;
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Fbar") == 0)
	  {
      ChkMat.Is_Locking_Control_Fbar = true;
      New_Material.Locking_Control_Fbar = Activate_Options(Parameter_pars[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parameter_pars[0],"Fbar-alpha") == 0)
	  {
      ChkMat.Is_alpha_Fbar = true;
			New_Material.alpha_Fbar = atof(Parameter_pars[1]);

			if((New_Material.alpha_Fbar < 0.0) 
      || (New_Material.alpha_Fbar > 1.0))
			{
				sprintf(Error_message,"The range for Fbar-alpha is [0,1]");
	   		standard_error(Error_message); 
			}
	  }    
    /**************************************************/
    else if(Parser_status > 0)
    {
      sprintf(Error_message,"%s %s","Undefined",Parameter_pars[0]);
      standard_error(); 
    }
  
  }

  strcpy(New_Material.Type,Material_Model);
  
  check_Compressible_Newtonian_Fluid_Material(New_Material,ChkMat,Material_Idx);

  /* Return outputs */
  return New_Material;
}

/***************************************************************************/

static Check_Material Initialise_Check_Material()
{
  Check_Material ChkMat;

  ChkMat.Is_rho = false; // Reference fensity
  ChkMat.Is_Compressibility = false; // Bulk stiffness
  ChkMat.Is_ReferencePressure = false;
  ChkMat.Is_Viscosity = false;
  ChkMat.Is_n_Macdonald_model = false;
  ChkMat.Is_Locking_Control_Fbar = false;
  ChkMat.Is_alpha_Fbar = false;

  return ChkMat;
}


/**********************************************************************/

static void check_Compressible_Newtonian_Fluid_Material(Material Mat_particle, Check_Material ChkMat, int Idx)
{
  if(ChkMat.Is_rho 
  && ChkMat.Is_Compressibility
  && ChkMat.Is_ReferencePressure
  && ChkMat.Is_Viscosity
  && ChkMat.Is_n_Macdonald_model)
  {
    printf("\t -> %s \n","Compressible Newtonian Fluid material");
    printf("\t \t -> %s : %i \n","Idx",Idx);
    printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
    printf("\t \t -> %s : %f \n","Compressibility",Mat_particle.Compressibility);
    printf("\t \t -> %s : %f \n","Reference Pressure",Mat_particle.ReferencePressure);
    printf("\t \t -> %s : %f \n","Viscosity",Mat_particle.Viscosity);
    printf("\t \t -> %s : %f \n","Macdonald-parameter",Mat_particle.n_Macdonald_model);

    if(ChkMat.Is_Locking_Control_Fbar)
    {
      printf("\t \t -> %s : %s \n","F-bar","Enabled");

      if(ChkMat.Is_alpha_Fbar)
      {
        printf("\t \t -> %s : %f \n","alpha F-bar",Mat_particle.alpha_Fbar);
      }
    }
    else
    {
      printf("\t \t -> %s : %s \n","F-bar","Disabled");
    }

  }
  else
  {
    fprintf(stderr,"%s : %s \n",
      "Error in GramsMaterials()",
      "Some parameter is missed for Compressible Newtonian Fluid material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    fputs(ChkMat.Is_Compressibility ? "Compressibility : true \n" : "Compressibility : false \n", stdout);
    fputs(ChkMat.Is_ReferencePressure ? "Reference Pressure : true \n" : "Reference Pressure : false \n", stdout);
    fputs(ChkMat.Is_Viscosity ? "Viscosity : true \n" : "Viscosity : false \n", stdout);
    fputs(ChkMat.Is_n_Macdonald_model ? "Macdonald-parameter : true \n" : "Macdonald-parameter : false \n", stdout);
    exit(EXIT_FAILURE);
  }
}

/***************************************************************************/

static bool Activate_Options(char * status_text)
{
  bool status;

  if(strcmp(status_text,"true") == 0)
  {
    return true;
  }
  else if(strcmp(status_text,"false") == 0)
  {
    return false;
  }
  else
  {
    sprintf(Error_message,"The status was %s. Please, use : true/false",status_text);
    standard_error(); 
  }

  return status;
}

/**********************************************************************/

static void standard_error()
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in Define-Material()",Error_message);
    exit(EXIT_FAILURE);
}

/***************************************************************************/