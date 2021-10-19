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

  bool Is_rho; // Reference fensity
  bool Is_E; // Young modulus
  bool Is_nu; // Poisson cefficient
  bool Is_Yield_Function_Frictional; // Kind of yield function
  bool Is_m; // 
  bool Is_c0; //
  bool Is_alpha_Hardening_Borja; // Hardening parameter
  bool Is_a1_Hardening_Borja; // Hardening parameter
  bool Is_a2_Hardening_Borja; // Hardening parameter
  bool Is_a3_Hardening_Borja; // Hardening parameter
  bool Is_atmospheric_pressure; // Reference pressure
  bool Is_friction_angle; // Friction angle
  bool Is_dilatancy_angle; // Dilatancy angle
  bool Is_Locking_Control_Fbar; // Locking control
  bool Is_alpha_Fbar; // Tunning paramer for the F-bar

} Check_Material;

static void standard_error();
static bool Activate_Options(char *);
static Check_Material Initialise_Check_Material();
static void check_Frictional_Material(Material,Check_Material,int);

/**********************************************************************/

Material Define_Frictional(
  FILE * Simulation_file,
  char * Material_Model,
  int Material_Idx)
{

  int Ndim = NumberDimensions;
  /* Define outputs */
  Material Frictional_Material;

  /* Variables for reading purposes */
  char Parameter_line[MAXC] = {0};
  char * Parameter_pars[MAXW] = {NULL};
  int Parser_status;

  /* Check variables for sintax */
  bool Is_Open = false;
  bool Is_Close = false;
  Check_Material ChkMat = Initialise_Check_Material();

  /* Auxiliar parameters for granular materials */
  double rad_friction_angle;
  double rad_dilatancy_angle;

  /* Default parameters */
  Frictional_Material.Locking_Control_Fbar = false;
  Frictional_Material.alpha_Fbar = 0.0;

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
      Frictional_Material.rho = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"E") == 0)
    {
      ChkMat.Is_E = true;
      Frictional_Material.E = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"nu") == 0)
    {
      ChkMat.Is_nu = true;
      Frictional_Material.nu = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Yield-Function") == 0)
    {

		  if((strcmp(Parameter_pars[1],"Matsuoka-Nakai") == 0) 
      || (strcmp(Parameter_pars[1],"Lade-Duncan") == 0) 
      || (strcmp(Parameter_pars[1],"Modified-Lade-Duncan") == 0))
      {
        strcpy(Frictional_Material.Yield_Function_Frictional,Parameter_pars[1]);
        ChkMat.Is_Yield_Function_Frictional = true;
      }
      else 
      {
        fprintf(stderr,"%s : %s \n",
        "Error in GramsMaterials()","Unrecognized yield function");
        exit(EXIT_FAILURE);	
      }

    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"m") == 0)
    {
      ChkMat.Is_m = true;
      Frictional_Material.m_Frictional = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"c0") == 0)
    {
      ChkMat.Is_c0 = true;
      Frictional_Material.c0_Frictional = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"alpha-Hardening-Borja") == 0)
    {
      ChkMat.Is_alpha_Hardening_Borja = true;
      Frictional_Material.alpha_Hardening_Borja = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"a1-Hardening-Borja") == 0)
    {
      ChkMat.Is_a1_Hardening_Borja = true;
      Frictional_Material.a_Hardening_Borja[0] = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"a2-Hardening-Borja") == 0)
    {
      ChkMat.Is_a2_Hardening_Borja = true;
      Frictional_Material.a_Hardening_Borja[1] = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"a3-Hardening-Borja") == 0)
    {
      ChkMat.Is_a3_Hardening_Borja = true;
      Frictional_Material.a_Hardening_Borja[2] = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Atmospheric-pressure") == 0)
    {
      ChkMat.Is_atmospheric_pressure = true;
      Frictional_Material.atmospheric_pressure = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Friction-angle") == 0)
    {
      ChkMat.Is_friction_angle = true;
      Frictional_Material.phi_Frictional = (PI__MatrixLib__/180)*atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Dilatancy-angle") == 0)
    {
      ChkMat.Is_dilatancy_angle = true;
      Frictional_Material.psi_Frictional = (PI__MatrixLib__/180)*atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Fbar") == 0)
	  {
      ChkMat.Is_Locking_Control_Fbar = true;
      Frictional_Material.Locking_Control_Fbar = Activate_Options(Parameter_pars[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parameter_pars[0],"Fbar-alpha") == 0)
	  {
      ChkMat.Is_alpha_Fbar = true;
			Frictional_Material.alpha_Fbar = atof(Parameter_pars[1]);

			if((Frictional_Material.alpha_Fbar < 0.0) 
      || (Frictional_Material.alpha_Fbar > 1.0))
			{
				sprintf(Error_message,"The range for Fbar-alpha is [0,1]");
	   		standard_error(Error_message); 
			}
	  }
	  /**************************************************/
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

  strcpy(Frictional_Material.Type,Material_Model);
  
  check_Frictional_Material(Frictional_Material,ChkMat,Material_Idx);

  TOL_Radial_Returning = 1E-12;
  Max_Iterations_Radial_Returning = 30;

  /* Return outputs */
  return Frictional_Material;
}

/***************************************************************************/

static Check_Material Initialise_Check_Material()
{
  Check_Material ChkMat;

  ChkMat.Is_rho = false; // Reference fensity
  ChkMat.Is_E = false; // Young modulus
  ChkMat.Is_nu = false; // Poisson cefficient
  ChkMat.Is_Yield_Function_Frictional = false; // Kind of yield function
  ChkMat.Is_m = false; // 
  ChkMat.Is_c0 = false; //
  ChkMat.Is_alpha_Hardening_Borja = false; // Hardening parameter
  ChkMat.Is_a1_Hardening_Borja = false; // Hardening parameter
  ChkMat.Is_a2_Hardening_Borja = false; // Hardening parameter
  ChkMat.Is_a3_Hardening_Borja = false; // Hardening parameter
  ChkMat.Is_atmospheric_pressure = false; // Reference pressure
  ChkMat.Is_friction_angle = false; // Friction angle
  ChkMat.Is_dilatancy_angle = false; // Dilatancy angle
  ChkMat.Is_Locking_Control_Fbar = false; // Locking control
  ChkMat.Is_alpha_Fbar = false; // Tunning parameter for the F-bar

  return ChkMat;
}

/**********************************************************************/

static void check_Frictional_Material(Material Mat_particle, Check_Material ChkMat, int Idx)
{
	if(ChkMat.Is_rho 
  && ChkMat.Is_E 
  && ChkMat.Is_nu 
  && ChkMat.Is_Yield_Function_Frictional
  && ChkMat.Is_atmospheric_pressure
  && ChkMat.Is_alpha_Hardening_Borja
  && ChkMat.Is_a1_Hardening_Borja
  && ChkMat.Is_a2_Hardening_Borja
  && ChkMat.Is_a3_Hardening_Borja
  && ChkMat.Is_friction_angle
  && ChkMat.Is_dilatancy_angle)
	{

		printf("\t -> %s \n","Frictional material");
		printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
		printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
		printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
		printf("\t \t -> %s : %s \n","Yield function",Mat_particle.Yield_Function_Frictional);
    printf("\t \t -> %s : %f \n","Atmospheric-pressure",Mat_particle.atmospheric_pressure);
    printf("\t \t -> %s : %f \n","alpha-Hardening-Borja",Mat_particle.alpha_Hardening_Borja);
    printf("\t \t -> %s : %f \n","a1-Hardening-Borja",Mat_particle.a_Hardening_Borja[0]);	
    printf("\t \t -> %s : %f \n","a2-Hardening-Borja",Mat_particle.a_Hardening_Borja[1]);	
    printf("\t \t -> %s : %f \n","a3-Hardening-Borja",Mat_particle.a_Hardening_Borja[2]);	
    printf("\t \t -> %s : %f \n","Friction-angle",Mat_particle.phi_Frictional*(180/PI__MatrixLib__));	
    printf("\t \t -> %s : %f \n","Dilatancy-angle",Mat_particle.psi_Frictional*(180/PI__MatrixLib__));	
    
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
        

    if(ChkMat.Is_m 
    && ChkMat.Is_c0)
    {
      printf("\t \t -> %s : %f \n","m",Mat_particle.m_Frictional);
      printf("\t \t -> %s : %f \n","c0",Mat_particle.c0_Frictional);
    }


	}
	else
	{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()",
			"Some parameter is missed for Frictional material");
		fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
		fputs(ChkMat.Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
		fputs(ChkMat.Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
		fputs(ChkMat.Is_Yield_Function_Frictional ? "Yield Function : true \n" : "Yield Function : false \n", stdout);
    fputs(ChkMat.Is_m ? "m : true \n" : "m : false \n", stdout);
    fputs(ChkMat.Is_c0 ? "c0 : true \n" : "c0 : false \n", stdout);
    fputs(ChkMat.Is_atmospheric_pressure ? "Atmospheric-pressure : true \n" : "Atmospheric-pressure : false \n", stdout);
		fputs(ChkMat.Is_alpha_Hardening_Borja ? "alpha-Hardening-Borja : true \n" : "alpha-Hardening-Borja : false \n", stdout);
		fputs(ChkMat.Is_a1_Hardening_Borja ? "a1-Hardening-Borja : true \n" : "a1-Hardening-Borja : false \n", stdout);
		fputs(ChkMat.Is_a2_Hardening_Borja ? "a2-Hardening-Borja : true \n" : "a2-Hardening-Borja : false \n", stdout);
		fputs(ChkMat.Is_a3_Hardening_Borja ? "a3-Hardening-Borja : true \n" : "a3-Hardening-Borja : false \n", stdout);
    fputs(ChkMat.Is_friction_angle ? "Friction-angle : true \n" : "Friction-angle: false \n", stdout);
    fputs(ChkMat.Is_dilatancy_angle ? "Dilatancy-angle : true \n" : "Dilatancy-angle: false \n", stdout);
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