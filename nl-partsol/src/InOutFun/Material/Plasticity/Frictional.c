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
  bool Is_alpha_Hardening_Borja; // Hardening parameter
  bool Is_a1_Hardening_Borja; // Hardening parameter
  bool Is_a2_Hardening_Borja; // Hardening parameter
  bool Is_a3_Hardening_Borja; // Hardening parameter
  bool Is_friction_angle; // Friction angle
  bool Is_dilatancy_angle; // Dilatancy angle

} Check_Material;

static void standard_error();
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

  TOL_Radial_Returning = 1E-10;
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
  ChkMat.Is_alpha_Hardening_Borja = false; // Hardening parameter
  ChkMat.Is_a1_Hardening_Borja = false; // Hardening parameter
  ChkMat.Is_a2_Hardening_Borja = false; // Hardening parameter
  ChkMat.Is_a3_Hardening_Borja = false; // Hardening parameter
  ChkMat.Is_friction_angle = false; // Friction angle
  ChkMat.Is_dilatancy_angle = false; // Dilatancy angle

  return ChkMat;
}

/**********************************************************************/

static void check_Frictional_Material(Material Mat_particle, Check_Material ChkMat, int Idx)
{
	if(ChkMat.Is_rho 
    && ChkMat.Is_E 
    && ChkMat.Is_nu 
    && ChkMat.Is_Yield_Function_Frictional)
	{
		printf("\t -> %s \n","Frictional material");
		printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
		printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
		printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
		printf("\t \t -> %s : %s \n","Yield function",Mat_particle.Yield_Function_Frictional);

		if(strcmp(Mat_particle.Yield_Function_Frictional,"Matsuoka-Nakai") == 0)
		{

			if(Mat_particle.Hardening_Borja)
			{
				if(ChkMat.Is_alpha_Hardening_Borja 
                && ChkMat.Is_a1_Hardening_Borja 
                && ChkMat.Is_a2_Hardening_Borja 
                && ChkMat.Is_a3_Hardening_Borja)
				{
					printf("\t \t -> %s : %f \n","alpha-Hardening-Borja",Mat_particle.alpha_Hardening_Borja);
					printf("\t \t -> %s : %f \n","a1-Hardening-Borja",Mat_particle.a_Hardening_Borja[0]);	
					printf("\t \t -> %s : %f \n","a2-Hardening-Borja",Mat_particle.a_Hardening_Borja[1]);	
					printf("\t \t -> %s : %f \n","a3-Hardening-Borja",Mat_particle.a_Hardening_Borja[2]);	
				}
				else
				{
					fprintf(stderr,"%s : %s \n","Error in GramsMaterials()",
					"Some parameter is missed for Borja Hardening");
					fputs(ChkMat.Is_alpha_Hardening_Borja ? "alpha-Hardening-Borja : true \n" : "alpha-Hardening-Borja : false \n", stdout);
					fputs(ChkMat.Is_a1_Hardening_Borja ? "a1-Hardening-Borja : true \n" : "a1-Hardening-Borja : false \n", stdout);
					fputs(ChkMat.Is_a2_Hardening_Borja ? "a2-Hardening-Borja : true \n" : "a2-Hardening-Borja : false \n", stdout);
					fputs(ChkMat.Is_a3_Hardening_Borja ? "a3-Hardening-Borja : true \n" : "a3-Hardening-Borja : false \n", stdout);
				}

			}

		}
		else if(strcmp(Mat_particle.Yield_Function_Frictional,"Lade-Duncan") == 0)
		{


				if(ChkMat.Is_alpha_Hardening_Borja 
                && ChkMat.Is_a1_Hardening_Borja 
                && ChkMat.Is_a2_Hardening_Borja 
                && ChkMat.Is_a3_Hardening_Borja)
				{
					printf("\t \t -> %s : %f \n","alpha-Hardening-Borja",Mat_particle.alpha_Hardening_Borja);
					printf("\t \t -> %s : %f \n","a1-Hardening-Borja",Mat_particle.a_Hardening_Borja[0]);	
					printf("\t \t -> %s : %f \n","a2-Hardening-Borja",Mat_particle.a_Hardening_Borja[1]);	
					printf("\t \t -> %s : %f \n","a3-Hardening-Borja",Mat_particle.a_Hardening_Borja[2]);	
				}
				else
				{
					fprintf(stderr,"%s : %s \n","Error in GramsMaterials()",
					"Some parameter is missed for Borja Hardening");
					fputs(ChkMat.Is_alpha_Hardening_Borja ? "alpha-Hardening-Borja : true \n" : "alpha-Hardening-Borja : false \n", stdout);
					fputs(ChkMat.Is_a1_Hardening_Borja ? "a1-Hardening-Borja : true \n" : "a1-Hardening-Borja : false \n", stdout);
					fputs(ChkMat.Is_a2_Hardening_Borja ? "a2-Hardening-Borja : true \n" : "a2-Hardening-Borja : false \n", stdout);
					fputs(ChkMat.Is_a3_Hardening_Borja ? "a3-Hardening-Borja : true \n" : "a3-Hardening-Borja : false \n", stdout);
				}

		}
		else if(strcmp(Mat_particle.Yield_Function_Frictional,"Modified-Lade-Duncan") == 0)
		{

			if(ChkMat.Is_alpha_Hardening_Borja 
            && ChkMat.Is_a1_Hardening_Borja 
            && ChkMat.Is_a2_Hardening_Borja 
            && ChkMat.Is_a3_Hardening_Borja)
			{
				printf("\t \t -> %s : %f \n","alpha-Hardening-Borja",Mat_particle.alpha_Hardening_Borja);
				printf("\t \t -> %s : %f \n","a1-Hardening-Borja",Mat_particle.a_Hardening_Borja[0]);	
				printf("\t \t -> %s : %f \n","a2-Hardening-Borja",Mat_particle.a_Hardening_Borja[1]);	
				printf("\t \t -> %s : %f \n","a3-Hardening-Borja",Mat_particle.a_Hardening_Borja[2]);	
			}
			else
			{
				fprintf(stderr,"%s : %s \n","Error in GramsMaterials()",
				"Some parameter is missed for Borja Hardening");
				fputs(ChkMat.Is_alpha_Hardening_Borja ? "alpha-Hardening-Borja : true \n" : "alpha-Hardening-Borja : false \n", stdout);
				fputs(ChkMat.Is_a1_Hardening_Borja ? "a1-Hardening-Borja : true \n" : "a1-Hardening-Borja : false \n", stdout);
				fputs(ChkMat.Is_a2_Hardening_Borja ? "a2-Hardening-Borja : true \n" : "a2-Hardening-Borja : false \n", stdout);
				fputs(ChkMat.Is_a3_Hardening_Borja ? "a3-Hardening-Borja : true \n" : "a3-Hardening-Borja : false \n", stdout);
                exit(EXIT_FAILURE);
			}

		}
		else
		{
			fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()","Unrecognized yield function");
			exit(EXIT_FAILURE);	
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
		exit(EXIT_FAILURE);
	}
}

/***************************************************************************/

static void standard_error()
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in Define-Material()",Error_message);
    exit(EXIT_FAILURE);
}

/***************************************************************************/