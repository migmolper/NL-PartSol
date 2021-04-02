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

  bool Is_rho; // Reference fensity
  bool Is_E; // Young modulus
  bool Is_nu; // Poisson cefficient
  bool Is_Compressibility; // Bulk stiffness
  bool Is_E_p0; // Initial plastic strain
  bool Is_yield_stress; // Initial Yield stress
  bool Is_H; // Hardening parameter (modulus)
  bool Is_Hexp; // Hardening parameter (exponent)
  bool Is_cohesion; // Initial cohesion (Equiv yield_stress)
  bool Is_friction_angle; // Friction angle
  bool Is_dilatancy_angle; // Dilatancy angle

} Check_Material;

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
static Material  Define_Material(FILE *, Param_Index_and_Model);
static bool Activate_Options(char *, char *);
static Check_Material Initialise_Check_Material();
static void check_Solid_Rigid_Material(Material,Check_Material,int);
static void check_Linear_Elastic_Material(Material,Check_Material,int);
static void check_Saint_Venant_Kirchhoff_Material(Material,Check_Material,int);
static void check_Neo_Hookean_Wriggers_Material(Material,Check_Material,int);
static void check_Von_Mises_Material(Material,Check_Material,int);
static void check_Von_Mises_Material(Material,Check_Material,int);
static void check_Drucker_Prager_Material(Material,Check_Material,int);
static void check_Compressible_Newtonian_Fluid_Material(Material,Check_Material,int);
static void standard_error();
static void standard_output(char *);
static FILE * Open_and_Check_simulation_file(char *);

/**********************************************************************/

Material * Read_Materials__InOutFun__(char * SimulationFile, int NumberMaterials)
/*

Define-Material(idx=0,Model=Drucker-Prager-Plane-Strain)
{
	E=1
	nu=0.3
	Reference-Plastic-Strain=
	Yield_stress=
	Hardening_modulus=
	Hardening_exponent=
	Cohesion=
	Friction_angle=
	Dilatancy_angle=
}

*/
{
	/* Simulation file */
	FILE * Sim_dat;

  /* Variables for reading purposes */
  char line[MAXC] = {0};
  char * kwords[MAXW] = {NULL};
  int nkwords;

	/* Index for the materials */
	int idx = 0;
  
  /* Auxiliar parameter */
  Param_Index_and_Model Index_and_Model;

	/* Allocate table with the material */
  Material * List_Materials = (Material *)malloc(NumberMaterials*sizeof(Material));
  if(List_Materials == NULL)
  {
  	sprintf(Error_message,"%s","Memory error for table of material");
		standard_error();
  }

  /* Open and check file */
  Sim_dat = Open_and_Check_simulation_file(SimulationFile);

  while(fgets(line, sizeof line, Sim_dat) != NULL)
  {

  		/* Read the line with the delimiter_1 */
    	nkwords = parse (kwords, line, delimiters_1);
    	if (nkwords < 0)
    	{
        sprintf(Error_message,"%s","Parser failed");
    		standard_error();
    	}

    	/* Read Initial-nodal-values */
    	if ((nkwords > 0) && (strcmp(kwords[0],"Define-Material") == 0 ))
    	{

    		/* Read index and model */
    		Index_and_Model = Read_Index_and_Model(kwords[1],kwords[2]);

    		/* Read material properties  and asign to the material */
    		List_Materials[idx] = Define_Material(Sim_dat,Index_and_Model);

    		idx++;
    	}

	}

  /* Close  file */
  fclose(Sim_dat);

	return List_Materials;
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

  /* Fill index of the material */
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

static Material Define_Material(FILE * Simulation_file, 
                                Param_Index_and_Model Index_and_Model)
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

  /* Auxiliar parameters for granular materials */
  double rad_friction_angle;
  double rad_dilatancy_angle;

  /* Default values */
  New_Material.thickness = 1;

  /* Default options */
  New_Material.Hardening_Ortiz = false;
  New_Material.Eigenerosion = false;
  New_Material.Eigensoftening = false;

  while(fgets(Parameter_line, sizeof(Parameter_line), Simulation_file) != NULL)
  {
    /* Parse line */    
    Parser_status = parse(Parameter_pars,Parameter_line,delimiters_2);

    if((strcmp(Parameter_pars[0],"{") == 0) && (Parser_status == 1))
    {
      Is_Open = true;
    }
    else if(strcmp(Parameter_pars[0],"rho") == 0)
    {
      ChkMat.Is_rho = true;
      New_Material.rho = atof(Parameter_pars[1]);
    }
    else if(strcmp(Parameter_pars[0],"E") == 0)
    {
      ChkMat.Is_E = true;
      New_Material.E = atof(Parameter_pars[1]);
    }
    else if(strcmp(Parameter_pars[0],"nu") == 0)
    {
      ChkMat.Is_nu = true;
      New_Material.nu = atof(Parameter_pars[1]);
    }
    else if(strcmp(Parameter_pars[0],"Compressibility") == 0)
    {
      ChkMat.Is_Compressibility = true;
      New_Material.Compressibility = atof(Parameter_pars[1]);
    }
    else if(strcmp(Parameter_pars[0],"Reference-Plastic-Strain") == 0)
    {
      ChkMat.Is_E_p0 = true;
      New_Material.E_plastic_reference = atof(Parameter_pars[1]);
    } 
    else if(strcmp(Parameter_pars[0],"Yield-stress") == 0)
    {
      ChkMat.Is_yield_stress = true;
      New_Material.yield_stress_0 = atof(Parameter_pars[1]);
    }
    else if(strcmp(Parameter_pars[0],"Hardening-modulus") == 0)
    {
      ChkMat.Is_H = true;
      New_Material.hardening_modulus = atof(Parameter_pars[1]);
    }
    else if(strcmp(Parameter_pars[0],"Hardening-exponent") == 0)
    {
      ChkMat.Is_Hexp = true;
      New_Material.hardening_exp = atof(Parameter_pars[1]);
    }
    else if(strcmp(Parameter_pars[0],"Cohesion") == 0)
    {
      ChkMat.Is_cohesion = true;
      New_Material.cohesion_reference = atof(Parameter_pars[1]);
    }
    else if(strcmp(Parameter_pars[0],"Friction-angle") == 0)
    {
      ChkMat.Is_friction_angle = true;
      New_Material.friction_angle = atof(Parameter_pars[1]);
      rad_friction_angle  = (PI__MatrixLib__/180)*New_Material.friction_angle;
    }
    else if(strcmp(Parameter_pars[0],"Dilatancy-angle") == 0)
    {
      ChkMat.Is_dilatancy_angle = true;
      New_Material.dilatancy_angle = atof(Parameter_pars[1]);
      rad_dilatancy_angle = (PI__MatrixLib__/180)*New_Material.dilatancy_angle;
    }
    else if(strcmp(Parameter_pars[0],"Hardening-Ortiz") == 0)
    {
      New_Material.Hardening_Ortiz = Activate_Options("Hardening-Ortiz", Parameter_pars[1]);
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

  strcpy(New_Material.Type,Index_and_Model.Model);

  if(strcmp(New_Material.Type,"Solid-Rigid") == 0)
  {
    check_Solid_Rigid_Material(New_Material,ChkMat,Index_and_Model.Idx);
  }
  else if(strcmp(New_Material.Type,"LE") == 0)
  {
    check_Linear_Elastic_Material(New_Material,ChkMat,Index_and_Model.Idx);
    New_Material.Cel = sqrt(New_Material.E/New_Material.rho);
  }
  else if(strcmp(New_Material.Type,"Saint-Venant-Kirchhoff") == 0)
  { 
    check_Saint_Venant_Kirchhoff_Material(New_Material,ChkMat,Index_and_Model.Idx);
    New_Material.Cel = sqrt(New_Material.E/New_Material.rho);
  }
  else if(strcmp(New_Material.Type,"Neo-Hookean-Wriggers") == 0)
  {
    check_Neo_Hookean_Wriggers_Material(New_Material,ChkMat,Index_and_Model.Idx);
    New_Material.Cel = sqrt(New_Material.E/New_Material.rho);
  }
  else if(strcmp(New_Material.Type,"Von-Mises") == 0)
  { 
    check_Von_Mises_Material(New_Material,ChkMat,Index_and_Model.Idx); 
    New_Material.Cel = sqrt(New_Material.E/New_Material.rho);
    New_Material.E_plastic_reference = New_Material.yield_stress_0/New_Material.hardening_modulus;
    TOL_Radial_Returning = 1E-10;
    Max_Iterations_Radial_Returning = 300;
  }
  else if(strcmp(New_Material.Type,"Drucker-Prager-Plane-Strain") == 0)
  { 
    check_Drucker_Prager_Material(New_Material,ChkMat,Index_and_Model.Idx);
    New_Material.Cel = sqrt(New_Material.E/New_Material.rho);
    if(!New_Material.Hardening_Ortiz)
    {
      New_Material.E_plastic_reference = New_Material.cohesion_reference/New_Material.hardening_modulus; 
      New_Material.hardening_exp = 1;
    }
    New_Material.alpha_F_Drucker_Prager = sqrt(2/3.)*tan(rad_friction_angle)/sqrt(3+4*DSQR(tan(rad_friction_angle)));
    New_Material.alpha_Q_Drucker_Prager = sqrt(2/3.)*tan(rad_dilatancy_angle)/sqrt(3+4*DSQR(tan(rad_dilatancy_angle)));
    New_Material.beta_Drucker_Prager    = sqrt(2/3.)*3/sqrt(3+4*DSQR(tan(rad_friction_angle)));
    TOL_Radial_Returning = 1E-10;
    Max_Iterations_Radial_Returning = 30;
  }
  else if(strcmp(New_Material.Type,"Drucker-Prager-Outer-cone") == 0)
  { 
    check_Drucker_Prager_Material(New_Material,ChkMat,Index_and_Model.Idx);
    New_Material.Cel = sqrt(New_Material.E/New_Material.rho);
    if(!New_Material.Hardening_Ortiz)
    {
      New_Material.E_plastic_reference = New_Material.cohesion_reference/New_Material.hardening_modulus; 
      New_Material.hardening_exp = 1;
    }
    New_Material.alpha_F_Drucker_Prager = sqrt(2/3.)*2*sin(rad_friction_angle)/(3-sin(rad_friction_angle));
    New_Material.alpha_Q_Drucker_Prager = sqrt(2/3.)*2*sin(rad_dilatancy_angle)/(3-sin(rad_dilatancy_angle));
    New_Material.beta_Drucker_Prager    = sqrt(2/3.)*6*cos(rad_friction_angle)/(3-sin(rad_friction_angle));
    TOL_Radial_Returning = 1E-10;
    Max_Iterations_Radial_Returning = 30;
  }
  else if(strcmp(New_Material.Type,"Compressible-Newtonian-Fluid") == 0)
  {
    check_Compressible_Newtonian_Fluid_Material(New_Material,ChkMat,Index_and_Model.Idx);
    New_Material.Cel = sqrt(New_Material.Compressibility/New_Material.rho);
  }
  else
  {
    sprintf(Error_message,"%s","Unrecognized kind of material");
    standard_error();
  }

  /* Return outputs */
  return New_Material;
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

static Check_Material Initialise_Check_Material()
{
  Check_Material ChkMat;

  ChkMat.Is_rho = false;
  ChkMat.Is_E = false;
  ChkMat.Is_nu = false;
  ChkMat.Is_E_p0 = false;
  ChkMat.Is_yield_stress = false;
  ChkMat.Is_H = false;
  ChkMat.Is_Hexp = false;
  ChkMat.Is_cohesion = false;
  ChkMat.Is_friction_angle = false;
  ChkMat.Is_dilatancy_angle = false;
  ChkMat.Is_Compressibility = false;

  return ChkMat;
}

/***************************************************************************/
static void check_Solid_Rigid_Material(Material Mat_particle, Check_Material ChkMat, int Idx)
{
  if(ChkMat.Is_rho)
  {
    printf("\t -> %s \n","Solid rigid material");
    printf("\t \t -> %s : %i \n","Idx",Idx);
    printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
  }
  else
  {
    fprintf(stderr,"%s : %s \n",
      "Error in GramsMaterials()",
      "Some parameter is missed for Solid Rigid material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    exit(EXIT_FAILURE);
  }
}

/**********************************************************************/

static void check_Linear_Elastic_Material(Material Mat_particle, Check_Material ChkMat, int Idx)
{
  if(ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu)
  {
    printf("\t -> %s \n","Linear elastic material");
    printf("\t \t -> %s : %i \n","Idx",Idx);
    printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
    printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
    printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
  }
  else
  {
    fprintf(stderr,"%s : %s \n",
      "Error in GramsMaterials()",
      "Some parameter is missed for Linear Elastic material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    fputs(ChkMat.Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
    fputs(ChkMat.Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
    exit(EXIT_FAILURE);
  } 

}

/**********************************************************************/

static void check_Saint_Venant_Kirchhoff_Material(Material Mat_particle, Check_Material ChkMat, int Idx)
{
  if(ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu)
  {
    printf("\t -> %s \n","Saint-Venant-Kirchhoff material");
    printf("\t \t -> %s : %i \n","Idx",Idx);
    printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
    printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
    printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
  }
  else
  {
    fprintf(stderr,"%s : %s \n",
      "Error in GramsMaterials()",
      "Some parameter is missed for Saint-Venant-Kirchhoff material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    fputs(ChkMat.Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
    fputs(ChkMat.Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
    exit(EXIT_FAILURE);
  }

}

/**********************************************************************/

static void check_Neo_Hookean_Wriggers_Material(Material Mat_particle, Check_Material ChkMat, int Idx)
{
  if(ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu)
  {
    printf("\t -> %s \n","Neo-Hookean material (Wriggers model)");
    printf("\t \t -> %s : %i \n","Idx",Idx);
    printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
    printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
    printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
  }
  else
  {
    fprintf(stderr,"%s : %s \n",
      "Error in GramsMaterials()",
      "Some parameter is missed for Neo-Hookean material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    fputs(ChkMat.Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
    fputs(ChkMat.Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
    exit(EXIT_FAILURE);
  }

}

/**********************************************************************/

static void check_Von_Mises_Material(Material Mat_particle, Check_Material ChkMat, int Idx)
{
  if(ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu &&
   ChkMat.Is_yield_stress && ChkMat.Is_H)
  {
    printf("\t -> %s \n","Von-Mises material");
    printf("\t \t -> %s : %i \n","Idx",Idx);
    printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
    printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
    printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
    printf("\t \t -> %s : %f \n","Yield stress",Mat_particle.yield_stress_0);
    printf("\t \t -> %s : %f \n","Hardening modulus",Mat_particle.hardening_modulus);
  }
  else
  {
    fprintf(stderr,"%s : %s \n",
      "Error in GramsMaterials()",
      "Some parameter is missed for Von-Mises material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    fputs(ChkMat.Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
    fputs(ChkMat.Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
    fputs(ChkMat.Is_yield_stress  ? "Yield stress : true \n" : "Yield stress : false \n", stdout);
    fputs(ChkMat.Is_H  ? "Hardening modulus : true \n" : "Hardening modulus : false \n", stdout);
    exit(EXIT_FAILURE);
  }

}

/**********************************************************************/

static void check_Drucker_Prager_Material(Material Mat_particle, Check_Material ChkMat, int Idx)
{
  if(ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu && ChkMat.Is_cohesion && 
     ChkMat.Is_friction_angle && ChkMat.Is_dilatancy_angle)
  {
    printf("\t -> %s \n","Drucker-Prager material");
    printf("\t \t -> %s : %i \n","Idx",Idx);
    printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
    printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
    printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
    printf("\t \t -> %s : %f \n","Cohesion",Mat_particle.cohesion_reference);
    printf("\t \t -> %s : %f \n","Friction angle",Mat_particle.friction_angle);
    printf("\t \t -> %s : %f \n","Dilatancy angle",Mat_particle.dilatancy_angle);
  }
  else
  {
    fprintf(stderr,"%s : %s \n",
      "Error in GramsMaterials()",
      "Some parameter is missed for Drucker-Prager material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    fputs(ChkMat.Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
    fputs(ChkMat.Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
    fputs(ChkMat.Is_cohesion  ? "Cohesion : true \n" : "Cohesion : false \n", stdout);
    fputs(ChkMat.Is_friction_angle  ? "Friction angle : true \n" : "Friction angle : false \n", stdout);
    fputs(ChkMat.Is_dilatancy_angle  ? "Dilatancy angle : true \n" : "Dilatancy angle : false \n", stdout);
    exit(EXIT_FAILURE);
  }

  if(Mat_particle.Hardening_Ortiz)
  {
    if(ChkMat.Is_Hexp && ChkMat.Is_E_p0)
    {
      printf("\t \t -> %s : %f \n","Reference plastic strain",Mat_particle.E_plastic_reference);
      printf("\t \t -> %s : %f \n","Hardening exponent",Mat_particle.hardening_exp);
    }
    else 
    {
      fprintf(stderr,"%s : %s \n",
      "Error in GramsMaterials()",
      "Some parameter is missed for Hardening-Ortiz in Drucker-Prager material");
      fputs(ChkMat.Is_E_p0 ? "Reference plastic strain : true \n" : "Reference plastic strain : false \n", stdout);
      fputs(ChkMat.Is_Hexp  ? "Hardening exponent : true \n" : "Hardening exponent : false \n", stdout);
      exit(EXIT_FAILURE);      
    }
  }
  else
  {
    if(ChkMat.Is_H)
    {
      printf("\t \t -> %s : %f \n","Hardening modulus",Mat_particle.hardening_modulus);
    }
    else
    {
      fprintf(stderr,"%s : %s \n",
      "Error in GramsMaterials()",
      "Some parameter is missed for Drucker-Prager material");
      fputs(ChkMat.Is_H  ? "Hardening modulus : true \n" : "Hardening modulus : false \n", stdout);
      exit(EXIT_FAILURE);
    }
  }

}


/**********************************************************************/

static void check_Compressible_Newtonian_Fluid_Material(Material Mat_particle, Check_Material ChkMat, int Idx)
{
  if(ChkMat.Is_rho && ChkMat.Is_Compressibility)
  {
    printf("\t -> %s \n","Compressible Newtonian Fluid material");
    printf("\t \t -> %s : %i \n","Idx",Idx);
    printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
    printf("\t \t -> %s : %f \n","Compressibility",Mat_particle.Compressibility);
  }
  else
  {
    fprintf(stderr,"%s : %s \n",
      "Error in GramsMaterials()",
      "Some parameter is missed for Compressible Newtonian Fluid material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    fputs(ChkMat.Is_Compressibility ? "Compressibility : true \n" : "Compressibility : false \n", stdout);
    exit(EXIT_FAILURE);
  }
}

/**********************************************************************/

static void standard_error()
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in Define-Material()",Error_message);
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
