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
  bool Is_Cel; // Material Celerity
  bool Is_Compressibility; // Bulk stiffness
  bool Is_Plastic_solver; // FE or BE solver
  bool Is_yield_stress; // Initial Yield stress
  bool Is_friction_angle; // Friction angle
  bool Is_dilatancy_angle; // Dilatancy angle

  bool Is_Hardening_modulus;
  bool Is_Hardening;
  bool Is_Hardening_Hughes;
  bool Is_Parameter_Hardening_Hughes;
  bool Is_Hardening_Cervera;
  bool Is_Hardening_Ortiz;
  bool Is_Exponent_Hardening_Ortiz;
  bool Is_Reference_Plastic_Strain_Ortiz;

  bool Is_Viscous_regularization;
  bool Is_fluidity_param;

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

Material * Read_Materials__InOutFun__(
  char * SimulationFile,
  int NumberMaterials)
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

static Param_Index_and_Model Read_Index_and_Model(
  char * String_Index,
  char * String_Model)
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

static Material Define_Material(
  FILE * Simulation_file,
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

  /* Default options */
  New_Material.Hardening_Ortiz = false;
  New_Material.Eigenerosion = false;
  New_Material.Eigensoftening = false;
  New_Material.Hardening_Hughes = false;
  New_Material.Hardening_Cervera = false;
  New_Material.Exponent_Hardening_Ortiz = false;
  New_Material.Viscous_regularization = false;

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
    else if(strcmp(Parameter_pars[0],"E") == 0)
    {
      ChkMat.Is_E = true;
      New_Material.E = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"nu") == 0)
    {
      ChkMat.Is_nu = true;
      New_Material.nu = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Cel") == 0)
    {
      ChkMat.Is_Cel = true;
      New_Material.Cel = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Compressibility") == 0)
    {
      ChkMat.Is_Compressibility = true;
      New_Material.Compressibility = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Plastic-Solver") == 0)
    {
      ChkMat.Is_Plastic_solver = true;
      strcpy(New_Material.Plastic_Solver,Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Yield-stress") == 0)
    {
      ChkMat.Is_yield_stress = true;
      New_Material.yield_stress_0 = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Hardening-Criteria") == 0)
    {
      ChkMat.Is_Hardening = true;

      if (strcmp(Parameter_pars[1],"Hughes") == 0)
      {
        ChkMat.Is_Hardening_Hughes = true;
        New_Material.Hardening_Hughes = true;
      }
      else if (strcmp(Parameter_pars[1],"Cervera") == 0)
      {
        ChkMat.Is_Hardening_Cervera = true;
        New_Material.Hardening_Cervera = true;
      }
      else if (strcmp(Parameter_pars[1],"Ortiz") == 0)
      {
        ChkMat.Is_Hardening_Ortiz = true;
        New_Material.Hardening_Ortiz = true;
      }
      else
      {
        sprintf(Error_message,"%s","Options for Hardening-Criteria -> Hughes/Cervera/Ortiz");
      standard_error(Error_message);
      }
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Hardening-Modulus") == 0)
    {
      ChkMat.Is_Hardening_modulus = true;
      New_Material.Hardening_modulus = atof(Parameter_pars[1]);
    }
    /**************************************************/

    else if(strcmp(Parameter_pars[0],"Parameter-Hardening-Hughes") == 0)
    {
      ChkMat.Is_Parameter_Hardening_Hughes = true;
      New_Material.Parameter_Hardening_Hughes = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Exponent-Hardening-Ortiz") == 0)
    {
      ChkMat.Is_Exponent_Hardening_Ortiz = true;
      New_Material.Exponent_Hardening_Ortiz = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Reference-Plastic-Strain_Ortiz") == 0)
    {
      ChkMat.Is_Reference_Plastic_Strain_Ortiz = true;
      New_Material.Reference_Plastic_Strain_Ortiz = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Viscous-regularization") == 0)
    {
      ChkMat.Is_Viscous_regularization = true;

      if (strcmp(Parameter_pars[1],"true") == 0)
      {
        New_Material.Viscous_regularization = true;
      }
      else if(strcmp(Parameter_pars[1],"false") == 0)
      {
        New_Material.Viscous_regularization = false;
      }
      else
      {
        sprintf(Error_message,"%s","Options for Viscous-regularization -> true/false");
        standard_error(Error_message);
      }     
    }    
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Fluidity-Parameter") == 0)
    {
      ChkMat.Is_fluidity_param = true;
      New_Material.fluidity_param = atof(Parameter_pars[1]);
    }    
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Friction-angle") == 0)
    {
      ChkMat.Is_friction_angle = true;
      New_Material.friction_angle = atof(Parameter_pars[1]);
      rad_friction_angle  = (PI__MatrixLib__/180)*New_Material.friction_angle;
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Dilatancy-angle") == 0)
    {
      ChkMat.Is_dilatancy_angle = true;
      New_Material.dilatancy_angle = atof(Parameter_pars[1]);
      rad_dilatancy_angle = (PI__MatrixLib__/180)*New_Material.dilatancy_angle;
    }
    /**************************************************/
    else if(strcmp(Parameter_pars[0],"Hardening-Ortiz") == 0)
    {
      New_Material.Hardening_Ortiz = Activate_Options("Hardening-Ortiz", Parameter_pars[1]);
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
    TOL_Radial_Returning = 1E-10;
    Max_Iterations_Radial_Returning = 300;
  }
  else if(strcmp(New_Material.Type,"Drucker-Prager-Plane-Strain") == 0)
  { 
    check_Drucker_Prager_Material(New_Material,ChkMat,Index_and_Model.Idx);
    New_Material.Cel = sqrt(New_Material.E/New_Material.rho);
    TOL_Radial_Returning = 1E-10;
    Max_Iterations_Radial_Returning = 30;
  }
  else if(strcmp(New_Material.Type,"Drucker-Prager-Outer-cone") == 0)
  { 
    check_Drucker_Prager_Material(New_Material,ChkMat,Index_and_Model.Idx);
    New_Material.Cel = sqrt(New_Material.E/New_Material.rho);
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

  ChkMat.Is_rho = false; // Reference fensity
  ChkMat.Is_E = false; // Young modulus
  ChkMat.Is_nu = false; // Poisson cefficient
  ChkMat.Is_Cel = false; // Celerity
  ChkMat.Is_Compressibility = false; // Bulk stiffness
  ChkMat.Is_Plastic_solver = false; // FE or BE solver
  ChkMat.Is_yield_stress = false; // Initial Yield stress
  ChkMat.Is_friction_angle = false; // Friction angle
  ChkMat.Is_dilatancy_angle = false; // Dilatancy angle

  ChkMat.Is_Hardening_modulus = false;
  ChkMat.Is_Hardening = false;
  ChkMat.Is_Hardening_Hughes = false;
  ChkMat.Is_Parameter_Hardening_Hughes = false;
  ChkMat.Is_Hardening_Cervera = false;
  ChkMat.Is_Hardening_Ortiz = false;
  ChkMat.Is_Exponent_Hardening_Ortiz = false;
  ChkMat.Is_Reference_Plastic_Strain_Ortiz = false;

  ChkMat.Is_Viscous_regularization = false;
  ChkMat.Is_fluidity_param = false;

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
  if(ChkMat.Is_rho && ChkMat.Is_Cel && ChkMat.Is_E && 
   ChkMat.Is_nu && ChkMat.Is_yield_stress && ChkMat.Is_Plastic_solver)
  {
    printf("\t -> %s \n","Von-Mises material");
    printf("\t \t -> %s : %f \n","Celerity",Mat_particle.Cel);
    printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
    printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
    printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
    printf("\t \t -> %s : %f \n","Yield stress",Mat_particle.yield_stress_0);
    printf("\t \t -> %s : %s \n","Plastic solver",Mat_particle.Plastic_Solver);
    
    if(ChkMat.Is_Hardening)
    {
      if(ChkMat.Is_Hardening_Hughes)
      {
        if(ChkMat.Is_Hardening_modulus && ChkMat.Is_Parameter_Hardening_Hughes)
        {
          printf("\t \t -> %s : %f \n","Hardening-Modulus",Mat_particle.Hardening_modulus);
          printf("\t \t -> %s : %f \n","Parameter-Hardening-Hughes",Mat_particle.Parameter_Hardening_Hughes); 
        }
        else
        {
          fprintf(stderr,"%s : %s \n",
          "Error in GramsMaterials()",
          "Some parameter is missed for Von-Mises material (Hughes Hardening)");
          fputs(ChkMat.Is_Hardening_modulus  ? "Hardening-Modulus : true \n" : "Hardening-Modulus : false \n", stdout);
          fputs(ChkMat.Is_Parameter_Hardening_Hughes  ? "Parameter-Hardening-Hughes : true \n" : "Parameter-Hardening-Hughes : false \n", stdout);
          exit(EXIT_FAILURE);
        }
  
      }
      else if(ChkMat.Is_Hardening_Cervera)
      {

        if(strcmp(Mat_particle.Plastic_Solver,"Forward-Euler") == 0)
        {
          fprintf(stderr,"%s : %s \n",
            "Error in GramsMaterials()",
            "Switch to Backward-Euler for non-linear Hardening laws)");
          exit(EXIT_FAILURE); 
        }

        if(ChkMat.Is_Hardening_modulus)
        {
          printf("\t \t -> %s : %f \n","Hardening-Modulus",Mat_particle.Hardening_modulus);
        }
        else
        {
          fprintf(stderr,"%s : %s \n",
          "Error in GramsMaterials()",
          "Some parameter is missed for Von-Mises material (Cervera Hardening)");
          fputs(ChkMat.Is_Hardening_modulus  ? "Hardening-Modulus : true \n" : "Hardening-Modulus : false \n", stdout);
          exit(EXIT_FAILURE);
        }
      }
      else if(ChkMat.Is_Hardening_Ortiz)
      {

        if(strcmp(Mat_particle.Plastic_Solver,"Forward-Euler") == 0)
        {
          fprintf(stderr,"%s : %s \n",
            "Error in GramsMaterials()",
            "Switch to Backward-Euler for non-linear Hardening laws)");
          exit(EXIT_FAILURE); 
        }

        if(ChkMat.Is_Exponent_Hardening_Ortiz && ChkMat.Is_Reference_Plastic_Strain_Ortiz)
        {
          printf("\t \t -> %s : %f \n","Exponent-Hardening-Ortiz",Mat_particle.Exponent_Hardening_Ortiz);
          printf("\t \t -> %s : %f \n","Reference-Plastic-Strain_Ortiz",Mat_particle.Reference_Plastic_Strain_Ortiz); 
        }
        else
        {
          fprintf(stderr,"%s : %s \n",
          "Error in GramsMaterials()",
          "Some parameter is missed for Von-Mises material (Ortiz Hardening)");
          fputs(ChkMat.Is_Exponent_Hardening_Ortiz  ? "Exponent-Hardening-Ortiz : true \n" : "Exponent-Hardening-Ortiz : false \n", stdout);
          fputs(ChkMat.Is_Reference_Plastic_Strain_Ortiz  ? "Reference-Plastic-Strain_Ortiz : true \n" : "Reference-Plastic-Strain_Ortiz : false \n", stdout);
          exit(EXIT_FAILURE);
        }

      }
    }

    if(ChkMat.Is_Viscous_regularization)
    {
      if(ChkMat.Is_fluidity_param)
      {
        printf("\t \t -> %s : %f \n","Fluidity parameter",Mat_particle.fluidity_param);
      }
      else
      {
        fprintf(stderr,"%s : %s \n","Error in GramsMaterials()",
          "Some parameter is missed for Von-Mises material (Viscoplastic regularization)");
        fputs(ChkMat.Is_fluidity_param ? "Fluidity parameter : true \n" : "Fluidity parameter : false \n", stdout);
        exit(EXIT_FAILURE);
      }

    }


  }
  else
  {
    fprintf(stderr,"%s : %s \n",
      "Error in GramsMaterials()",
      "Some parameter is missed for Von-Mises material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    fputs(ChkMat.Is_Cel ? "Celerity : true \n" : "Celerity : false \n", stdout);
    fputs(ChkMat.Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
    fputs(ChkMat.Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
    fputs(ChkMat.Is_yield_stress  ? "Yield stress : true \n" : "Yield stress : false \n", stdout);
    exit(EXIT_FAILURE);
  }
}

/**********************************************************************/

static void check_Drucker_Prager_Material(Material Mat_particle, Check_Material ChkMat, int Idx)
{
  if(ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu && 
     ChkMat.Is_friction_angle && ChkMat.Is_dilatancy_angle)
  {
    printf("\t -> %s \n","Drucker-Prager material");
    printf("\t \t -> %s : %i \n","Idx",Idx);
    printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
    printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
    printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
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
    fputs(ChkMat.Is_friction_angle  ? "Friction angle : true \n" : "Friction angle : false \n", stdout);
    fputs(ChkMat.Is_dilatancy_angle  ? "Dilatancy angle : true \n" : "Dilatancy angle : false \n", stdout);
    exit(EXIT_FAILURE);
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
