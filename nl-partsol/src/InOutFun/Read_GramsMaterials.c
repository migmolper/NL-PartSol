#include "nl-partsol.h"

/*
  Call global variables
*/
double Thickness_Plain_Stress;
double TOL_Radial_Returning;
int Max_Iterations_Radial_Returning;

/*
	Global variables for Read_GramsMateriasl.c
*/
bool Is_Id = false;
bool Is_Cel = false;
bool Is_rho = false;
bool Is_E = false;
bool Is_nu = false;
bool Is_Ceps = false;
bool Is_Gf = false;
bool Is_ft = false;
bool Is_E_p0 = false;
bool Is_heps = false;
bool Is_Wc = false;
bool Is_yield_stress = false;

bool Is_H = false;
bool Is_Hexp = false;

bool Is_isotropic_hardening_modulus = false;
bool Is_isotropic_hardening_theta = false;
bool Is_kinematic_hardening_modulus = false;
bool Is_kinematic_hardening_beta = false;

bool Is_fluidity_param = false;

bool Is_cohesion = false;
bool Is_friction_angle = false;
bool Is_dilatancy_angle = false;
bool Is_Compressibility = false;
bool Is_ReferencePressure = false;
bool Is_Viscosity = false;
bool Is_n_Macdonald_model = false;
bool Is_Locking_Control_Fbar = false;

/*
  Auxiliar functions 
*/
static void check_Solid_Rigid_Material(Material);
static void check_Linear_Elastic_Material(Material);
static void check_Saint_Venant_Kirchhoff_Material(Material);
static void check_Neo_Hookean_Wriggers_Material(Material);
static void check_Newtonian_Fluid_Compressible_Material(Material);
static void check_Von_Mises_Material(Material);
static void check_Von_Mises_Perzyna_Material(Material);
static void check_Drucker_Prager_Material(Material);
static void check_Eigenerosion(Material);
static void check_Eigensoftening(Material);
static void standard_error(char *);

/**********************************************************************/

Material * GramsMaterials(char * Name_File, Particle GP_Mesh, int GPxElement)
/*
GramsMaterials (Particles=route.txt) {
               Id=0
               Type=LE
	       Cel=1
	       rho=20
	       E=6.e9
	       nu=0.2
	       Fracture=TRUE
	       Ceps=1.5
	       Gf=0.00001
}  
*/
{
  /* Asign material library to an auxiliar variable */
  Material Mat_GP;
  Material * Mat_Table;
  int Ndim = NumberDimensions;
  /* Simulation file */
  FILE * Sim_dat;

  /* Error message */
  char Error_message[MAXC] = {0};

  /* Index of the material */
  int Aux_Mat_id;
  char * Parse_Mat_id[MAXW] = {NULL};

  /* Material properties */
  char Line_Material_Prop[MAXC] = {0};
  char * Parse_Mat_Prop[MAXW] = {NULL};

  /* Parse file name with the list of nodes */
  char Route_Nodes[MAXC] = {0};
  char FileNodesRoute[MAXC];

  /* Array */
  int Num_Nodes;
  ChainPtr Chain_Nodes = NULL;
  int * Array_Nodes;

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  /* Auxiliar variable for status */
  char * STATUS_LINE;
  int CountMaterials = 0;

  /* Auxiliar parameters for granular materials */
  double rad_friction_angle;
  double rad_dilatancy_angle;
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL)
  {
  	sprintf(Error_message,"%s %s","Incorrect lecture of",Name_File);
	standard_error(Error_message);
  }

  /* Generate route */
  generate_route(Route_Nodes,Name_File);
   
  /* Allocate table with the material */
  Mat_Table = (Material *)malloc(GP_Mesh.NumberMaterials*sizeof(Material));
  if(Mat_Table == NULL)
  {
  	sprintf(Error_message,"%s","Memory error for table of material");
	standard_error(Error_message);
  }
  
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");
    if (nkwords < 0)
    {
    	sprintf(Error_message,"%s \n","Parser failed");
		standard_error(Error_message);
    }

    if ((nkwords>0) &&
	(strcmp(kwords[0],"GramsMaterials") == 0 )){

      /* Count the number of materials */
      ++CountMaterials;

      /* Read the index of the material */
      Aux_Mat_id = parse (Parse_Mat_id, kwords[1],"(=)");
      if( (Aux_Mat_id != 2) || (strcmp(Parse_Mat_id[0],"Particles") != 0))
      {
      	sprintf(Error_message,"%s","Use this format -> (Particles=route.txt)");
		standard_error(Error_message);
      }

      /* Read file with the nodes */
      sprintf(FileNodesRoute,"%s%s",Route_Nodes,Parse_Mat_id[1]);
      printf("\t -> %s : %s \n","Material points",FileNodesRoute);

      /* Get an array with the nodes */
      Chain_Nodes = File2Chain(FileNodesRoute);
      Num_Nodes = lenght__SetLib__(Chain_Nodes);
      Array_Nodes = set_to_memory__SetLib__(Chain_Nodes,Num_Nodes);
      free__SetLib__(&Chain_Nodes);

      /* Id of the material */
      Mat_GP.Id = -1;
      /* Celerity */
      Mat_GP.Cel = NAN;
      /* Density */      
      Mat_GP.rho = NAN;
      /* Linear elastic parameters */
      Mat_GP.E = NAN;
      Mat_GP.nu = NAN;
      /* Fracture module */
      Mat_GP.Eigenerosion = false;
      Mat_GP.Eigensoftening = false;
      /* Parameters for Eigenerosion */
      Mat_GP.Ceps = NAN;
      Mat_GP.Gf = NAN;
      /* Parameters for Eigensoftening */
      Mat_GP.ft = NAN;
      Mat_GP.heps = NAN;
      Mat_GP.Wc = NAN;
      /* Parameters for plastic simulations */
      Mat_GP.yield_stress_0 = NAN;
      Mat_GP.cohesion_reference = NAN;
      Mat_GP.friction_angle = NAN;
      Mat_GP.dilatancy_angle = NAN;
      Mat_GP.hardening_modulus = NAN;
	  Mat_GP.hardening_exp = NAN;

	  /* Parameters for isotropic/kinematic hardening */
	  Mat_GP.isotropic_hardening_modulus = NAN;
	  Mat_GP.isotropic_hardening_theta = 1.0;
	  Mat_GP.kinematic_hardening_modulus = NAN;
	  Mat_GP.kinematic_hardening_beta = 0.0;

	  /* Fluidity parameters for the viscoplasticity */
	  Mat_GP.fluidity_param = NAN;

      /* Look for the curly brace { */
      if(strcmp(kwords[2],"{") == 0){

	/* Initial line */
	STATUS_LINE = fgets(Line_Material_Prop, sizeof(Line_Material_Prop), Sim_dat);
	if(STATUS_LINE == NULL)
	{
		sprintf(Error_message,"%s","Unspected EOF");
		standard_error(Error_message);
	}
	Aux_Mat_id = parse(Parse_Mat_Prop,Line_Material_Prop," =\t\n");
	if(strcmp(Parse_Mat_Prop[0],"}") == 0)
	{
		sprintf(Error_message,"%s","The material was not defined");
		standard_error(Error_message);
	}
	while(STATUS_LINE != NULL){

	  if(Aux_Mat_id != 2)
	  {
	  	sprintf(Error_message,"%s","The input format -> Propertie=value");
		standard_error(Error_message);
	  }
	  /**************************************************/
	  if(strcmp(Parse_Mat_Prop[0],"Id") == 0)
	  {
	    Is_Id = true;
	    Mat_GP.Id = atoi(Parse_Mat_Prop[1]);
	    if(Mat_GP.Id >= GP_Mesh.NumberMaterials)
	    {
	    	sprintf(Error_message,"%s %i !!! \n","Id should go from 0 to",GP_Mesh.NumberMaterials-1);
			standard_error(Error_message);
	    }
	    for(int i = 0 ; i<Num_Nodes ; i++)
	    {
	      for(int j = 0 ; j<GPxElement ; j++)
	      {
	    	GP_Mesh.MatIdx[Array_Nodes[i]*GPxElement+j] = Mat_GP.Id;
	      }
	    }
	  }
	  /**************************************************/
 	  else if(strcmp(Parse_Mat_Prop[0],"Type") == 0)
 	  {
	    strcpy(Mat_GP.Type,Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Cel") == 0)
	  {
	    Is_Cel = true;
	    Mat_GP.Cel = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"rho") == 0)
	  {
	    Is_rho = true;
	    Mat_GP.rho = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"E") == 0)
	  {
	    Is_E = true;
	    Mat_GP.E = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"nu") == 0)
	  {
	    Is_nu = true;
	    Mat_GP.nu = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/	    
	  else if(strcmp(Parse_Mat_Prop[0],"Fracture") == 0)
	  {
	    if (strcmp(Parse_Mat_Prop[1],"Eigenerosion") == 0)
	    {
	      Mat_GP.Eigenerosion = true;
	    }
	    else if (strcmp(Parse_Mat_Prop[1],"Eigensoftening") == 0)
	    {
	      Mat_GP.Eigensoftening=true;
	    }
	    else
	    {
	    	sprintf(Error_message,"%s","Options -> Eigenerosion/Eigensoftening");
			standard_error(Error_message);
	    }
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Fbar") == 0)
	  {
	   	if(strcmp(Parse_Mat_Prop[1],"true") == 0)
	   	{
	   		Is_Locking_Control_Fbar = true;
	   		Mat_GP.Locking_Control_Fbar = true;
	   	}
	   	else if(strcmp(Parse_Mat_Prop[1],"false") == 0)
	   	{
		    Is_Locking_Control_Fbar = false;
	   		Mat_GP.Locking_Control_Fbar = false;
	   	}
	   	else
	   	{
	   		sprintf(Error_message,"The input was %s. Please, use : true/false",Parse_Mat_Prop[1]);
	   		standard_error(Error_message); 
	   	}
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Ceps") == 0)
	  {
	    Is_Ceps = true;
	    Mat_GP.Ceps = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Gf") == 0)
	  {
	    Is_Gf = true;
	    Mat_GP.Gf = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"ft") == 0)
	  {
	    Is_ft = true;
	    Mat_GP.ft = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"heps") == 0)
	  {
	    Is_heps = true;
	    Mat_GP.heps = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Reference-Plastic-Strain") == 0)
	  {
	    Is_E_p0 = true;
	    Mat_GP.E_plastic_reference = atof(Parse_Mat_Prop[1]);
	  } 
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Wc") == 0)
	  {
	    Is_Wc = true;
	    Mat_GP.Wc = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Yield_stress") == 0)
	  {
	    Is_yield_stress = true;
	    Mat_GP.yield_stress_0 = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Hardening_modulus") == 0)
	  {
	    Is_H = true;
	    Mat_GP.hardening_modulus = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/

	  else if(strcmp(Parse_Mat_Prop[0],"Isotropic-Hardening-Modulus") == 0)
	  {
	    Is_isotropic_hardening_modulus = true;
	    Mat_GP.isotropic_hardening_modulus = atof(Parse_Mat_Prop[1]);
	  }

	  /**************************************************/

	  else if(strcmp(Parse_Mat_Prop[0],"Isotropic-Hardening-Theta") == 0)
	  {
	    Is_isotropic_hardening_theta = true;
	    Mat_GP.isotropic_hardening_theta = atof(Parse_Mat_Prop[1]);
	  }

	  /**************************************************/

	  else if(strcmp(Parse_Mat_Prop[0],"Kinematic-Hardening-Modulus") == 0)
	  {
	    Is_kinematic_hardening_modulus = true;
	    Mat_GP.kinematic_hardening_modulus = atof(Parse_Mat_Prop[1]);
	  }

	  /**************************************************/

	  else if(strcmp(Parse_Mat_Prop[0],"Kinematic-Hardening-Beta") == 0)
	  {
	    Is_kinematic_hardening_beta = true;
	    Mat_GP.kinematic_hardening_beta = atof(Parse_Mat_Prop[1]);
	  }

	 /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Fluidity-Parameter") == 0)
	  {
	  	Is_fluidity_param = true;
	  	Mat_GP.fluidity_param = atof(Parse_Mat_Prop[1]);
	  }

	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Hardening_exponent") == 0)
	  {
	    Is_Hexp = true;
	    Mat_GP.hardening_exp = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Cohesion") == 0)
	  {
	    Is_cohesion = true;
	    Mat_GP.cohesion_reference = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Friction_angle") == 0)
	  {
	    Is_friction_angle = true;
	    Mat_GP.friction_angle = atof(Parse_Mat_Prop[1]);
	    rad_friction_angle  = (PI__MatrixLib__/180)*Mat_GP.friction_angle;
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Dilatancy_angle") == 0)
	  {
	    Is_dilatancy_angle = true;
	    Mat_GP.dilatancy_angle = atof(Parse_Mat_Prop[1]);
	    rad_dilatancy_angle = (PI__MatrixLib__/180)*Mat_GP.dilatancy_angle;
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Compressibility") == 0)
	  {
	  	Is_Compressibility = true;
	  	Mat_GP.Compressibility = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Reference-Pressure") == 0)
	  {
	  	Is_ReferencePressure = true;
	  	Mat_GP.ReferencePressure = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Viscosity") == 0)
	  {
	  	Is_Viscosity = true;
	  	Mat_GP.Viscosity = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Macdonald-parameter") == 0)
	  {
	  	Is_n_Macdonald_model = true;
	  	Mat_GP.n_Macdonald_model = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else
	  {
	  	sprintf(Error_message,"%s %s %s","the propertie",Parse_Mat_Prop[0],"is not defined");
		standard_error(Error_message);
	  }
	  /* Read next line and check */
	  STATUS_LINE = fgets(Line_Material_Prop, sizeof(Line_Material_Prop), Sim_dat);
	  Aux_Mat_id = parse(Parse_Mat_Prop,Line_Material_Prop," =\t\n");
	  if(strcmp(Parse_Mat_Prop[0],"}") == 0)
	  {
	    break;
	  }
	}
	if(STATUS_LINE == NULL)
	{
		sprintf(Error_message,"%s","you forget to put a }");
		standard_error(Error_message);
	}

	if(strcmp(Parse_Mat_Prop[0],"}") == 0)
	{
	  /**************************************************/
	  if(Is_Id)
	  {
	    printf("\t \t -> %s : %i \n","Index of the material",Mat_GP.Id);
	  }
	  /**************************************************/
	  /* Solid rigid material */
	  if(strcmp(Mat_GP.Type,"Solid-Rigid") == 0)
	  {
	  	check_Solid_Rigid_Material(Mat_GP);
		Mat_GP.Cel = 0;
	  }
	  /* Parameters for a linear elastic material */
	  else if(strcmp(Mat_GP.Type,"LE") == 0)
	  {
	  	check_Linear_Elastic_Material(Mat_GP);
	  }
	  /* Paramters for a Saint Venant Kirchhoff material */
	  else if(strcmp(Mat_GP.Type,"Saint-Venant-Kirchhoff") == 0)
	  { 
		check_Saint_Venant_Kirchhoff_Material(Mat_GP);
	  }
	  /* Parameters for a Neohookean (Wriggers) material */
	  else if(strcmp(Mat_GP.Type,"Neo-Hookean-Wriggers") == 0)
	  {
		check_Neo_Hookean_Wriggers_Material(Mat_GP);
	  }
	  /* Parameters for a Newtonian Compressible fluid */
	  else if(strcmp(Mat_GP.Type,"Newtonian-Fluid-Compressible") == 0)
	  {
	  	check_Newtonian_Fluid_Compressible_Material(Mat_GP);
	  }
	  /* Parameters for a Von Mises Yield criterium */
	  else if(strcmp(Mat_GP.Type,"Von-Mises") == 0)
	  { 
		check_Von_Mises_Material(Mat_GP);	
	  	TOL_Radial_Returning = 1E-10;
		Max_Iterations_Radial_Returning = 300;
	  }
	  /* Parameters for a Von Mises Yield criterium + Perzyna visplasticity */
	  else if(strcmp(Mat_GP.Type,"Von-Mises-Perzyna") == 0)
	  { 
		check_Von_Mises_Perzyna_Material(Mat_GP);	
	  	TOL_Radial_Returning = 1E-10;
		Max_Iterations_Radial_Returning = 300;
	  }
	  /* Parameters for a Drucker-Prager Yield criterium */
	  else if(strcmp(Mat_GP.Type,"Drucker-Prager-Plane-Strain") == 0)
	  { 
	  	check_Drucker_Prager_Material(Mat_GP);

		/*	Plane strain yield surface */
		Mat_GP.alpha_F_Drucker_Prager = sqrt(2/3.)*tan(rad_friction_angle)/sqrt(3+4*DSQR(tan(rad_friction_angle)));
		Mat_GP.alpha_Q_Drucker_Prager = sqrt(2/3.)*tan(rad_dilatancy_angle)/sqrt(3+4*DSQR(tan(rad_dilatancy_angle)));
		Mat_GP.beta_Drucker_Prager    = sqrt(2/3.)*3/sqrt(3+4*DSQR(tan(rad_friction_angle)));
		TOL_Radial_Returning = 1E-10;
		Max_Iterations_Radial_Returning = 30;
	  }
	  else if(strcmp(Mat_GP.Type,"Drucker-Prager-Outer-cone") == 0)
	  { 
	  	check_Drucker_Prager_Material(Mat_GP);

		/*	Outer cone yield surface */
		Mat_GP.alpha_F_Drucker_Prager = sqrt(2/3.)*2*sin(rad_friction_angle)/(3-sin(rad_friction_angle));
		Mat_GP.alpha_Q_Drucker_Prager = sqrt(2/3.)*2*sin(rad_dilatancy_angle)/(3-sin(rad_dilatancy_angle));
		Mat_GP.beta_Drucker_Prager    = sqrt(2/3.)*6*cos(rad_friction_angle)/(3-sin(rad_friction_angle));
		TOL_Radial_Returning = 1E-10;
		Max_Iterations_Radial_Returning = 30;
	  }
	  else
	  {
	  	sprintf(Error_message,"%s","Unrecognized kind of material");
		standard_error(Error_message);
	  }
	  /**************************************************/
	  if(Mat_GP.Eigenerosion)
	  { 
		check_Eigenerosion(Mat_GP);
	  }
	  /**************************************************/ 
	  if(Mat_GP.Eigensoftening)
	  { 
		check_Eigensoftening(Mat_GP);
	  }
	  /**************************************************/
	  
	  /* Transfere information */
	  Mat_Table[Mat_GP.Id] = Mat_GP;
	  
	  /* break; */
	}
      }
      else
      {
      	sprintf(Error_message,"%s","Use this format -> (Particles=route.txt) { ");
		standard_error(Error_message);
      }
    }

    /* /\* Free array nodes *\/ */
    /* free(Array_Nodes); */
  }
    
  /* Check the number of materials */
  if(CountMaterials != GP_Mesh.NumberMaterials)
  {
  	sprintf(Error_message,"%s %i %s %i %s \n","Spected",GP_Mesh.NumberMaterials, "materials, but",CountMaterials,"where defined");
	standard_error(Error_message);
  }

  /* Close .dat file */
  /* Final message */
  fclose(Sim_dat);

  printf("\t * End of read data file !!! \n");

  return Mat_Table;
}

/**********************************************************************/

static void check_Solid_Rigid_Material(Material Mat_particle)
{
	if(Is_rho)
	{
		printf("\t -> %s \n","Solid rigid material");
		printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
	}
	else
	{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()",
			"Some parameter is missed for Solid Rigid material");
		fputs(Is_rho ? "Density : true \n" : "Density : false \n", stdout);
		exit(EXIT_FAILURE);
	}
}

/**********************************************************************/

static void check_Linear_Elastic_Material(Material Mat_particle)
{
	if(Is_rho && Is_Cel && Is_E && Is_nu)
	{
		printf("\t -> %s \n","Linear elastic material");
		printf("\t \t -> %s : %f \n","Celerity",Mat_particle.Cel);
		printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
		printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
		printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
	}
	else
	{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()",
			"Some parameter is missed for Linear Elastic material");
		fputs(Is_rho ? "Density : true \n" : "Density : false \n", stdout);
		fputs(Is_Cel ? "Celerity : true \n" : "Celerity : false \n", stdout);
		fputs(Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
		fputs(Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
		exit(EXIT_FAILURE);
	}	
}

/**********************************************************************/

static void check_Saint_Venant_Kirchhoff_Material(Material Mat_particle)
{
	if(Is_rho && Is_Cel && Is_E && Is_nu)
	{
		printf("\t -> %s \n","Saint-Venant-Kirchhoff material");
		printf("\t \t -> %s : %f \n","Celerity",Mat_particle.Cel);
		printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
		printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
		printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
	}
	else
	{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()",
			"Some parameter is missed for Saint-Venant-Kirchhoff material");
		fputs(Is_rho ? "Density : true \n" : "Density : false \n", stdout);
		fputs(Is_Cel ? "Celerity : true \n" : "Celerity : false \n", stdout);
		fputs(Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
		fputs(Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
		exit(EXIT_FAILURE);
	}
}

/**********************************************************************/

static void check_Neo_Hookean_Wriggers_Material(Material Mat_particle)
{
	if(Is_rho && Is_Cel && Is_E && Is_nu)
	{
		printf("\t -> %s \n","Neo-Hookean material (Wriggers model)");
		printf("\t \t -> %s : %f \n","Celerity",Mat_particle.Cel);
		printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
		printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
		printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
	}
	else
	{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()",
			"Some parameter is missed for Neo-Hookean material");
		fputs(Is_rho ? "Density : true \n" : "Density : false \n", stdout);
		fputs(Is_Cel ? "Celerity : true \n" : "Celerity : false \n", stdout);
		fputs(Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
		fputs(Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
		exit(EXIT_FAILURE);
	}
}

/**********************************************************************/

static void check_Newtonian_Fluid_Compressible_Material(Material Mat_particle)
{	

	if(Is_rho && Is_Cel && Is_Compressibility && Is_ReferencePressure && Is_Viscosity && Is_n_Macdonald_model)
	{
		printf("\t -> %s \n","Newtonian-Compressible fluid");
		printf("\t \t -> %s : %f \n","Celerity",Mat_particle.Cel);
		printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
		printf("\t \t -> %s : %f \n","Compressibility",Mat_particle.Compressibility);
		printf("\t \t -> %s : %f \n","Reference Pressure",Mat_particle.ReferencePressure);
		printf("\t \t -> %s : %f \n","Viscosity",Mat_particle.Viscosity);
		printf("\t \t -> %s : %f \n","Macdonald-parameter",Mat_particle.n_Macdonald_model);
	}
	else
	{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()",
			"Some parameter is missed for Newtonian-Compressible fluid");
		fputs(Is_rho ? "Density : true \n" : "Density : false \n", stdout);
		fputs(Is_Compressibility ? "Compressibility : true \n" : "Compressibility : false \n", stdout);
		fputs(Is_ReferencePressure ? "Reference Pressure : true \n" : "Reference Pressure : false \n", stdout);
		fputs(Is_Viscosity ? "Viscosity : true \n" : "Viscosity : false \n", stdout);
		fputs(Is_n_Macdonald_model ? "Macdonald-parameter : true \n" : "Macdonald-parameter : false \n", stdout);
		exit(EXIT_FAILURE);
	}

}

/**********************************************************************/

static void check_Von_Mises_Material(Material Mat_particle)
{
	if(Is_rho && Is_Cel && Is_E && 
	 Is_nu && Is_yield_stress &&
	 Is_isotropic_hardening_modulus &&
	 Is_isotropic_hardening_theta &&
	 Is_kinematic_hardening_modulus &&
	 Is_kinematic_hardening_beta)
	{
		printf("\t -> %s \n","Von-Mises material");
		printf("\t \t -> %s : %f \n","Celerity",Mat_particle.Cel);
		printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
		printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
		printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
		printf("\t \t -> %s : %f \n","Yield stress",Mat_particle.yield_stress_0);
		printf("\t \t -> %s : %f \n","Isotropic hardening modulus",Mat_particle.isotropic_hardening_modulus);
		printf("\t \t -> %s : %f \n","Isotropic hardening theta",Mat_particle.isotropic_hardening_theta);		
		printf("\t \t -> %s : %f \n","Kinematic hardening modulus",Mat_particle.kinematic_hardening_modulus);
		printf("\t \t -> %s : %f \n","Kinematic hardening beta",Mat_particle.kinematic_hardening_beta);
	}
	else
	{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()",
			"Some parameter is missed for Von-Mises material");
		fputs(Is_rho ? "Density : true \n" : "Density : false \n", stdout);
		fputs(Is_Cel ? "Celerity : true \n" : "Celerity : false \n", stdout);
		fputs(Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
		fputs(Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
		fputs(Is_yield_stress  ? "Yield stress : true \n" : "Yield stress : false \n", stdout);
		fputs(Is_isotropic_hardening_modulus  ? "Isotropic hardening modulus : true \n" : "Isotropic hardening modulus : false \n", stdout);
		fputs(Is_isotropic_hardening_theta  ? "Isotropic hardening theta : true \n" : "Isotropic hardening theta : false \n", stdout);
		fputs(Is_kinematic_hardening_modulus  ? "Kinematic hardening modulus : true \n" : "Kinematic hardening modulus : false \n", stdout);
		fputs(Is_kinematic_hardening_beta  ? "Kinematic hardening beta : true \n" : "Kinematic hardening beta : false \n", stdout);
		exit(EXIT_FAILURE);
	}
}

/**********************************************************************/

static void check_Von_Mises_Perzyna_Material(Material Mat_particle)
{
	if(Is_rho && Is_Cel && Is_E && 
	 Is_nu && Is_yield_stress &&
	 Is_isotropic_hardening_modulus &&
	 Is_isotropic_hardening_theta &&
	 Is_kinematic_hardening_modulus &&
	 Is_kinematic_hardening_beta &&
	 Is_fluidity_param)
	{
		printf("\t -> %s \n","Von-Mises material");
		printf("\t \t -> %s : %f \n","Celerity",Mat_particle.Cel);
		printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
		printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
		printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
		printf("\t \t -> %s : %f \n","Yield stress",Mat_particle.yield_stress_0);
		printf("\t \t -> %s : %f \n","Isotropic hardening modulus",Mat_particle.isotropic_hardening_modulus);
		printf("\t \t -> %s : %f \n","Isotropic hardening theta",Mat_particle.isotropic_hardening_theta);		
		printf("\t \t -> %s : %f \n","Kinematic hardening modulus",Mat_particle.kinematic_hardening_modulus);
		printf("\t \t -> %s : %f \n","Kinematic hardening beta",Mat_particle.kinematic_hardening_beta);
		printf("\t \t -> %s : %f \n","Fluidity parameter",Mat_particle.fluidity_param);
	}
	else
	{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()",
			"Some parameter is missed for Von-Mises material");
		fputs(Is_rho ? "Density : true \n" : "Density : false \n", stdout);
		fputs(Is_Cel ? "Celerity : true \n" : "Celerity : false \n", stdout);
		fputs(Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
		fputs(Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
		fputs(Is_yield_stress  ? "Yield stress : true \n" : "Yield stress : false \n", stdout);
		fputs(Is_isotropic_hardening_modulus  ? "Isotropic hardening modulus : true \n" : "Isotropic hardening modulus : false \n", stdout);
		fputs(Is_isotropic_hardening_theta  ? "Isotropic hardening theta : true \n" : "Isotropic hardening theta : false \n", stdout);
		fputs(Is_kinematic_hardening_modulus  ? "Kinematic hardening modulus : true \n" : "Kinematic hardening modulus : false \n", stdout);
		fputs(Is_kinematic_hardening_beta  ? "Kinematic hardening beta : true \n" : "Kinematic hardening beta : false \n", stdout);
		fputs(Is_fluidity_param  ? "Fluidity parameter : true \n" : "Fluidity parameter : false \n", stdout);
		exit(EXIT_FAILURE);
	}
}

/**********************************************************************/

static void check_Drucker_Prager_Material(Material Mat_particle)
{
	if(Is_rho && Is_Cel && Is_E && Is_nu && Is_cohesion && Is_Hexp && Is_E_p0 && Is_friction_angle && Is_dilatancy_angle)
	{
		printf("\t -> %s \n","Drucker-Prager material");
		printf("\t \t -> %s : %f \n","Celerity",Mat_particle.Cel);
		printf("\t \t -> %s : %f \n","Density",Mat_particle.rho);
		printf("\t \t -> %s : %f \n","Elastic modulus",Mat_particle.E);
		printf("\t \t -> %s : %f \n","Poisson modulus",Mat_particle.nu);
		printf("\t \t -> %s : %f \n","Cohesion",Mat_particle.cohesion_reference);
		printf("\t \t -> %s : %f \n","Reference plastic strain",Mat_particle.E_plastic_reference);
		printf("\t \t -> %s : %f \n","Hardening exponent",Mat_particle.hardening_exp);
		printf("\t \t -> %s : %f \n","Friction angle",Mat_particle.friction_angle);
		printf("\t \t -> %s : %f \n","Dilatancy angle",Mat_particle.dilatancy_angle);
	}
	else
	{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()",
			"Some parameter is missed for Drucker-Prager material");
		fputs(Is_rho ? "Density : true \n" : "Density : false \n", stdout);
		fputs(Is_Cel ? "Celerity : true \n" : "Celerity : false \n", stdout);
		fputs(Is_E   ? "Elastic modulus : true \n" : "Elastic modulus : false \n", stdout);
		fputs(Is_nu  ? "Poisson modulus : true \n" : "Poisson modulus : false \n", stdout);
		fputs(Is_cohesion  ? "Cohesion : true \n" : "Cohesion : false \n", stdout);
		fputs(Is_E_p0 ? "Reference plastic strain : true \n" : "Reference plastic strain : false \n", stdout);
		fputs(Is_Hexp  ? "Hardening exponent : true \n" : "Hardening exponent : false \n", stdout);
		fputs(Is_friction_angle  ? "Friction angle : true \n" : "Friction angle : false \n", stdout);
		fputs(Is_dilatancy_angle  ? "Dilatancy angle : true \n" : "Dilatancy angle : false \n", stdout);
		exit(EXIT_FAILURE);
	}
}

/**********************************************************************/

static void check_Eigenerosion(Material Mat_particle)
{
	if(Is_Ceps && Is_Gf)
	{
		printf("\t \t -> %s : %s \n","Fracture criterion","Eigenerosion");
		printf("\t \t -> %s : %f \n","Normalizing constant",Mat_particle.Ceps);
		printf("\t \t -> %s : %f \n","Failure energy",Mat_particle.Gf);
	}
	else
	{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()",
			"Some parameters are missed for Eigenerosion");
		fputs(Is_Ceps ? "Normalizing constan : true \n" : "Normalizing constan : false \n", stdout);
		fputs(Is_Gf ? "Failure energy : true \n" : "Failure energy : false \n", stdout);
		exit(EXIT_FAILURE);
	}
}

/**********************************************************************/

static void check_Eigensoftening(Material Mat_particle)
{
	if(Is_Ceps && Is_ft && Is_heps && Is_Wc)
	{
		printf("\t \t -> %s : %s \n","Fracture criterion","Eigensoftening");
		printf("\t \t -> %s : %f \n","Normalizing constant",Mat_particle.Ceps);
		printf("\t \t -> %s : %f \n","Tensile strengt",Mat_particle.ft);
		printf("\t \t -> %s : %f \n","Bandwidth Bazant",Mat_particle.heps);
		printf("\t \t -> %s : %f \n","Critical opening",Mat_particle.Wc);
	}
	else
	{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsMaterials()",
			"Some parameters are missed for Eigensoftening");
		fputs(Is_Ceps ? "Normalizing constan : true \n" : "Normalizing constan : false \n", stdout);
		fputs(Is_ft ? "Tensile strengt : true \n" : "Tensile strengt : false \n", stdout);
		fputs(Is_heps ? "Bandwidth Bazant : true \n" : "Bandwidth Bazant : false \n", stdout);
		fputs(Is_Wc ? "Critical opening : true \n" : "Critical opening : false \n", stdout);
		exit(EXIT_FAILURE);
	}
}

/**********************************************************************/

static void standard_error(char * Error_message)
{
	fprintf(stderr,"%s : %s !!! \n",
	   "Error in GramsMaterials()",Error_message);
    exit(EXIT_FAILURE);
}

/**********************************************************************/
