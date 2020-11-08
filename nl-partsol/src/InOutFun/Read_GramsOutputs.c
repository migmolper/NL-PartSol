#include "nl-partsol.h"

/*
  Call global variables
*/

int ResultsTimeStep;
int NumTimeStep;
char OutputDir[MAXC];

bool Out_global_coordinates = false;
bool Out_element_coordinates = false;
bool Out_mass = false;
bool Out_density = false;
bool Out_damage = false;
bool Out_nodal_idx = false;
bool Out_material_idx = false;
bool Out_velocity = false;
bool Out_acceleration = false;
bool Out_displacement = false;
bool Out_stress = false;
bool Out_eigenvalues_stress = false;
bool Out_volumetric_stress = false;
bool Out_strain = false;
bool Out_eigenvalues_strain = false;
bool Out_deformation_gradient = false;
bool Out_energy = false;

/*
  Auxiliar functions 
*/
static bool Is_Output_Activate(char *);

/**********************************************************************/

void GramsOutputs(char * Name_File)
/*
  Example : 
  GramsOutputs (i=100) {
  DIR=test/Sulsky_MPM	
  }
*/
{
  /* Simulation file */
  FILE * Sim_dat;

  /* Temporal integator */
  int Aux_Out_id;
  char * Parse_Out_id[MAXW] = {NULL};

  /* Temporal integrator properties */
  char Line_Out_Prop[MAXC] = {0};
  char * Parse_Out_Prop[MAXW] = {NULL};

  /* Parse file for the route */
  bool Is_OutputDir = false;
  char Route_Outs[MAXC] = {0};

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  /* Auxiliar variable for status */
  char * STATUS_LINE;

  /* Initial message */  
  puts("*************************************************");
  printf(" \t %s : \n\t %s \n",
	 "* Read Outputs properties ",
	 Name_File);
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	   "Error in GramsOutputs()",
	   "Incorrect lecture of",
	   Name_File);
    exit(EXIT_FAILURE);
  }

  /* Generate route */
  generate_route(Route_Outs,Name_File);
    
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");
    if (nkwords < 0){
      fprintf(stderr,"%s : %s \n",
	     "Error in GramsOutputs()",
	     "Parser failed");
      exit(EXIT_FAILURE);
    }

    if ((nkwords > 0) &&
	(strcmp(kwords[0],"GramsOutputs") == 0 )){

      /* Read temporal integrator scheme */
      Aux_Out_id = parse (Parse_Out_id, kwords[1],"(=)");
      if( (Aux_Out_id != 2) ||
	  (strcmp(Parse_Out_id[0],"i") != 0)){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsOutputs()",
	       "Use this format -> (i=int) !!!");
	exit(EXIT_FAILURE);
      }
      ResultsTimeStep = atoi(Parse_Out_id[1]);
      if(ResultsTimeStep > NumTimeStep){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsOutputs()",
	       "The result time step should be less than the total time steps !!!");
	exit(EXIT_FAILURE);	
      }
      printf("\t -> %s : %i \n","Output values each",ResultsTimeStep);

      /* Look for the curly brace { */
      if(strcmp(kwords[2],"{") == 0){
	/* Initial line */
	STATUS_LINE = fgets(Line_Out_Prop,
			    sizeof(Line_Out_Prop),
			    Sim_dat);
	if(STATUS_LINE == NULL){
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsOutputs()",
		  "Unspected EOF !!!");
	  exit(EXIT_FAILURE);	
	}
	Aux_Out_id = parse(Parse_Out_Prop,Line_Out_Prop," =\t\n");
	if(strcmp(Parse_Out_Prop[0],"}") == 0){
	  /* Check output dir */
	  if(!Is_OutputDir){
	    fprintf(stderr,"%s : %s \n",
		    "Error in GramsOutputs()",
		    "Non output dir defined !!!");
	    exit(EXIT_FAILURE);
	  }
	  break;
	}
	while(STATUS_LINE != NULL){
	  
	  if(Aux_Out_id != 2){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsOutputs()",
		   "Use this format -> Propertie = value !!!");
	    exit(EXIT_FAILURE);
	  }

 	  if(strcmp(Parse_Out_Prop[0],"DIR") == 0)
 	  {
	    sprintf(OutputDir,"%s%s",Route_Outs,Parse_Out_Prop[1]);
	    printf("\t -> %s : %s \n","Output directory",OutputDir);
	    Is_OutputDir = true;
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-global-coordinates") == 0)
	  {
		Out_global_coordinates = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-element-coordinates") == 0)
	  {
		Out_element_coordinates = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
  	  else if(strcmp(Parse_Out_Prop[0],"Out-mass") == 0)
	  {
		Out_mass = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-density") == 0)
	  {
		Out_density = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-damage") == 0)
	  {
		Out_damage = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-nodal-idx") == 0)
	  {
		Out_nodal_idx = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-material-idx") == 0)
	  {
		Out_material_idx = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-velocity") == 0)
	  {
		Out_velocity = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-acceleration") == 0)
	  {
		Out_acceleration = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-displacement") == 0)
	  {
		Out_displacement = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-stress") == 0)
	  {
		Out_stress = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-eigenvalues-stress") == 0)
	  {
		Out_eigenvalues_stress = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-volumetric-stress") == 0)
	  {
		Out_volumetric_stress = Is_Output_Activate(Parse_Out_Prop[1]);
	  }	  
	  else if(strcmp(Parse_Out_Prop[0],"Out-strain") == 0)
	  {
		Out_strain = Is_Output_Activate(Parse_Out_Prop[1]);
	  }
	  else if(strcmp(Parse_Out_Prop[0],"Out-eigenvalues-strain") == 0)
	  {
		Out_eigenvalues_strain = Is_Output_Activate(Parse_Out_Prop[1]);
	  }	
	  else if(strcmp(Parse_Out_Prop[0],"Out-deformation-gradient") == 0)
	  {
		Out_deformation_gradient = Is_Output_Activate(Parse_Out_Prop[1]);
	  }		  
	  else if(strcmp(Parse_Out_Prop[0],"Out-energy") == 0)
	  {
		Out_energy = Is_Output_Activate(Parse_Out_Prop[1]);
	  }	
	  else
	  {
	    fprintf(stderr,"%s : %s %s \n",
		   "Error in GramsOutputs()",
		   "Undefined",Parse_Out_Prop[0]);
	    exit(EXIT_FAILURE);
	  }
	  /* Read next line and check */
	  STATUS_LINE = fgets(Line_Out_Prop,
			      sizeof(Line_Out_Prop),
			      Sim_dat);
	  Aux_Out_id = parse(Parse_Out_Prop,Line_Out_Prop," =\t\n");
	  if(strcmp(Parse_Out_Prop[0],"}") == 0){
	    break;
	  }
	}
	if(STATUS_LINE == NULL){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsOutputs()",
	       "you forget to put a } !!!");
	exit(EXIT_FAILURE);	  
	}

	if(strcmp(Parse_Out_Prop[0],"}") == 0){
	  /* Check output dir */
	  if(!Is_OutputDir){
	    fprintf(stderr,"%s : %s \n",
		    "Error in GramsOutputs()",
		    "Non output dir defined !!!");
	    exit(EXIT_FAILURE);
	  }
	  break;
	}
      }
      else{
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsOutputs()",
	       "Use this format -> GramsOutputs (Type=string) { !!!");
	exit(EXIT_FAILURE);
      }
    }
  }

  /* Close .dat file */
  /* Final message */
  fclose(Sim_dat);

}

/***************************************************************************/

static bool Is_Output_Activate(char * status)
{
	if(strcmp(status,"true") == 0)
	{
		return true;
	}
	else if(strcmp(status,"false") == 0)
	{
		return false;
	}
	else
	{
		fprintf(stderr,"%s : %s \n",
		     "Error in GramsOutputs()",
		     "Options for status output -> true/false");
	     exit(EXIT_FAILURE);
	}
}

/***************************************************************************/

/***************************************************************************/
