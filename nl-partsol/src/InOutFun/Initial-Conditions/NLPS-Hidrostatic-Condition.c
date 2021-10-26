#include "nl-partsol.h"
#include <sys/stat.h>


/*
	Auxiliar structures
*/
typedef struct {

	ChainPtr Particles_List;

	Tensor Origin;

	Tensor Direction;

	double Gravity;

	int MatIndx;

} Parameters;

/*
  Auxiliar functions 
*/

#ifdef _WIN32
static char * delimiters = " =\t\r\n"; 
#else
static char * delimiters = " =\t\n"; 
#endif

static Parameters Read_Hidrostatic_Parameters__InOutFun__(FILE *,int *);
static void assign_hidrostatic_condition(Parameters,Particle,int,int *);

/*****************************************************************/

void Hidrostatic_condition_particles__InOutFun__(
	char * Name_File, 
	Particle MPM_Mesh, 
	int GPxElement,
	int * STATUS)
/*!
 * Hydrostatic-condition
 * {
 * 	Particles = ListInit.txt
 * 	Origin = {0 ; 10}
 * 	Direction = {0 ; -1}
 * 	Gravity = 10
 * 	MatIndx = 0
 * } 
 * */
{
	// Simulation file
	FILE * Simulation_file;

	// Error variables
	int INFO_Read_file = 0;
	int INFO_Read_Parser = 0;
	int INFO_Read_Parameters = 0;
	int INFO_Assign_Parameters = 0;

	char Route_NODES[MAXC];
	Parameters hidrostatic_parameters;

	// Special variables for line-reading
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  // Initial message
  puts("*************************************************");
  printf(" \t %s \n", "* Read hidrostatic condition");

	// Open and check file
  Simulation_file = Open_and_Check_File__InOutFun__(Name_File,&INFO_Read_file);

  if(INFO_Read_file)
	{
		fprintf(stderr,"Error in : %s\n",
		"Hidrostatic_condition_particles__InOutFun__");
		(*STATUS) = 1;
		return;
	}

	// Read the file line by line
  while(fgets(line, sizeof line, Simulation_file) != NULL)
  {

  	// Read the line with the space as separators
    nkwords = Parse_string__InOutFun__(kwords, line," \n\t",&INFO_Read_Parser);
    
    if(INFO_Read_Parser)
    {
   		fprintf(stderr,"Error in : %s\n",
				"Hidrostatic_condition_particles__InOutFun__");
      (*STATUS) = 1;
      return;
    }

		// Read Hydrostatic-condition
    if ((nkwords > 0) && (strcmp(kwords[0],"Hydrostatic-condition") == 0 ))
    {

   		hidrostatic_parameters = Read_Hidrostatic_Parameters__InOutFun__(Simulation_file,&INFO_Read_Parameters);

   		if(INFO_Read_Parameters)
   		{
   			fprintf(stderr,"Error in : %s\n",
					"Hidrostatic_condition_particles__InOutFun__");
      	(*STATUS) = 1;
      	return;
   		}

		  assign_hidrostatic_condition(hidrostatic_parameters,MPM_Mesh,GPxElement,&INFO_Assign_Parameters);

		  if(INFO_Assign_Parameters)
			{
				fprintf(stderr,"Error in : %s\n",
					"Hidrostatic_condition_particles__InOutFun__");
      	(*STATUS) = 1;
      	return;
      }

	   	free__SetLib__(&hidrostatic_parameters.Particles_List);
	   	free__TensorLib__(hidrostatic_parameters.Origin);
	   	free__TensorLib__(hidrostatic_parameters.Direction);

    }

  }

  // Close .dat file
  fclose(Simulation_file);

}

/***************************************************************************/

static Parameters Read_Hidrostatic_Parameters__InOutFun__(
	FILE * Simulation_file,
	int * STATUS)
/*!
 *
 * Hydrostatic-condition
 * {
 * Particles = ListInit.txt
 * Origin = {0, 10}
 * Direction = {0,-1}
 * Gravity = 10
 * MatIndx = 0
 * }
 *  
 * */
{

  int Ndim = NumberDimensions;
  int INFO_Particles = 0;
  int INFO_Origin = 0;
	int INFO_Direction = 0;

	Parameters hidrostatic_parameters;

  /* Variables for reading purposes */
  char Parameter_line[MAXC] = {0};
  char * Parameter_pars[MAXW] = {NULL};
  int Parser_status;

  /* Check variables for sintax */
  bool Is_Open = false;
  bool Is_Particle_File = false;
  bool Is_Origin = false;
  bool Is_Direction = false;
  bool Is_Gravity = false;
  bool Is_MatIndx = false;
  bool Is_Close = false;

  while(fgets(Parameter_line, sizeof(Parameter_line), Simulation_file) != NULL)
  {

    /* Parse line */    
    Parser_status = parse(Parameter_pars,Parameter_line,delimiters);

    if((strcmp(Parameter_pars[0],"{") == 0) && (Parser_status == 1))
    {
      Is_Open = true;
    }
    else if(strcmp(Parameter_pars[0],"Particles") == 0)
    {
    	Is_Particle_File = true;
      hidrostatic_parameters.Particles_List = File_to_Chain__InOutFun__(Parameter_pars[1], &INFO_Particles);
    }
    else if(strcmp(Parameter_pars[0],"Origin") == 0)
    {
    	Is_Origin = true;
			hidrostatic_parameters.Origin = Read_Vector__InOutFun__(Parameter_pars[1],&INFO_Origin);

      if(INFO_Origin)
      {
      	fprintf(stderr,"Error in : %s. %s\n",
				"Read_Hidrostatic_Parameters__InOutFun__","Incorrect number of dimensions for the origin vector");
       	(*STATUS) = 1;
      	return hidrostatic_parameters;
      }

    }
    else if(strcmp(Parameter_pars[0],"Direction") == 0)
    {
    	Is_Direction = true;
      hidrostatic_parameters.Direction = Read_Vector__InOutFun__(Parameter_pars[1],&INFO_Direction);

      if(INFO_Direction)
      {
      	fprintf(stderr,"Error in : %s. %s\n",
				"Read_Hidrostatic_Parameters__InOutFun__","Incorrect number of dimensions for the direction vector");
       	(*STATUS) = 1;
      	return hidrostatic_parameters;
      }

    }
    else if(strcmp(Parameter_pars[0],"Gravity") == 0)
    {
    	Is_Gravity = true;
      hidrostatic_parameters.Gravity = atof(Parameter_pars[1]);
    } 
    else if(strcmp(Parameter_pars[0],"MatIndx") == 0)
    {
    	Is_MatIndx = true;
      hidrostatic_parameters.MatIndx = atoi(Parameter_pars[1]);
    }
    else if((strcmp(Parameter_pars[0],"}") == 0) && (Parser_status == 1))
    {
    	Is_Close = true;
    	break;
    }
    else if(Parser_status > 0)
    {
			fprintf(stderr,"Error in : %s. %s %s\n",
				"Read_Hidrostatic_Parameters__InOutFun__","Undefined paramter",Parameter_pars[0]);
      (*STATUS) = 1;
      return hidrostatic_parameters;
    }

  }

  if(!Is_Open || !Is_Close)
  {
  	fprintf(stderr,"Error in : %s. %s\n","Read_Hidrostatic_Parameters__InOutFun__","Unbalanced braces");
  	(*STATUS) = 1;
  	return hidrostatic_parameters;
  }

  if(!Is_Particle_File || 
  	!Is_Origin || 
  	!Is_Direction || 
    !Is_Gravity || 
    !Is_MatIndx)
  {
    fprintf(stderr,"%s : %s \n","Error in Read_Hidrostatic_Parameters__InOutFun__","Some parameters are missed");
    fputs(Is_Particle_File ? "Particle File : true \n" : "Particle File : false \n", stdout);
    fputs(Is_Origin  ? "Origin : true \n" : "Origin : false \n", stdout);
    fputs(Is_Direction ? "Direction : true \n" : "Direction : false \n", stdout);
    fputs(Is_Gravity ? "Gravity : true \n" : "Gravity : false \n", stdout);
    fputs(Is_MatIndx ? "MatIndx : true \n" : "MatIndx : false \n", stdout);
		(*STATUS) = 1;
    return hidrostatic_parameters;
  }

	return hidrostatic_parameters;
}

/***************************************************************************/

static void assign_hidrostatic_condition(
	Parameters hidrostatic_parameters,
	Particle MPM_Mesh,
	int GPxElement,
	int * STATUS)
{
	int Ndim = NumberDimensions;
	int p;
	double g = hidrostatic_parameters.Gravity;
	double rho = MPM_Mesh.Mat[hidrostatic_parameters.MatIndx].rho;
	double P0 = MPM_Mesh.Mat[hidrostatic_parameters.MatIndx].ReferencePressure;
	double gamma = g*rho;
	double distance;
	double Pressure;
	Tensor Direction = hidrostatic_parameters.Direction;
	Tensor X_0 = hidrostatic_parameters.Origin;
	Tensor X_p;
	Tensor X_0p;
	ChainPtr Particles_List = hidrostatic_parameters.Particles_List;

	while(Particles_List != NULL)
	{
		for(int i = 0 ; i<GPxElement ; i++)
		{
			p = (Particles_List->Idx)*GPxElement+i;

			if((p < 0) || (p >= MPM_Mesh.NumGP))
			{
				fprintf(stderr,"%s : %s %i %s \n",
					"assign_hidrostatic_condition",
					"The particle",p,"does not exists");
				(*STATUS) = 1;
				return;
			}

			X_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.x_GC.nM[p],1);

			X_0p = subtraction__TensorLib__(X_p, X_0);

			distance = inner_product__TensorLib__(X_0p, Direction);

			Pressure = gamma*distance + P0;

		  for(int j = 0 ; j<Ndim ; j++)
		  {
		    MPM_Mesh.Phi.Stress.nM[p][j*Ndim + j] += Pressure;
		  }

		  free__TensorLib__(X_0p);

		}

		Particles_List = Particles_List->next;
	}

}

/***************************************************************************/

