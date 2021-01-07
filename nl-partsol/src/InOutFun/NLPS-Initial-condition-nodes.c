#include "nl-partsol.h"
#include <sys/stat.h>

/*
  Auxiliar functions 
*/
static void standard_error(char *);
static void standard_output(char *);
static bool Check_List_Nodes(char *);
static FILE * Open_and_Check_simulation_file(char *);
static void Read_list_nodes_file(char *, char *, char *);
static Tensor Read_initial_values(FILE *);
static Tensor Read_initial_vectorial_field(char *);
static void interpolate_initial_values_particles(GaussPoint, Mesh, char *, Tensor);

/***************************************************************************/

void Initial_condition_nodes__InOutFun__(char * Name_File, GaussPoint MPM_Mesh, Mesh FEM_Mesh)
/*
  Initial-nodal-values (NODES=ListInit.txt) 
  {
  	VELOCITY=[5.0,0.0,0.0]
  }
*/
{
	/* Simulation file */
	FILE * Sim_dat;

	char Route_NODES[MAXC];
	Tensor Field;

	/* Special variables for line-reading */
  	char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  	char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  	int nkwords; /* Number of element in the line , just for check */

  	/* Initial message */  
  	puts("*************************************************");
  	printf(" \t %s : \n\t %s \n", "* Read initial nodal values ", Name_File);

	/* Open and check file */
  	Sim_dat = Open_and_Check_simulation_file(Name_File);

	/* Read the file line by line */
  	while( fgets(line, sizeof line, Sim_dat) != NULL )
  	{

    	/* Read the line with the space as separators */
    	nkwords = parse (kwords, line," \n\t");
    	if (nkwords < 0)
    	{
    		standard_error("Parser failed");
    	}

		/* Read Initial-nodal-values */
    	if ((nkwords > 0) && (strcmp(kwords[0],"Initial-nodal-values") == 0 ))
    	{

    		/* Read route with the nodes for the initial condition */
    		Read_list_nodes_file(Route_NODES, kwords[1], Name_File);

	    	/* Read csv parameters */
   			Field = Read_initial_values(Sim_dat);

	   		/* Fill csv parameters and period */
	   		interpolate_initial_values_particles(MPM_Mesh, FEM_Mesh, Route_NODES, Field);

	   		/* Free memory */
	   		free__TensorLib__(Field);

    	}

  	}

  /* Close .dat file */
  /* Final message */
  fclose(Sim_dat);

}

/***************************************************************************/

static void standard_error(char * Error_message)
{
	fprintf(stderr,"%s : %s !!! \n",
	   "Error in Initial-nodal-values()",Error_message);
    exit(EXIT_FAILURE);
}

/***************************************************************************/

static void standard_output(char * Status_message)
{
	fprintf(stdout,"%s \n",Status_message);
}

/***************************************************************************/

static FILE * Open_and_Check_simulation_file(char * Name_File)
{
  FILE * Simulation_file = fopen(Name_File,"r");  
  char Error_message[MAXW];
  
  if (Simulation_file==NULL)
  {
  	sprintf(Error_message,"%s %s","Incorrect lecture of",Name_File);
	standard_error(Error_message); 
  }  

  return Simulation_file;
}

/***************************************************************************/

static void Read_list_nodes_file(char * Route_NODES, char * List_nodes_string, char * Name_File)
{
	char Error_message[MAXW] = {0};
	char Route_Path[MAXC] = {0};
	
	int Parser_status;
	char * Parameter_pars[MAXW] = {NULL};

	bool Is_NODES = false;

	Parser_status = parse (Parameter_pars, List_nodes_string,"(=)");
	 
    /* Check list of nodes */
    if((strcmp(Parameter_pars[0],"NODES") == 0) && (Parser_status== 2))
      {
      	generate_route(Route_Path,Name_File);
	  	sprintf(Route_NODES,"%s%s",Route_Path,Parameter_pars[1]);
      	Is_NODES = Check_List_Nodes(Route_NODES);
      }
    else
      {
	  	sprintf(Error_message,"%s %s","Undefined",Parameter_pars[0]);
	  	standard_error(Error_message); 
      }

     if(!Is_NODES)
     {
     	sprintf(Error_message,"%s","File with the list of nodes was not defined");
	  	standard_error(Error_message);
     }

}

/***************************************************************************/

static bool Check_List_Nodes(char * PATH_Name)
{
	struct stat info;
	stat(PATH_Name,&info);
	char Error_message[MAXW];
	bool status_check;

	if(S_ISREG(info.st_mode))
	{
		printf("\t -> %s : %s \n","List of nodes",PATH_Name);
		status_check = true;
	}
	else
	{
		sprintf(Error_message,"\t -> %s : %s %s \n","List of nodes",PATH_Name,"does not exists");
		standard_error(Error_message);
	} 

	return status_check;
}

/***************************************************************************/

static Tensor Read_initial_values(FILE * Simulation_file)
{

	Tensor Field;

	/* Log messages */
	char Status_message[MAXW];
	char Error_message[MAXW];

	/* Variables for reading purposes */
	char Parameter_line[MAXC] = {0};
  	char * Parameter_pars[MAXW] = {NULL};
  	int Parser_status;

  	/* Check variables for sintax */
  	bool Is_VELOCITY = false;
  	bool Is_Open = false;
  	bool Is_Close = false;

  	while(fgets(Parameter_line, sizeof(Parameter_line), Simulation_file) != NULL)
  	{
  		/* Parse line */  	
		Parser_status = parse(Parameter_pars,Parameter_line," =\t\n");

		if((strcmp(Parameter_pars[0],"{") == 0) && (Parser_status == 1))
		{
			Is_Open = true;
		}
		else if((strcmp(Parameter_pars[0],"VELOCITY") == 0) && (Parser_status == 2)) 
 	  	{
 	  		Field = Read_initial_vectorial_field(Parameter_pars[1]);
 	  		sprintf(Status_message,"\t -> %s","The initial velocity is :");
			standard_output(Status_message);
			print__TensorLib__(Field);
 	  		Is_VELOCITY = true;
	  	}
	  	else if((strcmp(Parameter_pars[0],"}") == 0) && (Parser_status == 1))
	  	{
	  		Is_Close = true;
	    	break;
	  	}	  
		else if(Parser_status > 0)
	  	{
	  		sprintf(Error_message,"%s %s","Undefined",Parameter_pars[0]);
	  		standard_error(Error_message); 
	  	}
	
	}

	if(!Is_Open && !Is_Close)
	{
	  	sprintf(Error_message,"%s","Unbalanced curls {}");
	  	standard_error(Error_message); 
	}

	if(!Is_VELOCITY)
	{
		sprintf(Error_message,"%s","You forget to define the initial condition");
	  	standard_error(Error_message);
	}

	return Field;
}

/***************************************************************************/

static Tensor Read_initial_vectorial_field(char * String_Values)
{
	int Ndim = NumberDimensions;
	char Error_message[MAXW];

	Tensor Field;

	char * Parameter_pars[MAXW] = {NULL};
  	int Parser_status;

	Parser_status = parse(Parameter_pars,String_Values,"[,]");

	if(Parser_status == Ndim)
	{
		Field = alloc__TensorLib__(1);

		for(int i = 0 ; i<Ndim ; i++)
		{
			Field.n[i] = atof(Parameter_pars[i]);
		}
	}
	else
	{
		sprintf(Error_message,"%s %i-D %s %i-D %s","You put a",Parser_status,"vector in a",Ndim,"problem");
	  	standard_error(Error_message);
	}

	return Field;
}

/***************************************************************************/

static void interpolate_initial_values_particles(GaussPoint MPM_Mesh, Mesh FEM_Mesh, char * Route_NODES, Tensor Field)
{

	/* Auxiliar variables to read the chain of nodes */
	ChainPtr Chain_Nodes = NULL;
	int Lenght_Path;
	int * Idx_Path;
	int * Nodes_mesh = (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh,sizeof(int));
	int Idx;

	/* Auxiliar variables for the particles */
	int Ndim = NumberDimensions;
	int Np = MPM_Mesh.NumGP;
	Element Nodes_p; /* Element for each particle */
	int Ap;
	Matrix ShapeFunction_p; /* Value of the shape-function in the particle */
	double ShapeFunction_pI; /* Nodal value for the particle */
	double Field_pI;

	Chain_Nodes = File2Chain(Route_NODES);
    Lenght_Path = lenght__SetLib__(Chain_Nodes);
    Idx_Path = set_to_memory__SetLib__(Chain_Nodes,Lenght_Path);
    free__SetLib__(&Chain_Nodes);

	/* Nodes with initial condition == 1 */
    for(int i = 0 ; i<Lenght_Path ; i++)
    {
    	Idx = Idx_Path[i];
    	Nodes_mesh[Idx] = 1;
    }

    /* iterate over the particles */
  	for(int p = 0 ; p<Np ; p++)
    {
    	/* Define element of the particle */
      	Nodes_p = nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      	/*
		Evaluate the shape function and gradient in the coordinates of the particle 
      	*/
      	ShapeFunction_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

      	/* Iterate over the nodes of the particle */
      	for(int A = 0; A<Nodes_p.NumberNodes; A++)
		{
			Ap = Nodes_p.Connectivity[A];	

			/*
	    		Evaluate the GP function in the node 
	  		*/
	  		ShapeFunction_pI = ShapeFunction_p.nV[A];

	  		/*
	    		Update initial value for each particle
	  		*/
	  		for(int i = 0 ; i<Ndim ; i++)
	    	{
	      		Field_pI = ShapeFunction_pI*Nodes_mesh[Ap]*Field.n[i];	      
	      		MPM_Mesh.Phi.vel.nM[p][i] += Field_pI;
	    	}	 
	    }
	}

    free(Idx_Path);
    free(Nodes_mesh);
}

/***************************************************************************/
