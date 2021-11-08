#include "nl-partsol.h"
#include <sys/stat.h>

/*
  Call global variables
*/
Event * Out_nodal_path_csv;
int Number_Out_nodal_path_csv;
int NumTimeStep;
char OutputDir[MAXC];

/*
	Auxiliar structures
*/
typedef struct {

	char DIR_Name[MAXC];
	char PATH_Name[MAXC];
	bool Out_velocity;
  	bool Out_acceleration;
  	bool Out_displacement;
  	bool Out_forces;
  	bool Out_reactions;
  	bool Out_residual; 

} Parameters;

typedef struct {

	int i_start;
	int i_step;
	int i_end;

} Intervals;

/*
  Auxiliar functions 
*/
static void standard_error(char *);
static int Number_nodal_csv_events(char *);
static FILE * Open_and_Check_simulation_file(char *);
static bool Check_Output_directory(char *);
static bool Check_Path(char *);
static Intervals read_CSV_Intervals(char *);
static Parameters read_CSV_Parameters(FILE *, char *);
static bool Is_Output_Activate(char *, char *);
static Event fill_CSV_Parameters(Intervals,Parameters);

/**********************************************************************/

void NLPS_Out_nodal_path_csv__InOutFun__(char * Name_File)
/*
  Example : 
  Out-nodal-path-csv (i_ini=0;i_step=10;i_end=100) 
  {
  	DIR=test/Sulsky_MPM	
  	PATH="Nodes.txt"
	VAR=velocity
  }
*/
{
  /* Simulation file */
  FILE * Sim_dat;

  char Error_message[MAXW];

  int i_CSV = 0;

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */


  Parameters CSV_Parameters;
  Intervals CSV_Intervals;

  int Num_NCE = Number_nodal_csv_events(Name_File);
  Out_nodal_path_csv = (Event *)malloc(Num_NCE*sizeof(Event));
  Number_Out_nodal_path_csv = Num_NCE;

  if(Num_NCE > 0)
  {
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

		/* Read Out-nodal-path-csv */
	    if ((nkwords > 0) && (strcmp(kwords[0],"Out-nodal-path-csv") == 0 ))
	    {

	    	/* Read output period */
	    	CSV_Intervals = read_CSV_Intervals(kwords[1]);
		    /* Read csv parameters */
	   		CSV_Parameters = read_CSV_Parameters(Sim_dat,Name_File);
		   	/* Fill csv parameters and period */
		   	Out_nodal_path_csv[i_CSV] = fill_CSV_Parameters(CSV_Intervals,CSV_Parameters);

		   	/* Update counter */
		   	i_CSV++;

	    }

	  }

	   /* Close .dat file */
  	/* Final message */
  	fclose(Sim_dat);

  }
  else
  {
  	printf("\t * %s \n","No nodal path results were defined");
  }

}


/***************************************************************************/

static void standard_error(char * Error_message)
{
	fprintf(stderr,"%s : %s !!! \n",
	   "Error in Out-nodal-path-csv()",Error_message);
    exit(EXIT_FAILURE);
}

/***************************************************************************/

static int Number_nodal_csv_events(char * Name_File)
{
	int NNCE = 0;

	/* Special variables for line-reading */
  	char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  	char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  	int nkwords; /* Number of element in the line , just for check */

	/* Simulation file */
  	FILE * Sim_dat = Open_and_Check_simulation_file(Name_File);
	
	/* Read the file line by line */
  	while( fgets(line, sizeof line, Sim_dat) != NULL )
  	{
  		nkwords = parse (kwords, line," \n\t");

  		if ((nkwords > 0) && (strcmp(kwords[0],"Out-nodal-path-csv") == 0 ))
    	{
    		NNCE++;
    	}
  	}

	fclose(Sim_dat);

	return NNCE;
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
		standard_error(Error_message); 
	} 

	return status_check;
}

/***************************************************************************/

static bool Check_Path(char * PATH_Name)
{
	struct stat info;
	stat(PATH_Name,&info);
	char Error_message[MAXW];
	bool status_check;

	if(S_ISREG(info.st_mode))
	{
		printf("\t -> %s : %s \n","Path file",PATH_Name);
		status_check = true;
	}
	else
	{
		sprintf(Error_message,"\t -> %s : %s %s \n","Path file",PATH_Name,"does not exists");
		standard_error(Error_message);
	} 

	return status_check;
}

/***************************************************************************/

static Intervals read_CSV_Intervals(char * Interval_message)
{
	Intervals CSV_Intervals;

	char Error_message[MAXW];
	
	int Interval_status_1;
	char * Aux_Parse_1[MAXW] = {NULL};
	int Interval_status_2;
	char * Aux_Parse_2[MAXW] = {NULL};

  	/* Set parameters to default */
  	CSV_Intervals.i_start = 0;
  	CSV_Intervals.i_step = 1;
  	CSV_Intervals.i_end = NumTimeStep;

	Interval_status_1 = parse (Aux_Parse_1, Interval_message,"(;)");

	/* Check format */
	if(Interval_status_1 > 3)
	  {
	  	standard_error("You have exceded the maximum number of parameters");
      }
	

	for(int i = 0 ; i<Interval_status_1 ; i++)
	{
		Interval_status_2 = parse (Aux_Parse_2, Aux_Parse_1[i],"="); 

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
	      	standard_error(Error_message);
    	}
	}
	 
    /* Check interval output */
    if(CSV_Intervals.i_end < CSV_Intervals.i_step)
      {
      	standard_error("The result interval step should be less than final time");
      }
    if(CSV_Intervals.i_end < CSV_Intervals.i_start)
      {
      	standard_error("The initial result time step should be less than the final result time");
      }
    if(CSV_Intervals.i_end > NumTimeStep)
      {
      	standard_error("The final result time step should be less than the final time");
      }

	/* Print some info */
    printf("\t -> %s : %i \n","i_start",CSV_Intervals.i_start);
    printf("\t -> %s : %i \n","i_step",CSV_Intervals.i_step);
    printf("\t -> %s : %i \n","i_end",CSV_Intervals.i_end);

	return CSV_Intervals;
}

/***************************************************************************/

static Parameters read_CSV_Parameters(FILE * Simulation_file, char * Name_File)
{
	Parameters Output_csv;

	/* Default output paramters	*/
	Output_csv.Out_velocity = false;
  	Output_csv.Out_acceleration = false;
  	Output_csv.Out_displacement = false;
  	Output_csv.Out_forces = false;
  	Output_csv.Out_reactions = false;
  	Output_csv.Out_residual = false;

	char Error_message[MAXW];

	/* Variables for reading purposes */
	char Line_Out_Prop[MAXC] = {0};
  	char * Parse_Out_Prop[MAXW] = {NULL};
  	int Aux_Out_id;
  	char Route_Outs[MAXC] = {0};
  	char Route_Path[MAXC] = {0};

  	/* Check variables for sintax */
  	bool Is_DIR = false;
  	bool Is_PATH = false;
  	bool Is_Open = false;
  	bool Is_Close = false;

  	while(fgets(Line_Out_Prop, sizeof(Line_Out_Prop), Simulation_file) != NULL)
  	{
  		/* Parse line */  	
		Aux_Out_id = parse(Parse_Out_Prop,Line_Out_Prop," =\t\n");

		if((strcmp(Parse_Out_Prop[0],"{") == 0) && (Aux_Out_id == 1))
		{
			Is_Open = true;
		}
		else if((strcmp(Parse_Out_Prop[0],"DIR") == 0) && (Aux_Out_id == 2)) 
 	  	{
 	  		generate_route(Route_Outs,Name_File);
	    	sprintf(Output_csv.DIR_Name,"%s%s",Route_Outs,Parse_Out_Prop[1]);
	    	Is_DIR = Check_Output_directory(Output_csv.DIR_Name);
	  	}
		else if((strcmp(Parse_Out_Prop[0],"PATH") == 0) && (Aux_Out_id == 2))
	  	{
	  		generate_route(Route_Path,Name_File);
	  		sprintf(Output_csv.PATH_Name,"%s%s",Route_Path,Parse_Out_Prop[1]);
	  		Is_PATH = Check_Path(Output_csv.PATH_Name);
	  	}
	  	else if(strcmp(Parse_Out_Prop[0],"Out-velocity") == 0)
	  	{
			Output_csv.Out_velocity = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
	  	}
	  	else if(strcmp(Parse_Out_Prop[0],"Out-acceleration") == 0)
	  	{
			Output_csv.Out_acceleration = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
	  	}
	  	else if(strcmp(Parse_Out_Prop[0],"Out-displacement") == 0)
	  	{
			Output_csv.Out_displacement = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
	  	}
	  	else if(strcmp(Parse_Out_Prop[0],"Out-forces") == 0)
	  	{
			Output_csv.Out_forces = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
	  	}
	  	else if(strcmp(Parse_Out_Prop[0],"Out-reactions") == 0)
	  	{
			Output_csv.Out_reactions = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
	  	}
	  	else if(strcmp(Parse_Out_Prop[0],"Out-residual") == 0)
	  	{
			Output_csv.Out_residual = Is_Output_Activate(Parse_Out_Prop[0],Parse_Out_Prop[1]);
	  	}
	  	else if((strcmp(Parse_Out_Prop[0],"}") == 0) && (Aux_Out_id == 1))
	  	{
	  		Is_Close = true;
	    	break;
	  	}	  
		else if(Aux_Out_id > 0)
	  	{
	  		sprintf(Error_message,"%s %s","Undefined",Parse_Out_Prop[0]);
	  		standard_error(Error_message); 
	  	}
	
	}

	if(Is_Open && Is_Close)
	{
		if(!Is_DIR || !Is_PATH)
		{
	  		sprintf(Error_message,"%s","You forgot to define something");
	  		standard_error(Error_message); 
		}
	}
	else
	{
	  	sprintf(Error_message,"%s","Unbalanced curls {}");
	  	standard_error(Error_message); 
	}

	/* Check syntax	*/

	return Output_csv;
}

/***************************************************************************/

static bool Is_Output_Activate(char * output_field, char * status_text)
{
	bool status;
	char Error_message[MAXW];

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
		standard_error(Error_message); 
	}

	return status;
}

/***************************************************************************/

static Event fill_CSV_Parameters(Intervals CSV_Intervals,Parameters CSV_Parameters)
{
	Event CSV_Event;

	/* Set outputs to default */
	CSV_Event.Out_csv_particles_path_Damage = false;
  	CSV_Event.Out_csv_particles_path_Velocity = false;
  	CSV_Event.Out_csv_particles_path_Acceleration = false;
  	CSV_Event.Out_csv_particles_path_Displacement = false;
  	CSV_Event.Out_csv_particles_path_Stress = false;
  	CSV_Event.Out_csv_particles_path_Strain = false;
  	CSV_Event.Out_csv_particles_path_Deformation_gradient = false;
  	  	
	CSV_Event.Out_csv_nodes_path_Velocity = false;
  	CSV_Event.Out_csv_nodes_path_Acceleration = false;
  	CSV_Event.Out_csv_nodes_path_D_Displacement = false;
  	CSV_Event.Out_csv_nodes_path_Forces = false;
  	CSV_Event.Out_csv_nodes_path_Reactions = false;
  	CSV_Event.Out_csv_nodes_path_Residual = false; 

	/* Auxiliar variable to read the chain of nodes */
	ChainPtr Chain_Nodes = NULL;

	/* Read outputs intervals */
	CSV_Event.i_start = CSV_Intervals.i_start;
	CSV_Event.i_step  = CSV_Intervals.i_step;
	CSV_Event.i_end   = CSV_Intervals.i_end;

	/* Write name of the ouput directory */
	strcpy(CSV_Event.Directory,CSV_Parameters.DIR_Name);

	/* Write nodal path */
	Chain_Nodes = File2Chain(CSV_Parameters.PATH_Name);
    CSV_Event.Lenght_Path = lenght__SetLib__(Chain_Nodes);
    CSV_Event.Idx_Path = set_to_memory__SetLib__(Chain_Nodes,CSV_Event.Lenght_Path);
    free__SetLib__(&Chain_Nodes);

    /* Select ouput variable */ 
	if(CSV_Parameters.Out_velocity)
	{
		CSV_Event.Out_csv_nodes_path_Velocity = true;
	}
  	else if(CSV_Parameters.Out_acceleration)
	{
	  	CSV_Event.Out_csv_nodes_path_Acceleration = true;
	}
	else if(CSV_Parameters.Out_displacement)
	{
	  	CSV_Event.Out_csv_nodes_path_D_Displacement = true;
	}
	else if(CSV_Parameters.Out_forces)
	{
	  	CSV_Event.Out_csv_nodes_path_Forces = true;
	}
	else if(CSV_Parameters.Out_reactions)
	{
	  	CSV_Event.Out_csv_nodes_path_Reactions = true;
	}
	else if(CSV_Parameters.Out_residual)
	{
	  	CSV_Event.Out_csv_nodes_path_Residual = true;
	}

	return CSV_Event;
}

/***************************************************************************/

