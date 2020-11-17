#include "nl-partsol.h"

/*
  Call global variables
*/

/*
  Auxiliar functions 
*/
static void standard_error(char *);
static void print_nodes_scalar_variable_to_csv(FILE *, Matrix, Mask, Event);
static void print_nodes_vectorial_variable_to_csv(FILE *, Matrix, Mask, Event);
static void print_nodes_tensorial_variable_to_csv(FILE *, Matrix, Mask, Event);
static void print_particles_scalar_variable_to_csv(FILE *, Matrix, Event);
static void print_particles_vectorial_variable_to_csv(FILE *, Matrix, Event);
static void print_particles_tensorial_variable_to_csv(FILE *, Matrix, Event);

/*****************************************************************/

void path_nodes_analysis_csv__InOutFun__(Matrix Nodal_Field, Mask ActiveNodes, Event Output_Commands,
										 int TimeStep, double DeltaTimeStep)
{
	int Ndim    = NumberDimensions;
	int NumCols = Nodal_Field.N_cols;
	int k_start = Output_Commands.i_start;
	int k_end   = Output_Commands.i_end;
	int k_step  = Output_Commands.i_step;
	char * Error_message;
	char * Name_File = Output_Commands.File;
	char * Name_file_Ts;
	FILE * csv_file;

	if(k_start <= TimeStep && TimeStep <= k_end && TimeStep % k_step == 0)
	{

		sprintf(Name_file_Ts,"%s/%s_%i.csv",OutputDir,Name_File,TimeStep);
		csv_file = fopen(Name_file_Ts,"w");

		/* Scalar variable */
		if(NumCols == 1)
		{
			print_nodes_scalar_variable_to_csv(csv_file,Nodal_Field,ActiveNodes,Output_Commands);
		}
		/* Vectorial variable */
		else if(NumCols == Ndim)
		{
			print_nodes_vectorial_variable_to_csv(csv_file,Nodal_Field,ActiveNodes,Output_Commands);
		}
		/* Tensorial variable */
		else if(NumCols == Ndim*Ndim)
		{
			print_nodes_tensorial_variable_to_csv(csv_file,Nodal_Field,ActiveNodes,Output_Commands);
		}

		/* Close the file */
		fclose(csv_file);
	}
}

/*****************************************************************/

void path_particles_analysis_csv__InOutFun__(Matrix Particle_Field, Event Output_Commands,
											 int TimeStep, double DeltaTimeStep)
{
	int Ndim    = NumberDimensions;
	int NumCols = Particle_Field.N_cols;
	int k_start = Output_Commands.i_start;
	int k_end   = Output_Commands.i_end;
	int k_step  = Output_Commands.i_step;
	char * Error_message;
	char * Name_File = Output_Commands.File;
	char * Name_file_Ts;
	FILE * csv_file;

	if(k_start <= TimeStep && TimeStep <= k_end && TimeStep % k_step == 0)
	{

		sprintf(Name_file_Ts,"%s/%s_%i.csv",OutputDir,Name_File,TimeStep);
		csv_file = fopen(Name_file_Ts,"w");

		/* Scalar variable */
		if(NumCols == 1)
		{
			print_particles_scalar_variable_to_csv(csv_file,Particle_Field,Output_Commands);
		}
		/* Vectorial variable */
		else if(NumCols == Ndim)
		{
			print_particles_vectorial_variable_to_csv(csv_file,Particle_Field,Output_Commands);
		}
		/* Tensorial variable */
		else if(NumCols == Ndim*Ndim)
		{
			print_particles_tensorial_variable_to_csv(csv_file,Particle_Field,Output_Commands);
		}

		/* Close the file */
		fclose(csv_file);
	}
}

/*****************************************************************/

static void standard_error(char * Error_message)
{
	fprintf(stderr,"%s : %s !!! \n",
	   "Error in GramsOutputs()",Error_message);
    exit(EXIT_FAILURE);
}

/*****************************************************************/

static void print_nodes_scalar_variable_to_csv(FILE * csv_file, Matrix Variable, Mask ActiveNodes, Event Output_Commands)
{
	int Idx;
	int Idx_mask;

	for(int i = 0 ; i< Output_Commands.Lenght_Path ; i++)
	{
		Idx = Output_Commands.Idx_Path[i];
		Idx_mask = ActiveNodes.Nodes2Mask[Idx];
		
		if(Idx_mask == - 1)
		{
        	fprintf(csv_file,"0.0 \n");
    	}
    	else
    	{
          	fprintf(csv_file,"%lf \n",Variable.nV[Idx_mask]);  		
    	}
	}

}

/*****************************************************************/

static void print_nodes_vectorial_variable_to_csv(FILE * csv_file, Matrix Variable, Mask ActiveNodes, Event Output_Commands)
{
	int Idx;
	int Idx_mask;
	int Ndim = NumberDimensions;
	
	for(int i = 0 ; i<Output_Commands.Lenght_Path ; i++)
	{
		Idx = Output_Commands.Idx_Path[i];
		Idx_mask = ActiveNodes.Nodes2Mask[Idx];

		if(Idx_mask == - 1)
		{
			if(Ndim == 2)
			{
				fprintf(csv_file,"0.0;0.0 \n");
			}
			else if(Ndim == 3)
			{
				fprintf(csv_file,"0.0;0.0;0.0 \n");
			}
		}
		else
		{
			if(Ndim == 2)
			{
				fprintf(csv_file,"%lf;%lf \n",
						Variable.nM[Idx_mask][0],
						Variable.nM[Idx_mask][1]);
			}
			else if(Ndim == 3)
			{
				fprintf(csv_file,"%lf;%lf;%lf \n",
						Variable.nM[Idx_mask][0],
						Variable.nM[Idx_mask][1],
						Variable.nM[Idx_mask][2]);
			}			
		}
	}

}


/*****************************************************************/

static void print_nodes_tensorial_variable_to_csv(FILE * csv_file, Matrix Variable, Mask ActiveNodes, Event Output_Commands)
{
	int Idx;
	int Idx_mask;
	int Ndim = NumberDimensions;
	
	for(int i = 0 ; i<Output_Commands.Lenght_Path ; i++)
	{
		Idx = Output_Commands.Idx_Path[i];
		Idx_mask = ActiveNodes.Nodes2Mask[Idx];

		if(Idx_mask == - 1)
		{
			if(Ndim == 2)
			{
				fprintf(csv_file,
						"0.0;0.0;0.0;0.0 \n");
			}
			else if(Ndim == 3)
			{
				fprintf(csv_file,
						"0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0 \n");
			}
		}
		else
		{
			if(Ndim == 2)
			{
				fprintf(csv_file,
						"%lf;%lf;%lf;%lf \n",
						Variable.nM[Idx_mask][0],
						Variable.nM[Idx_mask][1],
						Variable.nM[Idx_mask][2],
						Variable.nM[Idx_mask][3]);
			}
			else if(Ndim == 3)
			{
				fprintf(csv_file,
						"%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf \n",
						Variable.nM[Idx_mask][0],
						Variable.nM[Idx_mask][1],
						Variable.nM[Idx_mask][2],
						Variable.nM[Idx_mask][3],
						Variable.nM[Idx_mask][4],
						Variable.nM[Idx_mask][5],
						Variable.nM[Idx_mask][6],
						Variable.nM[Idx_mask][7],
						Variable.nM[Idx_mask][8]);
			}
		}
	}

}

/*****************************************************************/

static void print_particles_scalar_variable_to_csv(FILE * csv_file, Matrix Variable, Event Output_Commands)
{
	int Idx;

	for(int i = 0 ; i< Output_Commands.Lenght_Path ; i++)
	{
		Idx = Output_Commands.Idx_Path[i];

        fprintf(csv_file,"%lf \n",Variable.nV[Idx]);
	}

}

/*****************************************************************/

static void print_particles_vectorial_variable_to_csv(FILE * csv_file, Matrix Variable, Event Output_Commands)
{
	int Idx;
	int Ndim = NumberDimensions;
	
	for(int i = 0 ; i<Output_Commands.Lenght_Path ; i++)
	{
		Idx = Output_Commands.Idx_Path[i];

		if(Ndim == 2)
			{
				fprintf(csv_file,"%lf;%lf \n",
						Variable.nM[Idx][0],
						Variable.nM[Idx][1]);
			}
		else if(Ndim == 3)
			{
				fprintf(csv_file,"%lf;%lf;%lf \n",
						Variable.nM[Idx][0],
						Variable.nM[Idx][1],
						Variable.nM[Idx][2]);
			}
	}

}


/*****************************************************************/

static void print_particles_tensorial_variable_to_csv(FILE * csv_file, Matrix Variable, Event Output_Commands)
{
	int Idx;
	int Ndim = NumberDimensions;
	
	for(int i = 0 ; i<Output_Commands.Lenght_Path ; i++)
	{
		Idx = Output_Commands.Idx_Path[i];

		if(Ndim == 2)
			{
				fprintf(csv_file,"%lf;%lf;%lf;%lf \n",
						Variable.nM[Idx][0],
						Variable.nM[Idx][1],
						Variable.nM[Idx][2],
						Variable.nM[Idx][3]);
			}
		else if(Ndim == 3)
			{
				fprintf(csv_file,"%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf \n",
						Variable.nM[Idx][0],
						Variable.nM[Idx][1],
						Variable.nM[Idx][2],
						Variable.nM[Idx][3],
						Variable.nM[Idx][4],
						Variable.nM[Idx][5],
						Variable.nM[Idx][6],
						Variable.nM[Idx][7],
						Variable.nM[Idx][8]);
			}
	}

}

/*****************************************************************/

