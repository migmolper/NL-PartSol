
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

/***************************************************************************/


void ReadDatFile(char * Name_File)
/*
  Read data from the .DAT file and initialize the variables
  
  Inputs
  - Name_file : Name of the file
  - 
  Outputs
  - Conectivity matrix
  - Coordenates of the nodes
  - phi_n : Values of the variables in the n step, 
        initializerd with the Initial conditions
  - DeltaT : Time-step
  - A_el : Area of the element
  - type_elem : Type of the element (1)
  - 
  - N_nodes : Number of nodes
  - N_elem : Number of elements
  - N_steps : Number of time steps
*/
{

  /* Simulation file */
  FILE * Sim_dat;

  /* Auxiliar variable for reading the lines in the files */
  char line[MAXC] = {0};
  
  /* Number of element in the line , just for check */
  int nwords;
  int nSimParameter;
  int nKindAnalysis;
  char * words[MAXW] = {NULL};
  char * SimParameter[MAXW] = {NULL};
  char * KindAnalysis[MAXW] = {NULL};
  int Element_i,Nodes_i;

  /* Initialize parser to read files */
  ParserDictionary Dict = InitParserDictionary();
  char * delim_spa = Dict.sep[6];
  char * delim_equ = Dict.sep[1];
  char * delim_perc = Dict.sep[7];

  printf("Begin of read data file : %s \n",Name_File);
  
  /* Open and check .dat file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    puts("Error during the lecture of .dat file");
    exit(0);
  }

  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the white space as separators */
    nwords = parse (words, line, delim_spa);
    if(nwords>=1){
      /* In the line, read the words */
      for(int i = 0; i<nwords ; i++){
	/* Use the equal (=) separator */
	nSimParameter = parse (SimParameter, words[i], delim_equ);
	if(nSimParameter > 1){

	  if ( strcmp(SimParameter[0],"G") == 0 ){
	    g = atof(SimParameter[1]);
	    printf("Set gravity to : %f \n",g);
	  }

	  if( strcmp(SimParameter[0],"KIND_ANALYSIS") == 0 ){
	    nKindAnalysis = parse (KindAnalysis, SimParameter[1], delim_perc);
	    if(nKindAnalysis == 3){
	      printf("Kind of analysis : \n");

	      if( strcmp(KindAnalysis[0],"FEM") == 0 ){
		printf("\t This is a finite element method simulation \n");
	      }
	      if( strcmp(KindAnalysis[1],"SIGMA_V") == 0 ){
		printf("\t The stress tensor and the velocity will be the analysis fields \n");
	      }
	      if( strcmp(KindAnalysis[2],"2STG\n") == 0 ){
		printf("\t The temporal discretization will be done with Two-step Taylor-Galerkin \n");
	      }
       
	    }
	  }

	  /* Time parameters */
	  if( strcmp(SimParameter[0],"TIME_STEP") == 0 ){
	    DeltaTimeStep = atof(SimParameter[1]);
	    printf("Set increment of time step to : %f \n",DeltaTimeStep);
	  }
	  if( strcmp(SimParameter[0],"NUM_STEP") == 0 ){
	    NumTimeStep = atoi(SimParameter[1]);
	    printf("Set number of time steps to : %i \n",NumTimeStep);
	  }

	  /* Names of files */
	  if( strcmp(SimParameter[0],"MESH_FILE") == 0 ){
	    MeshFileName = SimParameter[1];
	    printf("Set name of the mesh file : %s \n",MeshFileName);
	  }
	  if( strcmp(SimParameter[0],"COND_INIT") == 0 ){
	    InitCondFileName = SimParameter[1];
	    printf("Set name of the initial conditions file : %s \n",InitCondFileName);
	  }
	  if( strcmp(SimParameter[0],"BOUND_COND") == 0 ){
	    BounCondFileName = SimParameter[1];
	    printf("Set name of the boundary conditions file : %s \n",BounCondFileName);
	  }
	  
	} /* End if nSimParameter */     
      } /* End for nwords */
    } /* End if nwords */  
  } /* End while */   
      
  /* Close .dat file */
  printf("End of read : %s \n",Name_File);
  fclose(Sim_dat);
  

} /* void read_dat(char * Name_File) */



/***************************************************************************/
