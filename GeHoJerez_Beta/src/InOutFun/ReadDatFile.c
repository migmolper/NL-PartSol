
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
  char * KIND_ANALYSIS[MAXW] = {NULL};

  /* Initialize parser to read files */
  ParserDictionary Dict = InitParserDictionary();
  char * delim_spa = Dict.sep[6];
  char * delim_equ = Dict.sep[1];
  char * delim_perc = Dict.sep[7];

  /* Auxiliar variables for G */
  int AUX_G_DAT;
  char * G_DAT[MAXW] = {NULL};

  /* Initial message */
  printf("************************************************* \n");
  printf("Begin of read data file : %s !!! \n",Name_File);
  
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
	  /***********************************************************************/
	  if( strcmp(SimParameter[0],"KIND_ANALYSIS") == 0 ){
	    nKindAnalysis = parse (KIND_ANALYSIS, SimParameter[1], delim_perc);
	    if(nKindAnalysis == 4){
	      printf(" * Kind of analysis : \n");
	      /* First parameter of KIND_ANALYSIS : FEM/MPM */
	      if( strcmp(KIND_ANALYSIS[0],"FEM") == 0 ){
		KindAnalysis = KIND_ANALYSIS[0];
		printf("\t -> %s : Finite Element Method \n",KindAnalysis); 
	      }
	      if(strcmp(KIND_ANALYSIS[0],"MPM") == 0){
		KindAnalysis = KIND_ANALYSIS[0];
		printf("\t -> %s : Material Point Method \n",KindAnalysis); 
	      }
	      /*************************************************************/
	      /* Second parameter of KIND_ANALYSIS : U/SIGMA_V */
	      if( strcmp(KIND_ANALYSIS[1],"SIGMA_V") == 0 ){
		FieldsAnalysis = KIND_ANALYSIS[1];
		printf("\t -> %s : Stress-Velocity \n",FieldsAnalysis); 
	      }
	      if( strcmp(KIND_ANALYSIS[1],"U") == 0 ){
		FieldsAnalysis = KIND_ANALYSIS[1];
		printf("\t -> %s : Displacement \n",FieldsAnalysis);
	      }
	      /*************************************************************/
	      /* Third parameter of KIND_ANALYSIS : 1D/2D/3D */
	      if( strcmp(KIND_ANALYSIS[2],"1D") == 0 ){
		NumberDimensions = 1;
		printf("\t -> 1D \n");
	      }
	      if( strcmp(KIND_ANALYSIS[2],"2D") == 0 ){
		NumberDimensions = 2;
		printf("\t -> 2D \n");
	      }
	      if( strcmp(KIND_ANALYSIS[2],"3D") == 0 ){
		NumberDimensions = 3;
		printf("\t -> 3D \n");
	      }
	      /*************************************************************/
	      /* Third parameter of KIND_ANALYSIS : 2STG/VerletLF */
	      if( strcmp(KIND_ANALYSIS[3],"2STG") == 0 ){
		TimeIntegration = KIND_ANALYSIS[3];
		printf("\t -> %s : Two-Steps Taylor Galerkin \n",TimeIntegration); 
	      }
	      if( strcmp(KIND_ANALYSIS[3],"VerletLF") == 0 ){
		TimeIntegration = KIND_ANALYSIS[3];
		printf("\t -> %s : Leapfrog Verlet \n",TimeIntegration); 
	      }
	      /*************************************************************/
	    }
	    else{
	      printf("Error in ReadDatFile() : KIND_ANALYSIS !!! \n");
	      exit(0);
	    }
	  }
	  /***********************************************************************/
	  if ( strcmp(SimParameter[0],"G") == 0 ){
	    g = MatAlloc(NumberDimensions,1);
	    AUX_G_DAT = parse(G_DAT,SimParameter[1],"{,}\n");
	    if(AUX_G_DAT != NumberDimensions){
	      puts("Error in ReadDatFile() : Wrong number of dimensions of gravity");
	      exit(0);
	    }
	    printf(" * Set gravity to : {");
	    for(int j = 0; j<NumberDimensions ; j++){
	      g.nV[j] = atof(G_DAT[j]);
	      printf(" %f ",g.nV[j]);
	    }
	    printf("}\n");
	    
	  }
	  /***********************************************************************/
	  if( strcmp(SimParameter[0],"DENSITY") == 0 ){
	    Density = atof(SimParameter[1]);
	    printf(" * Set the value of the density : %f \n",Density);
	  }
	  /***********************************************************************/
	  if( strcmp(SimParameter[0],"ELASTIC_MODULUS") == 0 ){
	    ElasticModulus = atof(SimParameter[1]);
	    printf(" * Set the value of E : %f \n",ElasticModulus);
	  }
	  /***********************************************************************/
	  if( strcmp(SimParameter[0],"POISSON_MODULUS") == 0 ){
	    PoissonModulus = atof(SimParameter[1]);
	    printf(" * Set the value of mu : %f \n",PoissonModulus);
	  }
	  /***********************************************************************/
	  if( strcmp(SimParameter[0],"TIME_STEP") == 0 ){
	    DeltaTimeStep = atof(SimParameter[1]);
	    printf(" * Set increment of time step to : %f \n"
		   ,DeltaTimeStep);
	  }
	  /***********************************************************************/
	  if( strcmp(SimParameter[0],"NUM_STEP") == 0 ){
	    NumTimeStep = atoi(SimParameter[1]);
	    printf(" * Set number of time steps to : %i \n",
		   NumTimeStep);
	  }
	  /***********************************************************************/
	  if( strcmp(SimParameter[0],"FEM_FILE_NAME") == 0 ){
	    FEM_MeshFileName = SimParameter[1];
	    printf(" * Set name of the mesh file : \n");
	    printf("\t -> %s \n",FEM_MeshFileName);
	  }
	  if( strcmp(SimParameter[0],"MPM_FILE_NAME") == 0 ){
	    MPM_MeshFileName = SimParameter[1];
	    printf(" * Set name of the mesh file : \n");
	    printf("\t -> %s \n",MPM_MeshFileName);
	  }
	  /***********************************************************************/
	  if( strcmp(SimParameter[0],"COND_INIT") == 0 ){
	    InitCondFileName = SimParameter[1];
	    printf(" * Set name of the initial conditions file :  \n");
	    printf("\t -> %s \n",InitCondFileName);
	  }
	  /***********************************************************************/
	  if( strcmp(SimParameter[0],"BOUND_COND") == 0 ){
	    BounCondFileName = SimParameter[1];
	    printf(" * Set name of the boundary conditions file : \n");
	    printf("\t -> %s \n",BounCondFileName); 
	  }
	  /***********************************************************************/
	  if( strcmp(SimParameter[0],"OUTPUT_DIR") == 0 ){
	    OutputDir = SimParameter[1];
	    printf(" * Set route for the outputs  : \n");
	    printf("\t -> %s \n",OutputDir); 
	  }
	  /***********************************************************************/
	} /* End if nSimParameter */     
      } /* End for nwords */
    } /* End if nwords */
  } /* End while */

  /* Set the number of DOFs for each node depending of the kind of analysis */
  if( strcmp(FieldsAnalysis,"U") == 0 ){
    if(NumberDimensions == 1)
      NumberDOF = 1;
    if(NumberDimensions == 2)
      NumberDOF = 2;
    if(NumberDimensions == 3)
      NumberDOF = 3;
  }
  if( strcmp(FieldsAnalysis,"SIGMA_V") == 0 ){
    if(NumberDimensions == 1)
      NumberDOF = 2;
    if(NumberDimensions == 2)
      NumberDOF = 5;
    if(NumberDimensions == 3)
      NumberDOF = 9;
  }

  /* Close .dat file */
  /* Final message */
  printf("End of read data file !!! \n");
  fclose(Sim_dat);
  

} /* void read_dat(char * Name_File) */



/***************************************************************************/
