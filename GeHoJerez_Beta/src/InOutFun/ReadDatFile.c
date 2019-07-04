#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

/***************************************************************************/


void Read_GeneralParameters(char * Name_File)
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

  /* Initial message */
  printf("************************************************* \n");
  printf(" \t * Begin of read general parameters in : \n\t %s !!! \n",Name_File);
  
  /* Open and check .dat file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    puts("Error during the lecture of .dat file");
    exit(0);
  }

  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){
    /* Read the line with the white space as separators */
    nwords = parse (words, line," \n\t");
    if(nwords>=1){
      /* In the line, read the words */
      for(int i = 0; i<nwords ; i++){
	/* Use the equal (=) separator */
	nSimParameter = parse (SimParameter, words[i], "=\n");
	if(nSimParameter > 1){
	  /***********************************************************************/
	  if( strcmp(SimParameter[0],"KIND_ANALYSIS") == 0 ){
	    nKindAnalysis = parse (KIND_ANALYSIS, SimParameter[1],"%\n");

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
  printf("\t * End of read data file !!! \n");
  fclose(Sim_dat);
  

} /* void read_dat(char * Name_File) */


/***************************************************************************/

Boundaries Set_FEM_BCC(char * Name_File, Mesh FEM_Mesh)
/*
  Read the boundary conditions file :
  Inputs
  - Name_file : Name of the file
  BCC_TOP FIELD#V DIR={0,1} CURVE#{curve.txt}
  BCC_BOTTOM FIELD#V DIR={0,1} CURVE#{curve.txt}
  BCC_RIGHT FIELD#V DIR={0,1} CURVE#{curve.txt}
  BCC_LEFT FIELD#V DIR={0,1} CURVE#{curve.txt}
  Note : Only read those lines with a BCC in the init 
*/
{

  /* Define output */
  Boundaries FEM_BCC;
  
  /* Pointer to the FEM boundary conditions file */
  FILE * File_BCC;

  /* Number of boundary conditions */
  int BCC_NUM = 0;
  
  /* Auxiliar structure with the load curve */
  Curve Curve_BCC;
  
  /* Auxiliar variable for reading the lines in the files */
  char line[MAXC] = {0};
  
  /* Number of element in the line , just for check */
  int nkwords,nparam;
  char * kwords[MAXW] = {NULL};
  char * param[MAXW] = {NULL};

  /* Read the field */
  char * FIELD;


  /* Asign the nu,ber of boundaries */
  FEM_BCC.NumBounds = 4;

  /* Generate boundaries for the domain */
  FEM_BCC.BCC_i = (Load *)Allocate_Array(FEM_BCC.NumBounds,sizeof(Load));

  /*
    - Mover aqui la subrutina para generar los constornos
    - BCC LABEL=TOP FIELD#V DIR={0,1} CURVE#{curve.txt}
    - Guardar LABEL en FEM_BCC.BCC_i.Info
  */

  /* Initial message */  
  printf("************************************************* \n");
  printf(" \t * Begin of read boundary conditions in : \n\t %s \n",Name_File);
  
  /* Open and check .bcc file */
  File_BCC = fopen(Name_File,"r");  
  if (File_BCC==NULL){
    puts("Error during the lecture of .bcc file");
    exit(0);
  }

  /* Read the file line by line */
  while( fgets(line, sizeof line, File_BCC) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");

    /* When the parser find the keyword BCC_NUM : start reading BCC */
    if(strcmp(kwords[0],"BCC_NUM") == 0 ){

      BCC_NUM = atoi(kwords[1]);

      for(int i = 0 ; i<BCC_NUM ; i++){

	/* Read as line as BCC */
	fgets(line, sizeof line, File_BCC);
	/* Read the line with the space as separators */
	nkwords = parse (kwords, line," \n\t");
	if(nkwords != 3){
	  puts("Error in Read_FEM_BCC() : Wrong format of BCC !!!");
	  exit(0);
	}

	/* Set BCC in those node in the left of the domain */
	if ( ( strcmp(kwords[0],"BCC_TOP") == 0 ) ||
	     ( strcmp(kwords[0],"BCC_BOTTOM") == 0 ) ||
	     ( strcmp(kwords[0],"BCC_RIGHT") == 0 ) ||
	     ( strcmp(kwords[0],"BCC_LEFT") == 0) ){


	  /* Read the field to impose the boundary condition */
	  nparam = parse (param,kwords[1],"=\n");
	  if(strcmp(param[0],"FIELD") == 0){
	    FIELD = param[1];
	  }
	  else{
	    puts("Error in Read_FEM_BCC() : Wrong format !!!");
	    exit(0);
	  }

	  /* Read the curve to impose the boundary condition */
	  nparam = parse (param,kwords[2],"={}\n");
	  if(nparam != 2){
	    puts("Error in Read_FEM_BCC() : Wrong format !!!");
	    exit(0);
	  }
	  if(strcmp(param[0],"CURVE") == 0){
	    Curve_BCC = ReadCurve(param[1]);
	  }
	  else{
	    puts("Error in Read_FEM_BCC() : Wrong format !!!");
	    exit(0);
	  }

	  /* Apply this boundary conditions and copy information of the BCC */
	  if(strcmp(kwords[0],"BCC_TOP") == 0){
	    FEM_Mesh->TOP.Value = Curve_BCC;
	    strcpy(FEM_Mesh->TOP.Info,FIELD);
	  }
	  else if(strcmp(kwords[0],"BCC_BOTTOM") == 0){
	    FEM_Mesh->BOTTOM.Value = Curve_BCC;
	    strcpy(FEM_Mesh->BOTTOM.Info,FIELD);
	  }
	  else if(strcmp(kwords[0],"BCC_RIGHT") == 0){
	    FEM_Mesh->RIGHT.Value = Curve_BCC;
	    strcpy(FEM_Mesh->RIGHT.Info,FIELD);
	  }
	  else if(strcmp(kwords[0],"BCC_LEFT") == 0){
	    FEM_Mesh->LEFT.Value = Curve_BCC;
	    strcpy(FEM_Mesh->LEFT.Info,FIELD);
	  }
      
	} /* End of read BC */
      }
      
    } /* End of if(strcmp(kwords[0],"BCC_NUM") == 0 ) */
    
  }  /* End of while */
  printf("\t * End of read boundary conditions file !!! \n");
  fclose(File_BCC);
  
} /* BoundayConditions ReadBCC(char * Name_File) */

/**********************************************************************/

LoadCase Read_MPM_LoadCase_ExtForces(char * Name_File,GaussPoint GP_Mesh)
/*

  Read loads over the GP (.bcc)

  - Load format examples :
  - - External forces
  F_LOAD_NUM number
  F_LOAD_GP CURVE#{curve.txt} NUM_NODES#integer
  .
  . integer (NLIST)
  .

*/
{

  /* Define new load case for the contact forces */
  LoadCase GP_Loads;

  /* Simulation file */
  FILE * Sim_dat;

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */
  int nparam;
  char * param[MAXW] = {NULL};

  /* Read the curve file */
  int aux_CURVE;
  char * CURVE[MAXW] = {NULL};

  /* Number of load cases (by default is zero) */
  GP_Loads.NumLoads = 0;

  /* Initial message */  
  printf("************************************************* \n");
  printf(" \t * Begin of read contact forces in : \n\t %s \n",Name_File);
  
  /* Open and check .bcc file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    printf("Error in Read_MPM_LoadCase_ExtForces() during the lecture of : \n\t %s",Name_File);
    exit(0);
  }

  
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");  

    if (strcmp(kwords[0],"F_LOAD_NUM") == 0 ){

      /* Read the number of loads of the load case */
      GP_Loads.NumLoads = atoi(kwords[1]);
      /* Allocate a table of loads in the case NumLoads > 0 */
      if(GP_Loads.NumLoads > 0)
	GP_Loads.Load_i =
	  (Load *)Allocate_Array(GP_Loads.NumLoads,sizeof(Load));

      /* Fill the table of loads */
      for(int i = 0 ; i<GP_Loads.NumLoads ; i++){

	/* Read each line and split it in individual words */
	fgets(line, sizeof line, Sim_dat);
	nkwords = parse (kwords, line," \n\t");
	if(nkwords != 4){
	  puts("Error in Read_MPM_LoadCase_ExtForces() : Wrong format !!!");
	  exit(0);
	}

	/* Read the curve associated to the load */
	aux_CURVE = parse (CURVE,kwords[2],"={}\n");
	if(aux_CURVE != 2){
	  puts("Error in Read_MPM_LoadCase_ExtForces() : Wrong format !!!");
	  exit(0);
	}
	if(strcmp(CURVE[0],"CURVE") == 0){
	  GP_Loads.Load_i[i].Value = ReadCurve(CURVE[1]);
	}
	else{
	  puts("Error in Read_MPM_LoadCase_ExtForces() : Wrong format !!!");
	  exit(0);
	}

	/* Read the number of nodes associated to this initial condition */
	nparam = parse (param,kwords[3],"=\n");
	if(strcmp(param[0],"NUM_NODES") == 0){
	  if(strcmp(param[1],"ALL_NODES") == 0){
	    GP_Loads.Load_i[i].Nodes = 
	      (int *)Allocate_Array(GP_Mesh.NumGP,sizeof(int));
	    for(int j = 0 ; j<GP_Mesh.NumGP ; j++){
	      GP_Loads.Load_i[i].Nodes[j] = j;
	    }
	    GP_Loads.Load_i[i].NumNodes = GP_Mesh.NumGP;
	  }
	  else{
	    /* Allocate array with the size */
	    GP_Loads.Load_i[i].Nodes = 
	      (int *)Allocate_Array(atoi(param[1]),sizeof(int));
	    GP_Loads.Load_i[i].NumNodes = atoi(param[1]);
   
	    /* Fill the list of nodes */
	    for(int j = 0 ; j<GP_Loads.Load_i[i].NumNodes ; j++){
	      fgets(line, sizeof line, Sim_dat);
	      nparam = parse(param,line," \n");
	      if(nparam == 1){
		GP_Loads.Load_i[i].Nodes[j] = atoi(param[0]);
	      }
	      else{
		puts("Error in Read_MPM_LoadCase_ExtForces() : Check the list of nodes ");
		exit(0);
	      }
	    }
	  }
	}
	else{
	  puts("Error in Read_MPM_LoadCase_ExtForces() : Wrong format !!!");
	  exit(0);
	}
      }

    }
  }

  /* Close .dat file */
  /* Final message */
  printf("\t * End of read data file !!! \n");
  fclose(Sim_dat);

  return GP_Loads;
}


/**********************************************************************/

LoadCase Read_MPM_LoadCase_BodyForces(char * Name_File,GaussPoint GP_Mesh)
/*
  - - Gravity load
  B_LOAD_NUM number
  B_LOAD_GP CURVE#{curve.txt} NUM_NODES#ALL_NODES
*/
{

  /* Define new load case for the contact forces */
  LoadCase GP_Loads;

  /* Simulation file */
  FILE * Sim_dat;

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */
  int nparam;
  char * param[MAXW] = {NULL};

  /* Number of initial conditions (by default is zero) */
  GP_Loads.NumLoads = 0;

  /* Read the curve file */
  int aux_CURVE;
  char * CURVE[MAXW] = {NULL};
 
  /* Initial message */  
  printf("************************************************* \n");
  printf(" \t * Begin of read body forces in : \n\t %s \n",Name_File);
  
  /* Open and check .bcc file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    printf("Error in Read_MPM_LoadCase_BodyForces() during the lecture of : \n\t %s",Name_File);
    exit(0);
  }

  
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");  

    if (strcmp(kwords[0],"B_LOAD_NUM") == 0 ){

      /* Read the number of loads of the load case */
      GP_Loads.NumLoads = atoi(kwords[1]);
      /* Allocate a table of loads in the case NumLoads > 0 */
      if(GP_Loads.NumLoads > 0)
	GP_Loads.Load_i =
	  (Load *)Allocate_Array(GP_Loads.NumLoads,sizeof(Load));

      /* Fill the table of loads */
      for(int i = 0 ; i<GP_Loads.NumLoads ; i++){

	/* Read each line and split it in individual words */
	fgets(line, sizeof line, Sim_dat);
	nkwords = parse (kwords, line," \n\t");
	if(nkwords != 3){
	  puts("Error in Read_MPM_LoadCase_BodyForces() : Wrong format !!!");
	  exit(0);
	}

	/* Read the curve associated to the load */
	aux_CURVE = parse (CURVE,kwords[1],"={}\n");
	if(aux_CURVE != 2){
	  puts("Error in Read_MPM_LoadCase_BodyForces() : Wrong format !!!");
	  exit(0);	  
	}
	if(strcmp(CURVE[0],"CURVE") == 0){
	  GP_Loads.Load_i[i].Value = ReadCurve(CURVE[1]);
	}
	else{
	  puts("Error in Read_MPM_LoadCase_BodyForces() : Wrong format !!!");
	  exit(0);
	}

	/* Read the number of nodes associated to this initial condition */
	nparam = parse (param,kwords[2],"=\n");
	if(strcmp(param[0],"NUM_NODES") == 0){
	  if(strcmp(param[1],"ALL_NODES") == 0){
	    GP_Loads.Load_i[i].Nodes = 
	      (int *)Allocate_Array(GP_Mesh.NumGP,sizeof(int));
	    for(int j = 0 ; j<GP_Mesh.NumGP ; j++){
	      GP_Loads.Load_i[i].Nodes[j] = j;
	    }
	    GP_Loads.Load_i[i].NumNodes = GP_Mesh.NumGP;
	  }
	  else{
	    /* Allocate array with the size */
	    GP_Loads.Load_i[i].Nodes = 
	      (int *)Allocate_Array(atoi(param[1]),sizeof(int));
	    GP_Loads.Load_i[i].NumNodes = atoi(param[1]);
   
	    /* Fill the list of nodes */
	    for(int j = 0 ; j<GP_Loads.Load_i[i].NumNodes ; j++){
	      fgets(line, sizeof line, Sim_dat);
	      nparam = parse(param,line," \n");
	      if(nparam == 1){
		GP_Loads.Load_i[i].Nodes[j] = atoi(param[0]);
	      }
	      else{
		puts("Error in Read_MPM_LoadCase_BodyForces() : Check the list of nodes ");
		exit(0);
	      }
	    }
	  }
	}
	else{
	  puts("Error in Read_MPM_LoadCase_BodyForces() : Wrong format !!!");
	  exit(0);
	}
      }

    }
  }

  /* Close .dat file */
  /* Final message */
  printf("\t * End of read data file !!! \n");
  fclose(Sim_dat);

  return GP_Loads;  
}

/**********************************************************************/

void Read_MPM_InitVal(char * Name_File, GaussPoint GP_Mesh)
/*
  If the first word is INIT_GP read it : set an initial value
  
  - Initial values format :
  - - Initial velocities
  INIT_GP FIELD#V CURVE#{curve.txt} NUM_NODES#integer
  .
  . integer (NLIST)
  .
  - - ALL_NODES
  INIT_GP FIELD#V CURVE#{curve.txt} NUM_NODES#ALL_NODES
*/
{

  /* Simulation file */
  FILE * Sim_dat;

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */
  int nparam;
  char * param[MAXW] = {NULL};

  /* Number of initial conditions (by default is zero) */
  int INIT_NUM = 0;
  
  /* Special variables for the initial conditions parser */
  char FIELD[MAXC] = {0}; /* Name of the field of the initial condition */
  int DIM_FIELD; /* Number of dimensions of the field */
  int AUX_VAL; /* Output of the parser for the field value */
  char * READ_VAL[MAXW] = {NULL}; /* Variable to store the parser */
  double * FIELD_VAL; /* Value of the field */
  
  int NUM_NODES; /* Number of nodes whis this initial condition */
  int AUX_NODES; /* Output of the parser for the nodes value */
  char * READ_NODES[MAXW] = {NULL}; /* Variable to store the parser */

  /* Initial message */  
  printf("************************************************* \n");
  printf(" \t * Begin of read initial conditions in : \n\t %s \n",Name_File);
  
  /* Open and check .bcc file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    printf("Error in Read_MPM_InitVal() during the lecture of : %s",Name_File);
    exit(0);
  }

  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");

    /* When the parser find the keyword INIT_NUM : start reading initial conditions */
    if(strcmp(kwords[0],"INIT_NUM") == 0 ){

      /* Get the number of initial conditions in the problem */
      INIT_NUM = atoi(kwords[1]);

      /* Loop over the initial conditions */
      for(int i = 0 ; i<INIT_NUM ; i++){

	/* Read as line as initial conditions */
	fgets(line, sizeof line, Sim_dat);
	/* Read the line with the space as separators */
	nkwords = parse (kwords, line," \n\t");
	if(nkwords != 4){
	  puts("Error in Read_MPM_InitVal() : Wrong format !!!");
	  exit(0);
	}
      
	if (strcmp(kwords[0],"INIT_GP") == 0 ){

	  /* Read the field to impose the initial condition */
	  nparam = parse (param,kwords[1],"=\n");
	  if(strcmp(param[0],"FIELD") == 0){
	    strcpy(FIELD,param[1]);
	    if (strcmp(FIELD,"V") == 0){
	      DIM_FIELD = GP_Mesh.Phi.vel.N_cols;
	      FIELD_VAL = (double *)Allocate_Array(DIM_FIELD,sizeof(double));
	    }
	  }
	  else{
	    puts("Error in Read_MPM_InitVal() : Wrong format !!!");
	    exit(0);
	  }
      
	  /* Read the curve associated to the initial condition */
	  nparam = parse (param,kwords[2],"=\n");
	  if(strcmp(param[0],"VALUE") == 0){
	    /* Read file of the curve */
	    AUX_VAL = parse(READ_VAL,param[1],"{,}\n");
	    if(AUX_VAL == DIM_FIELD){
	      for(int j = 0 ; j<DIM_FIELD ; j++){
		FIELD_VAL[j] = atof(READ_VAL[j]);
	      }
	    }
	    else{
	      printf("Error in Read_MPM_InitVal() : Wrong number of dimensions for VALUE={,}");
	      exit(0);
	    }
	  }
	  else{
	    puts("Error in Read_MPM_InitVal() : Wrong format !!!");
	    exit(0);
	  }

	  /* Read the number of nodes associated to this initial condition */
	  nparam = parse (param,kwords[3],"=\n");
	  if(nparam != 2){
	    puts("Error in Read_MPM_InitVal() : Wrong format !!!");
	    exit(0);
	  }
	  if(strcmp(param[0],"NUM_NODES") == 0){
	    /* Read the number of nodes associated to this initial condition */
	    if(strcmp(param[1],"ALL_NODES") == 0){
	      for(int j = 0 ; j<GP_Mesh.NumGP ; j++){
		for(int k = 0 ; k<DIM_FIELD ; k++){
		  if (strcmp(FIELD,"V") == 0)
		    GP_Mesh.Phi.vel.nM[j][k] = FIELD_VAL[k];
		}		
	      }
	    }
	    else{
	      NUM_NODES = atoi(param[1]);
	      for(int j = 0 ; j<NUM_NODES ; j++){		
		fgets(line, sizeof line, Sim_dat);
		AUX_NODES = parse(READ_NODES,line," \n");
		if(AUX_NODES == 1){
		  for(int k = 0 ; k<DIM_FIELD ; k++){
		    if (strcmp(FIELD,"V") == 0)
		      GP_Mesh.Phi.vel.nM[atoi(READ_NODES[0])][k] = FIELD_VAL[k];
		  }
		  
		}
		else{
		  puts("Error in ReadLoads_GP() : Check the list of nodes ");
		  exit(0);
		}
	      }
	    }
	  }
	  else{
	    puts("Error in Read_MPM_InitVal() : Wrong format of INIT_GP !!!");
	    exit(0);
	  }

	  /* Asign initial conditions*/

      
	}
      }
      
    } /* End if(strcmp(kwords[0],"INIT_NUM") == 0 ) */
    
  } /* end While */
  
    /* Close .dat file */
  /* Final message */
  printf("\t * End of read data file !!! \n");
  fclose(Sim_dat);
  
}



/**********************************************************************/

/* void BCC_GP_Forces(GaussPoint MeshGP, Load * LoadsGP, int NumLoadsGP, int TimeStep) */
/* /\* */
/*   Forces defined in the Gauss Points : */
/*   Inputs */
/* *\/ */
/* { */
/*   /\* 0º Loop over the loads *\/ */
/*   for(int i = 0 ; i<NumLoadsGP ; i++){ */
  
/*     /\* 1º  Check the time step *\/ */
/*     if( (TimeStep < 0) || */
/* 	(TimeStep > LoadsGP[i].Value.Num)){ */
/*       puts("Error in BCC_GP_Forces() : The time step is out of the curve !!"); */
/*       exit(0); */
/*     } */
  
/*     /\* 2º Fill the matrix with the local forces *\/ */
/*     for(int j = 0 ; j<Loads[i].NumNodes ; j++){ */
      
/*       /\* 3º Check if this GP has a force applied *\/ */
/*       if( (Loads[i].Nodes[j] > GP_Mesh.NumGP) || */
/* 	  (Loads[i].Nodes[j] < 0)){ */
/* 	puts("Error in BCC_GP_Forces() : This GP does not exist !!"); */
/* 	exit(0); */
/*       } */
      
/*       /\* 4º Loop over the dimensions *\/ */
/*       for(int k = 0 ; k<NumberDimensions ; k++){ */
/* 	/\* 5º Apply the force in the node *\/ */
/* 	GP_Mesh.Phi.F.nM[Loads[i].Nodes[j]][k] += */
/* 	  Loads[i].Value.Fx[TimeStep][k]; */
/*       } */
      
/*     } */

/*   } */

/* } */

/**********************************************************************/
