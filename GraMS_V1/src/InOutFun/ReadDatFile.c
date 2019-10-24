#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../MathTools/MathTools.h"
#include "InOutFun.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))


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
	  if( strcmp(SimParameter[0],"RESULT_STEP") == 0 ){
	    ResultsTimeStep = atoi(SimParameter[1]);
	    printf(" * Set time steps to generate outputs : %i \n",
		   ResultsTimeStep);
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
  BCC LABEL=TOP FIELD=V DIR={0,1} CURVE={Curve.txt}
  BCC LABEL=BOTTOM FIELD=V DIR={0,1} CURVE={Curve.txt}
  BCC LABEL=RIGHT FIELD=V DIR={1,0} CURVE={Curve.txt}
  BCC LABEL=LEFT FIELD=V DIR={1,0} CURVE={Curve.txt}
  
  Note : Only read those lines with a BCC in the init 
*/
{

  /* Auxiliar variable with the number of nodes in a element */
  int NumNodesElem;
  
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

  /* Auxiliar table to store the boundaries */
  int * NodesBound_aux;
  int * CounterNodesBound;
  int NumNodesBound;

  /* Count the number of elements that share this node */
  ChainPtr Elem_Conn; /* Loop over the connectivity chain */
  int Repeat_Nod; /* Counter */
  /* Variables that fills the boundaries nodes */
  int aux_RIGHT = 0; 
  int aux_TOP = 0;
  int aux_LEFT = 0;
  int aux_BOTTOM = 0;
  /* Set to zero X and Y min and max values of the mesh */
  double MAX_X = 0;
  double MAX_Y = 0;
  double MIN_X = 0;
  double MIN_Y = 0;
  /* Boundaries labels */
  char * BoundLabels [4] = {"BOTTOM", "RIGHT", "TOP", "LEFT"};
  /* Index of the boundary to apply the BCC from (0 -> N_bounds -1) */
  int IndexBoundary = -999;

  /* Kind of domain */
  strcpy(FEM_BCC.Info,"SQUARE");      
  
  /* Asign the number of boundaries for a square domain */
  FEM_BCC.NumBounds = 4;

  /* Generate boundaries for the domain */
  FEM_BCC.BCC_i = (Load *)Allocate_Array(FEM_BCC.NumBounds,sizeof(Load));

  /* Fill some parameters */
  for(int i = 0 ; i<FEM_BCC.NumBounds  ; i++){
    /* Set to zero the number of nodes in the boundary of the mesh */
    FEM_BCC.BCC_i[i].NumNodes = 0;
    /* Set the labels for a square domain */
    strcpy(FEM_BCC.BCC_i[i].Info,BoundLabels[i]);
  }

 
  /*
    - BCC LABEL=TOP FIELD#V DIR={0,1} CURVE#{curve.txt}
  */
    
  /* Find the nodes in the boundary */
  switch(NumberDimensions){
    
  case 1: /******************** 1D mesh ********************/
    /*
     *-----*
     0     1
     */
    /* In a 1D mesh we only have two nodes in the boundary */
    FEM_BCC.BCC_i[0].NumNodes = 1;
    FEM_BCC.BCC_i[1].NumNodes = 1;
    /* Allocate the size of the array with the nodes */
    FEM_BCC.BCC_i[0].Nodes = (int *)Allocate_ArrayZ(1,sizeof(int));
    FEM_BCC.BCC_i[1].Nodes = (int *)Allocate_ArrayZ(1,sizeof(int));
    /* Set the nodes of the boundaries */
    FEM_BCC.BCC_i[0].Nodes[0] = 0;
    FEM_BCC.BCC_i[0].Nodes[0] = 1-FEM_Mesh.NumNodesMesh;
    
    break; /******************** 2D mesh ********************/
    
  case 2: /******************** 2D mesh ********************/
    /*  2
     *-----*
     |     |
     |3    |1
     |     | 
     *-----*
        0
    */

    /* FIX THIS ---> TEMPORARY */
    NumNodesElem = FEM_Mesh.NumNodesElem[0];
    
    if(NumNodesElem == 4){ /* Quadrilateral elements */
      /* 0º Allocate an array of zeros to assign a 1 to those nodes in the boundary */
      NodesBound_aux = (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh,sizeof(int));
      CounterNodesBound = (int *)Allocate_ArrayZ(4,sizeof(int));
      /* 1º Set to zero the number of nodes in the boundary */
      NumNodesBound = 0;
      /* 2º Iterate over the nodes to fin the nodes in the boundary */
      for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){

  	/* 3º Get the max values of the boundary */
  	MAX_X = MAXVAL(MAX_X,FEM_Mesh.Coordinates.nM[i][0]);
  	MAX_Y = MAXVAL(MAX_Y,FEM_Mesh.Coordinates.nM[i][1]);
  	MIN_X = MINVAL(MIN_X,FEM_Mesh.Coordinates.nM[i][0]);
  	MIN_Y = MINVAL(MIN_Y,FEM_Mesh.Coordinates.nM[i][1]);
	
  	/* 4º Set the counter to zero */
  	Repeat_Nod = 0;
	
  	/* 5º Loop over the connectivity mesh */
	for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){
	  Elem_Conn = FEM_Mesh.Connectivity[j];
	  while(Elem_Conn != NULL){
	    if((Elem_Conn->I) == i){
	      Repeat_Nod++;
	    }
	    Elem_Conn = Elem_Conn->next;
	  }
	}
	
  	/* 6º Add this element to the boundary */
  	if (Repeat_Nod < 4){
  	  NodesBound_aux[i] = 1;
  	  NumNodesBound++;
  	}
      }
      
      /* 7º Count the number of nodes in each boundarie */
      for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
  	if(NodesBound_aux[i] == 1){

	  if(FEM_Mesh.Coordinates.nM[i][1] == MIN_Y){
	    /* Count number of nodes in the bottom */
	    CounterNodesBound[0]++;
  	  }
  	  if(FEM_Mesh.Coordinates.nM[i][0] == MAX_X){
	    /* Count number of nodes in the right */
	    CounterNodesBound[1]++;
  	  }
  	  if(FEM_Mesh.Coordinates.nM[i][1] == MAX_Y){
	    /* Count the number of nodes in the top */
	    CounterNodesBound[2]++;
  	  }
  	  if(FEM_Mesh.Coordinates.nM[i][0] == MIN_X){
	    /* Count the number of nodes in the left */
	    CounterNodesBound[3]++;
  	  }
	    
  	}
      }
      
      /* Allocate the arrays with the boundary nodes */
      for(int i = 0 ; i<FEM_BCC.NumBounds ; i++){
	FEM_BCC.BCC_i[i].NumNodes = CounterNodesBound[i];
	FEM_BCC.BCC_i[i].Nodes =
	  (int *)Allocate_ArrayZ(CounterNodesBound[i],sizeof(int));
      }

      /* Free data */
      free(CounterNodesBound);
      

      /* Fill the arrays  */
      for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
  	if(NodesBound_aux[i] == 1){

	  if(FEM_Mesh.Coordinates.nM[i][1] == MIN_Y){
  	    FEM_BCC.BCC_i[0].Nodes[aux_BOTTOM] = i;
  	    aux_BOTTOM++;
  	  }
  	  if(FEM_Mesh.Coordinates.nM[i][0] == MAX_X){
  	    FEM_BCC.BCC_i[1].Nodes[aux_RIGHT] = i;
  	    aux_RIGHT++;
  	  }
  	  if(FEM_Mesh.Coordinates.nM[i][1] == MAX_Y){
  	    FEM_BCC.BCC_i[2].Nodes[aux_TOP] = i;
  	    aux_TOP++;
  	  }
  	  if(FEM_Mesh.Coordinates.nM[i][0] == MIN_X){
  	    FEM_BCC.BCC_i[3].Nodes[aux_LEFT] = i;
  	    aux_LEFT++;
  	  }

	    
  	}
      }

      
    } /* Quadrilateral elements */
    if(NumNodesElem == 3){ /* Triangular elements */
      puts("Error in Set_FEM_BCC() : Boundary nodes localization for T3 not implemented yet !");
    } /* Triangular elements */
    
    /* Free data */ 
    free(NodesBound_aux);

    break; /******************** 2D mesh ********************/
    
  case 3: /******************** 3D mesh ********************/
    printf("************************************************* \n");
    puts("Error in Set_FEM_BCC() : 3D cases not implemented yet !");
    printf("************************************************* \n");
    exit(0);
    break; /******************** 2D mesh ********************/
    
  default :
    printf("************************************************* \n");
    puts("Error in Set_FEM_BCC() : Wrong number of dimensions !");
    printf("************************************************* \n");
    exit(0);
  }


  /* Initial message */  
  printf(" \t -> Begin of read boundary conditions in : \n\t %s \n",Name_File);
  
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

	/* Set BCC in those node in the left of the domain */
	if( (nkwords == 5) && ( strcmp(kwords[0],"BCC") == 0)){

	  /* Read the label to impose the boundary condition */
	  nparam = parse (param,kwords[1],"=\n");
	  if(strcmp(param[0],"LABEL") == 0){
	    for(int j = 0; j<FEM_BCC.NumBounds ; j++){
	      if(strcmp(param[1],BoundLabels[j]) == 0)
		IndexBoundary = j;
	    }
	  }
	  else{
	    puts("Error in Read_FEM_BCC() : Wrong format !!!");
	    exit(0);
	  }

	
	  /* Read the field to impose the boundary condition */
	  nparam = parse (param,kwords[2],"=\n");
	  if( (strcmp(param[0],"FIELD") == 0 ) && (nparam == 2) ){
	    strcpy(FEM_BCC.BCC_i[IndexBoundary].Info,param[1]);
	  }
	  else{
	    puts("Error in Set_FEM_BCC() : Wrong format !!!");
	    exit(0);
	  }

	  /* Read the dirrection to impose the boundary condition */
	  nparam = parse (param,kwords[3],"={,}\n");
	  if( (strcmp(param[0],"DIR") == 0 ) && (nparam >= 2) ){
	    /* Number of dimensions of the BCC */
	    FEM_BCC.BCC_i[IndexBoundary].Dim = nparam - 1;
	    /* Direction of the BCC */
	    FEM_BCC.BCC_i[IndexBoundary].Dir =
	      (int *)Allocate_Array(FEM_BCC.BCC_i[IndexBoundary].Dim,sizeof(int));
	    /* Fill the direction of the BCC */
	    for(int j = 0 ; j<FEM_BCC.BCC_i[IndexBoundary].Dim ; j++){
	      FEM_BCC.BCC_i[IndexBoundary].Dir[j] = atoi(param[j+1]);
	    }	    
	  }
	  else{
	    puts("Error in Set_FEM_BCC() : Wrong format !!!");
	    exit(0);
	  }

	  /* Read the curve to impose the boundary condition */
	  nparam = parse (param,kwords[4],"={,}\n");
	  if( (strcmp(param[0],"CURVE") == 0) || (nparam >= 2) ){

	    /* Alocate the table of curves */
	    FEM_BCC.BCC_i[IndexBoundary].Value =
	      (Curve *)Allocate_Array(FEM_BCC.BCC_i[IndexBoundary].Dim,sizeof(Curve));
	    /* Fill the curve table */
	    for(int j = 0 ; j<FEM_BCC.BCC_i[IndexBoundary].Dim ; j++){
	      if(strcmp(param[j+1],"NULL") != 0)
		FEM_BCC.BCC_i[IndexBoundary].Value[j] = ReadCurve(param[j+1]);
	    }    
	  }
	  else{
	    puts("Error in Set_FEM_BCC() : Wrong format !!!");
	    exit(0);
	  }
      
	} /* End of read BC */
	else{
	  puts("Error in Set_FEM_BCC() : Wrong format of BCC !!!");
	  exit(0);
	}
	
      }
      
    } /* End of if(strcmp(kwords[0],"BCC_NUM") == 0 ) */
    
  }  /* End of while */
  printf("\t -> End of read boundary conditions file !!! \n");
  fclose(File_BCC);

  return FEM_BCC;
  
} /* BoundayConditions ReadBCC(char * Name_File) */

/**********************************************************************/

LoadCase Read_MPM_LoadCase_ExtForces(char * Name_File,GaussPoint GP_Mesh)
/*

  Read external loads over the GP (.bcc)

  - Load format examples :
  - - External forces
  F_LOAD_NUM number
  F_LOAD_GP DIR={int,int} CURVE={curve.txt,curve.txt} NUM_NODES#integer
  .
  . integer (NLIST)
  .
  Cases : 
  DIR={0,0} CURVE={NULL,NULL}
  DIR={0,1} CURVE={NULL,curve.txt}
  DIR={1,0} CURVE={curve.txt,NULL}
  DIR={1,1} CURVE={curve.txt,curve.txt}
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

  /* Read the direction */
  int aux_DIR;
  char * DIR[MAXW] = {NULL};
  
  /* Read the curve file */
  int aux_CURVE;
  char * CURVE[MAXW] = {NULL};

  /* Number of load cases (by default is zero) */
  GP_Loads.NumLoads = 0;

  /* Initial message */  
  printf("************************************************* \n");
  printf(" \t * Begin of read contact forces in : \n\t %s \n",
	 Name_File);
  
  /* Open and check .bcc file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    printf("Error in Read_MPM_LoadCase_ExtForces() during the lecture of : \n\t %s",
	   Name_File);
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
	if( (nkwords != 4) || (strcmp(kwords[0],"F_LOAD_GP") != 0) ){
	  puts("Error in Read_MPM_LoadCase_ExtForces() : Wrong format !!!");
	  exit(0);
	}

	/* Read the curve associated to the load */
	aux_DIR = parse (DIR,kwords[1],"={,}\n");
	if( (aux_DIR >= 2 ) && (strcmp(DIR[0],"DIR") == 0)){
	  /* Number of components of the load */
	  GP_Loads.Load_i[i].Dim = aux_DIR - 1;
	  /* Alocate the direction of the component */
	  GP_Loads.Load_i[i].Dir =
	    (int *)Allocate_Array(GP_Loads.Load_i[i].Dim,sizeof(int));
	  /* Fill the direction of the load */
	  for(int j = 0 ; j<GP_Loads.Load_i[i].Dim ; j++){
	    GP_Loads.Load_i[i].Dir[j] = atoi(DIR[j+1]);
	  }
	}
	else{
	  puts("Error in Read_MPM_LoadCase_ExtForces() : Wrong format DIR={,} !!!");
	  exit(0);
	}

	/* Read the curve associated to the load */
	aux_CURVE = parse (CURVE,kwords[2],"={,}\n");
	if( (aux_CURVE == aux_DIR) && (strcmp(CURVE[0],"CURVE") == 0) ){
	  /* Alocate the table of curves */
	  GP_Loads.Load_i[i].Value =
	    (Curve *)Allocate_Array(GP_Loads.Load_i[i].Dim,sizeof(Curve));
	  /* Fill the curve table */
	  for(int j = 0 ; j<GP_Loads.Load_i[i].Dim ; j++){
	    if(strcmp(CURVE[j+1],"NULL") != 0)
	      GP_Loads.Load_i[i].Value[j] = ReadCurve(CURVE[j+1]);
	  }
	}
	else{
	  puts("Error in Read_MPM_LoadCase_ExtForces() : Wrong format CURVE={,} !!!");
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
	  puts("Error in Read_MPM_LoadCase_ExtForces() : Wrong format NUM_NODES !!!");
	  exit(0);
	}
      }

    }
  }

  /* Close .dat file */
  /* Final message */
  printf("\t -> End of read data file !!! \n");
  fclose(Sim_dat);

  return GP_Loads;
}


/**********************************************************************/

LoadCase Read_MPM_LoadCase_BodyForces(char * Name_File,GaussPoint GP_Mesh)
/*
  Read body loads over the GP (.bcc)

  - Load format examples :
  - - Body forces
  B_LOAD_NUM number
  B_LOAD_GP DIR={int,int} CURVE={curve.txt,curve.txt} NUM_NODES=integer
  .
  . integer (NLIST)
  .
  Cases : 
  DIR={0,0} CURVE={NULL,NULL}
  DIR={0,1} CURVE={NULL,curve.txt}
  DIR={1,0} CURVE={curve.txt,NULL}
  DIR={1,1} CURVE={curve.txt,curve.txt}
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

  /* Read the direction */
  int aux_DIR;
  char * DIR[MAXW] = {NULL};
  
  /* Read the curve file */
  int aux_CURVE;
  char * CURVE[MAXW] = {NULL};

  /* Number of load cases (by default is zero) */
  GP_Loads.NumLoads = 0;

  /* Initial message */  
  printf("************************************************* \n");
  printf(" \t * Begin of read body forces in : \n\t %s \n",
	 Name_File);
  
  /* Open and check .bcc file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    printf("Error in Read_MPM_LoadCase_BodyForces() during the lecture of : \n\t %s",
	   Name_File);
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
	if( (nkwords != 4) || (strcmp(kwords[0],"B_LOAD_GP") != 0) ){
	  puts("Error in Read_MPM_LoadCase_BodyForces() : Wrong format !!!");
	  exit(0);
	}

	/* Read the curve associated to the load */
	aux_DIR = parse (DIR,kwords[1],"={,}\n");
	if( (aux_DIR >= 2 ) && (strcmp(DIR[0],"DIR") == 0)){
	  /* Number of components of the load */
	  GP_Loads.Load_i[i].Dim = aux_DIR - 1;
	  /* Alocate the direction of the component */
	  GP_Loads.Load_i[i].Dir =
	    (int *)Allocate_Array(GP_Loads.Load_i[i].Dim,sizeof(int));
	  /* Fill the direction of the load */
	  for(int j = 0 ; j<GP_Loads.Load_i[i].Dim ; j++){
	    GP_Loads.Load_i[i].Dir[j] = atoi(DIR[j+1]);
	  }
	}
	else{
	  puts("Error in Read_MPM_LoadCase_BodyForces() : Wrong format direction !!!");
	  exit(0);
	}

	/* Read the curve associated to the load */
	aux_CURVE = parse (CURVE,kwords[2],"={,}\n");
	if( (aux_CURVE == aux_DIR) && (strcmp(CURVE[0],"CURVE") == 0) ){
	  /* Alocate the table of curves */
	  GP_Loads.Load_i[i].Value =
	    (Curve *)Allocate_Array(GP_Loads.Load_i[i].Dim,sizeof(Curve));
	  /* Fill the curve table */
	  for(int j = 0 ; j<GP_Loads.Load_i[i].Dim ; j++){
	    if(strcmp(CURVE[j+1],"NULL") != 0)
	      GP_Loads.Load_i[i].Value[j] = ReadCurve(CURVE[j+1]);
	  }
	}
	else{
	  puts("Error in Read_MPM_LoadCase_BodyForces() : Wrong format curve !!!");
	  exit(0);
	}

	/* Read the number of nodes associated */
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
		puts("Error in Read_MPM_LoadCase_BodyForces() : Check the list of nodes ");
		exit(0);
	      }
	    }
	  }
	}
	else{
	  puts("Error in Read_MPM_LoadCase_BodyForces() : Wrong format num nodes !!!");
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
	      for(int j = 0 ; j<NUM_NODES; j++){		
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
