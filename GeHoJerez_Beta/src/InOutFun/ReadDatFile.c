
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
  printf("End of read data file !!! \n");
  fclose(Sim_dat);
  

} /* void read_dat(char * Name_File) */


/***************************************************************************/

void Read_FEM_BCC(char * Name_File, Mesh * FEM_Mesh)
/*
  Read the boundary conditions file :
  Inputs
  - Name_file : Name of the file
  BCC_TOP FIELD#V CURVE#{curve.txt}
  BCC_BOTTOM FIELD#V CURVE#{curve.txt}
  BCC_RIGHT FIELD#V CURVE#{curve.txt}
  BCC_LEFT FIELD#V CURVE#{curve.txt}
  Note : Only read those lines with a BCC in the init 
*/
{
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

  /* Read the curve file */
  int aux_CURVE;
  char * CURVE[MAXW] = {NULL};

  /* Read the field */
  char * FIELD;

  
  printf("************************************************* \n");
  printf("Begin of set boundary conditions !!! \n");
  printf(" * Begin of read boundary files : %s \n",Name_File);
  
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
	    puts("Error in Read_FEM_BCC() : Wrong format of BCC !!!");
	    exit(0);
	  }

	  /* Read the curve to impose the boundary condition */
	  nparam = parse (param,kwords[2],"=\n");
	  if(strcmp(param[0],"CURVE") == 0){
	    /* Read file of the curve */
	    aux_CURVE = parse(CURVE,param[1],"{}\n");
	    if(aux_CURVE != 1){
	      printf("Error in ReadBCC() : Wrong format for the CURVE#{File}");
	      exit(0);
	    }
	    Curve_BCC = ReadCurve(CURVE[0]); 
	  }
	  else{
	    puts("Error in Read_FEM_BCC() : Wrong format of BCC !!!");
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
  printf("End of read boundary conditions file !!! \n");
  fclose(File_BCC);
  
} /* BoundayConditions ReadBCC(char * Name_File) */

/**********************************************************************/

/* Load * Read_MPM_Loads(char * Name_File, GaussPoint * GP_Mesh) */
/* /\* */

/*   Read loads over the GP (.bcc) */

/*   - Load format examples : */
/*   - - External forces */
/*   LOAD_GP FIELD#F CURVE#{curve.txt} NUM_NODES#integer */
/*   . */
/*   . integer (NLIST) */
/*   . */

/*   - - Gravity load */
/*   LOAD_GP FIELD#G CURVE#{curve.txt} NUM_NODES#ALL_NODES */

/* *\/ */

/* { */

  
/*   /\* Read the file line by line *\/ */
/*   while( fgets(line, sizeof line, File_BCC) != NULL ){ */

/*     /\* Read the line with the space as separators *\/ */
/*     nkwords = parse (kwords, line," \n\t"); */
  
/*     /\* If the first word is SET_GP read it : Asigned value */

/*        - - External forces */
/*        SET_GP FIELD#F CURVE#{curve.txt} NUM_NODES#integer */
/*        . */
/*        . integer (NLIST) */
/*        . */

/*        - - Gravity load */
/*        SET_GP FIELD#G CURVE#{curve.txt} NUM_NODES#ALL_NODES */

/*     *\/ */
/*     if (strcmp(kwords[0],"LOAD_GP") == 0 ){ */

/*       /\* Read the field to impose the boundary condition *\/ */
/*       nparam = parse (param,kwords[1],"#\n"); */
/*       if(strcmp(param[0],"FIELD") == 0){ */
	
/*       } */
/*       else{ */
/* 	puts("Error in Asign_MPM_Values() : Wrong format of SET_GP !!!"); */
/* 	exit(0); */
/*       } */
	
/*       /\* Read the curve associated to the boundary condition *\/ */
/*       nparam = parse (param,kwords[2],"#\n"); */
/*       if(strcmp(param[0],"CURVE") == 0){ */

/*       } */
/*       else{ */
/* 	puts("Error in Asign_MPM_Values() : Wrong format of SET_GP !!!"); */
/* 	exit(0); */
/*       } */

/*       /\* Read the number of nodes associated to this initial condition *\/ */
/*       if(strcmp(param[0],"NUM_NODES") == 0){ */
/* 	if(strcmp(param[1],"ALL_NODES") == 0){ */
/* 	  BCC.Nodes = */
/* 	    (int *)Allocate_Array(GaussPoint->NumGP,sizeof(int)); */
/* 	} */
/* 	else{ */
/* 	  BCC.NumNodes = atoi(param[1]); */
/* 	  BCC.Nodes = */
/* 	    (int *)Allocate_Array(BCC.NumNodes,sizeof(int)); */
/* 	  /\* Fill the list of nodes *\/ */
/* 	  for(int j = 0 ; j<BCC.NumNodes ; j++){ */
/* 	    fgets(line_nodes, sizeof line_nodes, Sim_dat); */
/* 	    nparam = parse(param,line_nodes," \n"); */
/* 	    if(nparam == 1){ */
/* 	      BCC.Nodes[j] = atoi(param[0]); */
/* 	    } */
/* 	    else{ */
/* 	      puts("Error in ReadLoads_GP() : Check the list of nodes "); */
/* 	      exit(0); */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */

/*     } */

/*   } */


/*   return LoadsGP; */
/* } */

/**********************************************************************/

/* void Read_MPM_InitVal(char * Name_File, GaussPoint * GP_Mesh) */
/* /\* */
/*   - Initial values format : */
/*   - - Initial velocities */
/*   INIT_GP FIELD#V CURVE#{curve.txt} NUM_NODES#integer */
/*   . */
/*   . integer (NLIST) */
/*   . */

/*   - - ALL_NODES */
/*   INIT_GP FIELD#V CURVE#{curve.txt} NUM_NODES#ALL_NODES */
/* *\/ */
/* { */

/*   /\* Read the file line by line *\/ */
/*   while( fgets(line, sizeof line, File_BCC) != NULL ){ */

/*     /\* If the first word is INIT_GP read it : set an initial value  */

/*        - Initial values format : */
/*        - - Initial velocities */
/*        INIT_GP FIELD#V CURVE#{curve.txt} NUM_NODES#integer */
/*        . */
/*        . integer (NLIST) */
/*        . */
/*        - - ALL_NODES */
/*        INIT_GP FIELD#V CURVE#{curve.txt} NUM_NODES#ALL_NODES */

/*     *\/ */
/*     if (strcmp(kwords[0],"INIT_GP") == 0 ){ */

/*       /\* Read the field to impose the initial condition *\/ */
/*       nparam = parse (param,kwords[1],"#\n"); */
/*       if(strcmp(param[0],"FIELD") == 0){ */
	
/*       } */
/*       else{ */
/* 	puts("Error in Asign_MPM_Values() : Wrong format of INIT_GP !!!"); */
/* 	exit(0); */
/*       } */
      
/*       /\* Read the curve associated to the initial condition *\/ */
/*       nparam = parse (param,kwords[2],"#\n"); */
/*       if(strcmp(param[0],"CURVE") == 0){ */

/*       } */
/*       else{ */
/* 	puts("Error in Asign_MPM_Values() : Wrong format of INIT_GP !!!"); */
/* 	exit(0); */
/*       } */

/*       /\* Read the number of nodes associated to this initial condition *\/ */
/*       nparam = parse (param,kwords[2],"#\n"); */
/*       if(strcmp(param[0],"NUM_NODES") == 0){ */
	
/*       } */
/*       else{ */
/* 	puts("Error in Asign_MPM_Values() : Wrong format of INIT_GP !!!"); */
/* 	exit(0); */
/*       } */

      
/*     } */


/*   } */

/* } */


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
