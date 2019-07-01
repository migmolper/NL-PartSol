#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "../InOutFun/InOutFun.h"


/***************************************************************************/

void Read_FEM_BCC(char * Name_File, Mesh * FEM_Mesh)
/*
  Read the boundary conditions file :
  Inputs
  - Name_file : Name of the file
  BCC_TOP DIR#[0,1] CURVE#{curve.txt}
  BCC_BOTTOM DIR#[0,1] CURVE#{curve.txt}
  BCC_RIGHT DIR#[0,1] CURVE#{curve.txt}
  BCC_LEFT DIR#[0,1] CURVE#{curve.txt}
  Note : Only read those lines with a BCC in the init 
*/
{
  /* Pointer to the FEM boundary conditions file */
  FILE * File_BCC;
  
  /* Array with the direction of the boundary condition */
  int * Dir_BCC;
  Dir_BCC = (int *)Allocate_ArrayZ(NumberDimensions,sizeof(int));

  /* Auxiliar structure with the load curve */
  Curve Curve_BCC;
  
  /* Auxiliar variable for reading the lines in the files */
  char line[MAXC] = {0};
  
  /* Number of element in the line , just for check */
  int nkwords,nparam;
  char * kwords[MAXW] = {NULL};
  char * param[MAXW] = {NULL};

  /* Variable for the load parser */
  int aux_DIR;
  char * DIR[MAXW] = {NULL};

  int aux_CURVE;
  char * CURVE[MAXW] = {NULL};

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
    nkwords = parse (kwords, line," \n");

    /* Set BCC in those node in the left of the domain */
    if ( ( strcmp(kwords[0],"BCC_TOP") == 0 ) ||
	 ( strcmp(kwords[0],"BCC_BOTTOM") == 0 ) ||
	 ( strcmp(kwords[0],"BCC_RIGHT") == 0 ) ||
	 ( strcmp(kwords[0],"BCC_LEFT") == 0) ){

      /* Loop over the words in each line */
      for(int i  = 0 ; i<nkwords ; i++){ 
	/* Parse/descompose the words parameters with # */
	nparam = parse (param,kwords[i],"#\n");
	/* Read only keywords with an asignement # */
	if(nparam>1){ 
	  /* Parse the direction of the BCC */
	  if(strcmp(param[0],"DIR") == 0){ 
	    /* Read the DOFS to impose the BCC :*/
	    aux_DIR = parse(DIR,param[1],"[,]\n");
	    if(aux_DIR != NumberDimensions){
	      printf("Error in ReadBCC() : Invalid direction in BCC !!! ");
	      exit(0);
	    }
	    /* Fill the array with the direction */
	    for(int j = 0 ; j<aux_DIR ; j++){
	      Dir_BCC[j] = atoi(DIR[j]);
	    }
	  }
	  /* Parse the curve assigned to this BC */
	  if(strcmp(param[0],"CURVE") == 0){
	    /* Read file of the curve */
	    aux_CURVE = parse(CURVE,param[1],"{}\n");
	    if(aux_CURVE != 1){
	      printf("Error in ReadBCC() : Wrong format for the CURVE#{File}");
	      exit(0);
	    }
	    Curve_BCC = ReadCurve(CURVE[0]); 
	  }	  	  
	} 
      } /* Read word by word */

      /* Apply this boundary conditions and copy information of the BCC */
      if(strcmp(kwords[0],"BCC_TOP") == 0){
	/* Fill the array with the direction */
	for(int j = 0 ; j<NumberDimensions ; j++){
	  FEM_Mesh->TOP.Dir[j] = Dir_BCC[j];
	}
	FEM_Mesh->TOP.Value = Curve_BCC;
	strcpy(FEM_Mesh->TOP.Info,"V");
      }
      else if(strcmp(kwords[0],"BCC_BOTTOM") == 0){
	/* Fill the array with the direction */
	for(int j = 0 ; j<NumberDimensions ; j++){
	  FEM_Mesh->BOTTOM.Dir[j] = Dir_BCC[j];
	}
	FEM_Mesh->BOTTOM.Value = Curve_BCC;
	strcpy(FEM_Mesh->BOTTOM.Info,"V");
      }
      else if(strcmp(kwords[0],"BCC_RIGHT") == 0){
	/* Fill the array with the direction */
	for(int j = 0 ; j<NumberDimensions ; j++){
	  FEM_Mesh->RIGHT.Dir[j] = Dir_BCC[j];
	}
	FEM_Mesh->RIGHT.Value = Curve_BCC;
	strcpy(FEM_Mesh->RIGHT.Info,"V");
      }
      else if(strcmp(kwords[0],"BCC_LEFT") == 0){
	/* Fill the array with the direction */
	for(int j = 0 ; j<NumberDimensions ; j++){
	  FEM_Mesh->LEFT.Dir[j] = Dir_BCC[j];
	}
	FEM_Mesh->LEFT.Value = Curve_BCC;
	strcpy(FEM_Mesh->LEFT.Info,"V");
      }
      
    } /* End of read BC */
    
  }  /* End of read file */
  printf("End of read boundary conditions file !!! \n");
  fclose(File_BCC);

  /* Free auxiliar direction array */
  free(Dir_BCC);
  
} /* BoundayConditions ReadBCC(char * Name_File) */

/**********************************************************************/

void Asign_MPM_Values(char * Name_File, GaussPoint GP_Mesh)
/*
  Read boundary conditions, initial conditions, ... and
  assign the values to the GP

  - Load format examples :
  - - External forces
  SET_GP FIELD#F DIR#[0,1] CURVE#{curve.txt} NUM_NODES#integer
  .
  . integer (NLIST)
  .

  - - Gravity load
  SET_GP FIELD#G DIR#[0,1] CURVE#{curve.txt} ALL_NODES

  - Initial values format :
  - - Initial velocities
  INIT_GP FIELD#V DIR#[0,1] CURVE#{curve.txt} NUM_NODES#integer
  .
  . integer (NLIST)
  .

  - - ALL_NODES
  INIT_GP FIELD#V DIR#[0,1] CURVE#{curve.txt} NUM_NODES#integer
  LOAD_NODES ALL_NODES

*/

{


  /* Precribed loads in some material point */

  if (strcmp(kwords[0],"LOAD_GP") == 0 ){

    for(int i  = 1 ; i<nkwords ; i++){
      nparam = parse (param,kwords[i],"#\n");
      if(nparam == 2){

	if(strcmp(param[0],"DIR") == 0){

	  /* Read the DOFS to impose the BCC :*/
	  aux_DIR = parse(DIR,param[1],"[,]\n");
	  if(aux_DIR != NumberDimensions){
	    printf("Error in ReadBCC() : Invalid direction in LOAD_GP !!! ");
	    exit(0);
	  }
	  BCC.Dir = (int *)Allocate_ArrayZ(aux_DIR,sizeof(int));
	  for(int j = 0 ; j<aux_DIR ; j++){
	    BCC.Dir[j] = atoi(DIR[j]);
	  }
	    
	}
	/* Parse the curve assigned to this BC */
	if(strcmp(param[0],"CURVE") == 0){
	  /* Read file of the curve */
	  aux_CURVE = parse(CURVE,param[1],"{}\n");
	  if(aux_CURVE != 1){
	    printf("Error in ReadBCC() : Wrong format for the CURVE#{File}");
	    exit(0);
	  }
	  BCC.Value = ReadCurve(CURVE[0]);
	}
	if(strcmp(param[0],"NUM_NODES") == 0){ // integer
	  BCC.NumNodes = atoi(param[1]);
	  BCC.Nodes =
	    (int *)Allocate_Array(BCC.NumNodes,sizeof(int));
	}
	  
      }
    }
  }
  if ( strcmp(kwords[0],"LIST_NODES") == 0 ){

    /* Fill the list of nodes */
    for(int i = 0 ; i<BCC.NumNodes ; i++){
      fgets(line_nodes, sizeof line_nodes, Sim_dat);
      nparam = parse(param,line_nodes," \n");
      if(nparam == 1){
	BCC.Nodes[i] = atoi(param[0]);
      }
      else{
	puts("Error in ReadLoads_GP() : Check the list of nodes ");
	exit(0);
      }
    }
  }

}

/* void BCC_GP_Forces(GaussPoint GP_Mesh, BoundaryConditions Loads, int TimeStep) */
/* /\* */
/*   Forces defined in the Gauss Points : */
/*   Inputs */
/* *\/ */
/* { */
/*   /\* 0º  Check the time step *\/ */
/*   if( (TimeStep < 0) || */
/*       (TimeStep > Loads.Value.Num)){ */
/*     puts("Error in BCC_GP_Forces() : The time step is out of the curve !!"); */
/*     exit(0); */
/*   } */
  
/*   /\* 1º Fill the matrix with the local forces *\/ */
/*   for(int i = 0 ; i<Loads.NumNodes ; i++){ */
/*     /\* 2º Check if this GP has a force applied *\/ */
/*     if( (Loads.Nodes[i] > GP_Mesh.NumGP) || */
/* 	(Loads.Nodes[i] < 0)){ */
/*       puts("Error in BCC_GP_Forces() : This GP does not exist !!"); */
/*       exit(0); */
/*     } */
/*     /\* 3º Loop over the dimensions *\/ */
/*     for(int j = 0 ; j<NumberDimensions ; j++){ */
/*       /\* 5º Assign the BC to the nodes *\/ */
/*       if((Loads.Dir[j] == 1) || */
/* 	 (Loads.Dir[j] == -1)){ */
/* 	/\* 6º Apply the force in the node *\/ */
/* 	GP_Mesh.Phi.F.nM[Loads.Nodes[i]][j] = */
/* 	  (double)Loads.Dir[j]*Loads.Value.Fx[TimeStep]; */
/*       }       */
/*     }    */
/*   } */

/* } */

/**********************************************************************/

BoundaryConditions SetBCC(int * Dir, int NumNodes, int * Nodes,
			  char * Info, char * Curve_File){
  
  /* Define boundary conditions */
  BoundaryConditions BCC;

  /* Apply this boundary conditions */
  BCC.Dir = Dir;
  BCC.Nodes = Nodes;
  BCC.NumNodes = NumNodes;

  /* Asign curve of values */
  BCC.Value = ReadCurve(Curve_File);

  /* Copy information of the BCC */
  strcpy(BCC.Info,Info);

  /* Return boundary condition */
  return BCC;
}


/**********************************************************************/

void BCC_Nod_Momentum(Mesh FEM_Mesh, Matrix Nodal_MOMENTUM, int TimeStep)
/*
  Apply the momentum boundary conditions over the nodes 
*/
{

  /* 1º Define auxilar variables */
  int Id_BCC; /* Index of the node where we apply the BCC */
  

  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.TOP.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }
  /* 2º Loop over the nodes of the TOP boundary */
  for(int i = 0 ; i<FEM_Mesh.TOP.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = FEM_Mesh.TOP.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      if((FEM_Mesh.TOP.Dir[j] == 1) ||
	 (FEM_Mesh.TOP.Dir[j] == -1)){
	Nodal_MOMENTUM.nM[j][Id_BCC] =
	  (double)FEM_Mesh.TOP.Dir[j]*FEM_Mesh.TOP.Value.Fx[TimeStep];
      }
    }          
  }

  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.BOTTOM.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }
  /* 2º Loop over the nodes of the BOTTOM boundary */
  for(int i = 0 ; i<FEM_Mesh.BOTTOM.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = FEM_Mesh.BOTTOM.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      if((FEM_Mesh.BOTTOM.Dir[j] == 1) ||
	 (FEM_Mesh.BOTTOM.Dir[j] == -1)){
	Nodal_MOMENTUM.nM[j][Id_BCC] =
	  (double)FEM_Mesh.BOTTOM.Dir[j]*FEM_Mesh.BOTTOM.Value.Fx[TimeStep];
      }
    }          
  }

  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.RIGHT.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }
  /* 2º Loop over the nodes of the RIGHT boundary */
  for(int i = 0 ; i<FEM_Mesh.RIGHT.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = FEM_Mesh.RIGHT.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      if((FEM_Mesh.RIGHT.Dir[j] == 1) ||
	 (FEM_Mesh.RIGHT.Dir[j] == -1)){
	Nodal_MOMENTUM.nM[j][Id_BCC] =
	  (double)FEM_Mesh.RIGHT.Dir[j]*FEM_Mesh.RIGHT.Value.Fx[TimeStep];
      }
    }          
  }

  
  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.LEFT.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }  
  /* 2º Loop over the nodes of the LEFT boundary */
  for(int i = 0 ; i<FEM_Mesh.LEFT.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = FEM_Mesh.LEFT.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      if((FEM_Mesh.LEFT.Dir[j] == 1) ||
	 (FEM_Mesh.LEFT.Dir[j] == -1)){
	Nodal_MOMENTUM.nM[j][Id_BCC] =
	  (double)FEM_Mesh.LEFT.Dir[j]*FEM_Mesh.LEFT.Value.Fx[TimeStep];
      }
    }          
  } 
  

}

/*********************************************************************/



/*********************************************************************/

/* void BC_Nod_Springs(GaussPoint GP_Mesh, BoundaryConditions Loads, int TimeStep) */
/* /\* */
/*   Add to a surface springs to simulate the winkler */
/* *\/ */
/* { */
  
/* } */
