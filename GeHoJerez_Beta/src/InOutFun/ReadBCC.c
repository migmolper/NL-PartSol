#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

/***************************************************************************/

BoundaryConditions ReadBCC(char * Name_File, Mesh FEM_Mesh)
/*
  Read the boundary conditions file :
  Inputs
  - Name_file : Name of the file
  BCC_BOTTOM V#[0:1]={NAN,curve.txt}

  - Load format : 
  LOAD_GP DIM#integer CURVE#curve.txt NUM_NODES#integer
  LOAD_NODES 
  .
  . integer (NLIST)
  .
  END LIST_NODES

  Note : Only read those lines with a BCC in the init 
*/
{

  /* Create a array with the number of steps and fill it with  V#[1]={-1} V#[6]={0} */
  BoundaryConditions BCC;
  
  /* Simulation file */
  FILE * Sim_dat;

  /* Auxiliar variable for reading the lines in the files */
  char line[MAXC] = {0};
  char line_nodes[MAXC] = {0};
  
  /* Number of element in the line , just for check */
  int nkwords,nparam;
  char * kwords[MAXW] = {NULL};
  char * param[MAXW] = {NULL};

  /* /\* Variable for the time parser *\/ */
  /* int auxT; */
  /* char * T_range[MAXW] = {NULL}; */
  /* Variable for the velocity parser */
  int auxV,auxV_nod,auxV_val;
  char * V[MAXW] = {NULL};
  char * V_nod[MAXW] = {NULL};
  char * V_val[MAXW] = {NULL};
  /* /\* Variable for the stress parser *\/ */
  /* int auxSIGMA,auxSIGMA_nod,auxSIGMA_val; */
  /* char * SIGMA[MAXW] = {NULL}; */
  /* char * SIGMA_nod[MAXW] = {NULL}; */
  /* char * SIGMA_val[MAXW] = {NULL}; */

  printf("************************************************* \n");
  printf("Begin of set boundary conditions !!! \n");
  printf(" * Begin of read boundary files : %s \n",Name_File);
  
  /* Open and check .bcc file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    puts("Error during the lecture of .bcc file");
    exit(0);
  }

  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n");

    
    /* Set BCC in those node in the top of the domain */
    if ( strcmp(kwords[0],"BCC_TOP") == 0 ){

      /* Screen output */
      printf(" * This boundary condition is in the TOP \n");
      
      for(int i  = 0 ; i<nkwords ; i++){ /* Loop over the words */
	/* Parse the keywords parameters */
	nparam = parse (param,kwords[i],"#\n");
	if(nparam>1){ /* Read only keywords with an asignement */

	  /* Apply this boundary conditions */
	  BCC.Nodes = FEM_Mesh.TOP;
	  BCC.NumNodes = FEM_Mesh.NumTOP;
	  
	  /* Parse the Velocity BCC */
	  if(strcmp(param[0],"V") == 0){

	    /* Screen output */
	    printf(" * This boundary condition is over the velocity field \n");
	    
	    auxV = parse(V,param[1],"=\n");
	    if(auxV != 2){
	      puts("Error in ReadBCC() : Check the input format of the velocity !!!");
	      exit(0);
	    }

	    /* Read the DOFS to impose the BCC :*/
	    auxV_nod = parse(V_nod,V[0],"[:]\n");
	    for(int j = 0 ; j<auxV_nod ; j++){
	      if(atoi(V_nod[j]) == 1){ /* This DOF is impossed */
		BCC.Dim = atoi(V_nod[j]);
	      }
	    }

	    /* Screen output */
	    printf(" * This boundary condition is over the %i direction \n",BCC.Dim);
	    
	    /* Read curve with the value to impose */
	    auxV_val = parse(V_val,V[1],"{,}\n");
	    if(auxV_val != NumberDimensions){
	      puts("Error in ReadBCC() : Check the number of dim of the velocity !!!");
	      exit(0);
	    }	    
	    /* Init reading curve */
	    printf(" * Begin of read load curve file !!! \n");
	    printf(" \t -> WORKING ... \n");
	    BCC.Value = ReadCurve(V_val[BCC.Dim]);
	    printf(" \t -> DONE !!! \n");
	    printf(" * End of read load curve file !!! \n");

	    /* Copy information of the BCC */
	    strcpy(BCC.Info,"V_TOP");
	    
	  }
	  
	} /* End parse Velocity BCC */ 
      } /* Read word by word */
    } /* End of read BC in the top */

    /* Set BCC in those node in the bottom of the domain */
    if ( strcmp(kwords[0],"BCC_BOTTOM") == 0 ){

      /* Screen output */
      printf(" * This boundary condition is in the BOTTOM \n");
      
      for(int i  = 0 ; i<nkwords ; i++){ /* Loop over the words */
	/* Parse the keywords parameters */
	nparam = parse (param,kwords[i],"#\n");
	if(nparam>1){ /* Read only keywords with an asignement */

	  /* Apply this boundary conditions */
	  BCC.Nodes = FEM_Mesh.BOTTOM;
	  BCC.NumNodes = FEM_Mesh.NumBOTTOM;
	  	  
	  /* Parse the Velocity BCC */
	  if(strcmp(param[0],"V") == 0){

	    /* Screen output */
	    printf(" * This boundary condition is over the velocity field \n");
	    
	    auxV = parse(V,param[1],"=\n");
	    if(auxV != 2){
	      puts("Error in ReadBCC() : Check the input format of the velocity !!!");
	      exit(0);
	    }

	    /* Read the DOFS to impose the BCC :*/
	    auxV_nod = parse(V_nod,V[0],"[:]\n");
	    for(int j = 0 ; j<auxV_nod ; j++){
	      if(atoi(V_nod[j]) == 1){ /* This DOF is impossed */
		BCC.Dim = atoi(V_nod[j]);
	      }
	    }

	    /* Screen output */
	    printf(" * This boundary condition is over the %i direction \n",BCC.Dim);
	    
	    /* Read curve with the value to impose */
	    auxV_val = parse(V_val,V[1],"{,}\n");
	    if(auxV_val != NumberDimensions){
	      puts("Error in ReadBCC() : Check the number of dim of the velocity !!!");
	      exit(0);
	    }	    
	    /* Init reading curve */
	    printf(" * Begin of read load curve file !!! \n");
	    printf(" \t -> WORKING ... \n");
	    BCC.Value = ReadCurve(V_val[BCC.Dim]);
	    printf(" \t -> DONE !!! \n");
	    printf(" * End of read load curve file !!! \n");

	    /* Copy information of the BCC */
	    strcpy(BCC.Info,"V_BOTTOM");
	    
	  }
	  
	} /* End parse Velocity BCC */ 
      } /* Read word by word */
    } /* End of read BC in the bottom */

    /* Set BCC in those node in the right of the domain */
    if ( strcmp(kwords[0],"BCC_RIGHT") == 0 ){

      /* Screen output */
      printf(" * This boundary condition is in the RIGHT \n");
      
      for(int i  = 0 ; i<nkwords ; i++){ /* Loop over the words */
	/* Parse the keywords parameters */
	nparam = parse (param,kwords[i],"#\n");
	if(nparam>1){ /* Read only keywords with an asignement */

	  /* Apply this boundary conditions */
	  BCC.Nodes = FEM_Mesh.RIGHT;
	  BCC.NumNodes = FEM_Mesh.NumRIGHT;
 	  	   
	  /* Parse the Velocity BCC */
	  if(strcmp(param[0],"V") == 0){

	    /* Screen output */
	    printf(" * This boundary condition is over the velocity field \n");
	    
	    auxV = parse(V,param[1],"=\n");
	    if(auxV != 2){
	      puts("Error in ReadBCC() : Check the input format of the velocity !!!");
	      exit(0);
	    }

	    /* Read the DOFS to impose the BCC :*/
	    auxV_nod = parse(V_nod,V[0],"[:]\n");
	    for(int j = 0 ; j<auxV_nod ; j++){
	      if(atoi(V_nod[j]) == 1){ /* This DOF is impossed */
		BCC.Dim = atoi(V_nod[j]);
	      }
	    }

	    /* Screen output */
	    printf(" * This boundary condition is over the %i direction \n",BCC.Dim);
	    
	    /* Read curve with the value to impose */
	    auxV_val = parse(V_val,V[1],"{,}\n");
	    if(auxV_val != NumberDimensions){
	      puts("Error in ReadBCC() : Check the number of dim of the velocity !!!");
	      exit(0);
	    }	    
	    /* Init reading curve */
	    printf(" * Begin of read load curve file !!! \n");
	    printf(" \t -> WORKING ... \n");
	    BCC.Value = ReadCurve(V_val[BCC.Dim]);
	    printf(" \t -> DONE !!! \n");
	    printf(" * End of read load curve file !!! \n");

	    /* Copy information of the BCC */
	    strcpy(BCC.Info,"V_RIGHT");
	    
	  }
	  
	} /* End parse Velocity BCC */ 
      } /* Read word by word */
    } /* End of read BC in the right */

    /* Set BCC in those node in the left of the domain */
    if ( strcmp(kwords[0],"BCC_LEFT") == 0 ){

      /* Screen output */
      printf(" * This boundary condition is in the LEFT \n");
     
      for(int i  = 0 ; i<nkwords ; i++){ /* Loop over the words */
	/* Parse the keywords parameters */
	nparam = parse (param,kwords[i],"#\n");
	if(nparam>1){ /* Read only keywords with an asignement */

	  /* Apply this boundary conditions */
	  BCC.Nodes = FEM_Mesh.LEFT;
	  BCC.NumNodes = FEM_Mesh.NumLEFT;
	  
	  /* Parse the Velocity BCC */
	  if(strcmp(param[0],"V") == 0){

	    /* Screen output */
	    printf(" * This boundary condition is over the velocity field \n");
	    
	    auxV = parse(V,param[1],"=\n");
	    if(auxV != 2){
	      puts("Error in ReadBCC() : Check the input format of the velocity !!!");
	      exit(0);
	    }

	    /* Read the DOFS to impose the BCC :*/
	    auxV_nod = parse(V_nod,V[0],"[:]\n");
	    for(int j = 0 ; j<auxV_nod ; j++){
	      if(atoi(V_nod[j]) == 1){ /* This DOF is impossed */
		BCC.Dim = atoi(V_nod[j]);
	      }
	    }

	    /* Screen output */
	    printf(" * This boundary condition is over the %i direction \n",BCC.Dim);
	    
	    /* Read curve with the value to impose */
	    auxV_val = parse(V_val,V[1],"{,}\n");
	    if(auxV_val != NumberDimensions){
	      puts("Error in ReadBCC() : Check the number of dim of the velocity !!!");
	      exit(0);
	    }	    
	    /* Init reading curve */
	    printf(" * Begin of read load curve file !!! \n");
	    printf(" \t -> WORKING ... \n");
	    BCC.Value = ReadCurve(V_val[BCC.Dim]);
	    printf(" \t -> DONE !!! \n");
	    printf(" * End of read load curve file !!! \n");

	    /* Copy information of the BCC */
	    strcpy(BCC.Info,"V_LEFT");
	    
	  }
	  
	} /* End parse Velocity BCC */ 
      } /* Read word by word */
    } /* End of read BC in the left */

    /* Precribed loads in some material point */

    if (strcmp(kwords[0],"LOAD_GP") == 0 ){

      for(int i  = 1 ; i<nkwords ; i++){
	nparam = parse (param,kwords[i],"#\n");	
	if(nparam == 2){

	  if(strcmp(param[0],"DIM") == 0){ 
	    BCC.Dim = atoi(param[1]); 
	  }	  
	  if(strcmp(param[0],"CURVE") == 0){ 
	    BCC.Value = ReadCurve(param[1]);
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
	nparam = parse (param,line_nodes," \n");
	if(nparam == 1){
	  BCC.Nodes[i] = atoi(param[0]);
	}
	else{
	  puts("Error in ReadLoads_GP() : Check the list of nodes ");
	  exit(0);
	}
      }
    }

    
  }  /* End of read file */
  printf("End of read boundary conditions file !!! \n");
  fclose(Sim_dat);

  /* Return output */
  return BCC;  
  
} /* BoundayConditions ReadBCC(char * Name_File) */
