
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

/***************************************************************************/

void ReadBCC(char * Name_File, Mesh FEM_Mesh)
/*
  Read the boundary conditions file :
  Inputs
  - Name_file : Name of the file
  FORMAT example 1D : BCC T#[1:100] SIGMA#[1]={0} V#[6]={0}
  FORMAT example 2D : BCC T#[1:100] SIGMA#[1]={0,0,0} V#[6]={0,0}
  FORMAT example 3D : BCC T#[1:100] SIGMA#[1]={0,0,0,0,0,0} V#[6]={0,0,0}

  Note : Only read those lines with a BCC in the init 
*/
{

  /* create a array with the number of steps and fill it with  V#[1]={-1} V#[6]={0} */
  
  /* Simulation file */
  FILE * Sim_dat;

  /* Auxiliar variable for reading the lines in the files */
  char line[MAXC] = {0};
  
  /* Number of element in the line , just for check */
  int nkwords,nparam;
  char * kwords[MAXW] = {NULL};
  char * param[MAXW] = {NULL};

  /* Variable for the kind of analysis parser */
  int nKindAnalysis;
  char * KindAnalysisParam[MAXW] = {NULL};
  /* Variable for the time parser */
  int auxT;
  char * T_range[MAXW] = {NULL};
  /* Variable for the velocity parser */
  int auxV,auxV_nod,auxV_val;
  char * V[MAXW] = {NULL};
  char * V_nod[MAXW] = {NULL};
  char * V_val[MAXW] = {NULL};
  /* Variable for the velocity parser */
  int auxSIGMA,auxSIGMA_nod,auxSIGMA_val;
  char * SIGMA[MAXW] = {NULL};
  char * SIGMA_nod[MAXW] = {NULL};
  char * SIGMA_val[MAXW] = {NULL};
  /* Variables for the TOP, BOTTOM, RIGHT, LEFT */
  double TOP_val = NAN;
  double BOTTOM_val = NAN;
  double RIGHT_val = NAN;
  double LEFT_val = NAN;
  /* Apply the BC */
  int BC_nod;

  printf("************************************************* \n");
  printf("Begin of set boundary conditions !!! \n");
  printf(" * Begin of read boundary files : %s \n",Name_File);
  printf(" * Boundary conditions values : \n");
  
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
    
    /* General BCC */
    if ( strcmp(kwords[0],"BCC") == 0 ){
      for(int i  = 0 ; i<nkwords ; i++){ /* Loop over the words */
	/* Parse the keywords parameters */
	nparam = parse (param,kwords[i],"#\n");
	if(nparam>1){ /* Read only keywords with an asignement */

	  /* Parse the time range (init:end)*/
	  if( strcmp(param[0],"T") == 0 ){
	    auxT = parse (T_range,param[1],"[:]\n");
	    if(auxT==2){
	      printf("\t -> Range : [%i -> %i] \n",atoi(T_range[0]),atoi(T_range[1]));
	    }
	    else if(auxT==1){
	      printf("\t -> Instant : %i \n",atoi(T_range[0]));
	    }
	  } /* End parse Time */

	  /* Parse the Velocity BCC */
	  if(strcmp(param[0],"V") == 0){
	    auxV = parse(V,param[1],"=\n");
	    /* Read the nodes to impose the BCC :*/
	    auxV_nod = parse(V_nod,V[0],"[:]\n");
	    /* Read the value to impose */
	    auxV_val = parse(V_val,V[1],"={,}\n");
	    for(int i = 0 ; i<auxV_nod ; i++){
	      printf("\t V[%i] = ",atoi(V_nod[i]));
	      printf("{");
	      for(int j  = 0 ; j<auxV_val ; j++ ){
		printf(" %f ",atof(V_val[j]));
	      }
	      printf("}\n");
	    }
	  } /* End parse Velocity BCC */

	  /* Parse the Stress BCC */
	  if(strcmp(param[0],"SIGMA") == 0){
	    auxSIGMA = parse(SIGMA,param[1],"=\n");
	    /* Read the nodes to impose the BCC :*/
	    auxSIGMA_nod = parse(SIGMA_nod,SIGMA[0],"[,]\n");
	    /* Read the value to impose */
	    auxSIGMA_val = parse(SIGMA_val,SIGMA[1],"={,}\n");
	    for(int i = 0 ; i<auxSIGMA_nod ; i++){
	      printf("\t SIGMA[%i] = ",atoi(SIGMA_nod[i]));
	      printf("{");
	      for(int j  = 0 ; j<auxSIGMA_val ; j++ ){
		printf(" %f ",atof(SIGMA_val[j]));
	      }
	      printf("}\n");
	    }	    
	  } /* End parse Stress BCC */
	  
	  /* Other parsers */

	  /* Fill output array */
	  /* if(auxT==1){ /\* Fill a instant *\/ */
	  /*   int T_i = atoi(T_range[0]);	     */
	  /* } */
	  /* else if(auxT == 2){ /\* Fill a time range *\/ */
	  /*   for(int T_i = atoi(T_range[0]) ; i<= atoi(T_range[1]) ; i++){ */
	  /*     BCC[T_i][auxSIGMA_nod] = ; */
	  /*   } */
	  /* } */
	  
	} /* Read word by word */
      } /* Read only keywords with an asignement */
      printf("\n");
    } /* Read general BCC */

    /* Set BCC in those node in the top of the domain */
    if ( strcmp(kwords[0],"BCC_TOP") == 0 ){
      for(int i  = 0 ; i<nkwords ; i++){ /* Loop over the words */
	/* Parse the keywords parameters */
	nparam = parse (param,kwords[i],"#\n");
	if(nparam>1){ /* Read only keywords with an asignement */

	  /* Parse the time range (init:end)*/
	  if( strcmp(param[0],"T") == 0 ){
	    auxT = parse (T_range,param[1],"[:]\n");
	    if(auxT==2){
	      printf("\t -> Range : [%i -> %i] \n",atoi(T_range[0]),atoi(T_range[1]));
	    }
	    else if(auxT==1){
	      printf("\t -> Instant : %i \n",atoi(T_range[0]));
	    }
	    else if(auxT==0){
	      printf("\t -> This boundary condition will be applied all the simulation \n");
	    }
	  } /* End parse Time */

	  /* Parse the position of the boundary */
	  if( strcmp(param[0],"TOP") == 0 ){
	    auxT = parse (T_range,param[1],"[=]\n");
	    TOP_val = atof(T_range[1]);
	    if(auxT==2){	      
	      printf("\t -> Boundary : %s = %f \n",T_range[0],TOP_val);	      
	    }
	  } /* End parse position of the boundary */	  
	  
	  /* Parse the Velocity BCC */
	  if(strcmp(param[0],"U") == 0){
	    auxV = parse(V,param[1],"=\n");
	    /* Read the DOFS to impose the BCC :*/
	    auxV_nod = parse(V_nod,V[0],"[:]\n");
	    /* Read the value to impose */
	    auxV_val = parse(V_val,V[1],"{,}\n");
	    for(int j = 0 ; j<auxV_nod ; j++){
	      if(atoi(V_nod[j]) == 1){ /* This DOF is impossed */
	  	printf("\t V[%i] = %f \n",j,atof(V_val[j]));
	  	/* Apply this boundary conditions */
	  	for(int k = 0 ; k<FEM_Mesh.NumNodesBound ; k++){
	  	  BC_nod = FEM_Mesh.NodesBound[k][0];
	  	    if(FEM_Mesh.Coordinates.nM[BC_nod][1] >= TOP_val){
		      FEM_Mesh.NodesBound[k][j+1] = 1;
	  	      FEM_Mesh.ValueBC.nM[k][j] = atof(V_val[j]);
	  	    }
	  	} /* End of apply this boundary conditions */
	      }
	    }
	  } /* End parse Velocity BCC */
	  
	} /* Read word by word */
      } /* Read only keywords with an asignement */
      printf("\n");
      
    } /* End of read BC in the top */

    /* Set BCC in those node in the bottom of the domain */
    if ( strcmp(kwords[0],"BCC_BOTTOM") == 0 ){
      for(int i  = 0 ; i<nkwords ; i++){ /* Loop over the words */
	/* Parse the keywords parameters */
	nparam = parse (param,kwords[i],"#\n");
	if(nparam>1){ /* Read only keywords with an asignement */

	  /* Parse the time range (init:end)*/
	  if( strcmp(param[0],"T") == 0 ){
	    auxT = parse (T_range,param[1],"[:]\n");
	    if(auxT==2){
	      printf("\t -> Range : [%i -> %i] \n",atoi(T_range[0]),atoi(T_range[1]));
	    }
	    else if(auxT==1){
	      printf("\t -> Instant : %i \n",atoi(T_range[0]));
	    }
	    else if(auxT==0){
	      printf("\t -> This boundary condition will be applied all the simulation \n");
	    }
	  } /* End parse Time */

	  /* Parse the position of the boundary */
	  if( strcmp(param[0],"BOTTOM") == 0 ){
	    auxT = parse (T_range,param[1],"[=]\n");
	    BOTTOM_val = atof(T_range[1]);
	    if(auxT==2){	      
	      printf("\t -> Boundary : %s = %f \n",T_range[0],BOTTOM_val);	      
	    }
	  } /* End parse position of the boundary */

	  
	  /* Parse the Velocity BCC */
	  if(strcmp(param[0],"U") == 0){
	    auxV = parse(V,param[1],"=\n");
	    /* Read the DOFS to impose the BCC :*/
	    auxV_nod = parse(V_nod,V[0],"[:]\n");
	    /* Read the value to impose */
	    auxV_val = parse(V_val,V[1],"{,}\n");
	    for(int j = 0 ; j<auxV_nod ; j++){
	      if(atoi(V_nod[j]) == 1){ /* This DOF is impossed */
	  	printf("\t V[%i] = %f \n",j,atof(V_val[j]));
	  	/* Apply this boundary conditions */
	  	for(int k = 0 ; k<FEM_Mesh.NumNodesBound ; k++){
	  	  BC_nod = FEM_Mesh.NodesBound[k][0];
	  	    if(FEM_Mesh.Coordinates.nM[BC_nod][1] <= BOTTOM_val){
		      FEM_Mesh.NodesBound[k][j+1] = 1;
	  	      FEM_Mesh.ValueBC.nM[k][j] = atof(V_val[j]);
	  	    }
	  	} /* End of apply this boundary conditions */
	      }
	    }
	  } /* End parse Velocity BCC */

	} /* Read word by word */
      } /* Read only keywords with an asignement */
      printf("\n");
    } /* End of read BC in the bottom */

    /* Set BCC in those node in the right of the domain */
    if ( strcmp(kwords[0],"BCC_RIGHT") == 0 ){
      for(int i  = 0 ; i<nkwords ; i++){ /* Loop over the words */
    	/* Parse the keywords parameters */
    	nparam = parse (param,kwords[i],"#\n");
    	if(nparam>1){ /* Read only keywords with an asignement */

    	  /* Parse the time range (init:end)*/
    	  if( strcmp(param[0],"T") == 0 ){
    	    auxT = parse (T_range,param[1],"[:]\n");
    	    if(auxT==2){
    	      printf("\t -> Range : [%i -> %i] \n",atoi(T_range[0]),atoi(T_range[1]));
    	    }
    	    else if(auxT==1){
    	      printf("\t -> Instant : %i \n",atoi(T_range[0]));
    	    }
    	    else if(auxT==0){
    	      printf("\t -> This boundary condition will be applied all the simulation \n");
    	    }
    	  } /* End parse Time */

    	  /* Parse the position of the boundary */
    	  if( strcmp(param[0],"RIGHT") == 0 ){
    	    auxT = parse (T_range,param[1],"[=]\n");
    	    RIGHT_val = atof(T_range[1]);
    	    if(auxT==2){
    	      printf("\t -> Boundary : %s = %f \n",T_range[0],RIGHT_val);
    	    }
    	  } /* End parse position of the boundary */

	  /* Parse the Velocity BCC */
	  if(strcmp(param[0],"U") == 0){
	    auxV = parse(V,param[1],"=\n");
	    /* Read the DOFS to impose the BCC :*/
	    auxV_nod = parse(V_nod,V[0],"[:]\n");
	    /* Read the value to impose */
	    auxV_val = parse(V_val,V[1],"{,}\n");
	    for(int j = 0 ; j<auxV_nod ; j++){
	      if(atoi(V_nod[j]) == 1){ /* This DOF is impossed */
	  	printf("\t V[%i] = %f \n",j,atof(V_val[j]));
	  	/* Apply this boundary conditions */
	  	for(int k = 0 ; k<FEM_Mesh.NumNodesBound ; k++){
	  	  BC_nod = FEM_Mesh.NodesBound[k][0];
	  	    if(FEM_Mesh.Coordinates.nM[BC_nod][0] >= RIGHT_val){
		      FEM_Mesh.NodesBound[k][j+1] = 1;
	  	      FEM_Mesh.ValueBC.nM[k][j] = atof(V_val[j]);
	  	    }
	  	} /* End of apply this boundary conditions */
	      }
	    }
	  } /* End parse Velocity BCC */

    	} /* Read word by word */
      } /* Read only keywords with an asignement */
      printf("\n");
    } /* End of read BC in the right */

    /* Set BCC in those node in the left of the domain */
    if ( strcmp(kwords[0],"BCC_LEFT") == 0 ){
      for(int i  = 0 ; i<nkwords ; i++){ /* Loop over the words */
    	/* Parse the keywords parameters */
    	nparam = parse (param,kwords[i],"#\n");
    	if(nparam>1){ /* Read only keywords with an asignement */

    	  /* Parse the time range (init:end)*/
    	  if( strcmp(param[0],"T") == 0 ){
    	    auxT = parse (T_range,param[1],"[:]\n");
    	    if(auxT==2){
    	      printf("\t -> Range : [%i -> %i] \n",atoi(T_range[0]),atoi(T_range[1]));
    	    }
    	    else if(auxT==1){
    	      printf("\t -> Instant : %i \n",atoi(T_range[0]));
    	    }
    	    else if(auxT==0){
    	      printf("\t -> This boundary condition will be applied all the simulation \n");
    	    }
    	  } /* End parse Time */

    	  /* Parse the position of the boundary */
    	  if( strcmp(param[0],"LEFT") == 0 ){
    	    auxT = parse (T_range,param[1],"[=]\n");
    	    LEFT_val = atof(T_range[1]);
    	    if(auxT==2){
    	      printf("\t -> Boundary : %s = %f \n",T_range[0],LEFT_val);
    	    }
    	  } /* End parse position of the boundary */

	  /* Parse the Velocity BCC */
	  if(strcmp(param[0],"U") == 0){
	    auxV = parse(V,param[1],"=\n");
	    /* Read the DOFS to impose the BCC :*/
	    auxV_nod = parse(V_nod,V[0],"[:]\n");
	    /* Read the value to impose */
	    auxV_val = parse(V_val,V[1],"{,}\n");
	    for(int j = 0 ; j<auxV_nod ; j++){
	      if(atoi(V_nod[j]) == 1){ /* This DOF is impossed */
	  	printf("\t V[%i] = %f \n",j,atof(V_val[j]));
	  	/* Apply this boundary conditions */
	  	for(int k = 0 ; k<FEM_Mesh.NumNodesBound ; k++){
	  	  BC_nod = FEM_Mesh.NodesBound[k][0];
	  	    if(FEM_Mesh.Coordinates.nM[BC_nod][0] <= LEFT_val){
		      FEM_Mesh.NodesBound[k][j+1] = 1;
	  	      FEM_Mesh.ValueBC.nM[k][j] = atof(V_val[j]);
	  	    }
	  	} /* End of apply this boundary conditions */
	      }
	    }
	  } /* End parse Velocity BCC */

    	} /* Read word by word */
      } /* Read only keywords with an asignement */
      printf("\n");
    } /* End of read BC in the left */

    
  }  /* End of read file */
  printf("End of read boundary conditions file !!! \n");
  fclose(Sim_dat);
} /* void ReadBCC(char * Name_File) */
