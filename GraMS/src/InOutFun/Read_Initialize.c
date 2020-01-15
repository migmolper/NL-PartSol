#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

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
  puts("*************************************************");
  printf(" \t * %s : \n\t %s \n",
	 "Begin of read initial conditions in",
	 Name_File);
  
  /* Open and check .bcc file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    printf("%s : %s",
	   "Error in Read_MPM_InitVal() during the lecture of",
	   Name_File);
    exit(0);
  }

  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");

    /* When the parser find the keyword INIT_NUM : 
       start reading initial conditions */
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
	  printf("%s : %s \n",
		 "Error in Read_MPM_InitVal()",
		 "Wrong format !!!");
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
	    printf("%s : %s \n",
		 "Error in Read_MPM_InitVal()",
		 "Wrong format !!!");
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
	      printf("%s : %s \n",
		     "Error in Read_MPM_InitVal()",
		     "Wrong number of dimensions for VALUE={,}");
	      exit(0);
	    }
	  }
	  else{
	    printf("%s : %s \n",
		 "Error in Read_MPM_InitVal()",
		 "Wrong format !!!");
	    exit(0);
	  }

	  /* Read the number of nodes associated to this initial condition */
	  nparam = parse (param,kwords[3],"=\n");
	  if(nparam != 2){
	    printf("%s : %s \n",
		   "Error in Read_MPM_InitVal()",
		   "Wrong format !!!");
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
		  printf("%s : %s \n",
			 "Error in ReadLoads_GP()",
			 "Check the list of nodes");
		  exit(0);
		}
	      }
	    }
	  }
	  else{
	    printf("%s : %s \n",
		   "Error in Read_MPM_InitVal()",
		   "Wrong format of INIT_GP !!!");
	    exit(0);
	  }

	  /* Asign initial conditions*/

      
	}
      }
      
    } /* End if(strcmp(kwords[0],"INIT_NUM") == 0 ) */
    
  } /* end While */
  
    /* Close .dat file */
  /* Final message */
  printf("\t %s \n",
	 "* End of read data file !!!");
  fclose(Sim_dat);
  
}

/**********************************************************************/
