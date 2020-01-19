#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))


/**********************************************************************/

Load * GramsContactForces(char * Name_File,GaussPoint GP_Mesh)
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
  Load * F = (Load *)Allocate_Array(GP_Mesh.NumberContactForces,
				    sizeof(Load));

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

  /* Initial message */  
  puts("*************************************************");
  printf(" \t * %s : \n\t %s \n",
	 "Begin of read contact forces in",
	 Name_File);
  
  /* Open and check .bcc file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    printf("%s : %s \n\t %s",
	   "Error in Read_MPM_LoadCase_ExtForces()",
	   "during the lecture of",
	   Name_File);
    exit(0);
  }

  
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");  

    if (strcmp(kwords[0],"F_LOAD_NUM") == 0 ){

      /* Fill the table of loads */
      for(int i = 0 ; i<GP_Mesh.NumberContactForces ; i++){

	/* Read each line and split it in individual words */
	fgets(line, sizeof line, Sim_dat);
	nkwords = parse (kwords, line," \n\t");
	if( (nkwords != 4) || (strcmp(kwords[0],"F_LOAD_GP") != 0) ){
	  printf("%s : %s \n",
		 "Error in Read_MPM_LoadCase_ExtForces()",
		 "Wrong format !!!");
	  exit(0);
	}

	/* Read the curve associated to the load */
	aux_DIR = parse (DIR,kwords[1],"={,}\n");
	if( (aux_DIR >= 2 ) && (strcmp(DIR[0],"DIR") == 0)){
	  /* Number of components of the load */
	  F[i].Dim = aux_DIR - 1;
	  /* Alocate the direction of the component */
	  F[i].Dir = (int *)Allocate_Array(F[i].Dim,sizeof(int));
	  /* Fill the direction of the load */
	  for(int j = 0 ; j<F[i].Dim ; j++){
	    F[i].Dir[j] = atoi(DIR[j+1]);
	  }
	}
	else{
	  printf("%s : %s \n",
		 "Error in Read_MPM_LoadCase_ExtForces()",
		 "Wrong format DIR={,} !!!");
	  exit(0);
	}

	/* Read the curve associated to the load */
	aux_CURVE = parse (CURVE,kwords[2],"={,}\n");
	if( (aux_CURVE == aux_DIR) && (strcmp(CURVE[0],"CURVE") == 0) ){
	  /* Alocate the table of curves */
	  F[i].Value = (Curve *)Allocate_Array(F[i].Dim,sizeof(Curve));
	  /* Fill the curve table */
	  for(int j = 0 ; j<F[i].Dim ; j++){
	    if(strcmp(CURVE[j+1],"NULL") != 0)
	      F[i].Value[j] = ReadCurve(CURVE[j+1]);
	  }
	}
	else{
	  printf("%s : %s \n",
		 "Error in Read_MPM_LoadCase_ExtForces()",
		 "Wrong format CURVE={,} !!!");
	  exit(0);
	}

	/* Read the number of nodes associated to this initial condition */
	nparam = parse (param,kwords[3],"=\n");
	if(strcmp(param[0],"NUM_NODES") == 0){
	  if(strcmp(param[1],"ALL_NODES") == 0){
	    F[i].Nodes = 
	      (int *)Allocate_Array(GP_Mesh.NumGP,sizeof(int));
	    for(int j = 0 ; j<GP_Mesh.NumGP ; j++){
	      F[i].Nodes[j] = j;
	    }
	    F[i].NumNodes = GP_Mesh.NumGP;
	  }
	  else{
	    /* Allocate array with the size */
	    F[i].Nodes = (int *)Allocate_Array(atoi(param[1]),sizeof(int));
	    F[i].NumNodes = atoi(param[1]);
   
	    /* Fill the list of nodes */
	    for(int j = 0 ; j<F[i].NumNodes ; j++){
	      fgets(line, sizeof line, Sim_dat);
	      nparam = parse(param,line," \n");
	      if(nparam == 1){
		F[i].Nodes[j] = atoi(param[0]);
	      }
	      else{
		printf("%s : %s \n",
		       "Error in Read_MPM_LoadCase_ExtForces()",
		       "Check the list of nodes ");
		exit(0);
	      }
	    }
	  }
	}
	else{
	  printf("%s : %s \n",
		 "Error in Read_MPM_LoadCase_ExtForces()",
		 "Wrong format NUM_NODES !!!");
	  exit(0);
	}
      }

    }
  }

  /* Close .dat file */
  /* Final message */
  printf("\t %s \n",
	 "-> End of read data file !!!");
  fclose(Sim_dat);

  return F;
}

