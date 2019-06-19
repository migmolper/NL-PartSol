#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"


Curve ReadCurve(char * Name_File)
/*
  Read the external forces file :
  Inputs
  - Name_file : Name of the file
  FORMAT example custom curve: 
  DAT_CURVE SCALE#double NUM#integer
  CUSTOM_CURVE
  .
  . double (NVALUE)
  .

  FORMAT example constant curve:
  DAT_CURVE SCALE#double NUM#integer
  CONSTANT_CURVE

*/
{
  /* Define load curve */
  Curve DatCurve;

  /* Scale parameter */
  double Scale;
  
  /* Define simulation file */
  FILE * Sim_dat;

  /* Auxiliar variable for reading the lines in the files */
  char line[MAXC] = {0};
  char line_values[MAXC] = {0};
  
  /* Number of element in the line , just for check */
  int nkwords,nparam;
  char * kwords[MAXW] = {NULL};
  char * param[MAXW] = {NULL};

  /* Open and check .load file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    puts("Error during the lecture of .load file");
    exit(0);
  }

  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n");

    /* Read data of the curve */
    if ( strcmp(kwords[0],"DAT_CURVE") == 0 ){
      for(int i  = 1 ; i<nkwords ; i++){ /* Loop over the words */
	nparam = parse (param,kwords[i],"#\n");
	if(nparam == 2){
	  /* Scale parameter for the curve */
	  if(strcmp(param[0],"SCALE") == 0){
	   Scale = atof(param[1]);
	  }
	  /* Number of values, it should be the same 
	     as the number of timesteps */
	  if(strcmp(param[0],"NUM") == 0){
	    DatCurve.Num = atoi(param[1]);
	  }
	}	
      }
    }

    if ( strcmp(kwords[0],"CUSTOM_CURVE") == 0 ){
      /* Allocate the curve */
      DatCurve.Fx = (double *)Allocate_Array(DatCurve.Num,sizeof(double));
      
      /* Fill the values of the curve */
      for(int i = 0 ; i<DatCurve.Num ; i++){
	fgets(line_values, sizeof line_values, Sim_dat);
	nparam = parse (param,line_values," \n");
	if(nparam == 1){
	  DatCurve.Fx[i] = Scale*atof(param[0]);
	}
	else{
	  puts("Error in ReadLoads_GP() : Check the values of the curve ");
	  exit(0);
	}
      }    
    }

    
    if ( strcmp(kwords[0],"CONSTANT_CURVE") == 0 ){
      /* Allocate the curve */
      DatCurve.Fx = (double *)Allocate_Array(DatCurve.Num,sizeof(double));

      /* Fill the values of the curve */
      for(int i = 0 ; i<DatCurve.Num ; i++){
	DatCurve.Fx[i] = Scale;
      }
    }

    
  }/* End of read file */
  fclose(Sim_dat);

  return DatCurve;
}
