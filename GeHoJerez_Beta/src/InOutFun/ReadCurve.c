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

  FORMAT example heaviside curve:
  DAT_CURVE SCALE#double NUM#integer
  HEAVISIDE_CURVE Tc#integer

  FORMAT example delta curve:
  DAT_CURVE SCALE#double NUM#integer
  DELTA_CURVE Tc#integer

  FORMAT example hat curve:
  DAT_CURVE SCALE#double NUM#integer
  HAT_CURVE T0#integer T1#integer

  FORMAT example ramp curve:
  DAT_CURVE SCALE#double NUM#integer
  RAMP_CURVE

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

  /* Axiliar variable for defining default curves */
  int Tc, T0, T1;

  /* Open and check .load file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    puts("Error during the lecture of .load file");
    exit(0);
  }

  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n");

    /* Read general data of the curve */
    /*********************************************************************/
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
    /*********************************************************************/

    /* CUSTOM CURVE */
    /*********************************************************************/
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
	  puts("Error in ReadCurve() : Check the values of the curve ");
	  exit(0);
	}
      }    
    }
    /*********************************************************************/

    /* CONSTANT CURVE */
    /*********************************************************************/
    if ( strcmp(kwords[0],"CONSTANT_CURVE") == 0 ){
      /* Allocate the curve */
      DatCurve.Fx = (double *)Allocate_Array(DatCurve.Num,sizeof(double));

      /* Fill the values of the curve */
      for(int i = 0 ; i<DatCurve.Num ; i++){
	DatCurve.Fx[i] = Scale;
      }
    }
    /*********************************************************************/

    /* RAMP CURVE */
    /*********************************************************************/
    if ( strcmp(kwords[0],"RAMP_CURVE") == 0 ){
      /* Allocate the curve */
      DatCurve.Fx = (double *)Allocate_Array(DatCurve.Num,sizeof(double));

      /* Fill the values of the curve */
      for(int i = 0 ; i<DatCurve.Num ; i++){
	DatCurve.Fx[i] = Scale*(double)i/DatCurve.Num;
      }
    }
    /*********************************************************************/

    /* HEAVISIDE CURVE */
    /*********************************************************************/
    if ( strcmp(kwords[0],"HEAVISIDE_CURVE") == 0 ){
      /* Allocate the curve */
      DatCurve.Fx = (double *)Allocate_Array(DatCurve.Num,sizeof(double));

      nparam = parse (param,kwords[1],"#\n");
      if( (strcmp(param[0],"Tc") == 0) || (nparam != 2) ){
	puts("Error in ReadCurve() : Wrong parameters of the Heaviside curve !!!");
	exit(0);
      }

      /* Critic time of the Heaviside */
      Tc = atoi(param[1]);
      if(Tc > DatCurve.Num){
	puts("Error in ReadCurve() : Tc > DatCurve.Num in the Heaviside curve !!!");
	exit(0);
      }
      if(Tc < 0){
	puts("Error in ReadCurve() : Tc < 0 in the Heaviside curve !!!");
	exit(0);
      }
      
      /* Fill the values of the curve */
      for(int i = 0 ; i< DatCurve.Num; i++){
	if(i<=Tc)
	  DatCurve.Fx[i] = 0;
	else
	  DatCurve.Fx[i] = Scale;
      }
    }
    /*********************************************************************/

    /* DELTA CURVE */
    /*********************************************************************/
    if ( strcmp(kwords[0],"DELTA_CURVE") == 0 ){
      /* Allocate the curve */
      DatCurve.Fx = (double *)Allocate_Array(DatCurve.Num,sizeof(double));

      nparam = parse (param,kwords[1],"#\n");
      if( (strcmp(param[0],"Tc") == 0) || (nparam != 2) ){
	puts("Error in ReadCurve() : Wrong parameters of the Delta curve !!!");
	exit(0);
      }

      /* Critic time of the Delta */
      Tc = atoi(param[1]);
      if(Tc > DatCurve.Num){
	puts("Error in ReadCurve() : Tc > DatCurve.Num in the Delta curve !!!");
	exit(0);
      }
      if(Tc < 0){
	puts("Error in ReadCurve() : Tc < 0 in the delta curve !!!");
	exit(0);
      }
      
      /* Fill the values of the curve */
      for(int i = 0 ; i< DatCurve.Num; i++){
	if(i == Tc)
	  DatCurve.Fx[i] = Scale;
	else
	  DatCurve.Fx[i] = 0;
      }
    }
    /*********************************************************************/

    /* HAT CURVE */
    /*********************************************************************/
    if ( strcmp(kwords[0],"HAT_CURVE") == 0 ){
      /* Allocate the curve */
      DatCurve.Fx = (double *)Allocate_Array(DatCurve.Num,sizeof(double));

      /* Read */
      nparam = parse (param,kwords[1],"#\n");
      if( (strcmp(param[0],"T0") == 0) || (nparam != 2) ){
	puts("Error in ReadCurve() : Wrong parameters of the Hat curve !!!");
	exit(0);
      }
      T0 = atoi(param[1]);

      nparam = parse (param,kwords[1],"#\n");
      if( (strcmp(param[0],"T1") == 0) || (nparam != 2) ){
	puts("Error in ReadCurve() : Wrong parameters of the Hat curve !!!");
	exit(0);
      }
      T1 = atoi(param[1]);


      /* Check the time T0 and T1 of the hat */
      if(T0 > T1){
	puts("Error in ReadCurve() : T0 > T1 in the Hat curve !!!");
	exit(0);
      }
      if(T1 > DatCurve.Num){
	puts("Error in ReadCurve() : T1 > DatCurve.Num in the Hat curve !!!");
	exit(0);
      }
      if(T0 < 0){
	puts("Error in ReadCurve() : T0 < 0 in the Hat curve !!!");
	exit(0);
      }
      
      /* Fill the values of the curve */
      for(int i = 0 ; i< DatCurve.Num; i++){
	if(i < T0)
	  DatCurve.Fx[i] = 0;
	else if ( (i >= T0) || (i <= T1) )
	  DatCurve.Fx[i] = Scale;
	else
	  DatCurve.Fx[i] = 0;
      }
    }
    /*********************************************************************/



    
  }/* End of read file */
  fclose(Sim_dat);

  return DatCurve;
}
