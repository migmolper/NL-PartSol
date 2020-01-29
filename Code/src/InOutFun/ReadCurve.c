#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"


Curve ReadCurve(char * Name_File)
/*
  Read the external forces file :
  Inputs
  - Name_file : Name of the file

  FORMAT example custom curve + constant: 
  DAT_CURVE NUM#3
  HEAVISIDE_CURVE SCALE#2.5 Tc#integer
  CONSTANT_CURVE SCALE#0
  
  Curve parameters :
  CONSTANT_CURVE SCALE#double
  RAMP CURVE SCALE#double 
  HEAVISIDE_CURVE SCALE#double Tc#integer
  DELTA_CURVE SCALE#double Tc#integer
  HAT_CURVE SCALE#double T0#integer T1#integer

  Future curves :
  CUSTOM_CURVE

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
  
  /* Number of element in the line , just for check */
  int nkwords,nparam;
  char * kwords[MAXW] = {NULL};
  char * param[MAXW] = {NULL};

  /* Axiliar variable for defining default curves */
  int Tc, T0, T1;

  /* Open and check .load file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    printf("Error in ReadCurve() during the lecture of %s \n",Name_File);
    exit(0);
  }

  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");

    /* Read general data of the curve */
    /*********************************************************************/
    if ( strcmp(kwords[0],"DAT_CURVE") == 0 ){
      for(int i  = 1 ; i<nkwords ; i++){ /* Loop over the words */
	nparam = parse (param,kwords[i],"#\n");
	if(nparam == 2){
	  /* Number of values, it should be the same 
	     as the number of timesteps */
	  if(strcmp(param[0],"NUM") == 0){
	    DatCurve.Num = atoi(param[1]);
	  }
	}	
      }
      /* ALLOCATE FX MATRIX */
      DatCurve.Fx = (double *)Allocate_ArrayZ(DatCurve.Num,sizeof(double));
    }
    /*********************************************************************/
        
    /* CONSTANT_CURVE SCALE#double  */
    /*********************************************************************/
    if ( strcmp(kwords[0],"CONSTANT_CURVE") == 0 ){

      /* Scale parameter of the curve */
      nparam = parse (param,kwords[1],"#\n");
      if(strcmp(param[0],"SCALE") == 0){
	Scale = atof(param[1]);
      }
      else{
	puts("Error in ReadCurve() : Wrong parameters of the constant curve !!!");
	exit(0);
      }

      /* Fill the values of the curve */
      for(int i = 0 ; i<DatCurve.Num ; i++){
	DatCurve.Fx[i] += Scale;
      }
    }
    /*********************************************************************/

    /* RAMP CURVE SCALE#double */
    /*********************************************************************/
    if ( strcmp(kwords[0],"RAMP_CURVE") == 0 ){

      /* Scale parameter of the curve */
      nparam = parse (param,kwords[1],"#\n");
      if(strcmp(param[0],"SCALE") == 0){
	Scale = atof(param[1]);
      }
      else{
	puts("Error in ReadCurve() : Wrong parameters of the ramp curve !!!");
	exit(0);
      }

      /* Fill the values of the curve */
      for(int i = 0 ; i<DatCurve.Num ; i++){
	DatCurve.Fx[i] = Scale*(double)i/DatCurve.Num;
      }
    }
    /*********************************************************************/

    /* HEAVISIDE_CURVE SCALE#double Tc#integer */
    /*********************************************************************/
    if ( strcmp(kwords[0],"HEAVISIDE_CURVE") == 0 ){

      /* Scale parameter of the curve */
      nparam = parse (param,kwords[1],"#\n");
      if(strcmp(param[0],"SCALE") == 0){
	Scale = atof(param[1]);
      }
      else{
	puts("Error in ReadCurve() : Wrong parameters of the Heaviside curve !!!");
	exit(0);
      }

      /* Critical time of the curve */
      nparam = parse (param,kwords[2],"#\n");
      if(strcmp(param[0],"Tc") == 0){
	Tc = atoi(param[1]);
	/* Check the critical time */
	if(Tc > DatCurve.Num){
	  puts("Error in ReadCurve() : Tc > DatCurve.Num in the Heaviside curve !!!");
	  exit(0);
	}
	if(Tc < 0){
	  puts("Error in ReadCurve() : Tc < 0 in the Heaviside curve !!!");
	  exit(0);
	}
      }
      else{
	puts("Error in ReadCurve() : Wrong parameters of the Heaviside curve !!!");
	exit(0);
      }
           
      /* Fill the values of the curve */
      for(int i = 0 ; i< DatCurve.Num; i++){
	if(i<=Tc)
	  DatCurve.Fx[i] += 0;
	else
	  DatCurve.Fx[i] += Scale;
      }
    }
    /*********************************************************************/

    /* DELTA_CURVE SCALE#double Tc#integer */
    /*********************************************************************/
    if ( strcmp(kwords[0],"DELTA_CURVE") == 0 ){

      /* Scale parameter of the curve */
      nparam = parse (param,kwords[1],"#\n");
      if(strcmp(param[0],"SCALE") == 0){
	Scale = atof(param[1]);
      }
      else{
	puts("Error in ReadCurve() : Wrong parameters of the delta curve !!!");
	exit(0);
      }

      /* Critical time of the delta */
      nparam = parse (param,kwords[2],"#\n");
      if(strcmp(param[0],"Tc") == 0){
	Tc = atoi(param[1]);
	/* Check the critical time */
	if(Tc > DatCurve.Num){
	  puts("Error in ReadCurve() : Tc > DatCurve.Num in the delta curve !!!");
	  exit(0);
	}
	if(Tc < 0){
	  puts("Error in ReadCurve() : Tc < 0 in the delta curve !!!");
	  exit(0);
	}
      }
      else{
	puts("Error in ReadCurve() : Wrong parameters of the delta curve !!!");
	exit(0);
      }
      
      /* Fill the values of the curve */
      for(int i = 0 ; i< DatCurve.Num; i++){
	if(i == Tc)
	  DatCurve.Fx[i] += Scale;
	else
	  DatCurve.Fx[i] += 0;
      }
    }
    /*********************************************************************/

    /* HAT_CURVE SCALE#double T0#integer T1#integer */
    /*********************************************************************/
    if ( strcmp(kwords[0],"HAT_CURVE") == 0 ){

      /* Scale parameter of the curve */
      nparam = parse (param,kwords[1],"#\n");
      if(strcmp(param[0],"SCALE") == 0){
	Scale = atof(param[1]);
      }
      else{
	puts("Error in ReadCurve() : Wrong parameters of the hat curve !!!");
	exit(0);
      }

      /* Read critical times */
      nparam = parse (param,kwords[2],"#\n");
      if( (strcmp(param[0],"T0") == 0) && (nparam == 2) ){
	T0 = atoi(param[1]);
      }
      else{
	puts("Error in ReadCurve() : Wrong parameters of the Hat curve !!!");
	exit(0);
      }
      
      nparam = parse (param,kwords[3],"#\n");
      if( (strcmp(param[0],"T1") == 0) && (nparam == 2) ){
	T1 = atoi(param[1]);
      }
      else{
	puts("Error in ReadCurve() : Wrong parameters of the Hat curve !!!");
	exit(0);
      }      


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
	  DatCurve.Fx[i] += 0;
	else if ( (i >= T0) || (i <= T1) )
	  DatCurve.Fx[i] += Scale;
	else
	  DatCurve.Fx[i] += 0;
      }
    }
    /*********************************************************************/

    
  }/* End of read file */
  fclose(Sim_dat);

  return DatCurve;
}
