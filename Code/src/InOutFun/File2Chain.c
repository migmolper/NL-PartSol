#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

ChainPtr File2Chain(char * Name_File){
  
  /* File pointer */
  FILE * Sim_dat;

  /* File lines */
  int Numers_Line;
  char Line_File[MAXC] = {0}; 
  char * Parse_File[MAXW] = {NULL};
  
  /* */
  ChainPtr File_Chain = NULL;
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in File2Chain()",
	    "Incorrect lecture of",
	    Name_File);
    exit(0);
  }

  /* Read file line by line */
  while(fgets(Line_File,sizeof(Line_File),Sim_dat) != NULL){
    Numers_Line = parse(Parse_File,Line_File," \n\t");
    PushNodeTop(&File_Chain,atoi(Parse_File[0]));    
  }

  /* Close file */
  fclose(Sim_dat);

  return File_Chain;
  
}
