#include "grams.h"

/***************************************************************************/

Matrix Read_CSV(char * Name_File, int NumData){

  /* Auxiliar variables */
  FILE * CSV_file;
  char line[MAXC] = {0};
  char * Aux_line; /* Check Variable */
  int Num_line = 0; /* Number of line */
  char * Field[MAXW] = {NULL};
  char Info[MAXC] = {0};
  int NumFields;
  int int_aux;

  /* Outputs */
  Matrix Mat_Out;

  /* Simulation file */
  CSV_file  = fopen(Name_File,"r");

  /* If the file is corrupted, put a wrong message */
  if( CSV_file == NULL ){
    perror(Name_File);
  }

  /* Read the first line with the header */
  Aux_line = fgets(line, sizeof line, CSV_file);
  if(Aux_line != NULL){
    strcpy(Info,line);
    NumFields = parse (Field, line, ";\n");
    Num_line++;
  }
  else{
    printf("Error in Read_CSV() : during line %i !!! \n",Num_line);
    exit(0);
  }
    
  /* Allocate data  */
  Mat_Out = MatAlloc(NumFields,NumData);
  strcpy(Mat_Out.Info,Info);  

  /* Read data */
  if(NumData == 1){ /* Multiple fields with only one line */
    Aux_line = fgets(line, sizeof line, CSV_file);
    if(Aux_line != NULL){
      int_aux = parse (Field, line, ";\n");
      if(int_aux != NumFields){
	fprintf(stderr,"%s : %s -> %i ! \n",
		"Error in Read_CSV()",
		"Incorrect number of fields",
		int_aux);
	exit(0);
      }
      Num_line++;
    }
    else{
      fprintf(stderr,"%s : %s %i !! \n",
	      "Error in Read_CSV()",
	      "during line",Num_line);
      exit(0);
    }
    for(int i = 0 ; i<NumFields ; i++){      
      Mat_Out.nV[i] = atof(Field[i]);
    }
  }
  else if(NumFields == 1){ /* Only one Field with multiple data */
    for(int i = 0 ; i<NumData ; i++){
      Aux_line = fgets(line, sizeof line, CSV_file);
      if(Aux_line != NULL){
	int_aux = parse (Field, line, ";\n");
	if(int_aux != 1){
	  fprintf(stderr,"%s : %s %i !! \n",
		  "Error in Read_CSV()",
		  "Incorrect number of fields in line",i);
	  exit(0);
	}
	Num_line++;
      }
      else{
	fprintf(stderr,"%s : %s %i !!! \n",
		"Error in Read_CSV()",
		"during line",Num_line);
	exit(0);
      }
      Mat_Out.nV[i] = atof(Field[0]);
    }
  }
  else if((NumData != 1) && (NumFields != 1)){
    for(int i = 0 ; i<NumData ; i++){ /*  */
      Aux_line = fgets(line, sizeof line, CSV_file);
      if(Aux_line != NULL){
	int_aux = parse (Field, line, ";\n");
	if(int_aux != NumFields){
	  fprintf(stderr,"%s : %s : %i !! \n",
		  "Error in Read_CSV()",
		  "Incorrect number of fields",int_aux);
	  fprintf(stderr,"%s %i \n",
		  "Error detected in line",i);
	  exit(0);
	}
	Num_line++;
      }
      else{
	fprintf(stderr,"%s : %s %i !!! \n",
		"Error in Read_CSV()",
		"during line",Num_line);
	exit(0);
      }
      for(int j = 0 ; j<NumFields ; j++){
	Mat_Out.nM[j][i] = atof(Field[j]);	
	/* Note that the information is stored in a matrix
	   with as rows as fields and as columns as data lines,
	   this is because the matrix is stored in memory as rows arrays 
	 and for future computations this storage system is more eficient */
      }
    }
  }

  /* Once we have read the file, we close it */
  fclose(CSV_file);
  
  return Mat_Out;
}
