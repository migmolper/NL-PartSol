#include "nl-partsol.h"


Matrix Read_CSV(char * Name_File, int NumData)

{

  FILE * CSV_file;
  char line[MAXC] = {0};

  /*!
    Check Variable 
  */
  char * Aux_line; 

  /*!
    Number of line 
  */
  int Num_line = 0; 

  char * Field[MAXW] = {NULL};
  char Info[MAXC] = {0};
  int NumFields;
  int int_aux;

  /*!
    Outputs 
  */
  Matrix Mat_Out;

  /*!
    Simulation file 
  */
  CSV_file  = fopen(Name_File,"r");

  /*!
    If the file is corrupted, put a wrong message 
  */
  if( CSV_file == NULL ){
    perror(Name_File);
  }

  /*!
    Read the first line with the header 
  */
  Aux_line = fgets(line, sizeof line, CSV_file);
  if(Aux_line != NULL)
    {
      strcpy(Info,line);
      NumFields = parse (Field, line, ";\n");
      Num_line++;
    }
  else
    {
      fprintf(stderr,"%s : %s %i !!! \n",
	      "Error in Read_CSV()","during line",Num_line);
      exit(EXIT_FAILURE);
    }
    
  /*!
    Allocate data  
  */
  Mat_Out = alloc__MatrixLib__(NumFields,NumData);
  strcpy(Mat_Out.Info,Info);  

  /*!
    Read data with multiple fields with only one line
  */
  if(NumData == 1)
    {
      Aux_line = fgets(line, sizeof line, CSV_file);
      if(Aux_line != NULL)
	{
	  int_aux = parse (Field, line, ";\n");
	  if(int_aux != NumFields)
	    {
	      fprintf(stderr,"%s : %s -> %i ! \n",
		      "Error in Read_CSV()","Incorrect number of fields",int_aux);
	      exit(EXIT_FAILURE);
	    }
	  Num_line++;
	}
      else
	{
	  fprintf(stderr,"%s : %s %i !! \n",
		  "Error in Read_CSV()","during line",Num_line);
	  exit(EXIT_FAILURE);
	}
      for(int i = 0 ; i<NumFields ; i++)
	{      
	  Mat_Out.nV[i] = atof(Field[i]);
	}
    }

  /*!
    Read data with only one Field with multiple data
  */
  else if(NumFields == 1)
    {
      for(int i = 0 ; i<NumData ; i++)
	{
	  Aux_line = fgets(line, sizeof line, CSV_file);
	  if(Aux_line != NULL)
	    {
	      int_aux = parse (Field, line, ";\n");
	      if(int_aux != 1)
		{
		  fprintf(stderr,"%s : %s %i !! \n",
			  "Error in Read_CSV()",
			  "Incorrect number of fields in line",i);
		  exit(EXIT_FAILURE);
		}
	      Num_line++;
	    }
	  else
	    {
	      fprintf(stderr,"%s : %s %i !!! \n",
		      "Error in Read_CSV()","during line",Num_line);
	      exit(EXIT_FAILURE);
	    }
	  Mat_Out.nV[i] = atof(Field[0]);
	}
    }

  /*!
    Read data with multiple Field with multiple data
  */
  else if((NumData != 1) && (NumFields != 1)){
    for(int i = 0 ; i<NumData ; i++)
      { 
	Aux_line = fgets(line, sizeof line, CSV_file);
	if(Aux_line != NULL)
	  {
	    int_aux = parse (Field, line, ";\n");
	    if(int_aux != NumFields)
	      {
		fprintf(stderr,"%s : %s : %i !! \n",
			"Error in Read_CSV()",
			"Incorrect number of fields",int_aux);
		fprintf(stderr,"%s %i \n",
			"Error detected in line",i);
		exit(EXIT_FAILURE);
	      }
	    Num_line++;
	  }
	else
	  {
	    fprintf(stderr,"%s : %s %i !!! \n",
		    "Error in Read_CSV()",
		    "during line",Num_line);
	    exit(EXIT_FAILURE);
	  }
	for(int j = 0 ; j<NumFields ; j++)
	  {
	    Mat_Out.nM[j][i] = atof(Field[j]);	
	  }
      }
  }

  /*!
    Once we have read the file, we close it 
  */
  fclose(CSV_file);
  
  return Mat_Out;
}
