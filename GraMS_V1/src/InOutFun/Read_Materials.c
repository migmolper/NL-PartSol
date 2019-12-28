#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "../GRAMS/grams.h"

/**********************************************************************/

Material * InitializeMaterials(char * Name_File, GaussPoint GP_Mesh)
/*
GramsMaterials (id=0) {
               Type=LE
	       E=6.e9
	       mu=0.2
	       Fracture=TRUE
	       Ceps=1.5
	       Gf=0.00001
}  
*/
{
  /* Asign material library to an auxiliar variable */
  Material Mat_GP;

  /* Simulation file */
  FILE * Sim_dat;

  /* Index of the material */
  int Aux_Mat_id;
  char * Parse_Mat_id[MAXW] = {NULL};
  int Mat_id;

  /* Material properties */
  char Line_Material_Prop[MAXC] = {0};
  char * Parse_Mat_Prop[MAXW] = {NULL};

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  /* Auxiliar variable for status */
  char * STATUS_LINE;
  int CountMaterials = 0;

  /* Initial message */  
  puts("*************************************************");
  printf(" \t %s : \n\t %s \n",
	 "* Begin of materials properties in",
	 Name_File);
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    printf("%s : \n\t %s %s",
	   "Error in Read_MPM_Materials()",
	   "Incorrect lecture of",
	   Name_File);
    exit(0);
  }

  /* Allocate table with the material */
  GP_Mesh.Mat = malloc(GP_Mesh.NumMat*sizeof(Material));
  if(GP_Mesh.Mat == NULL){
    printf("%s : %s \n",
	   "Error in Read_MPM_Materials()",
	   "Memory error for Mat");
    exit(0);
  }
  
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");
    if (nkwords < 0){
      printf("%s : %s \n",
	     "Error in Read_MPM_Materials()",
	     "Parser failed");
      exit(0);
    }

    if (strcmp(kwords[0],"GramsMaterials") == 0 ){

      /* Count the number of materials */
      ++CountMaterials;

      /* Read the index of the material */
      Aux_Mat_id = parse (Parse_Mat_id, kwords[1],"(=)");
      Mat_id = atoi(Parse_Mat_id[1]);
      if( (Aux_Mat_id != 2) ||
	  (strcmp(Parse_Mat_id[0],"id") != 0) ||
	  (Mat_id<0 || Mat_id> GP_Mesh.NumMat-1) ){
	printf("%s : %s \n",
	       "Error in Read_MPM_Materials()",
	       "Use this format -> (id=Integer) !!!");
	exit(0);
      }

      /* Set to default all it properties */
      Mat_GP.Fracture=false;
      Mat_GP.E = NAN;
      Mat_GP.mu = NAN;
      Mat_GP.Ceps = NAN;
      Mat_GP.Gf = NAN;

      /* Look for the curly brace { */
      if(strcmp(kwords[2],"{") == 0){
	/* Initial line */
	STATUS_LINE = fgets(Line_Material_Prop,
			    sizeof(Line_Material_Prop),
			    Sim_dat);
	if(STATUS_LINE == NULL){
	    printf("%s : %s \n",
		   "Error in Read_MPM_Materials()",
		   "Unspected EOF !!!");
	    exit(0);	
	}
	Aux_Mat_id = parse(Parse_Mat_Prop,Line_Material_Prop," =\t\n");
	while(STATUS_LINE != NULL){

	  if(Aux_Mat_id != 2){
	    printf("%s : %s \n",
		   "Error in Read_MPM_Materials()",
		   "Use this format -> Propertie = double !!!");
	    exit(0);
	  }

 	  if(strcmp(Parse_Mat_Prop[0],"Type") == 0){
	    strcpy(Mat_GP.Type,Parse_Mat_Prop[1]);
	  }
	  else if(strcmp(Parse_Mat_Prop[0],"E") == 0){
	    Mat_GP.E = atof(Parse_Mat_Prop[1]);
	  }
	  else if(strcmp(Parse_Mat_Prop[0],"mu") == 0){
	    Mat_GP.mu = atof(Parse_Mat_Prop[1]);
	  }
	  else if(strcmp(Parse_Mat_Prop[0],"Fracture") == 0){
	    if (strcmp(Parse_Mat_Prop[1],"TRUE") == 0){
	      Mat_GP.Fracture=true;	      
	    }
	    else if (strcmp(Parse_Mat_Prop[1],"FALSE") == 0){
	      Mat_GP.Fracture=false;
	    }
	    else{
	      printf("%s : %s \n",
		     "Error in Read_MPM_Materials()",
		     "TRUE/FALSE not detected");
	      exit(0);
	    }
	  }	  
	  else if(strcmp(Parse_Mat_Prop[0],"Ceps") == 0){
	    Mat_GP.Ceps = atof(Parse_Mat_Prop[1]);
	  }
	  else if(strcmp(Parse_Mat_Prop[0],"Gf") == 0){
	    Mat_GP.Gf = atof(Parse_Mat_Prop[1]);
	  }
	  else{
	    printf("%s : %s %s \n",
		   "Error in Read_MPM_Materials()",
		   "Undefined",Parse_Mat_Prop[0]);
	    printf("%s \n","Check the list of properties");
	    exit(0);
	  }
	  /* Read next line and check */
	  STATUS_LINE = fgets(Line_Material_Prop,
			      sizeof(Line_Material_Prop),
			      Sim_dat);
	  Aux_Mat_id = parse(Parse_Mat_Prop,Line_Material_Prop," =\t\n");
	  if(strcmp(Parse_Mat_Prop[0],"}") == 0){
	    break;
	  }
	}
	if(STATUS_LINE == NULL){
	printf("%s : %s \n",
	       "Error in Read_MPM_Materials()",
	       "you forget to put a } !!!");
	exit(0);	  
	}

	if(strcmp(Parse_Mat_Prop[0],"}") == 0){

	  /* Check fracture properties */
	  if( (Mat_GP.Fracture == true) &&
	      (Mat_GP.Ceps != Mat_GP.Ceps)){
	    printf("%s : %s \n",
		   "Error in Read_MPM_Materials()",
		   "Ceps is required for fracture !!!");
	    exit(0);
	  }
	  if( (Mat_GP.Fracture == true) &&
	      (Mat_GP.Gf != Mat_GP.Gf)){
	    printf("%s : %s \n",
	   "Error in Read_MPM_Materials()",
		   "Gf is required for fracture !!!");
	    exit(0);
	  }
	  
	  /* Transfere information */
	  GP_Mesh.Mat[Mat_id] = Mat_GP;
	  
	  break;
	}
      }
      else{
	printf("%s : %s \n",
	       "Error in Read_MPM_Materials()",
	       "Use this format -> MATERIALS_DEF (ID=integer) { !!!");
	exit(0);
      }
    }
  }
    
  /* Check the number of materials */
  if(CountMaterials != GP_Mesh.NumMat){
    printf("%s : %s %i %s %i %s \n",
	   "Error in Read_MPM_Materials()",
	   "Spected",GP_Mesh.NumMat, "materials, but",
	   CountMaterials,"where defined !!!");
    exit(0);
  }

  /* Close .dat file */
  /* Final message */
  printf("\t * End of read data file !!! \n");
  fclose(Sim_dat);

  return GP_Mesh.Mat;
}

/**********************************************************************/

Material * Read_MPM_Materials(char * Name_File, GaussPoint GP_Mesh)
/*
  Read list of materials :
  - Examples :

  MATERIALS_NUM number
  MATERIALS_DEF (ID=integer) {
  E=6.e9;
  mu=0.2;
  Ceps=1.5;
  Gf=0.00001;
  Info="Elastic";
  } 
*/
{
  /* Asign material library to an auxiliar variable */
  Material Mat_GP;

  /* Simulation file */
  FILE * Sim_dat;

  /* Index of the material */
  int Aux_Mat_id;
  char * Parse_Mat_id[MAXW] = {NULL};
  int Mat_id;

  /* Material properties */
  char Line_Material_Prop[MAXC] = {0};
  char * Parse_Mat_Prop[MAXW] = {NULL};

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  /* Auxiliar variable for status */
  char * STATUS_LINE;
  int CountMaterials = 0;

  /* Initial message */  
  puts("*************************************************");
  printf(" \t %s : \n\t %s \n",
	 "* Begin of materials properties in",
	 Name_File);
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    printf("%s : \n\t %s %s",
	   "Error in Read_MPM_Materials()",
	   "Incorrect lecture of",
	   Name_File);
    exit(0);
  }

  /* Allocate table with the material */
  GP_Mesh.Mat = malloc(GP_Mesh.NumMat*sizeof(Material));
  if(GP_Mesh.Mat == NULL){
    printf("%s : %s \n",
	   "Error in Read_MPM_Materials()",
	   "Memory error for Mat");
    exit(0);
  }
  
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");
    if (nkwords < 0){
      printf("%s : %s \n",
	     "Error in Read_MPM_Materials()",
	     "Parser failed");
      exit(0);
    }

    if (strcmp(kwords[0],"MATERIALS_DEF") == 0 ){

      /* Count the number of materials */
      ++CountMaterials;

      /* Read the index of the material */
      Aux_Mat_id = parse (Parse_Mat_id, kwords[1],"(=)");
      Mat_id = atoi(Parse_Mat_id[1]);
      if( (Aux_Mat_id != 2) ||
	  (strcmp(Parse_Mat_id[0],"ID") != 0) ||
	  (Mat_id<0 || Mat_id> GP_Mesh.NumMat-1) ){
	printf("%s : %s \n",
	       "Error in Read_MPM_Materials()",
	       "Use this format -> (ID=Integer) !!!");
	exit(0);
      }

      /* Set to NAN all it properties */
      
      Mat_GP.E = NAN;
      Mat_GP.mu = NAN;
      Mat_GP.Ceps = NAN;
      Mat_GP.Gf = NAN;
      GP_Mesh.Mat[Mat_id] = Mat_GP;

      /* Look for the curly brace { */
      if(strcmp(kwords[2],"{") == 0){
	/* Initial line */
	STATUS_LINE = fgets(Line_Material_Prop,
			    sizeof(Line_Material_Prop),
			    Sim_dat);
	if(STATUS_LINE == NULL){
	    printf("%s : %s \n",
		   "Error in Read_MPM_Materials()",
		   "Unspected EOF !!!");
	    exit(0);	
	}
	Aux_Mat_id = parse(Parse_Mat_Prop,Line_Material_Prop," =\t\n");
	while(STATUS_LINE != NULL){

	  if(Aux_Mat_id != 2){
	    printf("%s : %s \n",
		   "Error in Read_MPM_Materials()",
		   "Use this format -> Propertie = double !!!");
	    exit(0);
	  }

	  if(strcmp(Parse_Mat_Prop[0],"E") == 0){
	    GP_Mesh.Mat[Mat_id].E = atof(Parse_Mat_Prop[1]);
	  }
	  else if(strcmp(Parse_Mat_Prop[0],"mu") == 0){
	    GP_Mesh.Mat[Mat_id].mu = atof(Parse_Mat_Prop[1]);
	  }
	  else if(strcmp(Parse_Mat_Prop[0],"Ceps") == 0){
	    GP_Mesh.Mat[Mat_id].Ceps = atof(Parse_Mat_Prop[1]);
	  }
	  else if(strcmp(Parse_Mat_Prop[0],"Gf") == 0){
	    GP_Mesh.Mat[Mat_id].Gf = atof(Parse_Mat_Prop[1]);
	  }
	  else if(strcmp(Parse_Mat_Prop[0],"Type") == 0){
	    strcpy(GP_Mesh.Mat[Mat_id].Type,Parse_Mat_Prop[1]);
	  }
	  else{
	    printf("%s : %s %s \n",
		   "Error in Read_MPM_Materials()",
		   "Undefined",Parse_Mat_Prop[0]);
	    printf("%s \n","Check the list of properties");
	    exit(0);
	  }
	  
	  /* Read next line and check */
	  STATUS_LINE = fgets(Line_Material_Prop,
			      sizeof(Line_Material_Prop),
			      Sim_dat);
	  Aux_Mat_id = parse(Parse_Mat_Prop,Line_Material_Prop," =\t\n");
	  if(strcmp(Parse_Mat_Prop[0],"}") == 0){
	    break;
	  }
	}
	if(STATUS_LINE == NULL){
	printf("%s : %s \n",
	       "Error in Read_MPM_Materials()",
	       "you forget to put a } !!!");
	exit(0);	  
	}

	if(strcmp(Parse_Mat_Prop[0],"}") == 0){
	  break;
	}
      }
      else{
	printf("%s : %s \n",
	       "Error in Read_MPM_Materials()",
	       "Use this format -> MATERIALS_DEF (ID=integer)  !!!");
	exit(0);
      }

    }
  }

  /* Check the number of materials */
  if(CountMaterials != GP_Mesh.NumMat){
    printf("%s : %s %i %s %i %s \n",
	   "Error in Read_MPM_Materials()",
	   "Spected",GP_Mesh.NumMat, "materials, but",
	   CountMaterials,"where defined !!!");
    exit(0);
  }

  /* Close .dat file */
  /* Final message */
  printf("\t * End of read data file !!! \n");
  fclose(Sim_dat);

  return GP_Mesh.Mat;
}

/**********************************************************************/

