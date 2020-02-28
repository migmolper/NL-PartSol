#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "grams.h"

/**********************************************************************/

Material * GramsMaterials(char * Name_File, GaussPoint GP_Mesh, int GPxElement)
/*
GramsMaterials (Particles=route.txt) {
               Id=0
               Type=LE
	       rho=20
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
  Material * Mat_Table;
  
  /* Simulation file */
  FILE * Sim_dat;

  /* Parser num chars */
  int Num_words_parse;

  /* Index of the material */
  int Aux_Mat_id;
  char * Parse_Mat_id[MAXW] = {NULL};

  /* Material properties */
  char Line_Material_Prop[MAXC] = {0};
  char * Parse_Mat_Prop[MAXW] = {NULL};

  /* Parse file name with the list of nodes */
  char * Name_File_Copy = malloc(strlen(Name_File)); 
  char * Name_Parse[MAXW] = {NULL};
  char Route_Nodes[MAXC] = {0};
  char FileNodesRoute[MAXC];

  /* Array */
  int Num_Nodes;
  ChainPtr Chain_Nodes = NULL;
  int * Array_Nodes;

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
    fprintf(stderr,"%s : \n\t %s %s",
	   "Error in GramsMaterials()",
	   "Incorrect lecture of",
	   Name_File);
    exit(0);
  }

  /* Generate route */
  strcpy(Name_File_Copy, Name_File);
  Num_words_parse = parse(Name_Parse,Name_File_Copy,"(/)");
  strcat(Route_Nodes,"./");
  for(int i = 0 ; i<Num_words_parse-1 ; i++){
    strcat(Route_Nodes, Name_Parse[i]);
    strcat(Route_Nodes,"/");
  }

  /* Allocate table with the material */
  Mat_Table = (Material * )malloc(GP_Mesh.NumberMaterials*sizeof(Material));
  if(Mat_Table == NULL){
    fprintf(stderr,"%s : %s \n",
	   "Error in GramsMaterials()",
	   "Memory error for Mat");
    exit(0);
  }
  
  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n\t");
    if (nkwords < 0){
      fprintf(stderr,"%s : %s \n",
	     "Error in GramsMaterials ()",
	     "Parser failed");
      exit(0);
    }

    if ((nkwords>0) &&
	(strcmp(kwords[0],"GramsMaterials") == 0 )){

      /* Count the number of materials */
      ++CountMaterials;

      /* Read the index of the material */
      Aux_Mat_id = parse (Parse_Mat_id, kwords[1],"(=)");
      if( (Aux_Mat_id != 2) ||
	  (strcmp(Parse_Mat_id[0],"Particles") != 0)){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsMaterials()",
	       "Use this format -> (Particles=route.txt) !!!");
	exit(0);
      }

      /* Read file with the nodes */
      sprintf(FileNodesRoute,"%s%s",Route_Nodes,Parse_Mat_id[1]);

      /* Get an array with the nodes */
      Chain_Nodes = File2Chain(FileNodesRoute);
      Num_Nodes = LenghtChain(Chain_Nodes);
      Array_Nodes = ChainToArray(Chain_Nodes,Num_Nodes);
      FreeChain(Chain_Nodes);

      /* Set to default all it properties */
      Mat_GP.Id=0;      
      Mat_GP.Fracture=false;
      Mat_GP.rho = NAN;
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
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsMaterials()",
		   "Unspected EOF !!!");
	    exit(0);	
	}
	Aux_Mat_id = parse(Parse_Mat_Prop,Line_Material_Prop," =\t\n");
	if(strcmp(Parse_Mat_Prop[0],"}") == 0){

	  /* Check fracture properties */
	  if( (Mat_GP.Fracture == true) &&
	      (Mat_GP.Ceps != Mat_GP.Ceps)){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsMaterials()",
		   "Ceps is required for fracture !!!");
	    exit(0);
	  }
	  if( (Mat_GP.Fracture == true) &&
	      (Mat_GP.Gf != Mat_GP.Gf)){
	    fprintf(stderr,"%s : %s \n",
	   "Error in GramsMaterials()",
		   "Gf is required for fracture !!!");
	    exit(0);
	  }
	  
	  /* Transfere information */
	  Mat_Table[Mat_GP.Id] = Mat_GP;
	  
	  break;
	}
	while(STATUS_LINE != NULL){

	  if(Aux_Mat_id != 2){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsMaterials()",
		   "Use this format -> Propertie = double !!!");
	    exit(0);
	  }

	  if(strcmp(Parse_Mat_Prop[0],"Id") == 0){
	    Mat_GP.Id = atoi(Parse_Mat_Prop[1]);
	    if(Mat_GP.Id >= GP_Mesh.NumberMaterials){
	      fprintf(stderr,"%s : %s %i !!! \n",
		      "Error in GramsMaterials()",
		      "Id should go from 0 to",
		      GP_Mesh.NumberMaterials-1);
	      exit(1);
	    }
	    printf("\t -> %s : %i \n",
		   "Index of the material",Mat_GP.Id);

	    for(int i = 0 ; i<Num_Nodes ; i++){
	      for(int j = 0 ; j<GPxElement ; j++){
		GP_Mesh.MatIdx[Array_Nodes[i]*GPxElement+j] = Mat_GP.Id;
	      }
	    }
	  }
 	  else if(strcmp(Parse_Mat_Prop[0],"Type") == 0){
	    strcpy(Mat_GP.Type,Parse_Mat_Prop[1]);
	    printf("\t -> %s : %s \n",
		   "Law",Parse_Mat_Prop[1]);
	  }
	  else if( strcmp(Parse_Mat_Prop[0],"rho") == 0 ){
	    Mat_GP.rho = atof(Parse_Mat_Prop[1]);
	    printf("\t -> %s : %f \n",
		   "Density",Mat_GP.rho);
	  }
	  else if(strcmp(Parse_Mat_Prop[0],"E") == 0){
	    Mat_GP.E = atof(Parse_Mat_Prop[1]);
	    printf("\t -> %s : %f \n",
		   "Elastic modulus",Mat_GP.E);
	  }
	  else if(strcmp(Parse_Mat_Prop[0],"mu") == 0){
	    Mat_GP.mu = atof(Parse_Mat_Prop[1]);
	    printf("\t -> %s : %f \n",
		   "Poisson modulus",Mat_GP.mu);
	  }
	  else if(strcmp(Parse_Mat_Prop[0],"Fracture") == 0){
	    if (strcmp(Parse_Mat_Prop[1],"TRUE") == 0){
	      Mat_GP.Fracture=true;
	      printf("\t -> %s : %s \n",
		     "Fracture","ON");
	    }
	    else if (strcmp(Parse_Mat_Prop[1],"FALSE") == 0){
	      Mat_GP.Fracture=false;
	      printf("\t -> %s : %s \n",
		     "Fracture","OFF");
	    }
	    else{
	      fprintf(stderr,"%s : %s \n",
		     "Error in GramsMaterials()",
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
	    fprintf(stderr,"%s : %s %s \n",
		   "Error in GramsMaterials()",
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
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsMaterials()",
	       "you forget to put a } !!!");
	exit(0);	  
	}

	if(strcmp(Parse_Mat_Prop[0],"}") == 0){

	  /* Check fracture properties */
	  if( (Mat_GP.Fracture == true) &&
	      (Mat_GP.Ceps != Mat_GP.Ceps)){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsMaterials()",
		   "Ceps is required for fracture !!!");
	    exit(0);
	  }
	  if( (Mat_GP.Fracture == true) &&
	      (Mat_GP.Gf != Mat_GP.Gf)){
	    fprintf(stderr,"%s : %s \n",
	   "Error in GramsMaterials()",
		   "Gf is required for fracture !!!");
	    exit(0);
	  }
	  
	  /* Transfere information */
	  Mat_Table[Mat_GP.Id] = Mat_GP;
	  
	  break;
	}
      }
      else{
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsMaterials()",
	       "Use this format -> (Particles=route.txt) { !!!");
	exit(0);
      }
    }

    /* Free array nodes */
    free(Array_Nodes);
  }
    
  /* Check the number of materials */
  if(CountMaterials != GP_Mesh.NumberMaterials){
    fprintf(stderr,"%s : %s %i %s %i %s \n",
	   "Error in GramsMaterials()",
	   "Spected",GP_Mesh.NumberMaterials, "materials, but",
	   CountMaterials,"where defined !!!");
    exit(0);
  }

  /* Close .dat file */
  /* Final message */
  printf("\t * End of read data file !!! \n");
  fclose(Sim_dat);

  return Mat_Table;
}

/**********************************************************************/
