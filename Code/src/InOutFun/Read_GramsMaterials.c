#include "grams.h"

/**********************************************************************/

Material * GramsMaterials(char * Name_File, GaussPoint GP_Mesh, int GPxElement)
/*
GramsMaterials (Particles=route.txt) {
               Id=0
               Type=LE
	       Cel=1
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
  int Ndim = NumberDimensions;
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


  bool Is_Id;
  bool Is_Cel;
  bool Is_rho;
  bool Is_E;
  bool Is_mu;
  bool Is_thickness;
  bool Is_Ceps;
  bool Is_Gf;
  bool Is_ft;
  bool Is_heps;
  bool Is_Wc;

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
  Mat_Table = (Material *)malloc(GP_Mesh.NumberMaterials*sizeof(Material));
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
      printf("\t -> %s : %s \n",
	     "Material points",FileNodesRoute);

      /* Get an array with the nodes */
      Chain_Nodes = File2Chain(FileNodesRoute);
      Num_Nodes = get_Lenght_Set(Chain_Nodes);
      Array_Nodes = Set_to_Pointer(Chain_Nodes,Num_Nodes);
      free_Set(Chain_Nodes);

      /* Id of the material */
      Is_Id = false;
      Mat_GP.Id=-1;
      /* Celerity */
      Is_Cel = false;
      Mat_GP.Cel = NAN;
      /* Density */      
      Is_rho = false;
      Mat_GP.rho = NAN;
      /* Linear elastic parameters */
      Is_E = false;
      Mat_GP.E = NAN;
      Is_mu = false;
      Mat_GP.mu = NAN;
      /* Thickness of the mateial (2D) */
      Is_thickness = false;
      Mat_GP.thickness = 1;
      /* Fracture module */
      Mat_GP.Eigenerosion = false;
      Mat_GP.Eigensoftening = false;
      /* Parameters for Eigenerosion */
      Is_Ceps = false;
      Mat_GP.Ceps = NAN;
      Is_Gf = false;
      Mat_GP.Gf = NAN;
      /* Parameters for Eigensoftening */
      Is_ft = false;
      Mat_GP.ft = NAN;
      Is_heps = false;
      Mat_GP.heps = NAN;
      Is_Wc = false;
      Mat_GP.Wc = NAN;    

      /* Look for the curly brace { */
      if(strcmp(kwords[2],"{") == 0){

	/* Initial line */
	STATUS_LINE = fgets(Line_Material_Prop,
			    sizeof(Line_Material_Prop),
			    Sim_dat);
	if(STATUS_LINE == NULL){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsMaterials()","Unspected EOF !!!");
	    exit(0);	
	}
	Aux_Mat_id = parse(Parse_Mat_Prop,Line_Material_Prop," =\t\n");
	if(strcmp(Parse_Mat_Prop[0],"}") == 0){
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsMaterials()","Undefined material");
	  exit(0);
	}
	while(STATUS_LINE != NULL){

	  if(Aux_Mat_id != 2){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsMaterials()","Format -> Propertie=value");
	    fprintf(stderr,"%s (%i) : %s , %s\n",
		    "Printed",Aux_Mat_id,Parse_Mat_Prop[0],Parse_Mat_Prop[1]);
	    exit(0);
	  }
	  /**************************************************/
	  if(strcmp(Parse_Mat_Prop[0],"Id") == 0){
	    Is_Id = true;
	    Mat_GP.Id = atoi(Parse_Mat_Prop[1]);
	    if(Mat_GP.Id >= GP_Mesh.NumberMaterials){
	      fprintf(stderr,"%s : %s %i !!! \n",
		      "Error in GramsMaterials()","Id should go from 0 to",
		      GP_Mesh.NumberMaterials-1);
	      exit(1);
	    }
	    for(int i = 0 ; i<Num_Nodes ; i++){
	      for(int j = 0 ; j<GPxElement ; j++){
	    	GP_Mesh.MatIdx[Array_Nodes[i]*GPxElement+j] = Mat_GP.Id;
	      }
	    }
	  }
	  /**************************************************/
 	  else if(strcmp(Parse_Mat_Prop[0],"Type") == 0){
	    strcpy(Mat_GP.Type,Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Cel") == 0){
	    Is_Cel = true;
	    Mat_GP.Cel = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"rho") == 0){
	    Is_rho = true;
	    Mat_GP.rho = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"E") == 0){
	    Is_E = true;
	    Mat_GP.E = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"mu") == 0){
	    Is_mu = true;
	    Mat_GP.mu = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/	    
	  else if(strcmp(Parse_Mat_Prop[0],"thickness") == 0){
	    Is_thickness = true;
	    Mat_GP.thickness = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/	    
	  else if(strcmp(Parse_Mat_Prop[0],"Fracture") == 0){
	    if (strcmp(Parse_Mat_Prop[1],"Eigenerosion") == 0){
	      Mat_GP.Eigenerosion = true;
	    }
	    else if (strcmp(Parse_Mat_Prop[1],"Eigensoftening") == 0){
	      Mat_GP.Eigensoftening=true;
	    }
	    else{
	      fprintf(stderr,"%s : %s \n",
		     "Error in GramsMaterials()",
		     "Options -> Eigenerosion/Eigensoftening");
	      exit(0);
	    }
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Ceps") == 0){
	    Is_Ceps = true;
	    Mat_GP.Ceps = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Gf") == 0){
	    Is_Gf = true;
	    Mat_GP.Gf = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"ft") == 0){
	    Is_ft = true;
	    Mat_GP.ft = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"heps") == 0){
	    Is_heps = true;
	    Mat_GP.heps = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
	  else if(strcmp(Parse_Mat_Prop[0],"Wc") == 0){
	    Is_Wc = true;
	    Mat_GP.Wc = atof(Parse_Mat_Prop[1]);
	  }
	  /**************************************************/
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
	  /**************************************************/
	  if(Is_Id){
	    printf("\t \t -> %s : %i \n","Index of the material",Mat_GP.Id);
	  }
	  /**************************************************/
	  if(strcmp(Mat_GP.Type,"SR") == 0){ /* Solid rigid material */
	    if(Is_rho){
	      Mat_GP.Cel = 0;
	      printf("\t -> %s \n","Solid rigid material");
	      printf("\t \t -> %s : %f \n","Density",Mat_GP.rho);
	    }
	    else{
	      fprintf(stderr,"%s : %s \n",
		      "Error in GramsMaterials()",
		      "Density is required for LE materials");
	      exit(0);
	    }
	  }
	  else if(strcmp(Mat_GP.Type,"LE") == 0){ /* Linear elastic parameters */
	    if(Is_rho & Is_Cel && Is_E && Is_mu){
	      printf("\t -> %s \n","Linear elastic material");
	      printf("\t \t -> %s : %f \n","Celerity",Mat_GP.Cel);
	      printf("\t \t -> %s : %f \n","Density",Mat_GP.rho);
	      printf("\t \t -> %s : %f \n","Elastic modulus",Mat_GP.E);
	      printf("\t \t -> %s : %f \n","Poisson modulus",Mat_GP.mu);
	    }
	    else{
	      fprintf(stderr,"%s : %s \n",
		      "Error in GramsMaterials()",
		      "Rho, Cel, E and mu required for LE materials");
	      exit(0);
	    }
	  }
	  else{
	    fprintf(stderr,"%s : %s \n",
		    "Error in GramsMaterials()","Unrecognized kind of material");
	    exit(0);
	  }

	  /**************************************************/
	  if(Is_thickness){
	    if(Ndim != 2){
	      fprintf(stderr,"%s : %s\n",
		      "Error in GramsMaterials()",
		      "thickness is only for plain strain cases");
	      exit(0);
	    }
	    printf("\t \t -> %s : %f \n","thickness",Mat_GP.thickness);
	  }
	  else{
	    printf("\t \t -> %s : %f \n","thickness",Mat_GP.thickness);
	  }
	  /**************************************************/
	  if(Mat_GP.Eigenerosion){ /* Check eigenerosion properties */	  
	    if(Is_Ceps && Is_Gf){
	      printf("\t \t -> %s : %s \n","Fracture","Eigenerosion");
	      printf("\t \t -> %s : %f \n","Normalizing constant",Mat_GP.Ceps);
	      printf("\t \t -> %s : %f \n","Failure energy",Mat_GP.Gf);
	    }
	    else{
	      fprintf(stderr,"%s : %s \n",
		      "Error in GramsMaterials()",
		      "Ceps and Gf required for Eigenerosion");
	      exit(0);
	    }
	  }
	  /**************************************************/ 
	  if(Mat_GP.Eigensoftening){ /* Check eigensoftening properties */
	    if(Is_Ceps && Is_ft && Is_heps && Is_Wc){
	      printf("\t \t -> %s : %s \n","Fracture","Eigensoftening");
	      printf("\t \t -> %s : %f \n","Normalizing constant",Mat_GP.Ceps);
	      printf("\t \t -> %s : %f \n","Tensile strengt",Mat_GP.ft);
	      printf("\t \t -> %s : %f \n","Bandwidth Bazant",Mat_GP.heps);
	      printf("\t \t -> %s : %f \n","Critical opening",Mat_GP.Wc);
	    }
	    else{
	      fprintf(stderr,"%s : %s\n",
		      "Error in GramsMaterials()",
		      "ft, heps and Wc required for Eigensoftening");
	      exit(0);
	    }
	  }
	  /**************************************************/
	  
	  /* Transfere information */
	  Mat_Table[Mat_GP.Id] = Mat_GP;
	  
	  /* break; */
	}
      }
      else{
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsMaterials()",
	       "Use this format -> (Particles=route.txt) { !!!");
	exit(0);
      }
    }

    /* /\* Free array nodes *\/ */
    /* free(Array_Nodes); */
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
