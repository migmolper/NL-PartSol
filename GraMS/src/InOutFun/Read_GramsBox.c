#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "../GRAMS/grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

/**********************************************************************/

Mesh GramsBox(char * Name_File)
/*
    
         top = 2
         *-----*
         |     |
left = 3 |     | right = 1
         |     | 
         *-----*
        bottom = 0

GramsBox (Type=GID,File=FEM_Mesh.msh) { 
  left = GramsBoundary {
        BcDirichlet U NULL
 	BcDirichlet V NULL	
  }
  right = GramsBoundary {
  	BcDirichlet U NULL
   	BcDirichlet V NULL	
  }
  top = GramsBoundary {
      	BcDirichlet U NULL
    	BcDirichlet V NULL
  }
  bottom = GramsBoundary {
  	BcDirichlet U NULL
    	BcDirichlet V NULL
  }
}
*/
{
  /* Asign material library to an auxiliar variable */
  Mesh FEM_Mesh;

  /* Simulation file */
  FILE * Sim_dat;

  int Num_Name;
  char * Name_File_Copy = malloc(strlen(Name_File)); 
  char * Name_Parse[MAXW] = {NULL};
  char Route_Mesh[MAXC] = {0};

  /* Index of the material */
  int Aux_Mesh_id;
  char * Parse_Mesh_id[MAXW] = {NULL};
  int Mesh_id;

  /* Parser num chars */
  int Num_words_line;
  
  /* Parse line of GramsBox */
  char Line_GramsBox[MAXC] = {0};
  char * Parse_GramsBox[MAXW] = {NULL};

  /* Parse lines of GramsBoundary */
  char Line_GramsBoundary[MAXC] = {0};
  char * Parse_GramsBoundary[MAXW] = {NULL};

  /* Parse boundary conditions lines */
  char Line_BCC[MAXC] = {0};
  char * Parse_BCC[MAXW] = {NULL};

  /* Special variables for line-reading */
  char line[MAXC] = {0}; /* Variable for reading the lines in the files */
  char * kwords[MAXW] = {NULL}; /* Variable to store the parser of a line */
  int nkwords; /* Number of element in the line , just for check */

  /* Auxiliar variable for status */
  char * STATUS_LINE;

  /* Auxiliar table to store the boundaries */
  int * NodesBound_aux;
  int * CounterNodesBound;
  int NumNodesBound = 0;
  
  /* Count the number of elements that share this node */
  ChainPtr Elem_Conn; /* Loop over the connectivity chain */
  int Repeat_Nod; /* Counter */
  /* Variables that fills the boundaries nodes */
  int aux_RIGHT = 0; 
  int aux_TOP = 0;
  int aux_LEFT = 0;
  int aux_BOTTOM = 0;
  /* Set to zero X and Y min and max values of the mesh */
  double MAX_X, MAX_Y, MIN_X, MIN_Y;
  
  /* Boundaries labels */
  char * BoundLabels [4] = {"BOTTOM", "RIGHT", "TOP", "LEFT"};
  

  /* Initial message */
  puts("*************************************************");
  puts(" Generate the MPM mesh");
  printf(" \t %s : \n\t %s \n",
	 "* Begin of read mesh from in",Name_File);
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	   "Error in GramsBox()",
	   "Incorrect lecture of",
	   Name_File);
    exit(0);
  }

  /* Read the file line by line */
  while( fgets(Line_GramsBox, sizeof(Line_GramsBox), Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    Num_words_line = parse (Parse_GramsBox, Line_GramsBox," \n\t");
    if (Num_words_line < 0){
      fprintf(stderr,"%s : %s \n",
	     "Error in GramsBox ()",
	     "Parser failed");
      exit(0);
    }

    if (strcmp(Parse_GramsBox[0],"GramsBox") == 0 ){

      /* Read the index of the material */
      Num_words_line = parse (Parse_Mesh_id, Parse_GramsBox[1],"(=,)");
      Mesh_id = atoi(Parse_Mesh_id[1]);
      if( (Num_words_line != 4) ||
	  (strcmp(Parse_Mesh_id[0],"Type") != 0) ||
	  (strcmp(Parse_Mesh_id[2],"File") != 0) ){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsBox()",
	       "Use this format -> (Type=str,File=str) !!!");
	exit(0);
      }
      
      if(strcmp(Parse_Mesh_id[1],"GID") == 0){
	strcpy(Name_File_Copy, Name_File);
	FEM_MeshFileName = Parse_Mesh_id[3];
	Num_words_line = parse(Name_Parse,Name_File_Copy,"(/)");
	for(int i = 0 ; i<Num_words_line-1 ; i++){
	  strcat(Route_Mesh, Name_Parse[i]);
	  strcat(Route_Mesh,"/");
	}
	strcat(Route_Mesh,FEM_MeshFileName);
	free(Name_File_Copy);	  
	FEM_Mesh = ReadGidMesh(Route_Mesh);
	exit(0);
      }
      else{
	fprintf(stderr,"%s : %s %s \n",
		"Error in GramsBox(Type=*, )",
		"Unrecognized mesh Type",
		Parse_Mesh_id[1]);
      }

      /* Generate nodal connectivity of the mesh :
       list of elements near to a node */
      GetNodalConnectivity(FEM_Mesh);

      /* Initialize GPs connectivity of each element */
      FEM_Mesh.GPsElements =
	(ChainPtr *)malloc(FEM_Mesh.NumElemMesh*sizeof(ChainPtr));
      if(FEM_Mesh.GPsElements == NULL){
	printf("%s : %s \n",
	       "GetNodalConnectivity",
	       "Memory error for GPsElements");
	exit(0);
      }
      for(int i = 0 ; i<FEM_Mesh.NumElemMesh ; i++){
	FEM_Mesh.GPsElements[i] = NULL;
      }  

      /* Define boundary conditions and continue reading the file */

      /* Kind of domain */
      strcpy(FEM_Mesh.Bounds.Info,"Box");      
  
      /* Asign the number of boundaries for a square domain */
      FEM_Mesh.Bounds.NumBounds = 4;
      FEM_Mesh.Bounds.BCC_i = (Load *)Allocate_Array(1,sizeof(Load));

      /* Fill some parameters */
      for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++){
	/* Set to zero the number of nodes in the boundary of the mesh */
	FEM_Mesh.Bounds.BCC_i[i].NumNodes = 0;
	/* Set the labels for a square domain */
	strcpy(FEM_Mesh.Bounds.BCC_i[i].Info,BoundLabels[i]);
      }

      /* Allocate an array of zeros to assign a 1 to those nodes in the boundary */
      NodesBound_aux =
	(int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh,sizeof(int));
      CounterNodesBound =
	(int *)Allocate_ArrayZ(4,sizeof(int));
  
      /* Iterate over the nodes to fin the nodes in the boundary */

      /* Get the max values of the boundary */
      MAX_X = FEM_Mesh.Coordinates.nM[0][0];
      MAX_Y = FEM_Mesh.Coordinates.nM[0][1];
      MIN_X = FEM_Mesh.Coordinates.nM[0][0];
      MIN_Y = FEM_Mesh.Coordinates.nM[0][1];
      
      for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){

	/* Get the max values of the boundary */
	MAX_X = MAXVAL(MAX_X,FEM_Mesh.Coordinates.nM[i][0]);
	MAX_Y = MAXVAL(MAX_Y,FEM_Mesh.Coordinates.nM[i][1]);
	MIN_X = MINVAL(MIN_X,FEM_Mesh.Coordinates.nM[i][0]);
	MIN_Y = MINVAL(MIN_Y,FEM_Mesh.Coordinates.nM[i][1]);
	
	/* Set the counter to zero */
	Repeat_Nod = 0;
	
	/* Loop over the connectivity mesh */
	for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){
	  Elem_Conn = FEM_Mesh.Connectivity[j];
	  while(Elem_Conn != NULL){
	    if((Elem_Conn->I) == i){
	      Repeat_Nod++;
	    }
	    Elem_Conn = Elem_Conn->next;
	  }
	}
	
	/* Add this element to the boundary */
	if (Repeat_Nod < 4){
	  NodesBound_aux[i] = 1;
	  NumNodesBound++;
	}
      }
      
      /* Count the number of nodes in each boundarie */
      for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
	if(NodesBound_aux[i] == 1){

	  if(FEM_Mesh.Coordinates.nM[i][1] == MIN_Y){
	    /* Count number of nodes in the bottom */
	    CounterNodesBound[0]++;
	  }
	  if(FEM_Mesh.Coordinates.nM[i][0] == MAX_X){
	    /* Count number of nodes in the right */
	    CounterNodesBound[1]++;
	  }
	  if(FEM_Mesh.Coordinates.nM[i][1] == MAX_Y){
	    /* Count the number of nodes in the top */
	    CounterNodesBound[2]++;
	  }
	  if(FEM_Mesh.Coordinates.nM[i][0] == MIN_X){
	    /* Count the number of nodes in the left */
	    CounterNodesBound[3]++;
	  }
	    
	}
      }
      
      /* Allocate the arrays with the boundary nodes */
      for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++){
	FEM_Mesh.Bounds.BCC_i[i].NumNodes = CounterNodesBound[i];
	FEM_Mesh.Bounds.BCC_i[i].Nodes =
	  (int *)Allocate_ArrayZ(CounterNodesBound[i],sizeof(int));
      }

      /* Free data */
      free(CounterNodesBound);
      
      /* Fill the arrays  */
      for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
	if(NodesBound_aux[i] == 1){

	  if(FEM_Mesh.Coordinates.nM[i][1] == MIN_Y){
	    FEM_Mesh.Bounds.BCC_i[0].Nodes[aux_BOTTOM] = i;
	    aux_BOTTOM++;
	  }
	  if(FEM_Mesh.Coordinates.nM[i][0] == MAX_X){
	    FEM_Mesh.Bounds.BCC_i[1].Nodes[aux_RIGHT] = i;
	    aux_RIGHT++;
	  }
	  if(FEM_Mesh.Coordinates.nM[i][1] == MAX_Y){
	    FEM_Mesh.Bounds.BCC_i[2].Nodes[aux_TOP] = i;
	    aux_TOP++;
	  }
	  if(FEM_Mesh.Coordinates.nM[i][0] == MIN_X){
	    FEM_Mesh.Bounds.BCC_i[3].Nodes[aux_LEFT] = i;
	    aux_LEFT++;
	  }	    
	}
      }
    
      /* Free data */ 
      free(NodesBound_aux);

      /* Keep reading */
      exit(0);
      
      /* Look for the curly brace { */
      if(strcmp(Parse_GramsBox[2],"{") == 0){

	/* Initial line */
	STATUS_LINE =
	  fgets(Line_GramsBoundary,sizeof(Line_GramsBoundary),Sim_dat);
	if(STATUS_LINE == NULL){
	    fprintf(stderr,"%s : %s \n",
		   "Error in GramsBox()",
		   "Unspected EOF !!!");
	    exit(0);	
	}
	Num_words_line =
	  parse(Parse_GramsBoundary,Line_GramsBoundary," =\t\n");
	if((strcmp(Parse_GramsBox[3],"}") == 0) ||
	   (strcmp(Parse_GramsBoundary[0],"}") == 0)){
	  /* No boundary conditions applied */
	  break;
	}
	while(STATUS_LINE != NULL){

	  if((Num_words_line == 4) &&
	     (strcmp(Parse_GramsBox[1],"GramsBoundary") == 0)){ 
	  }

 	  if((strcmp(Parse_GramsBoundary[0],"left") == 0) &&
	     (strcmp(Parse_GramsBoundary[1],"GramsBoundary") == 0) &&
	     (strcmp(Parse_GramsBoundary[2],"{") == 0)){

	    STATUS_LINE =
	      fgets(Line_BCC,sizeof(Line_BCC),Sim_dat);
	    if(STATUS_LINE == NULL){
	      fprintf(stderr,"%s : %s \n",
		      "Error in GramsBox()",
		      "Unspected EOF !!!");
	      exit(0);	
	    }
	    Aux_Mesh_id = parse(Parse_BCC,Line_BCC," =\t\n");

	    strcpy(FEM_Mesh.Bounds.BCC_i[3].Info,Parse_BCC[1]);

	    nparam = parse (param,kwords[3],"={,}\n");
	    if( (strcmp(param[0],"DIR") == 0 ) && (nparam >= 2) ){
	      /* Number of dimensions of the BCC */
	      FEM_BCC.BCC_i[IndexBoundary].Dim = nparam - 1;
	      /* Direction of the BCC */
	      FEM_BCC.BCC_i[IndexBoundary].Dir =
		(int *)Allocate_Array(FEM_BCC.BCC_i[IndexBoundary].Dim,sizeof(int));
	      /* Fill the direction of the BCC */
	      for(int j = 0 ; j<FEM_BCC.BCC_i[IndexBoundary].Dim ; j++){
		FEM_BCC.BCC_i[IndexBoundary].Dir[j] = atoi(param[j+1]);
	      }
	    }

	    /* Read the curve to impose the boundary condition */
	    nparam = parse (param,kwords[4],"={,}\n");
	    if( (strcmp(param[0],"CURVE") == 0) || (nparam >= 2) ){

	      /* Alocate the table of curves */
	      FEM_BCC.BCC_i[IndexBoundary].Value =
		(Curve *)Allocate_Array(FEM_BCC.BCC_i[IndexBoundary].Dim,
					sizeof(Curve));
	      /* Fill the curve table */
	      for(int j = 0 ; j<FEM_BCC.BCC_i[IndexBoundary].Dim ; j++){
		if(strcmp(param[j+1],"NULL") != 0)
		  FEM_BCC.BCC_i[IndexBoundary].Value[j] = ReadCurve(param[j+1]);
	      }
	    }
	    
	  }
	  else if((strcmp(Parse_GramsBoundary[0],"right") == 0) &&
		  (strcmp(Parse_GramsBoundary[1],"GramsBoundary") == 0) &&
		  (strcmp(Parse_GramsBoundary[2],"{") == 0)){
	    exit(0);
	  }
	  else if((strcmp(Parse_GramsBoundary[0],"top") == 0) &&
		  (strcmp(Parse_GramsBoundary[1],"GramsBoundary") == 0) &&
		  (strcmp(Parse_GramsBoundary[2],"{") == 0)){
	    exit(0);
	  }
	  else if((strcmp(Parse_GramsBoundary[0],"bottom") == 0) &&
		  (strcmp(Parse_GramsBoundary[1],"GramsBoundary") == 0) &&
		  (strcmp(Parse_GramsBoundary[2],"{") == 0)){
	    exit(0);
	  }
	  else{
	    fprintf(stderr,"%s : %s \n",
		    "Error in GramsBox()", "Use this format");
	    fprintf(stderr,
		    "\t %s \n \t %s \n \t %s \n \t %s \n \t %s \n \t %s \n \t %s",
		    "GramsBox (Type=srt,File=str) {",
		    "GfsBox {", 
		    "left = GramsBoundary { }",
		    "right = GramsBoundary { }",
		    "top = GramsBoundary { }",
		    "bottom = GramsBoundary { }",
		    "}");
	    exit(0);
	  }
	  /* Read next line and check */
	  STATUS_LINE = fgets(Line_GramsBox,
			      sizeof(Line_GramsBox),
			      Sim_dat);
	  Aux_Mesh_id = parse(Parse_GramsBox,Line_GramsBox," =\t\n");
	  if(strcmp(Parse_GramsBox[0],"}") == 0){
	    break;
	  }
	}
	if(STATUS_LINE == NULL){
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsBox()",
	       "you forget to put a } !!!");
	exit(0);	  
	}

	if(strcmp(Parse_GramsBox[0],"}") == 0){

	  /* Check boundary conditions properties */	  
	  break;
	  
	}
      }
      else{
	fprintf(stderr,"%s : %s \n",
	       "Error in GramsBox()",
	       "Use this format -> GramsBox (Type=srt,File=str)  { !!!");
	exit(0);
      }
    }
  }
    
  /* Close .dat file */
  /* Final message */
  printf("\t * End of read data file !!! \n");
  fclose(Sim_dat);

  return FEM_Mesh;
}
