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
  char * Name_File_Copy = malloc(strlen(Name_File)); 
  char * Name_Parse[MAXW] = {NULL};
  char Route_Mesh[MAXC] = {0};

  /* Index of the material */
  int Aux_Mesh_id;
  char * Parse_Mesh_id[MAXW] = {NULL};

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

  /* Auxiliar variable for status */
  char * STATUS_LINE;

  /* Auxiliar table to store the boundaries */
  int * ListNodesBound;
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
  int IndexBoundary;
  

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

  /* Read next line and check GramsBoundary  */
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
      printf("%s %s %s \n",
	     Parse_GramsBox[0],Parse_GramsBox[1],Parse_GramsBox[2]);
      Num_words_line = parse (Parse_Mesh_id, Parse_GramsBox[1],"(=,)");
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
      FEM_Mesh.Bounds.BCC_i = (Load *)Allocate_Array(4,sizeof(Load));

      /* Fill some parameters */
      for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++){
	/* Set to zero the number of nodes in the boundary of the mesh */
	FEM_Mesh.Bounds.BCC_i[i].NumNodes = 0;
	/* Set the labels for a square domain */
	strcpy(FEM_Mesh.Bounds.BCC_i[i].Info,BoundLabels[i]);
      }

      /* Allocate an array of zeros to assign a 1 to those nodes in the boundary */
      CounterNodesBound = (int *)Allocate_ArrayZ(4,sizeof(int));
      ListNodesBound = (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh,sizeof(int));
  
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
	  ListNodesBound[i] = 1;
	  NumNodesBound++;
	}
      }
      
      /* Count the number of nodes in each boundary */
      for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
	if(ListNodesBound[i] == 1){

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
	/* Number of nodes in each boundary */
	FEM_Mesh.Bounds.BCC_i[i].NumNodes = CounterNodesBound[i];
	/* List of nodes of each boundary */
	FEM_Mesh.Bounds.BCC_i[i].Nodes =
	  (int *)Allocate_ArrayZ(CounterNodesBound[i],sizeof(int));
      }

      /* Free data */
      free(CounterNodesBound);
      
      /* Fill the arrays with the nodes */
      for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
	if(ListNodesBound[i] == 1){

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
      free(ListNodesBound);

      /* Allocate variable for the BCCs */
      if(strcmp(Formulation,"Velocity") == 0){
	for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++){
	  /* Curve for each dimension */
	  FEM_Mesh.Bounds.BCC_i[i].Value =
	    (Curve *)Allocate_Array(NumberDOF,sizeof(Curve));
	  /* Name of the BCC */
	  strcpy(FEM_Mesh.Bounds.BCC_i[i].Info,"Velocity");
	  /* Number of dimensions of the BCC */
	  FEM_Mesh.Bounds.BCC_i[i].Dim = NumberDOF;
	  /* Direction of the BCC */
	  FEM_Mesh.Bounds.BCC_i[i].Dir =
	    (int *)Allocate_ArrayZ(NumberDOF,sizeof(int));
	}
      }

      /* Fill boundary conditions */
      if(strcmp(Parse_GramsBox[2],"{") == 0){
	
	if((Num_words_line == 4) &&
	   (strcmp(Parse_GramsBox[3],"}") != 0)){  
	  /* No boundary conditions applied */
	  break;
	}

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
	if(strcmp(Parse_GramsBoundary[0],"}") == 0){  
	  /* No boundary conditions applied */
	  exit(0);
	  break;
	}

	while(STATUS_LINE != NULL){

 	  if((strcmp(Parse_GramsBoundary[1],"GramsBoundary") == 0) &&
	     (strcmp(Parse_GramsBoundary[2],"{") == 0)){
    
	    /* Asign index for the boundary label */
	    if (strcmp(Parse_GramsBoundary[0],"left") == 0){
	      IndexBoundary = 3;
	    }
	    else if (strcmp(Parse_GramsBoundary[0],"top") == 0){
	      IndexBoundary = 2;
	    }
	    else if (strcmp(Parse_GramsBoundary[0],"right") == 0){
	      IndexBoundary = 1;
	    }
	    else if (strcmp(Parse_GramsBoundary[0],"bottom") == 0){
	      IndexBoundary = 0;
	    }
	    else {
	      fprintf(stderr,"%s : %s \n",
		      "Error in GramsBox()",
		      "Invalid boundary label (left, right, top, bottom) !!!");
	      exit(0);
	    }

	    STATUS_LINE =
	      fgets(Line_BCC,sizeof(Line_BCC),Sim_dat);
	    if(STATUS_LINE == NULL){
	      fprintf(stderr,"%s : %s \n",
		      "Error in GramsBox()",
		      "Unspected EOF !!!");
	      exit(0);	
	    }
	    Aux_Mesh_id = parse(Parse_BCC,Line_BCC," =\t\n");
	    if(Aux_Mesh_id != 3){
	      fprintf(stderr,"%s : %s \n",
		      "Error in GramsBox()",
		      "Invalid format for the BCCs !!!");
	      exit(0);	
	    }

	    while(STATUS_LINE != NULL){

	      /* Velocity boundary condition */
	      if(strcmp(Formulation,"Velocity") == 0){
		if(strcmp(Parse_BCC[1],"V.x") == 0){
		  /* Fill the direction of the BCC */
		  FEM_Mesh.Bounds.BCC_i[IndexBoundary].Dir[0] = 1;
		  /* Read the curve to impose the boundary condition */
		  if(strcmp(Parse_BCC[2],"NULL") != 0){
		    FEM_Mesh.Bounds.BCC_i[IndexBoundary].Value[0] =
		      ReadCurve(Parse_BCC[2]);
		  }
		  else{
		    FEM_Mesh.Bounds.BCC_i[IndexBoundary].Dir[0] = 0;
		  }
		}
		else if(strcmp(Parse_BCC[1],"V.y") == 0){
		  /* Fill the direction of the BCC */
		  FEM_Mesh.Bounds.BCC_i[IndexBoundary].Dir[1] = 1;
		  /* Read the curve to impose the boundary condition */
		  if(strcmp(Parse_BCC[2],"NULL") != 0){
		    FEM_Mesh.Bounds.BCC_i[IndexBoundary].Value[1] =
		      ReadCurve(Parse_BCC[2]);
		  }
		  else{
		    FEM_Mesh.Bounds.BCC_i[IndexBoundary].Dir[1] = 0;
		  }
		}
		else if(strcmp(Parse_BCC[1],"V.z") == 0){
		  /* Fill the direction of the BCC */
		  FEM_Mesh.Bounds.BCC_i[IndexBoundary].Dir[2] = 1;
		  /* Read the curve to impose the boundary condition */
		  if(strcmp(Parse_BCC[2],"NULL") != 0){
		    FEM_Mesh.Bounds.BCC_i[IndexBoundary].Value[2] =
		      ReadCurve(Parse_BCC[2]);
		  }
		  else{
		    FEM_Mesh.Bounds.BCC_i[IndexBoundary].Dir[2] = 0;
		  }
		}
		else{
		  fprintf(stderr,"%s : %s \n",
			  "Error in GramsBox()", "Velocity BCC -> U, V, W");
		  exit(0);
		}
	      }
	      else{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsBox()",
			"Non formulation defined");
		exit(0);
	      }

	      /* Read next line and check BCC */
	      /* Initalize parser */
	      for(int i = 0 ; i<MAXW ; i++){
		Parse_BCC[i] = "\0";
	      }
	      STATUS_LINE = fgets(Line_BCC,
				  sizeof(Line_BCC),
				  Sim_dat);
	      Aux_Mesh_id = parse(Parse_BCC,Line_BCC," =\t\n");
	      if(strcmp(Parse_BCC[0],"}") == 0){
		break;
	      }

	    }
	    
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

	  /* Read next line and check GramsBoundary */
	  /* Initalize parser */
	  for(int i = 0 ; i<MAXW ; i++){
	    Parse_GramsBoundary[i] = "\0";
	  }
	  STATUS_LINE = fgets(Line_GramsBoundary,
			      sizeof(Line_GramsBoundary),
			      Sim_dat);
	  Aux_Mesh_id = parse(Parse_GramsBoundary,Line_GramsBoundary," =\t\n");
	  if(strcmp(Parse_GramsBoundary[0],"}") == 0){
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
	       "No boundary conditions readed !!!");
	exit(0);
      }
    }

    /* Initalize parser */
    for(int i = 0 ; i<MAXW ; i++){
      Parse_GramsBox[i] = "\0";
    }
    
  }
    
  /* Close .dat file */
  /* Final message */
  printf("\t * End of read data file !!! \n");
  fclose(Sim_dat);

  return FEM_Mesh;
}
