#include "nl-partsol.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))


/**********************************************************************/

Load * GramsNeumannBC(char * Name_File, int NumNeumannBC, int GPxElement)
/*
  GramsNeumannBC (Nodes=ListNodes.txt) {
  V.x Load_x.txt
  V.y Load_x.txt
}
*/
{

  /* Number of dimensions */
  int Ndim = NumberDimensions;
  
  /* Define new load case for the contact forces */
  Load * F = (Load *)Allocate_Array(NumNeumannBC,sizeof(Load));

  /* Simulation file */
  FILE * Sim_dat;

  /* Number of words */
  int Num_words_line;

  /* Parse lines of GramsNeumannBC */
  char Line_GramsNeumannBC[MAXC] = {0};
  char * Parse_GramsNeumannBC[MAXW] = {NULL};
  char * Parse_Nodes[MAXW] = {NULL};

  /* Boundaries iterator */
  int IndexLoad = 0;

  /* Parse file name with the list of nodes */
  char Route_Nodes[MAXC] = {0};
  char FileNodesRoute[MAXC];
  
  /* Array */
  int Num_Nodes;
  ChainPtr Chain_Nodes = NULL;
  int * Array_Nodes;

  /* Parse GramsNeumannBC properties */
  char Line_Properties[MAXC] = {0};
  char * Parse_Properties[MAXW] = {NULL};
  char FileLoadRoute[MAXC];

  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsNeumannBC()",
	    "Incorrect lecture of",
	    Name_File);
    exit(0);
  }

  /* Generate route */
  generate_route(Route_Nodes,Name_File);
  
  /* Read GramsNeumannBC line  */
  while( fgets(Line_GramsNeumannBC, sizeof(Line_GramsNeumannBC), Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    Num_words_line = parse (Parse_GramsNeumannBC, Line_GramsNeumannBC," \n\t");
    if (Num_words_line < 0){
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsNeumannBC ()",
	      "Parser failed");
      exit(0);
    }

    /* Find GramsNeumannBC line */
    if ((Num_words_line > 0) &&
    	(strcmp(Parse_GramsNeumannBC[0],"GramsNeumannBC") == 0 ) &&
	((strcmp(Parse_GramsNeumannBC[2],"{") == 0))){

      /* File with the nodes */
      Num_words_line = parse (Parse_Nodes, Parse_GramsNeumannBC[1],"(=)");
      if( (Num_words_line != 2) ||
	  (strcmp(Parse_Nodes[0],"Nodes") != 0)){
	fprintf(stderr,"%s : %s \n",
		"Error in GramsNeumannBC()",
		"Use this format -> (Nodes=str) !!!");
	exit(0);
      }
      
      /* Read file with the nodes */
      sprintf(FileNodesRoute,"%s%s",Route_Nodes,Parse_Nodes[1]);

      /* Get an array with the nodes */
      Chain_Nodes = File2Chain(FileNodesRoute);
      Num_Nodes = get_Lenght_Set(Chain_Nodes);
      Array_Nodes = Set_to_Pointer(Chain_Nodes,Num_Nodes);
      free_Set(Chain_Nodes);

      /* Fill GPs */
      F[IndexLoad].NumNodes = Num_Nodes*GPxElement;      
      F[IndexLoad].Nodes =
	(int *)Allocate_ArrayZ(F[IndexLoad].NumNodes,sizeof(int));      
      for(int i = 0 ; i<Num_Nodes ; i++){
	for(int j = 0 ; j<GPxElement ; j++){
	  F[IndexLoad].Nodes[i*GPxElement+j] = Array_Nodes[i]*GPxElement+j;
	}
      }

      /* Free array nodes */
      free(Array_Nodes);

      /* Number of dimensions of the BCC */
      F[IndexLoad].Dim = Ndim;
      /* Direction of the BCC */
      F[IndexLoad].Dir =
	(int *)Allocate_ArrayZ(Ndim,sizeof(int));
      /* Curve for each dimension */
      F[IndexLoad].Value =
	(Curve *)Allocate_Array(Ndim,sizeof(Curve));
      /* Information of the BCC */
      if(strcmp(Formulation,"-V") == 0){
  	/* Name of the BCC */
  	strcpy(F[IndexLoad].Info,"Velocity");
      }

      /* Read properties line  */
      while( fgets(Line_Properties, sizeof(Line_Properties), Sim_dat) != NULL ){

	/* Read the line with the space as separators */
	Num_words_line = parse (Parse_Properties, Line_Properties," \n\t");
	if (Num_words_line < 0){
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsNeumannBC ()",
		  "Parser failed");
	  exit(0);
	}
	if((Num_words_line > 0) &&
	   (strcmp(Parse_Properties[0],"}") != 0 )){
	    if(strcmp(Formulation,"-V") == 0){ 
	      if(strcmp(Parse_Properties[0],"V.x") == 0){
		/* Fill the direction of the BCC */
		F[IndexLoad].Dir[0] = 1;
		/* Read the curve to impose the boundary condition */
		if(strcmp(Parse_Properties[1],"NULL") != 0){
		  sprintf(FileLoadRoute,"%s%s",Route_Nodes,Parse_Properties[1]);
		  F[IndexLoad].Value[0] = ReadCurve(FileLoadRoute);
		  puts("*************************************************");
		  printf(" \t %s (%i) : \n \t %s %s \n",
			 "* Neumann BC",
			 F[IndexLoad].NumNodes,
			 Parse_Properties[1],FileLoadRoute);
		}
		else{
		  F[IndexLoad].Dir[0] = 0;
		}
	      }
	      else if(strcmp(Parse_Properties[0],"V.y") == 0){
		/* Fill the direction of the BCC */
		F[IndexLoad].Dir[1] = 1;
		/* Read the curve to impose the boundary condition */
		if(strcmp(Parse_Properties[1],"NULL") != 0){
		  sprintf(FileLoadRoute,"%s%s",Route_Nodes,Parse_Properties[1]);
		  F[IndexLoad].Value[1] =
		    ReadCurve(FileLoadRoute);
		  puts("*************************************************");
		  printf(" \t %s (%i) : \n \t %s %s \n",
			 "* Neumann BC",
			 F[IndexLoad].NumNodes,
			 Parse_Properties[1],FileLoadRoute);
		}
		else{
		  F[IndexLoad].Dir[1] = 0;
		}
	      }
	      else if(strcmp(Parse_Properties[0],"V.z") == 0){
		/* Fill the direction of the BCC */
		F[IndexLoad].Dir[2] = 1;
		/* Read the curve to impose the boundary condition */
		if(strcmp(Parse_Properties[1],"NULL") != 0){
		  sprintf(FileLoadRoute,"%s%s",Route_Nodes,Parse_Properties[1]);
		  F[IndexLoad].Value[2] =
		    ReadCurve(FileLoadRoute);
		  puts("*************************************************");
		  printf(" \t %s (%i) : \n \t %s %s \n",
			 "* Neumann BC",
			 F[IndexLoad].NumNodes,
			 Parse_Properties[1],FileLoadRoute);
		}
		else{
		  F[IndexLoad].Dir[2] = 0;
		}
	      }
	      else{
		fprintf(stderr,"%s : %s \n",
			"Error in GramsNeumannBC()", "Velocity BCC -> U, V, W");
		exit(0);
	      }
	    }
	}
	else{
	  break;
	}
	
      }
      
      /* Increment the index */
      IndexLoad++;
    }

  }

  return F;
}
