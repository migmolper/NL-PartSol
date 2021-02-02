#include "nl-partsol.h"

Boundaries GramsBoundary(char * Name_File,int NumBounds){

  /* Boundaries */
  Boundaries Bounds;

  /* Simulation file */
  FILE * Sim_dat;

  /* Number of words */
  int Num_words_line;

  /* Parse lines of GramsBoundary */
  char Line_GramsBoundary[MAXC] = {0};
  char * Parse_GramsBoundary[MAXW] = {NULL};
  char * Parse_Nodes[MAXW] = {NULL};

  /* Boundaries iterator */
  int IndexBoundary = 0;

  /* Parse file name with the list of nodes */
  char Route_Nodes[MAXC] = {0};
  char FileNodesRoute[MAXC];
  ChainPtr Chain_Nodes = NULL;

  /* Parse GramsBoundary properties */
  char Line_Properties[MAXC] = {0};
  char * Parse_Properties[MAXW] = {NULL};
  char FileLoadRoute[MAXC];

  /* Assign the number of bounds */
  Bounds.NumBounds = NumBounds;

  /* Allocate boundaries */
  Bounds.BCC_i = (Load *)Allocate_Array(NumBounds,sizeof(Load));

  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsBoundary()",
	    "Incorrect lecture of",
	    Name_File);
    exit(EXIT_FAILURE);
  }

  /* Generate route */
  generate_route(Route_Nodes,Name_File);

  /* Read GramsBoundary line  */
  while( fgets(Line_GramsBoundary, sizeof(Line_GramsBoundary), Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    Num_words_line = parse (Parse_GramsBoundary, Line_GramsBoundary," \n\t");
    if (Num_words_line < 0){
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsBoundary ()",
	      "Parser failed");
      exit(EXIT_FAILURE);
    }

    /* Find GramsBoundary line */
    if ((Num_words_line > 0) && (strcmp(Parse_GramsBoundary[0],"GramsBoundary") == 0 ) &&
	((strcmp(Parse_GramsBoundary[2],"{") == 0))){

      /* File with the nodes */
      Num_words_line = parse (Parse_Nodes, Parse_GramsBoundary[1],"(=)");
      if( (Num_words_line != 2) ||
	  (strcmp(Parse_Nodes[0],"File") != 0)){
	fprintf(stderr,"%s : %s \n",
		"Error in GramsBoundary()",
		"Use this format -> (File=str) !!!");
	exit(EXIT_FAILURE);
      }
      
      /* Read file with the nodes */
      sprintf(FileNodesRoute,"%s%s",Route_Nodes,Parse_Nodes[1]);

      /* Get an array with the nodes */
      Chain_Nodes = File2Chain(FileNodesRoute);
      Bounds.BCC_i[IndexBoundary].NumNodes = lenght__SetLib__(Chain_Nodes);
      Bounds.BCC_i[IndexBoundary].Nodes =
	set_to_memory__SetLib__(Chain_Nodes,Bounds.BCC_i[IndexBoundary].NumNodes);
      free__SetLib__(&Chain_Nodes);

      /*
       	Number of dimensions of the BCC
      */
      Bounds.BCC_i[IndexBoundary].Dim = NumberDOF;
      /*
      	Direction of the BCC 
      */
      Bounds.BCC_i[IndexBoundary].Dir = (int *)Allocate_ArrayZ(NumberDOF,sizeof(int));
      /*
      	Curve for each dimension
      */
      Bounds.BCC_i[IndexBoundary].Value = (Curve *)Allocate_Array(NumberDOF,sizeof(Curve));

      /*
      	Information of the BCC
      */
      if(strcmp(Formulation,"-u") == 0)
      {
  		strcpy(Bounds.BCC_i[IndexBoundary].Info,"Displacements");
      }
      else if(strcmp(Formulation,"-upw") == 0)
      {
		strcpy(Bounds.BCC_i[IndexBoundary].Info,"Displacements-Pore-Water-Pressure");
      }

      /* Read properties line  */
      while( fgets(Line_Properties, sizeof(Line_Properties), Sim_dat) != NULL )
      {

		/* Read the line with the space as separators */
		Num_words_line = parse (Parse_Properties, Line_Properties," \n\t");
		if (Num_words_line < 0)
		{
		  fprintf(stderr,"%s : %s \n",
			  "Error in GramsBoundary ()",
			  "Parser failed");
		  exit(EXIT_FAILURE);
		}
		if((Num_words_line > 0) && (strcmp(Parse_Properties[0],"}") != 0 ))
		{
		  if(strcmp(Parse_Properties[0],"BcDirichlet") == 0 )
		  {
		    if(strcmp(Formulation,"-u") == 0)
		    { 
		      if(strcmp(Parse_Properties[1],"V.x") == 0)
		      {
				/* Fill the direction of the BCC */
				Bounds.BCC_i[IndexBoundary].Dir[0] = 1;
				/* Read the curve to impose the boundary condition */
				if(strcmp(Parse_Properties[2],"NULL") != 0)
				{
				  sprintf(FileLoadRoute,"%s%s",Route_Nodes,Parse_Properties[2]);
				  Bounds.BCC_i[IndexBoundary].Value[0] = ReadCurve(FileLoadRoute);
				  puts("*************************************************");
				  printf(" \t %s (%i) : \n \t %s %s \n",
					 "* BcDirichlet",
					 Bounds.BCC_i[IndexBoundary].NumNodes,
					 Parse_Properties[1],FileLoadRoute);
				}
				else
				{
				  Bounds.BCC_i[IndexBoundary].Dir[0] = 0;
				}
		      }
		      else if(strcmp(Parse_Properties[1],"V.y") == 0)
		      {
				/* Fill the direction of the BCC */
				Bounds.BCC_i[IndexBoundary].Dir[1] = 1;
				/* Read the curve to impose the boundary condition */
				if(strcmp(Parse_Properties[2],"NULL") != 0)
				{
				  sprintf(FileLoadRoute,"%s%s",Route_Nodes,Parse_Properties[2]);
				  Bounds.BCC_i[IndexBoundary].Value[1] =
				    ReadCurve(FileLoadRoute);
				  puts("*************************************************");
				  printf(" \t %s (%i) : \n \t %s %s \n",
					 "* BcDirichlet",
					 Bounds.BCC_i[IndexBoundary].NumNodes,
					 Parse_Properties[1],FileLoadRoute);
				}
				else
				{
				  Bounds.BCC_i[IndexBoundary].Dir[1] = 0;
				}
		      }
		      else if(strcmp(Parse_Properties[1],"V.z") == 0)
		      {
				/* Fill the direction of the BCC */
				Bounds.BCC_i[IndexBoundary].Dir[2] = 1;
				/* Read the curve to impose the boundary condition */
				if(strcmp(Parse_Properties[2],"NULL") != 0)
				{
				  sprintf(FileLoadRoute,"%s%s",Route_Nodes,Parse_Properties[2]);
				  Bounds.BCC_i[IndexBoundary].Value[2] =
				    ReadCurve(FileLoadRoute);
				  puts("*************************************************");
				  printf(" \t %s (%i) : \n \t %s %s \n",
					 "* BcDirichlet",
					 Bounds.BCC_i[IndexBoundary].NumNodes,
					 Parse_Properties[1],FileLoadRoute);
				}
				else
				{
				  Bounds.BCC_i[IndexBoundary].Dir[2] = 0;
				}
		      }
		      else
		      {
				fprintf(stderr,"%s : %s \n", "Error in GramsBox()", "Velocity BCC -> U, V, W");
				exit(EXIT_FAILURE);
		      }
		    }
		  }
		}
		else{
		  break;
		}
	
      }
      
      /* Increment the index */
      IndexBoundary++;
    }

  }


  return Bounds;
}
  

