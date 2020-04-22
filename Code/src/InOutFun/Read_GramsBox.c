#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "grams.h"

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
  char FileMeshRoute[MAXC];
   
  /* Parse line of GramsBox */
  int Num_GramsBox;
  char Line_GramsBox[MAXC] = {0};
  char * Parse_GramsBox[MAXW] = {NULL};
  int Num_GramsBox_Prop;
  char * Parse_Mesh_id[MAXW] = {NULL};

  /* Parse lines of GramsBoundary */
  int Num_GramsBoundary;
  char Line_GramsBoundary[MAXC] = {0};
  char * Parse_GramsBoundary[MAXW] = {NULL};
  int NumBounds;

  /* Auxiliar variable for status */
  char * STATUS_LINE;
    
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsBox()",
	    "Incorrect lecture of",
	    Name_File);
    exit(0);
  }

  /* Generate route */
  strcpy(Name_File_Copy, Name_File);
  Num_GramsBox = parse(Name_Parse,Name_File_Copy,"(/)");
  strcat(Route_Mesh,"./");
  for(int i = 0 ; i<Num_GramsBox-1 ; i++){
    strcat(Route_Mesh, Name_Parse[i]);
    strcat(Route_Mesh,"/");
  }

  /* Read GramsBox line  */
  while( fgets(Line_GramsBox, sizeof(Line_GramsBox), Sim_dat) != NULL ){
    
    /* Read the line with the space as separators */
    Num_GramsBox = parse (Parse_GramsBox, Line_GramsBox," \n\t");
    
    if (Num_GramsBox < 0){
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsBox ()",
	      "Parser failed");
      exit(0);
    }

    /* Find GramsBox line */
    if ((Num_GramsBox > 0) &&
    	(strcmp(Parse_GramsBox[0],"GramsBox") == 0 ) &&
	((strcmp(Parse_GramsBox[2],"{") == 0))){

      /* Bidimensional mesh */
      NumberDimensions = 2;

      /* Read the type of the mesh and the name of the file */
      Num_GramsBox_Prop = parse (Parse_Mesh_id, Parse_GramsBox[1],"(=,)");
      if( (Num_GramsBox_Prop != 4) ||
	  (strcmp(Parse_Mesh_id[0],"Type") != 0) ||
	  (strcmp(Parse_Mesh_id[2],"File") != 0) ){
	fprintf(stderr,"%s : %s \n",
		"Error in GramsBox()",
		"Use this format -> (Type=str,File=str) !!!");
	exit(0);
      }

      /* Read GID-type mesh */
      if(strcmp(Parse_Mesh_id[1],"GID") == 0){
	/* Read file with the mesh */
	sprintf(FileMeshRoute,"%s%s",Route_Mesh,Parse_Mesh_id[3]);

	puts("*************************************************");
	printf(" \t %s : \n \t %s \n",
	       "* Read GID mesh in",FileMeshRoute);
	FEM_Mesh = ReadGidMesh(FileMeshRoute);
      }
      else{
	fprintf(stderr,"%s : %s %s \n",
		"Error in GramsBox(Type=*, )",
		"Unrecognized mesh Type",
		Parse_Mesh_id[1]);
      }

      /* By default the number of BCCs is 0 */
      NumBounds = 0;

      /* Check format GramsBox (Type=,File=) { } */
      if((Num_GramsBox == 4) &&
	 (strcmp(Parse_GramsBox[3],"}") == 0)){
	puts("*************************************************");
	printf(" \t %s : \n \t %i \n",
	       "* Number of boundaries",0);
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
      Num_GramsBoundary = parse(Parse_GramsBoundary,Line_GramsBoundary," =\t\n");

      /* Check format GramsBox (Type=,File=) { 
	 }
      */
      if((Num_GramsBoundary > 0) &&
	 (strcmp(Parse_GramsBoundary[0],"}") == 0)){
	puts("*************************************************");
	printf(" \t %s : \n \t %i \n",
	       "* Number of boundaries",0);
	break;
      }

      /* Loop to count the number of boundaries */
      while(STATUS_LINE != NULL){
	if((Num_GramsBoundary > 0) &&
	   (strcmp(Parse_GramsBoundary[0],"GramsBoundary") == 0)){
	  NumBounds++;
	}
	/* Continue reading to the end*/
	STATUS_LINE = fgets(Line_GramsBoundary,
			    sizeof(Line_GramsBoundary),Sim_dat);
	Num_GramsBoundary = parse(Parse_GramsBoundary,
			       Line_GramsBoundary," =\t\n");	
      }
      
    }
  }

  /* Close file */
  fclose(Sim_dat);

  /**************************************************/
  /********** Set the minimum mesh size *************/
  /**************************************************/
  FEM_Mesh.DeltaX = GetMinElementSize(FEM_Mesh);
  puts("*************************************************");
  printf(" \t %s : \n \t %f \n",
	 "* Mesh size",FEM_Mesh.DeltaX);
  /**************************************************/	
  /*** Generate nodal connectivity of the mesh : ****/
  /******* list of elements near to a node ***********/
  /**************************************************/
  GetNodalConnectivity(FEM_Mesh);  
  /**************************************************/	
  /** Initialize particle connectivity of each node */
  /**************************************************/
  FEM_Mesh.NumParticles =
    (int *)Allocate_ArrayZ(FEM_Mesh.NumNodesMesh,sizeof(int));
  FEM_Mesh.I_particles =
    (ChainPtr *)malloc(FEM_Mesh.NumNodesMesh*sizeof(ChainPtr));
  if(FEM_Mesh.I_particles == NULL){
    printf("%s : %s \n",
	   "GetNodalConnectivity",
	   "Memory error for I_particles");
    exit(0);
  }
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.I_particles[i] = NULL;
  }
  /**************************************************/  
  /******** Allocate array with the boundaries ******/
  /**************************************************/
  if(NumBounds > 0){
    puts("*************************************************");
    printf(" \t %s : \n \t %i \n",
	   "* Number of boundaries",NumBounds);
    FEM_Mesh.Bounds.NumBounds = NumBounds;
    FEM_Mesh.Bounds = GramsBoundary(Name_File,NumBounds);
  }

  /* Free data */
  free(Name_File_Copy);

  return FEM_Mesh;
}
