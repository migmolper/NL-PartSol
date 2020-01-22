#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "../GRAMS/grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

/*********************************************************************/

GaussPoint GramsSolid2D(char * Name_File, Mesh FEM_Mesh)
/*
  
 */
{

  /* Simulation file */
  FILE * Sim_dat;
  
  /* Parser num chars */
  int Num_words_parse;

  /* Special variables GramsSolid2D */
  char Line_GramsSolid2D[MAXC] = {0}; 
  char * Parse_GramsSolid2D[MAXW] = {NULL};
  char * Parse_Mesh_id[MAXW] = {NULL};

  /* Parse file name with the list of nodes */
  char * Name_File_Copy = malloc(strlen(Name_File)); 
  char * Name_Parse[MAXW] = {NULL};
  char Route_Mesh[MAXC] = {0};
  
  /* Material point mesh (Gauss-Points) */
  Mesh MPM_GID_Mesh;
  GaussPoint MPM_Mesh;
  Matrix Poligon_Coordinates;
  ChainPtr Poligon_Connectivity;
  ChainPtr Vertex;
  int NumVertex;
  Matrix Poligon_Centroid;

  /* Set to false check variables */
  bool Is_GramsSolid2D = false;
  bool Is_GramsTime = false;
  bool Is_GramsShapeFun = false;
  bool Is_GramsMaterials = false; 
  bool Is_GramsInitials = false;
  bool Is_GramsBodyForces = false;
  bool Is_GramsContactForces = false;

  /* Initialize counters */
  int Counter_Materials = 0;
  int Counter_BodyForces = 0;
  int Counter_ContactForces = 0;
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsInitials()",
	    "Incorrect lecture of",
	    Name_File);
    exit(0);
  }
  
  /**************************************************/
  /********* Read and check GramsSolid2D ************/
  /**************************************************/
  while(fgets(Line_GramsSolid2D,sizeof(Line_GramsSolid2D),Sim_dat) != NULL){
    
    /* Read the line with the space as separators */
    Num_words_parse = parse(Parse_GramsSolid2D,Line_GramsSolid2D," \n\t");
    if (Num_words_parse < 0){
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsSolid2D()",
	      "Parser failed");
      exit(0);
    }

    if ((Num_words_parse > 0) &&
	(strcmp(Parse_GramsSolid2D[0],"GramsSolid2D") == 0)){
      Is_GramsSolid2D = true;
      Num_words_parse = parse(Parse_Mesh_id, Parse_GramsSolid2D[1],"(=)");
      strcpy(Name_File_Copy, Name_File);
      MPM_MeshFileName = Parse_Mesh_id[1];
      Num_words_parse = parse(Name_Parse,Name_File_Copy,"(/)");
      for(int i = 0 ; i<Num_words_parse-1 ; i++){
	strcat(Route_Mesh, Name_Parse[i]);
	strcat(Route_Mesh,"/");
      }
      strcat(Route_Mesh,MPM_MeshFileName);
      free(Name_File_Copy);	  
    }
    if ((Num_words_parse > 0) &&
	(strcmp(Parse_GramsSolid2D[0],"GramsTime") == 0)){
      Is_GramsTime = true;
    }
    if ((Num_words_parse > 0) &&
	(strcmp(Parse_GramsSolid2D[0],"GramsShapeFun") == 0)){
      Is_GramsShapeFun = true;
    }
    if ((Num_words_parse > 0) &&
	(strcmp(Parse_GramsSolid2D[0],"GramsMaterials") == 0)){
      Is_GramsMaterials = true;
      Counter_Materials++;
    }
    if ((Num_words_parse > 0) &&
	(strcmp(Parse_GramsSolid2D[0],"GramsInitials") == 0)){
      Is_GramsInitials = true;
    }
    if ((Num_words_parse > 0) &&
	(strcmp(Parse_GramsSolid2D[0],"GramsBodyForces") == 0)){
      Is_GramsBodyForces = true;
      Counter_BodyForces++;
    }
    if ((Num_words_parse > 0) &&
	(strcmp(Parse_GramsSolid2D[0],"GramsContactForces") == 0)){
      Is_GramsContactForces = true;
      Counter_ContactForces++;
    }    
  }

  /**************************************************/
  /***************** Define MPM Mesh ****************/
  /**************************************************/
  /* Define GP Mesh */
  if(Is_GramsSolid2D){
    
    /* Read GP mesh */
    MPM_GID_Mesh = ReadGidMesh(Route_Mesh);

    /* The number of Gauss-Points is the same as the number of elements
       in the input mesh, because we set a GP in the middle of each element */
    MPM_Mesh.NumGP = MPM_GID_Mesh.NumElemMesh;
 
    /* Allocate fields */
    MPM_Mesh.Element_id = /* Index of the Element */
      (int *)Allocate_ArrayZ(MPM_Mesh.NumGP,sizeof(int));  
    MPM_Mesh.NumberNodes = /* Number of tributary nodes for each GP */
      (int *)Allocate_ArrayZ(MPM_Mesh.NumGP,sizeof(int));
    MPM_Mesh.ListNodes =   /* Nodal connectivity of each GP */
      (ChainPtr *)malloc(MPM_Mesh.NumGP*sizeof(ChainPtr));
    if(MPM_Mesh.ListNodes == NULL){
      printf("%s : %s \n",
	     "Define_GP_Mesh",
	     "Memory error for ListNodes");
      exit(0);
    }
    for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
      MPM_Mesh.ListNodes[i] = NULL;  
    }
    MPM_Mesh.Beps =  /* GPs near to each GP */
      (ChainPtr *)malloc(MPM_Mesh.NumGP*sizeof(ChainPtr));
    if(MPM_Mesh.Beps == NULL){
      printf("%s : %s \n",
	     "Define_GP_Mesh",
	     "Memory error for Beps");
      exit(0);
    }
    for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
      MPM_Mesh.Beps[i] = NULL;  
    }
  
    MPM_Mesh.Phi.x_GC = /* Coordinates of the GP (Global/Local)*/
      MatAllocZ(MPM_Mesh.NumGP,3);
    strcpy(MPM_Mesh.Phi.x_GC.Info,"Global Coordinates");

    /**************************************************/
    /* Read parameters of the time-integration scheme */
    /**************************************************/
    if(Is_GramsTime){
      GramsTime(Name_File);
    }
    else{
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsSolid2D()",
	      "GramsTime no defined");
      exit(0);
    }
    /**************************************************/
    /*********** Read Material parameters *************/
    /**************************************************/
    if(Is_GramsMaterials){
      MPM_Mesh.NumberMaterials = Counter_Materials;
      MPM_Mesh.Mat = GramsMaterials(Name_File,MPM_Mesh);
      MPM_Mesh.MatIdx = (int *)malloc(MPM_Mesh.NumGP*sizeof(int));
      for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
	MPM_Mesh.MatIdx[i] = 0;
      }
    }
    else{
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsSolid2D()",
	      "GramsMaterials no defined");
      exit(0);
    }
    /**************************************************/
    /********* Read Shape functions parameters ********/
    /**************************************************/
    if(Is_GramsShapeFun){
      GramsShapeFun(Name_File);
      /* Lenght of the Voxel (Only GIMP) */
      if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
	MPM_Mesh.lp = MatAllocZ(MPM_Mesh.NumGP,NumberDimensions);
	strcpy(MPM_Mesh.lp.Info,"Voxel lenght GP");
      }
      /* Lagrange Multipliers / Beta (Only LME ) */
      if(strcmp(ShapeFunctionGP,"LME") == 0){
	MPM_Mesh.lambda = MatAllocZ(MPM_Mesh.NumGP,NumberDimensions);
	strcpy(MPM_Mesh.lambda.Info,"Lagrange Multiplier");
	MPM_Mesh.Beta = MatAllocZ(MPM_Mesh.NumGP,NumberDimensions);
	strcpy(MPM_Mesh.Beta.Info,"Beta parameter");
      }
    }
    else{
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsSolid2D()",
	      "GramsShapeFun no defined");
      exit(0);
    }
    /**************************************************/    
    /****** Allocate vectorial/tensorial fields *******/
    /**************************************************/
    switch(NumberDimensions){
    case 1 :
      /* Natural coordinates (Vectorial) */
      MPM_Mesh.Phi.x_EC = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.x_EC.Info,"Element Coordinates GP");
      /* Displacement field (Vectorial) */
      MPM_Mesh.Phi.dis = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.dis.Info,"Displacement field GP");
      /* Velocity field (Vectorial) */
      MPM_Mesh.Phi.vel = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.vel.Info,"Velocity field GP");
      /* Acceleration field (Vectorial) */
      MPM_Mesh.Phi.acc = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.acc.Info,"Acceleration field GP");
      /* Strain field (Tensor) */
      MPM_Mesh.Phi.Strain = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.Strain.Info,"Strain field GP");
      /* Stress field (Tensor) */
      MPM_Mesh.Phi.Stress = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.Stress.Info,"Stress field GP");
      /* Deformation Energy (Scalar) */
      MPM_Mesh.Phi.W = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.W.Info,"Deformation Energy GP");
      /* Damage parameter (fracture) */
      MPM_Mesh.Phi.ji = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.ji.Info,"Damage parameter GP");
      /* Mass */
      MPM_Mesh.Phi.mass = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.mass.Info,"Mass GP");
      /* Density */
      MPM_Mesh.Phi.rho = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.rho.Info,"Density GP");
      break;
    case 2 :
      /* Natural coordinates (Vectorial) */
      MPM_Mesh.Phi.x_EC = MatAllocZ(MPM_Mesh.NumGP,2);
      strcpy(MPM_Mesh.Phi.x_EC.Info,"Element Coordinates GP");    
      /* Displacement field (Vectorial) */
      MPM_Mesh.Phi.dis = MatAllocZ(MPM_Mesh.NumGP,2);
      strcpy(MPM_Mesh.Phi.dis.Info,"Displacement field GP");
      /* Velocity field (Vectorial) */
      MPM_Mesh.Phi.vel = MatAllocZ(MPM_Mesh.NumGP,2);
      strcpy(MPM_Mesh.Phi.vel.Info,"Velocity field GP");
      /* Acceleration field (Vectorial) */
      MPM_Mesh.Phi.acc = MatAllocZ(MPM_Mesh.NumGP,2);
      strcpy(MPM_Mesh.Phi.acc.Info,"Acceleration field GP");
      /* Strain field (Tensor) */
      MPM_Mesh.Phi.Strain = MatAllocZ(MPM_Mesh.NumGP,3);
      strcpy(MPM_Mesh.Phi.Strain.Info,"Strain field GP");
      /* Stress field (Tensor) */
      MPM_Mesh.Phi.Stress = MatAllocZ(MPM_Mesh.NumGP,3);
      strcpy(MPM_Mesh.Phi.Stress.Info,"Stress field GP");
      /* Deformation Energy (Scalar) */
      MPM_Mesh.Phi.W = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.W.Info,"Deformation Energy GP");
      /* Damage parameter (fracture) */
      MPM_Mesh.Phi.ji = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.ji.Info,"Damage parameter GP");
      /* Mass */
      MPM_Mesh.Phi.mass = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.mass.Info,"Mass GP");
      /* Density */
      MPM_Mesh.Phi.rho = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.rho.Info,"Density GP");
      break;
    case 3:
      /* Natural coordinates (Vectorial) */
      MPM_Mesh.Phi.x_EC = MatAllocZ(MPM_Mesh.NumGP,3);
      strcpy(MPM_Mesh.Phi.x_EC.Info,"Element Coordinates GP");   
      /* Displacement field (Vectorial) */
      MPM_Mesh.Phi.dis = MatAllocZ(MPM_Mesh.NumGP,3);
      strcpy(MPM_Mesh.Phi.dis.Info,"Displacement field GP");
      /* Velocity field (Vectorial) */
      MPM_Mesh.Phi.vel = MatAllocZ(MPM_Mesh.NumGP,3);
      strcpy(MPM_Mesh.Phi.vel.Info,"Velocity field GP");
      /* Acceleration field (Vectorial) */
      MPM_Mesh.Phi.acc = MatAllocZ(MPM_Mesh.NumGP,3);
      strcpy(MPM_Mesh.Phi.acc.Info,"Acceleration field GP");
      /* Strain field (Tensor) */
      MPM_Mesh.Phi.Strain = MatAllocZ(MPM_Mesh.NumGP,9);
      strcpy(MPM_Mesh.Phi.Strain.Info,"Strain field GP");
      /* Stress field (Tensor) */
      MPM_Mesh.Phi.Stress = MatAllocZ(MPM_Mesh.NumGP,9);
      strcpy(MPM_Mesh.Phi.Stress.Info,"Stress field GP");
      /* Deformation Energy (Scalar) */
      MPM_Mesh.Phi.W = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.W.Info,"Deformation Energy GP");
      /* Damage parameter (fracture) */
      MPM_Mesh.Phi.ji = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.ji.Info,"Damage parameter GP");
      /* Mass */
      MPM_Mesh.Phi.mass = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.mass.Info,"Mass GP");
      /* Density */
      MPM_Mesh.Phi.rho = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.Phi.rho.Info,"Density GP");
      break;
    default:
      puts("Error in Initialize_GP_Mesh() : Wrong number of dimensions");
      exit(0);
    }

    /* Fill geometrical properties of the GP mesh */
    for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
      /* Get the connectivity of the elements vertex */
      Poligon_Connectivity = MPM_GID_Mesh.Connectivity[i];    
      /* Generate a matrix with the poligon coordinates */
      NumVertex = MPM_GID_Mesh.NumNodesElem[i];
      Poligon_Coordinates = MatAllocZ(NumVertex,NumberDimensions); 
      /* Initialize chain interator */
      Vertex = Poligon_Connectivity;
      /* Loop in the chain to fill the poligon */
      for(int k = 0, I_Vertex = 0;
	  (k<NumVertex) || (Vertex != NULL);
	  k++, Vertex = Vertex->next){
	I_Vertex = Vertex->I;
	for(int l = 0 ; l<NumberDimensions ; l++){
	  Poligon_Coordinates.nM[k][l] = MPM_GID_Mesh.Coordinates.nM[I_Vertex][l];
	}
      }
      /* Free data */
      FreeChain(MPM_GID_Mesh.Connectivity[i]);
      /* Get the area (Poligon_Centroid.n) 
	 and the position of the centroid (Poligon_Centroid.nV) */
      Poligon_Centroid = Centroid_Poligon(Poligon_Coordinates);
      /* Free data */
      FreeMat(Poligon_Coordinates);    
      /* Assign the mass parameter */
      MPM_Mesh.Phi.mass.nV[i] = Poligon_Centroid.n*MPM_Mesh.Mat[0].rho;
      /* Set the initial density */
      MPM_Mesh.Phi.rho.nV[i] = MPM_Mesh.Mat[0].rho;
      /* Get the coordinates of the centre */
      MPM_Mesh.Phi.x_GC.nM[i][0] = Poligon_Centroid.nV[0];
      MPM_Mesh.Phi.x_GC.nM[i][1] = Poligon_Centroid.nV[1];
      MPM_Mesh.Phi.x_GC.nM[i][2] = 0.0;         
      /* Local coordinates of the element */
      MPM_Mesh.Element_id[i] = -999;
      MPM_Mesh.NumberNodes[i] = 4;
      /* Free data */
      FreeMat(Poligon_Centroid);    
    }
    /**************************************************/    
    /************** Read initial values ***************/
    /**************************************************/    
    if(Is_GramsInitials){
      GramsInitials(Name_File,MPM_Mesh);
    }
    else{
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsSolid2D()",
	      "GramsInitials no defined");
      exit(0);
    }
    /**************************************************/    
    /************** Read external forces **************/
    /**************************************************/    
    if(Is_GramsContactForces){
      MPM_Mesh.NumberContactForces = Counter_ContactForces;
      MPM_Mesh.F = GramsContactForces(Name_File,MPM_Mesh);
    }
    else{
      MPM_Mesh.NumberContactForces = Counter_ContactForces;
      puts("*************************************************");
      printf(" \t %s : \n\t %s \n",
	     "* No contact forces defined in",
	     Name_File);
    }
    /**************************************************/    
    /*************** Read body forces *****************/
    /**************************************************/    
    if(Is_GramsBodyForces){
      MPM_Mesh.NumberBodyForces = Counter_BodyForces;
      MPM_Mesh.B = GramsBodyForces(Name_File,MPM_Mesh); 
    }
    else{
      MPM_Mesh.NumberBodyForces = Counter_BodyForces;
      puts("*************************************************");
      printf(" \t %s : \n\t %s \n",
	 "* No body forces defined in",
	     Name_File);
    }
    /**************************************************/
    /******** INITIALIZE SHAPE FUNCTIONS **************/
    /**************************************************/ 
    puts("*************************************************");
    printf("\t * %s \n",
	   "Initialize shape functions ...");
    if(strcmp(ShapeFunctionGP,"MPMQ4") == 0){
      Q4_Initialize(MPM_Mesh, FEM_Mesh);
    }
    else if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
      GIMP_Initialize(MPM_Mesh,FEM_Mesh);
    }
    else if(strcmp(ShapeFunctionGP,"LME") == 0){
      LME_Initialize(MPM_Mesh,FEM_Mesh);
    }
    
    /**************************************************/    
    /************* Free the input data ****************/
    /**************************************************/    
    FreeMat(MPM_GID_Mesh.Coordinates);
    free(MPM_GID_Mesh.Connectivity);
    free(MPM_GID_Mesh.ActiveNode);
    free(MPM_GID_Mesh.NumNeighbour);
    free(MPM_GID_Mesh.NodeNeighbour);

  } 
  else{
    fprintf(stderr,"%s : %s \n",
	    "Error in GramsSolid2D()",
	    "GramsBodyForces no defined");
    exit(0);
  }
  
  return MPM_Mesh;
}

/***************************************************************************/


