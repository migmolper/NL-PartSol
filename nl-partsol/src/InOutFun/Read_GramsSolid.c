#include "grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

/*********************************************************************/

GaussPoint GramsSolid2D(char * Name_File, Mesh FEM_Mesh)
/*
 */
{
  int Ndim = NumberDimensions;
  int NumParticles;

  /* Simulation file */
  FILE * Sim_dat;
  
  /* Parser num chars */
  int Num_words_parse;

  /* Special variables GramsSolid2D */
  char Line_GramsSolid2D[MAXC] = {0}; 
  char * Parse_GramsSolid2D[MAXW] = {NULL};
  char * Parse_Mesh_id[MAXW] = {NULL};

  /* Parse file name with the list of nodes */
  char Route_Mesh[MAXC] = {0};
  
  /* Material point mesh (Gauss-Points) */
  Mesh MPM_GID_Mesh;
  GaussPoint MPM_Mesh;
  Matrix Poligon_Coordinates;
  ChainPtr Poligon_Connectivity;
  ChainPtr Vertex;
  int NumVertex;

  double Area_Element, A_p, th_p, m_p, rho_p;
  int GPxElement = 1;
  int i_p;

  /* Set to false check variables */
  bool Is_GramsSolid2D = false;
  bool Is_GramsShapeFun = false;
  bool Is_GramsMaterials = false;
  bool Is_GramsInitials = false;
  bool Is_GramsBodyForces = false;
  bool Is_GramsNeumannBC = false;

  /* Initialize counters */
  int Counter_Materials = 0;
  int Counter_BodyForces = 0;
  int Counter_GramsNeumannBC = 0;
  
  /* Open and check file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsInitials()","Incorrect lecture of",
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
      MPM_MeshFileName = Parse_Mesh_id[1];
      generate_route(Route_Mesh,Name_File);
      strcat(Route_Mesh,MPM_MeshFileName);      
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
	(strcmp(Parse_GramsSolid2D[0],"GramsNeumannBC") == 0)){
      Is_GramsNeumannBC = true;
      Counter_GramsNeumannBC++;
    }    
  }

  /**************************************************/
  /***************** Define MPM Mesh ****************/
  /**************************************************/
  /* Define GP Mesh */
  if(Is_GramsSolid2D){
    
    /* Read GP mesh */
    MPM_GID_Mesh = ReadGidMesh(Route_Mesh);

    /* Define the number of particles */
    NumParticles = GPxElement*MPM_GID_Mesh.NumElemMesh;
    MPM_Mesh.NumGP = NumParticles;

    /* Closest node to the particle */
    MPM_Mesh.I0 = (int *)Allocate_ArrayZ(NumParticles,sizeof(int));

    /* Number of tributary nodes for each particle */
    MPM_Mesh.NumberNodes = (int *)Allocate_ArrayZ(NumParticles,sizeof(int));

    /* Tributary nodes for each particle */
    MPM_Mesh.ListNodes = (ChainPtr *)malloc(NumParticles*sizeof(ChainPtr));
    if(MPM_Mesh.ListNodes == NULL){
      printf("%s : %s \n",
	     "Define_GP_Mesh","Memory error for ListNodes");
      exit(0);
    }
    for(int i = 0 ; i<NumParticles ; i++){
      MPM_Mesh.ListNodes[i] = NULL;  
    }

    /* List of particles close to each particle */
    MPM_Mesh.Beps = (ChainPtr *)malloc(NumParticles*sizeof(ChainPtr));
    if(MPM_Mesh.Beps == NULL){
      printf("%s : %s \n",
	     "Define_GP_Mesh","Memory error for Beps");
      exit(0);
    }
    for(int i = 0 ; i<NumParticles ; i++){
      MPM_Mesh.Beps[i] = NULL;  
    }

    /**************************************************/
    /********* Read Shape functions parameters ********/
    /**************************************************/
    if(Is_GramsShapeFun){
      GramsShapeFun(Name_File);
      /* Lenght of the Voxel (Only GIMP) */
      if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
	MPM_Mesh.lp = MatAllocZ(NumParticles,Ndim);
	strcpy(MPM_Mesh.lp.Info,"Voxel lenght GP");
      }
      /* Lagrange Multipliers / Beta (Only LME ) */
      if(strcmp(ShapeFunctionGP,"LME") == 0){
	MPM_Mesh.lambda = MatAllocZ(NumParticles,Ndim);
	strcpy(MPM_Mesh.lambda.Info,"Lagrange Multiplier");
	MPM_Mesh.Beta = MatAllocZ(NumParticles,Ndim);
	strcpy(MPM_Mesh.Beta.Info,"Beta parameter");
      }
    }
    else{
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsSolid2D()","GramsShapeFun no defined");
      exit(0);
    }
    
    /**************************************************/    
    /****** Allocate vectorial/tensorial fields *******/
    /**************************************************/
    MPM_Mesh.Phi = allocate_Fields(NumParticles);

    /**************************************************/
    /********** Generate the initial layout ***********/
    /**************************************************/
    GetInitialGaussPointPosition(MPM_Mesh.Phi.x_GC,MPM_GID_Mesh,GPxElement);

    /**************************************************/
    /*********** Read Material parameters *************/
    /**************************************************/
    if(Is_GramsMaterials){
      MPM_Mesh.NumberMaterials = Counter_Materials;
      MPM_Mesh.MatIdx = (int *)malloc(NumParticles*sizeof(int));
      MPM_Mesh.Mat = GramsMaterials(Name_File,MPM_Mesh,GPxElement);
    }
    else{
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsSolid2D()","GramsMaterials no defined");
      exit(0);
    }

    /**************************************************/
    /********* Read Shape functions parameters ********/
    /**************************************************/
    if(Is_GramsShapeFun){
      GramsShapeFun(Name_File);
      /* Lenght of the Voxel (Only GIMP) */
      if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
	MPM_Mesh.lp = MatAllocZ(NumParticles,NumberDimensions);
	strcpy(MPM_Mesh.lp.Info,"Voxel lenght GP");
      }
      /* Lagrange Multipliers / Beta (Only LME ) */
      if(strcmp(ShapeFunctionGP,"LME") == 0){
	MPM_Mesh.lambda = MatAllocZ(NumParticles,NumberDimensions);
	strcpy(MPM_Mesh.lambda.Info,"Lagrange Multiplier");
	MPM_Mesh.Beta = MatAllocZ(NumParticles,NumberDimensions);
	strcpy(MPM_Mesh.Beta.Info,"Beta parameter");
      }
    }
    else{
      fprintf(stderr,"%s : %s \n",
	      "Error in GramsSolid2D()","GramsShapeFun no defined");
      exit(0);
    }

    /* Fill geometrical properties of the GP mesh */
    for(int i = 0 ; i<MPM_GID_Mesh.NumElemMesh ; i++){
      /* Get the connectivity of the elements vertex */
      Poligon_Connectivity = MPM_GID_Mesh.Connectivity[i];    
      /* Generate a matrix with the poligon coordinates */
      NumVertex = MPM_GID_Mesh.NumNodesElem[i];
      Poligon_Coordinates = MatAllocZ(NumVertex,2); 
      /* Initialize chain interator */
      Vertex = Poligon_Connectivity;
      /* Loop in the chain to fill the poligon */
      for(int k = 0, I_Vertex = 0;
	  (k<NumVertex) || (Vertex != NULL);
	  k++, Vertex = Vertex->next){
	I_Vertex = Vertex->I;
	for(int l = 0 ; l<Ndim ; l++){
	  Poligon_Coordinates.nM[k][l] =
	    MPM_GID_Mesh.Coordinates.nM[I_Vertex][l];
	}
      }
      /* Free data */
      free_Set(MPM_GID_Mesh.Connectivity[i]);
      /* Get the area (Poligon_Centroid.n) 
	 and the position of the centroid (Poligon_Centroid.nV) */
      Area_Element = Area_Poligon(Poligon_Coordinates);
      /* Free data */
      FreeMat(Poligon_Coordinates);
      
      for(int j = 0 ; j<GPxElement ; j++){

	/* Get the index of the material point */
	i_p = i*GPxElement+j;

	/* Get material properties */
	A_p = Area_Element/GPxElement;
	rho_p = MPM_Mesh.Mat[MPM_Mesh.MatIdx[i_p]].rho;
	th_p = MPM_Mesh.Mat[MPM_Mesh.MatIdx[i_p]].thickness;
	m_p = th_p*A_p*rho_p;
	
      	/* Set the initial density */
      	MPM_Mesh.Phi.rho.nV[i_p] = rho_p;	
      	/* Assign the mass parameter */
      	MPM_Mesh.Phi.mass.nV[i_p] = m_p;
	
      	/* Local coordinates of the element */
      	MPM_Mesh.I0[i_p] = -999;
      	MPM_Mesh.NumberNodes[i_p] = 4;
      }
      
    }

    /**************************************************/    
    /************** Read initial values ***************/
    /**************************************************/    
    if(Is_GramsInitials){
      GramsInitials(Name_File,MPM_Mesh,GPxElement);
    }
    else{
      puts("*************************************************");
      puts(" No initial conditions defined ");
    }
    /**************************************************/    
    /************** Read external forces **************/
    /**************************************************/    
    if(Is_GramsNeumannBC){
      MPM_Mesh.NumNeumannBC = Counter_GramsNeumannBC;
      MPM_Mesh.F = GramsNeumannBC(Name_File,Counter_GramsNeumannBC,GPxElement);
    }
    else{
      MPM_Mesh.NumNeumannBC = Counter_GramsNeumannBC;
      puts("*************************************************");
      printf(" \t %s : \n\t %s \n",
	     "* No Neumann boundary conditions defined in",
	     Name_File);
    }
    /**************************************************/    
    /*************** Read body forces *****************/
    /**************************************************/    
    if(Is_GramsBodyForces){
      MPM_Mesh.NumberBodyForces = Counter_BodyForces;
      MPM_Mesh.B = GramsBodyForces(Name_File,Counter_BodyForces,GPxElement); 
    }
    else{
      MPM_Mesh.NumberBodyForces = Counter_BodyForces;
      puts("*************************************************");
      printf(" \t %s : \n\t %s \n",
	     "* No body forces defined in",Name_File);
    }
    /**************************************************/
    /******** INITIALIZE SHAPE FUNCTIONS **************/
    /**************************************************/
    puts("*************************************************");
    if(strcmp(ShapeFunctionGP,"MPMQ4") == 0){
      printf("\t * %s \n","Initialize MPMQ4 shape functions ...");
      Q4_Initialize(MPM_Mesh, FEM_Mesh);
    }
    else if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
      printf("\t * %s \n","Initialize uGIMP shape functions ...");      
      uGIMP_Initialize(MPM_Mesh,FEM_Mesh);
    }
    else if(strcmp(ShapeFunctionGP,"LME") == 0){
      printf("\t * %s \n","Initialize LME shape functions ...");
      LME_Initialize(MPM_Mesh,FEM_Mesh);
    }
    printf("\t %s \n","DONE !!");

   
    /**************************************************/    
    /************* Free the input data ****************/
    /**************************************************/    
    FreeMat(MPM_GID_Mesh.Coordinates);
    free(MPM_GID_Mesh.Connectivity);
    free(MPM_GID_Mesh.NumParticles);
    free(MPM_GID_Mesh.NumNeighbour);
    free(MPM_GID_Mesh.NodeNeighbour);

  } 
  else{
    fprintf(stderr,"%s : %s \n",
	    "Error in GramsSolid2D()","GramsBodyForces no defined");
    exit(0);
  }
  
  return MPM_Mesh;
}

/***************************************************************************/

