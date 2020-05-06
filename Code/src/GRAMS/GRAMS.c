#include "grams.h"

void TIC_toc(clock_t);
void tic_TOC(clock_t);
void globalfree(Mesh,GaussPoint);

int main(int argc, char * argv[])
/*
  GraMS V2.3
*/
{
  
  /*********************************************************************/
  /******************* Check command-line arguments ********************/
  /*********************************************************************/
  if(argc == 4){

    /* Assign Number of dimensions */
    if(strcmp(argv[1],"-2D") == 0){
      NumberDimensions = 2;
    }
    
    /* Read simulation file and kind of simulation */
    SimulationFile = argv[3];
    Formulation = argv[2];
    
    /* Assign degree of freedom */
    if((strcmp(Formulation,"-V") == 0) &&
       (NumberDimensions == 2)){
      NumberDOF = 2;
    }

    /*********************************************************************/
    /******************* CRONOGRAPH CALCULUS : INIT **********************/
    /*********************************************************************/
    clock_t start_t;
    TIC_toc(start_t);
  
    /*********************************************************************/
    /**************** DEFINE CALCULUS MESH and BCCs **********************/
    /*********************************************************************/
    puts("*************************************************");
    puts("Generating the background mesh ...");
    Mesh FEM_Mesh = GramsBox(SimulationFile);

    /*********************************************************************/
    /******************* TIME INTEGRATION SCHEME *************************/
    /*********************************************************************/
    puts("*************************************************");
    puts("Read time integration scheme ...");
    GramsTime(SimulationFile);

    /*********************************************************************/
    /******************* DEFINE GAUSS-POINT MESH *************************/
    /*********************************************************************/
    puts("*************************************************");
    puts("Generating MPM mesh ...");
    GaussPoint MPM_Mesh = GramsSolid2D(SimulationFile,FEM_Mesh);

    /*********************************************************************/
    /*********************** OUTPUT VARIABLES ****************************/
    /*********************************************************************/
    puts("*************************************************");
    puts("Read outputs ...");
    GramsOutputs(SimulationFile);
    
    /*********************************************************************/
    /********************** RUN THE MPM CALCULUS *************************/
    /*********************************************************************/
    puts("*************************************************");
    puts("Run simulation ...");
    if(strcmp(TimeIntegration,"FE") == 0 ){ /* Forward Euler */
      U_FE(FEM_Mesh, MPM_Mesh);
    }
    if(strcmp(TimeIntegration,"GA") == 0 ){ /* Generalized-alpha */
      U_GA(FEM_Mesh, MPM_Mesh);
    }
    if(strcmp(TimeIntegration,"PCE") == 0 ){ /* Explicit predictor-corrector */
      U_PCE(FEM_Mesh, MPM_Mesh);
    }

        
    /*********************************************************************/
    /************************* FREE ALL FIELDS ***************************/
    /*********************************************************************/
    puts("*************************************************");
    puts("Free memory ...");
    globalfree(FEM_Mesh, MPM_Mesh);
  
    /*********************************************************************/
    /******************** CRONOGRAPH CALCULUS : END **********************/
    /*********************************************************************/
    tic_TOC(start_t);
         
    /*********************************************************************/
    /************************ CLOSE THE PROGRAM **************************/
    /*********************************************************************/
    puts("Exiting of the program...");
    exit(EXIT_SUCCESS);
    
  }
  /* Gauss Point analysis */
  else if((argc == 3) &&
	  (strcmp(argv[1],"-GPA") == 0)){
    
    /* Read simulation file and kind of simulation */
    SimulationFile = argv[2];

    /*********************************************************************/
    /******************* CRONOGRAPH CALCULUS : INIT **********************/
    /*********************************************************************/
    clock_t start_t;
    TIC_toc(start_t);


    /*********************************************************************/
    /******************** CRONOGRAPH CALCULUS : END **********************/
    /*********************************************************************/
    tic_TOC(start_t);
             
    /*********************************************************************/
    /************************ CLOSE THE PROGRAM **************************/
    /*********************************************************************/
    puts("Exiting of the program...");
    exit(EXIT_SUCCESS);
  }
  /* Non-defined simulation */
  else{
    
    puts("Undefined kind of simulation");
    exit(EXIT_FAILURE);
    
  }
  
}



/*********************************************************************/

/* Time measure function declaration */

void TIC_toc(clock_t start_t){

  /* Start of tic-toc */
  start_t = clock();
  printf("Starting of the program, start_t = %ld\n", start_t);
  
}

/*********************************************************************/

void tic_TOC(clock_t start_t){

  /* Define internal parameters */
  clock_t end_t, total_t;

  /* End of tic-toc */
  end_t = clock();
  printf("End of the big loop, end_t = %ld\n", end_t);
  total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("Total time taken by CPU: %ld\n", total_t  );
  
}

/*********************************************************************/

void globalfree(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
/*
  Function to free the reamaining memory
*/
{
  /* Free malloc in FEM_Mesh */
  FreeMat(FEM_Mesh.Coordinates);
  free(FEM_Mesh.NumNodesElem);
  free_SetTable(FEM_Mesh.Connectivity,FEM_Mesh.NumElemMesh);
  free(FEM_Mesh.NumNeighbour);
  free_SetTable(FEM_Mesh.NodeNeighbour,FEM_Mesh.NumNodesMesh);
  free(FEM_Mesh.SizeNodalLocality);
  free_SetTable(FEM_Mesh.NodalLocality,FEM_Mesh.NumNodesMesh);
  free(FEM_Mesh.NumParticles);
  free_SetTable(FEM_Mesh.I_particles,FEM_Mesh.NumNodesMesh);
  
  /* FEM_Mesh.Bounds */

  /* Free malloc in MPM_Mesh */
  free(MPM_Mesh.I0);
  free(MPM_Mesh.NumberNodes);
  free_SetTable(MPM_Mesh.ListNodes,MPM_Mesh.NumGP);
  free_SetTable(MPM_Mesh.Beps,MPM_Mesh.NumGP);
  free_Fields(MPM_Mesh.Phi);
  /* free_Fields(MPM_Mesh.Phi_n0);   */
  free(MPM_Mesh.MatIdx);

  /* Material * Mat; */

  /* Load * F; */
  /* Load * B; */

  if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
    FreeMat(MPM_Mesh.lp);
  }
  /* Lagrange Multipliers / Beta (Only LME ) */
  if(strcmp(ShapeFunctionGP,"LME") == 0){
    FreeMat(MPM_Mesh.lambda);
    FreeMat(MPM_Mesh.Beta);
  }
      
}
