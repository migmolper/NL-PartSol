#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"
#include "Utils.h"

int main(int argc, char * argv[])
/*
  Inputs parameters : Data file .gdf
*/
{
  
  /*********************************************************************/
  /******************* Check command-line arguments ********************/
  /*********************************************************************/
  if(argc == 4){
    SimulationFile = argv[3];
    Formulation = argv[2];
    if(strcmp(argv[1],"-2D") == 0){
      NumberDimensions = 2;
    }
  }

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
  if(strcmp(TimeIntegration,"FE") == 0 ){
    u_ForwardEuler(FEM_Mesh, MPM_Mesh);
  }
  if(strcmp(TimeIntegration,"GA") == 0 ){
    u_GeneralizedAlpha(FEM_Mesh, MPM_Mesh);
  }
  if(strcmp(TimeIntegration,"PCE") == 0 ){
    U_PCE(FEM_Mesh, MPM_Mesh);
  }
  
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

