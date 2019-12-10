#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "ToolsLib/GlobalVariables.h"
#include "InOutFun/InOutFun.h"
#include "ShapeFun/ShapeFun.h"
#include "MPM_Subroutines/MPM_Subroutines.h"
#include "Formulations/Formulations.h"

int main(int argc, char * argv[])
/*
  Inputs parameters : Data file .gdf
*/
{
  
  /*********************************************************************/
  /******************* Check command-line arguments ********************/
  /*********************************************************************/
  if(argc == 1){
    perror("Error in main(), insuficient number of input files !");
    exit(0);
  }

  /*********************************************************************/
  /******************* CRONOGRAPH CALCULUS : INIT **********************/
  /*********************************************************************/
  clock_t start_t;
  TIC_toc(start_t);
  
  /*********************************************************************/
  /******************* DEFINE GENERAL VARIABLES ************************/
  /*********************************************************************/
  /* Read the .gdf file */
  Read_GeneralParameters(argv[1]);

  /*********************************************************************/
  /********************* DEFINE CALCULUS MESH **************************/
  /*********************************************************************/
  Mesh FEM_Mesh = InitializeMesh(argv[1]);

  /*********************************************************************/
  /******************* DEFINE GAUSS-POINT MESH *************************/
  /*********************************************************************/
  GaussPoint GP_Mesh = InitializeGP(argv[1], FEM_Mesh);

  /*********************************************************************/
  /********************** RUN THE MPM CALCULUS *************************/
  /*********************************************************************/
  u_ForwardEuler(FEM_Mesh, GP_Mesh);

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

