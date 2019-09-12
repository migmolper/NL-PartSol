#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "ToolsLib/GlobalVariables.h"
#include "InOutFun/InOutFun.h"
#include "ElementsFunctions/ElementTools.h"
#include "MPM_Subroutines/MPM_Subroutines.h"
#include "GaussPointsFunctions/GaussPointsTools.h"
#include "Boundary_Conditions/BoundaryConditions.h"
#include "Constitutive/Constitutive.h"

int main(int argc, char *argv[])
/*
  Inputs parameters :
  * Mesh file
  * Data file
*/
{  

  /* Variables for the tic-toc */
  clock_t start_t, end_t, total_t;

  /* Start of tic-toc */
  start_t = clock();
  printf("Starting of the program, start_t = %ld\n", start_t);


  /* Variable decalaration for the program */
  
  /* Check command-line arguments */
  if(argc == 1){
    perror("Error in main(), insuficient number of input files !");
    exit(0);
  }

  /***********************************************************************/
  /******************** DEFINE GENERAL VARIABLES *************************/
  /***********************************************************************/
  /* Read the .dat file */
  Read_GeneralParameters(argv[1]);
  /* General variables */
  Matrix D_e; /* Constitutive matrix */
  Mesh FEM_Mesh; /* Define FEM mesh */
  GaussPoint GP_Mesh; /* Define MPM mesh */
  int TimeStep; /* Time step */
  /* Some auxiliar variables for the outputs */
  Matrix List_Fields;
  /* 2D linear elasticity */
  D_e = LinearElastic2D(PoissonModulus,ElasticModulus);
  /***********************************************************************/
  /***********************************************************************/


  /***********************************************************************/
  /******************* DEFINE FINITE ELEMENT MESH ************************/
  /***********************************************************************/
  puts("*************************************************");
  puts(" Generate the MPM mesh");
  puts(" \t Defining background FEM mesh ...");
  FEM_Mesh = ReadGidMesh(FEM_MeshFileName);
  puts(" \t DONE !!!");
  puts(" \t Searching neighbours elements for each node ...");
  FEM_Mesh.NodeNeighbour = GetNodalConnectivity(FEM_Mesh);
  puts(" \t DONE !!!");
  puts(" \t Reading FEM boundary conditions ...");
  FEM_Mesh.Bounds = Set_FEM_BCC(argv[1], FEM_Mesh);
  puts(" \t DONE !!! ");
  /***********************************************************************/
  /***********************************************************************/

  /***********************************************************************/
  /******************** DEFINE GAUSS-POINT MESH **************************/
  /***********************************************************************/
  puts("*************************************************");
  puts(" Generate the MPM mesh");
  puts(" \t Defining MPM mesh of GPs ...");
  GP_Mesh = Define_GP_Mesh(MPM_MeshFileName,Density,D_e);
  puts(" \t DONE !!!");
  puts(" \t Searching GPs in the FEM mesh ...");
  GlobalSearchGaussPoints(GP_Mesh,FEM_Mesh);
  puts(" \t DONE !!!");
  puts(" \t Reading MPM load cases ...");
  GP_Mesh.F = Read_MPM_LoadCase_ExtForces(argv[1],GP_Mesh);
  GP_Mesh.B = Read_MPM_LoadCase_BodyForces(argv[1],GP_Mesh);
  puts(" \t DONE !!!");
  puts(" \t Reading MPM initial conditions ...");
  Read_MPM_InitVal(argv[1],GP_Mesh);
  puts(" \t DONE !!!");  
  /***********************************************************************/
  /***********************************************************************/

  /***********************************************************************/
  /****** INITIALIZE AUXILIAR STRUCTURES TO STORE NODAL INFORMATION ******/
  /***********************************************************************/
  Matrix Nodal_MASS_MOMENTUM; /* Auxiliar variable for the mass and momentum */
  Matrix Nodal_MASS; /* Nodal mass */
  Nodal_MASS.N_rows = 1;
  Nodal_MASS.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_MASS.nM = NULL;
  strcpy(Nodal_MASS.Info,"MASS");
  Matrix Nodal_MOMENTUM; /* Nodal momentum */
  Nodal_MOMENTUM.N_rows = NumberDimensions;
  Nodal_MOMENTUM.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_MOMENTUM.nV = NULL;
  Nodal_MOMENTUM.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  strcpy(Nodal_MOMENTUM.Info,"MOMENTUM");
  Matrix Nodal_VELOCITY; /* Nodal velocities */
  Matrix Nodal_TOT_FORCES; /* Nodal total forces */
  /***********************************************************************/
  /***********************************************************************/


  /***********************************************************************/
  /********************** START THE MPM CALCULUS *************************/
  /***********************************************************************/

  for(TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    puts("*************************************************");
    printf("********************** STEP : %i \n",TimeStep);
    puts("*************************************************");
    
    /* First step : Output Gauss-Points values to Paraview */
    if(TimeStep % ResultsTimeStep == 0){
     puts("*************************************************");
     puts(" First step : Output Gauss-Points values to Paraview");
     puts(" \t WORKING ...");
      WriteVtk_MPM("MPM_VALUES",GP_Mesh,List_Fields,(int)TimeStep/ResultsTimeStep);
      printf(" \t DONE !!! \n");
    }
   
    /* Second step : Get the nodal mass and the momentum */
    puts("*************************************************");
    puts(" Second step : Get the nodal mass and the momentum");
    puts(" \t WORKING ...");
    Nodal_MASS_MOMENTUM = GetNodalMassMomentum(GP_Mesh,FEM_Mesh);
    Nodal_MASS.nV = Nodal_MASS_MOMENTUM.nM[0];
    Nodal_MOMENTUM.nM[0] = Nodal_MASS_MOMENTUM.nM[1];
    Nodal_MOMENTUM.nM[1] = Nodal_MASS_MOMENTUM.nM[2];
    printf(" \t DONE !!! \n");

    /* Third step : Set the essential boundary conditions (over p)*/
    puts("*************************************************");
    puts(" Third step : Set the essential BCC (over P)");
    puts(" \t WORKING ...");
    BCC_Nod_VALUE(FEM_Mesh,Nodal_MOMENTUM,TimeStep);
    puts(" DONE !!!");

    /* Four step : Update the particle stress state */
    puts("*************************************************");
    puts(" Four step : Update the particle stress state");
  
    /* a) Get the grid nodal velocity */
    puts(" \t a) Get the grid nodal velocity ... WORKING");
    Nodal_VELOCITY = GetNodalVelocity(FEM_Mesh,Nodal_MOMENTUM, Nodal_MASS);
    puts(" \t DONE !!!");
    /* b) Calculate the strain increment */
    puts(" \t b) Calculate the strain increment ... WORKING");
    UpdateGaussPointStrain(GP_Mesh,FEM_Mesh,Nodal_VELOCITY);
    puts(" \t DONE !!!");
    /* c) Update the particle stress state */
    puts(" \t c) Update the particle stress state ... WORKING");
    UpdateGaussPointStress(GP_Mesh);
    puts(" \t DONE !!!");

    /* Five step : Calculate total forces */
    puts("*************************************************");
    puts(" Five step : Calculate total forces forces");
    puts(" \t WORKING ...");
    /* BCC_GP_Forces(GP_Mesh, BCC_Loads, TimeStep); */
    Nodal_TOT_FORCES = GetNodalForces(GP_Mesh,FEM_Mesh,TimeStep);
    puts(" DONE !!!");    

    /* Six step : Integrate the grid nodal momentum equation */
    puts("*************************************************");
    puts(" Six step : Integrate the grid nodal momentum equation");
    puts(" \t WORKING ...");
    UpdateGridNodalMomentum(FEM_Mesh,Nodal_MOMENTUM,Nodal_TOT_FORCES);
    puts(" DONE !!!");

    /* Seven step : Update the particle velocity and position */
    puts("*************************************************");
    puts(" Seven step : Update the particle velocity and position");
    puts(" \t WORKING ...");
    UpdateVelocityAndPositionGP(GP_Mesh,FEM_Mesh,Nodal_MASS,
				Nodal_MOMENTUM,Nodal_TOT_FORCES);
    puts(" DONE !!!");

    
    /* Eight step : Search the GP in the mesh */
    puts("*************************************************");
    puts(" Eight step : Search the GP in the mesh");
    puts(" \t WORKING ...");
    LocalSearchGaussPoints(GP_Mesh,FEM_Mesh);
    puts(" DONE !!!");

    /* Nine step : Print nodal values */
    if(TimeStep % ResultsTimeStep == 0){
      puts("*************************************************");
      puts(" Nine step : Print nodal values");
      puts(" \t WORKING ...");
      WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_MOMENTUM,(int)TimeStep/ResultsTimeStep);
      puts(" DONE !!!");
    }
    
    /* Ten step : Store all the material properties in the particles
       so that the deformed grid can be discarted */
    puts("*************************************************");
    puts(" Ten step : Reset nodal values of the mesh");
    puts(" \t WORKING ...");
     /* Reset nodal values */
    FreeMat(Nodal_MASS_MOMENTUM);
    FreeMat(Nodal_VELOCITY);
    FreeMat(Nodal_TOT_FORCES);
    puts(" DONE !!!");

  } /* End of temporal integration */
  
  /***********************************************************************/
  /*********************** END THE MPM CALCULUS **************************/
  /***********************************************************************/


  /* End of tic-toc */
  end_t = clock();
  printf("End of the big loop, end_t = %ld\n", end_t);
  total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("Total time taken by CPU: %ld\n", total_t  );
  puts("Exiting of the program...");
     
      
  exit(EXIT_SUCCESS);
  
}

