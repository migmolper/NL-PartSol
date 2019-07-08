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
  printf("************************************************* \n");
  printf(" Generate the MPM mesh \n");
  printf(" \t Defining background FEM mesh ... \n");
  FEM_Mesh = ReadGidMesh(FEM_MeshFileName);
  printf(" \t DONE !!! \n");
  printf(" \t Searching neighbours elements for each node ... \n");
  FEM_Mesh.NodeNeighbour = GetNodalConnectivity(FEM_Mesh);
  printf(" \t DONE !!! \n");
  printf(" \t Reading FEM boundary conditions ... \n");
  FEM_Mesh.Bounds = Set_FEM_BCC(argv[1], FEM_Mesh);
  printf(" \t DONE !!! \n");
  /***********************************************************************/
  /***********************************************************************/

  /***********************************************************************/
  /******************** DEFINE GAUSS-POINT MESH **************************/
  /***********************************************************************/
  printf("************************************************* \n");
  printf(" Generate the MPM mesh \n");
  printf(" \t Defining MPM mesh of GPs ... \n");
  GP_Mesh = Define_GP_Mesh(MPM_MeshFileName,Density,D_e);
  printf(" \t DONE !!! \n");
  printf(" \t Searching GPs in the FEM mesh ... \n");
  GlobalSearchGaussPoints(GP_Mesh,FEM_Mesh);
  printf(" \t DONE !!! \n");
  printf(" \t Reading MPM load cases ... \n");
  GP_Mesh.F = Read_MPM_LoadCase_ExtForces(argv[1],GP_Mesh);
  GP_Mesh.B = Read_MPM_LoadCase_BodyForces(argv[1],GP_Mesh);
  printf(" \t DONE !!! \n");
  printf(" \t Reading MPM initial conditions ... \n");
  Read_MPM_InitVal(argv[1],GP_Mesh);
  printf(" \t DONE !!! \n");  
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

    printf("************************************************* \n");
    printf("********************** STEP : %i \n",TimeStep);
    printf("************************************************* \n");
    
    /* First step : Output Gauss-Points values to Paraview */
    printf("************************************************* \n");
    printf(" First step : Output Gauss-Points values to Paraview \n");
    printf(" \t WORKING ... \n");
    WriteVtk_MPM("MPM_VALUES",GP_Mesh,List_Fields,TimeStep);
    printf(" \t DONE !!! \n");
   
    /* Second step : Get the nodal mass and the momentum */
    printf("************************************************* \n");
    printf(" Second step : Get the nodal mass and the momentum \n");
    printf(" \t WORKING ... \n");
    Nodal_MASS_MOMENTUM = GetNodalMassMomentum(GP_Mesh,FEM_Mesh);
    Nodal_MASS.nV = Nodal_MASS_MOMENTUM.nM[0];
    Nodal_MOMENTUM.nM[0] = Nodal_MASS_MOMENTUM.nM[1];
    Nodal_MOMENTUM.nM[1] = Nodal_MASS_MOMENTUM.nM[2];
    printf(" \t DONE !!! \n");

    /* Third step : Set the essential boundary conditions (over p)*/
    printf("************************************************* \n");
    printf(" Third step : Set the essential BCC (over P) \n");
    printf(" \t WORKING ... \n");
    BCC_Nod_VALUE(FEM_Mesh,Nodal_MOMENTUM,TimeStep);
    printf(" DONE !!! \n");

    /* Four step : Update the particle stress state */
    printf("************************************************* \n");
    printf(" Four step : Update the particle stress state \n");
  
    /* a) Get the grid nodal velocity */
    printf(" \t a) Get the grid nodal velocity ... WORKING \n");
    Nodal_VELOCITY = GetNodalVelocity(FEM_Mesh,Nodal_MOMENTUM, Nodal_MASS);
    printf(" \t DONE !!! \n");
    /* b) Calculate the strain increment */
    printf(" \t b) Calculate the strain increment ... WORKING \n");
    UpdateGaussPointStrain(GP_Mesh,FEM_Mesh,Nodal_VELOCITY);
    printf(" \t DONE !!! \n");
    /* c) Update the particle stress state */
    printf(" \t c) Update the particle stress state ... WORKING \n");
    UpdateGaussPointStress(GP_Mesh);
    printf(" \t DONE !!! \n");

    /* Five step : Calculate total forces */
    printf("************************************************* \n");
    printf(" Five step : Calculate total forces forces\n");
    printf(" \t WORKING ... \n");
    /* BCC_GP_Forces(GP_Mesh, BCC_Loads, TimeStep); */
    Nodal_TOT_FORCES = GetNodalForces(GP_Mesh,FEM_Mesh,TimeStep);
    printf(" DONE !!! \n");    

    /* Six step : Integrate the grid nodal momentum equation */
    printf("************************************************* \n");
    printf(" Six step : Integrate the grid nodal momentum equation\n");
    printf(" \t WORKING ... \n");
    UpdateGridNodalMomentum(FEM_Mesh,Nodal_MOMENTUM,Nodal_TOT_FORCES);
    printf(" DONE !!! \n");

    /* Seven step : Update the particle velocity and position */
    printf("************************************************* \n");
    printf(" Seven step : Update the particle velocity and position\n");
    printf(" \t WORKING ... \n");
    UpdateVelocityAndPositionGP(GP_Mesh,FEM_Mesh,Nodal_MASS,
				Nodal_MOMENTUM,Nodal_TOT_FORCES);
    printf(" DONE !!! \n");

    
    /* Eight step : Search the GP in the mesh */
    printf("************************************************* \n");
    printf(" Eight step : Search the GP in the mesh \n");
    printf(" \t WORKING ... \n");
    LocalSearchGaussPoints(GP_Mesh,FEM_Mesh);
    printf(" DONE !!! \n");

    /* Nine step : Print nodal values */
    printf("************************************************* \n");
    printf(" Nine step : Print nodal values \n");
    printf(" \t WORKING ... \n");
    WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_MOMENTUM,TimeStep);
    printf(" DONE !!! \n");
    
    /* Ten step : Store all the material properties in the particles
       so that the deformed grid can be discarted */
    printf("************************************************* \n");
    printf(" Ten step : Reset nodal values of the mesh \n");
    printf(" \t WORKING ... \n");
     /* Reset nodal values */
    FreeMat(Nodal_MASS_MOMENTUM);
    FreeMat(Nodal_VELOCITY);
    FreeMat(Nodal_TOT_FORCES);
    printf(" DONE !!! \n");

  } /* End of temporal integration */
  
  /***********************************************************************/
  /*********************** END THE MPM CALCULUS **************************/
  /***********************************************************************/


  /* End of tic-toc */
  end_t = clock();
  printf("End of the big loop, end_t = %ld\n", end_t);
  total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("Total time taken by CPU: %ld\n", total_t  );
  printf("Exiting of the program...\n");
     
      
  exit(EXIT_SUCCESS);
  
}

