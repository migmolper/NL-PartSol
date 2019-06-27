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

  /* Read the .dat file */
  ReadDatFile(argv[1]);

  /* General variables */
  Matrix D_e; /* Constitutive matrix */
  Mesh FEM_Mesh; /* Define FEM mesh */
  GaussPoint GP_Mesh; /* Define MPM mesh */
  int TimeStep; /* Time step */

  /* Some auxiliar variables for the outputs */
  Matrix InputFields;
  InputFields.nM = NULL;
  Matrix List_Fields;

  /* 2D linear elasticity */
  D_e = LinearElastic2D(PoissonModulus,ElasticModulus);

  /* Read mesh data and initialize the element mesh */
  FEM_Mesh = ReadGidMesh(FEM_MeshFileName);
  /* Get the Connectivity of each node */
  FEM_Mesh.NodeNeighbour = GetNodalConnectivity(FEM_Mesh);

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
  /******** INITIALIZE AUXILIAR STRUCTURES FOR BOUNDARY CONDITIONS *******/
  /***********************************************************************/
  /* Define domain and its boundary conditions */
  /* BoundaryConditions BCC_Momentum_TOP; */
  BoundaryConditions BCC_Momentum_BOTTOM;
  /* BoundaryConditions BCC_Momentum_LEFT; */
  /* BoundaryConditions BCC_Momentum_RIGHT;  */

  /* Define loads */
  /* BoundaryConditions BCC_Loads; */
  
  /* Read the boundary conditions in the contours of the domain */
  /* BCC_Momentum_TOP = ReadBCC(BounCondFileName,FEM_Mesh); */
  BCC_Momentum_BOTTOM = ReadBCC(BounCondFileName,FEM_Mesh);
  /* BCC_Momentum_LEFT = ReadBCC(BounCondFileName,FEM_Mesh); */
  /* BCC_Momentum_RIGHT = ReadBCC(BounCondFileName,FEM_Mesh); */
  
  /* Read boundary conditions in the MPM nodes */
  /* BCC_Loads = ReadBCC(LoadsFileName,FEM_Mesh); */
  /***********************************************************************/
  /***********************************************************************/

  /***********************************************************************/
  /****************** INITIALIZE GAUSS-POINT MESH ************************/
  /***********************************************************************/
  /* Define Gauss point mesh and initialize it */
  GP_Mesh  = Initialize_GP_Mesh(MPM_MeshFileName,InputFields,Density,D_e);
  /***********************************************************************/
  /***********************************************************************/

  /* PrintMatrix(GP_Mesh.Phi.x_GC,GP_Mesh.NumGP,2); */
  /* PrintMatrix(FEM_Mesh.Coordinates,FEM_Mesh.NumNodesMesh,2); */

  printf("************************************************* \n");
  printf(" Begin of the global search of the GP over the mesh \n");
  printf(" \t WORKING ... \n");
  GlobalSearchGaussPoints(GP_Mesh,FEM_Mesh);
  printf(" DONE !!! \n");
  
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
    printf(" DONE !!! \n");
   
    /* Second step : Get the nodal mass and the momentum */
    printf("************************************************* \n");
    printf(" Second step : Get the nodal mass and the momentum \n");
    printf(" \t WORKING ... \n");
    Nodal_MASS_MOMENTUM = GetNodalMassMomentum(GP_Mesh,FEM_Mesh);
    Nodal_MASS.nV = Nodal_MASS_MOMENTUM.nM[0];
    Nodal_MOMENTUM.nM[0] = Nodal_MASS_MOMENTUM.nM[1];
    Nodal_MOMENTUM.nM[1] = Nodal_MASS_MOMENTUM.nM[2];
    printf(" DONE !!! \n");

    /* Third step : Set the essential boundary conditions (over p)*/
    printf("************************************************* \n");
    printf(" Third step : Set the essential BCC (over P) \n");
    printf(" \t WORKING ... \n");
    BCC_Nod_Momentum(FEM_Mesh,BCC_Momentum_BOTTOM,Nodal_MOMENTUM,TimeStep);
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
    Nodal_TOT_FORCES = GetNodalForces(GP_Mesh,FEM_Mesh);
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
