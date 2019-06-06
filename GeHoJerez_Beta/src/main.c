#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

  /* Some auxiliar variables */
  Matrix InputFields;
  InputFields.nM = NULL;
  Matrix List_Fields;

  /* 2D linear elasticity */
  D_e = LinearElastic2D(PoissonModulus,ElasticModulus);

  /* Read mesh data and initialize the element mesh */
  FEM_Mesh = ReadGidMesh(FEM_MeshFileName);

  /***********************************************************************/
  /****** INITIALIZE AUXILIAR STRUCTURES TO STORE NODAL INFORMATION ******/
  /***********************************************************************/
  Matrix Nod_Values;
  Matrix Nodal_MASS;
  Nodal_MASS.N_rows = 1;
  Nodal_MASS.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_MASS.nM = NULL;
  strcpy(Nodal_MASS.Info,"MASS");
  Matrix Nodal_VELOCITY;
  Nodal_VELOCITY.N_rows = NumberDimensions;
  Nodal_VELOCITY.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_VELOCITY.nV = NULL;
  strcpy(Nodal_VELOCITY.Info,"VELOCITY");
  Matrix Nodal_MOMENTUM;
  Nodal_MOMENTUM.N_rows = NumberDimensions;
  Nodal_MOMENTUM.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_MOMENTUM.nV = NULL;
  Nodal_MOMENTUM.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  strcpy(Nodal_MOMENTUM.Info,"MOMENTUM");
  Matrix Nodal_INT_FORCES; /* Nodal values of the internal forces */
  Nodal_INT_FORCES.N_rows = NumberDimensions;
  Nodal_INT_FORCES.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_INT_FORCES.nV = NULL;
  Nodal_INT_FORCES.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  strcpy(Nodal_INT_FORCES.Info,"Nodal Internal forces");
  Matrix Nodal_GRAVITY_FORCES; /* Nodal values of the gravity forces */
  Nodal_GRAVITY_FORCES.N_rows = NumberDimensions;
  Nodal_GRAVITY_FORCES.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_GRAVITY_FORCES.nV = NULL;
  Nodal_GRAVITY_FORCES.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  strcpy(Nodal_GRAVITY_FORCES.Info,"F_GRAV");
  Matrix Nodal_TOT_FORCES; /* Nodal total forces */
  Nodal_TOT_FORCES.nV = NULL;
  Nodal_TOT_FORCES.nM = (double **)malloc((unsigned)NumberDimensions*2*sizeof(double*));
  strcpy(Nodal_TOT_FORCES.Info,"F_TOT");
  /***********************************************************************/
  /***********************************************************************/

  /***********************************************************************/
  /****************** INITIALIZE GAUSS-POINT MESH ************************/
  /***********************************************************************/
  /* Read and imposse the boundary conditions */
  ReadBCC(BounCondFileName,FEM_Mesh);
  /* Define Gauss point mesh and initialize it */
  GP_Mesh  = Initialize_GP_Mesh(MPM_MeshFileName,InputFields,Density,D_e);
  /***********************************************************************/
  /***********************************************************************/


  /* Zero step : Search the GP in the mesh */
  printf("************************************************* \n");
  printf(" Zero step : Search the GP in the mesh \n");
  printf(" \t WORKING ... \n");
  LocateGP(GP_Mesh,FEM_Mesh,0);
  printf(" DONE !!! \n");
 
  /***********************************************************************/
  /********************** START THE MPM CALCULUS *************************/
  /***********************************************************************/

  for(TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    printf("************************************************* \n");
    printf("********************** STEP ********************* \n");    

    /* Zero step : Output Gauss-Points values to Paraview */
    printf("************************************************* \n");
    printf(" Zero step : Output Gauss-Points values to Paraview \n");
    printf(" \t WORKING ... \n");
    WriteVtk_MPM("MPM_VALUES",GP_Mesh,List_Fields,TimeStep);
    printf(" DONE !!! \n");
  
  
    /* First step : Get the nodal mass and the momentum */
    printf("************************************************* \n");
    printf(" First step : Get the nodal mass and the momentum \n");
    printf(" \t WORKING ... \n");
    char Get_values_1[MAXW] = "MASS;MOMENTUM";
    Nod_Values = GetNodalValuesFromGP(GP_Mesh,FEM_Mesh,Get_values_1);
    Nodal_MASS.nV = Nod_Values.nM[0];
    Nodal_MOMENTUM.nM[0] = Nod_Values.nM[1];
    Nodal_MOMENTUM.nM[1] = Nod_Values.nM[2];

    printf(" DONE !!! \n");

    /* Second step : Set the essential boundary conditions (over p)*/
    printf("************************************************* \n");
    printf(" Second step : Set the essential BCC (over P) \n");
    printf(" \t WORKING ... \n");
    ApplyBoundaryCondition_Nod(FEM_Mesh,Nodal_MOMENTUM,TimeStep);
    printf(" DONE !!! \n");

    /* Third step : Update the particle stress state */
    printf("************************************************* \n");
    printf(" Third step : Update the particle stress state \n");
  
    /* a) Get the grid nodal velocity */
    printf(" \t a) Get the grid nodal velocity ... WORKING \n");
    Nodal_VELOCITY = GetNodalVelocity(Nodal_MOMENTUM, Nodal_MASS);
    printf(" \t DONE !!! \n");
    /* b) Calculate the strain increment */
    printf(" \t b) Calculate the strain increment ... WORKING \n");
    UpdateGaussPointStrain(GP_Mesh,FEM_Mesh,Nodal_VELOCITY);
    printf(" \t DONE !!! \n");
    /* c) Update the particle density */
    printf(" \t c) Update the particle density ... WORKING \n");
    UpdateGaussPointDensity(GP_Mesh);
    printf(" \t DONE !!! \n");
    /* d) Update the particle stress state */
    printf(" \t d) Update the particle stress state ... WORKING \n");
    UpdateGaussPointStress(GP_Mesh,D_e);
    printf(" \t DONE !!! \n");

    /* Four step : Calculate the nodal internal, external forces */
    printf("************************************************* \n");
    printf(" Four step : Calculate internal and external forces\n");
    printf(" \t WORKING ... \n");
    char Get_values_2[MAXW] = "F_INT;F_GRAV";
    Nodal_TOT_FORCES = GetNodalValuesFromGP(GP_Mesh,FEM_Mesh,Get_values_2);
    Nodal_INT_FORCES.nM[0] = Nodal_TOT_FORCES.nM[0];
    Nodal_INT_FORCES.nM[1] = Nodal_TOT_FORCES.nM[1];
    Nodal_GRAVITY_FORCES.nM[0] = Nodal_TOT_FORCES.nM[2];
    Nodal_GRAVITY_FORCES.nM[1] = Nodal_TOT_FORCES.nM[3];
    printf(" DONE !!! \n");    

    /* Five step : Integrate the grid nodal momentum equation */
    printf("************************************************* \n");
    printf(" Five step : Integrate the grid nodal momentum equation\n");
    printf(" \t WORKING ... \n");
    UpdateGridNodalMomentum(FEM_Mesh,Nodal_MOMENTUM,Nodal_TOT_FORCES);
    printf(" DONE !!! \n");

    /* Six step : Update the particle velocity and position */
    printf("************************************************* \n");
    printf(" Six step : Update the particle velocity and position\n");
    printf(" \t WORKING ... \n");
    UpdateVelocityAndPositionGP(GP_Mesh,FEM_Mesh,Nodal_MASS,
				Nodal_MOMENTUM,Nodal_TOT_FORCES);
    printf(" DONE !!! \n");

    /* Seven step : Store all the material properties in the particles so that the deformed grid
       can be discarted */
    printf("************************************************* \n");
    printf(" Seven step : Reset nodal values of the mesh \n");
    printf(" \t WORKING ... \n");
    /* Print nodal values*/
    WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_MOMENTUM,TimeStep);
    /* Deallocate the array */
    free(Nod_Values.nM);
    free(Nodal_MASS.nV);
    free(Nodal_MOMENTUM.nM);
    free(Nodal_VELOCITY.nM);
    free(Nodal_INT_FORCES.nM);
    free(Nodal_GRAVITY_FORCES.nM);
    /* Allocate it again */
    Nodal_MOMENTUM.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
    Nodal_INT_FORCES.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
    Nodal_GRAVITY_FORCES.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
    Nodal_TOT_FORCES.nM = (double **)malloc((unsigned)NumberDimensions*2*sizeof(double*));  
    
    printf(" DONE !!! \n");

    /* Eight step : Search the GP in the mesh */
    printf("************************************************* \n");
    printf(" Eight step : Search the GP in the mesh \n");
    printf(" \t WORKING ... \n");
    LocateGP(GP_Mesh,FEM_Mesh,TimeStep);
    printf(" DONE !!! \n");


  } /* End of temporal integration */


  
  /***********************************************************************/
  /***********************************************************************/
   
  /* WriteVtk_FEM("Mesh",FEM_Mesh,NULL,0); */
  /* WriteVtk_FEM("Mesh_Equilibrium",FEM_Mesh,Nodal_TOT_FORCES,0); */
  
      
  exit(EXIT_SUCCESS);
  
}

