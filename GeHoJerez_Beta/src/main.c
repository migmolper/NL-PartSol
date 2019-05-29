#include <stdio.h>
#include <stdlib.h>
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

void main(int argc, char *argv[])
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

    Matrix D_e;
  Mesh FEM_Mesh;
  Matrix InputFields;
  InputFields.nM = NULL;
  GaussPoint GP_Mesh;
  char Get_values[MAXW];
  Matrix Nod_Values;
  Matrix Nodal_MASS;
  Matrix Nodal_VELOCITY;
  Nodal_VELOCITY.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  Matrix Nodal_MOMENTUM;
  Nodal_MOMENTUM.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  Matrix Nodal_F_TOT;
  Nodal_F_TOT.nM = (double **)malloc((unsigned)NumberDimensions*sizeof(double*));


  /**/
  D_e = LinearElastic2D(PoissonModulus,ElasticModulus);

  /* Read mesh data and initialize the element mesh */
  FEM_Mesh = ReadGidMesh(FEM_MeshFileName);

  /* Read and imposse the boundary conditions */
  ReadBCC(BounCondFileName,FEM_Mesh);

  /* Define Gauss point mesh and initialize it */
  GP_Mesh  = Initialize_GP_Mesh(MPM_MeshFileName,InputFields,Density,D_e);

  /* Search the GP in the mesh */
  LocateGP(GP_Mesh,FEM_Mesh,0);

  /* First step : Get the nodal mass and the momentum */
  printf("************************************************* \n");
  printf(" First step : Get the nodal mass and the momentum \n");
  printf(" \t Working ... \n");
  Get_values = "MASS;MOMENTUM";
  Nod_Values = GetNodalValuesFromGP(GP_Mesh,FEM_Mesh,Get_values);
  Nodal_MASS.nV = Nod_Values.nM[0];
  Nodal_MOMENTUM.nM[0] = Nod_Values.nM[1];
  Nodal_MOMENTUM.nM[1] = Nod_Values.nM[2];
  printf(" DONE !!! \n");

  /* Second step : Set the essential boundary conditions (over p)*/
  printf("************************************************* \n");
  printf(" Second step : Set the essential BCC (over P) \n");
  printf(" \t Working ... \n");
  ApplyBoundaryCondition_Nod(FEM_Mesh,Nodal_MOMENTUM,0);
  printf(" DONE !!! \n");

  /* Third step (USF only) : Update the particle stress state */
  printf("************************************************* \n");
  printf(" Third step : Update the particle stress state \n");
  
  /* a) Get the grid nodal velocity */
  printf(" \t a) Get the grid nodal velocity ... WORKING \n");
  Nodal_VELOCITY = GetNodalVelocity(Nodal_MOMENTUM, Nodal_MASS);
  /* b) Calculate the strain increment */
  
  /* c) Update the particle density */

  /* d) Update the particle stress state */

  /* Four step : Calculate the nodal internal, external and the total forces */

  /* Five step : Integrate the grid nodal momentum equation */

  /* Six step : Update the particle velocity and position */

  /* Seven step (MUSL only) : Recalculate the grid nodal momentum */

  /* Eight step (MUSL & USL only) : */

  /* Nine step : Store all the material properties in the particles so that the deformed grid
   can be discarted */
  
  /* Output for Paraview */
  Matrix List_Fields;
  WriteVtk_MPM("Initial_conditions",GP_Mesh,List_Fields,0);
  WriteVtk_FEM("Mesh",FEM_Mesh,Nod_Values,0);
    
  /* /\* Read the initial conditions fields as a CSV *\/ */
  /* Matrix InputFields = Read_CSV(InitCondFileName,9); */
  
  exit(EXIT_SUCCESS);
  
  
}


  /* /\********************************************************************\/ */
  /* /\* Print mesh properties *\/ */
  /* printf(" * Set Mesh properties : \n"); */
  /* printf("\t -> Number of nodes : %i \n",ElementMesh.NumNodesMesh); */
  /* printf("\t -> Number of boundary nodes : %i \n",ElementMesh.NumNodesBound); */
  /* printf("\t -> Number of elements : %i \n",ElementMesh.NumElemMesh); */
  /* printf("\t -> Order of the element : %s \n",ElementMesh.TypeElem); */
  /* printf("\t -> Number of nodes per element : %i \n",ElementMesh.NumNodesElem); */
  /* printf(" * Boundary nodes : \n"); */
  /* for(int i = 0 ; i<ElementMesh.NumNodesBound ; i++){ */
  /*   printf(" %i ",ElementMesh.NodesBound[i]); */
  /*   if( ((i+1)%4 == 0)) printf("\n"); */
  /* } */
  /* printf("\n"); */
  /* printf(" * Nodal coordinates : \n"); */
  /* for(int i = 0; i<ElementMesh.NumNodesMesh ; i++){ */
  /*   printf("\t [%f, %f, %f] \n", */
  /* 	   ElementMesh.Coordinates.nM[i][0], */
  /* 	   ElementMesh.Coordinates.nM[i][1], */
  /* 	   ElementMesh.Coordinates.nM[i][2]); */
  /* } */
  /* printf(" * Connectivity mesh : \n"); */
  /* for(int i =0; i<ElementMesh.NumElemMesh; i++){ */
  /*   printf("\t Element(%i) : [",i); */
  /*   for(int j = 0 ; j<ElementMesh.NumNodesElem; j++){ */
  /*     printf(" %i ",ElementMesh.Connectivity[i][j]); */
  /*   } */
  /*   printf("]\n"); */
  /* } */
  /* /\* Final messages *\/ */
  /* printf("End of set mesh properties !!! \n"); */
  /* /\********************************************************************\/ */
