#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "ToolsLib/GlobalVariables.h"
#include "InOutFun/InOutFun.h"
#include "ElementsFunctions/ElementTools.h"
#include "GaussPointsFunctions/GaussPointsTools.h"
#include "Constitutive/Constitutive.h"

void main(int argc, char *argv[])
/*
  Inputs parameters :
  * Mesh file
  * Data file
*/
{
  /* Check command-line arguments */
  if(argc == 1){
    perror("Error in main(), insuficient number of input files !");
    exit(0);
  }

  /* Read the .dat file */
  ReadDatFile(argv[1]);

  Matrix D_e = LinearElastic2D(PoissonModulus,ElasticModulus);

  /* Read mesh data and initialize the element mesh */
  Mesh FEM_Mesh = ReadGidMesh(FEM_MeshFileName);

  /* Read and imposse the boundary conditions */
  ReadBCC(BounCondFileName,FEM_Mesh);

  for(int i = 0 ; i<FEM_Mesh.NumNodesBound ; i++){
    for(int j = 1 ; j<=NumberDOF ; j++){
      if(FEM_Mesh.NodesBound[i][j] == 1){
	printf("Node %i -> Value %f, DOF : %i \n",
	       FEM_Mesh.NodesBound[i][0],
	       FEM_Mesh.ValueBC.nM[i][j-1],j-1);
      }
    }
  }

  /* Define Gauss point mesh and initialize it */
  Matrix InputFields;
  InputFields.nM = NULL;
  GaussPoint GP_Mesh  = Initialize_GP_Mesh(MPM_MeshFileName,InputFields,Density,D_e);

  /* Search the GP in the mesh */
  LocateGP(GP_Mesh,FEM_Mesh,0);

  /* Output for Paraview */
  Matrix List_Fields;
  WriteVtk_MPM("Initial_conditions",GP_Mesh,List_Fields,0);  
  WriteVtk_FEM("Mesh",FEM_Mesh,List_Fields,0);
    
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
