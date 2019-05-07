#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "ToolsLib/GlobalVariables.h"
#include "InOutFun/InOutFun.h"
#include "Constitutive/Constitutive.h"
#include "ElementsFunctions/ElementTools.h"
#include "GaussPointsFunctions/GaussPointsTools.h"
#include "Solvers/Solvers.h"

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

  /* Read mesh data and initialize the element mesh */
  Element ElementMesh = ReadGidMesh(MeshFileName);
  ElementMesh.NumberDOF = 2; /* Change this*/
    
  /* Read the initial conditions fields as a CSV */
  Matrix InputFields = Read_CSV(InitCondFileName, 5);
    
  /* Read the data file */  /* Physical parameters */
  double YoungModulus = 1;
  
  /* Get material Constitutive matrix */
  Matrix D = LinearElastic1D(YoungModulus);

  /* Define Gauss point mesh */
  GaussPoint GP_Mesh  = Initialize_GP_Mesh(InputFields,D,
					   ElementMesh);

  /* Localize all the Gauss-Points */
  UpdateElementLocationGP(GP_Mesh,ElementMesh);

  /* Get the mass matrix */
  Matrix M = Get_Geom_Mass_Matrix(GP_Mesh,ElementMesh);

  /* Get the Lumped-Mass matrix */
  Matrix M_l = Get_Lumped_Matrix(M);

  for(int i = 0; i<ElementMesh.NumNodesMesh ; i++){
    printf("Nodal coordinates : [%f, %f, %f] \n",
	   ElementMesh.Coordinates[i][0],
	   ElementMesh.Coordinates[i][1],
	   ElementMesh.Coordinates[i][2]);
  }

  printf("Connectivity mesh :\n");
  for(int i =0; i<ElementMesh.NumElemMesh; i++){
    printf("Element  %i : ",i+1);
    for(int j = 0 ; j<ElementMesh.NumNodesElem; j++){
      printf(" %i ",ElementMesh.Connectivity[i][j]);
    }
    printf(" -> Element active : %i \n",ElementMesh.ActiveElem[i]);
  }

  printf("D : %f \n",D.n);

  printf("Gauss points initial conditions : \n");
  printf("X_GP \t V_X \t SIGMA_X \n");
  for(int i = 0; i<5 ; i++){
    printf("%f \t %f \t %f \n",
	   InputFields.nM[0][i],
	   InputFields.nM[1][i],
	   InputFields.nM[2][i]);
  }

  Matrix X_g = MatAlloc(2,1);
  Matrix dNdX_Ref_GP;
  Matrix F_Ref;
  Matrix x_EC;
  int GP_Elem_id;

  for(int i = 0 ; i<GP_Mesh.NumGP ;i++){
    GP_Elem_id = GP_Mesh.Element_id[i]; 
    printf("The GP %i is in the element : %i \n",
	   i+1,GP_Elem_id+1);
    printf("\t X Glob coordinate : %f \n",
	   GP_Mesh.Phi.x_GC.nV[i]);
    printf("\t X Ref coordiante : %f \n",
	   GP_Mesh.Phi.x_EC.nV[i]);
    printf("\t Velocity : %f \n",
	   GP_Mesh.Phi.vel.nV[i]);
    printf("\t Stress : %f \n",
	   GP_Mesh.Phi.Stress.nV[i]);
    printf("\t Mass : %f \n",
	   GP_Mesh.Phi.mass.nV[i]);
  }

  printf("Geometrical mass matrix : \n");
  for(int j = 0 ; j<ElementMesh.NumNodesMesh ; j++){
    printf("[");
    for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){
      printf(" %f ",M.nM[i][j]);
    }
    printf("]\n");
  }

  printf("Inverse of the lumped version of the geometrical mass matrix : \n");
  printf("[");
  for(int i = 0 ; i<M_l.N_rows ; i++){
    printf(" %f ",1/M_l.nV[i]);
  }
  printf("]\n");
  
    
  exit(EXIT_SUCCESS);
  
}
