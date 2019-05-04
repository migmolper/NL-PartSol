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
  ReadDatFile(argv[3]);

  /* Read mesh data */
  ReadGidMesh(MeshFileName);
  
  /* Initialize the element mesh */
  char * TypeElement = "L2";
  Element ElementMesh = Initialize_Element(0,TypeElement);
  
  /* Read the initial conditions fields as a CSV */
  Matrix InputFields = Read_CSV(InitCondFileName, 5);
    
  /* Read the data file */  /* Physical parameters */
  double YoungModulus = 1;
  
  /* Get material Constitutive matrix */
  Matrix D = LinearElastic1D(YoungModulus);


  for(int i = 0; i<MeshProp.Nnodes ; i++){
    printf("Nodal coordinates : [%f, %f, %f] \n",
	   MeshProp.Coordinates[i][0],
	   MeshProp.Coordinates[i][1],
	   MeshProp.Coordinates[i][2]);
  }

  printf("Connectivity mesh :\n");
  for(int i =0; i<MeshProp.Nelem; i++){
    printf("Element  %i : ",i+1);
    for(int j = 0 ; j<ElemProp.Nnodes; j++){
      printf(" %i ",MeshProp.Connectivity[i][j]);
    }
    printf("\n");
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
  
  /* Define Gauss point mesh */
  GaussPoint GP_Mesh  = Initialize_GP_Mesh(InputFields,D);

  /* Localize all the Gauss-Points */
  UpdateElementLocationGP(GP_Mesh);
    
  exit(EXIT_SUCCESS);
  
}
