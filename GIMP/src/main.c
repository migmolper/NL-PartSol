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
#include "Num_Schemes/Num_Schemes.h"

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
  ElementMesh.NumberDOF = 2; /* Change this */
    
  /* Read the initial conditions fields as a CSV */
  Matrix InputFields = Read_CSV(InitCondFileName, 5);

  /* Read the boundary conditions */
  ReadBCC(BounCondFileName);
  
  /* Get material Constitutive matrix */
  Matrix D = LinearElastic1D(ElasticModulus);

  /* Define Gauss point mesh */
  GaussPoint GP_Mesh  = Initialize_GP_Mesh(InputFields,D,
					   ElementMesh);

  
  Matrix A = MatAllocZ(2,2);
  A.nM[0][0] = - D.n;
  A.nM[1][1] = - (double)1/Density;
  printf("************************************************* \n");
  printf("Coupling matrix A : \n");
  printf("\t [ %f %f ] \n",A.nM[0][0],A.nM[0][1]);
  printf("\t [ %f %f ] \n",A.nM[1][0],A.nM[1][1]);

  /* Run simulation  */
  printf("************************************************* \n");
  printf("Run simulation :  !!! \n");
  
  printf("************************************************* \n");
  printf("Time step : %i \n",0);
  Two_Steps_TG_Mc(ElementMesh,GP_Mesh,A,0);
  
  
    
  exit(EXIT_SUCCESS);
  
}
