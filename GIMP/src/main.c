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
  Matrix InputFields = Read_CSV(InitCondFileName, 10);

  /* Read the boundary conditions */
  ReadBCC(BounCondFileName);
  
  /* Get material Constitutive matrix */
  Matrix D = LinearElastic1D(ElasticModulus);

  /* Define Gauss point mesh */
  GaussPoint GP_Mesh  = Initialize_GP_Mesh(InputFields,D,
					   ElementMesh);

  
  Matrix A = MatAllocZ(2,2);
  A.nM[0][1] = - D.n;
  A.nM[1][0] = - (double)1/Density;
  printf("************************************************* \n");
  printf("Coupling matrix A : \n");
  printf("\t [ %f %f ] \n",A.nM[0][0],A.nM[0][1]);
  printf("\t [ %f %f ] \n",A.nM[1][0],A.nM[1][1]);

  /* Run simulation  */
  printf("************************************************* \n");
  printf("Run simulation :  !!! \n");

  for(int t_i = 0 ; t_i<10; t_i++){
    printf("************************************************* \n");

    GnuplotOutput1D(GP_Mesh.Phi.x_GC,
		    GP_Mesh.Phi.vel,
		    0.0, 5.0,
		    t_i,GP_Mesh.NumGP);

    printf("Time step : %i \n",t_i);
    Two_Steps_TG_Mc(ElementMesh,GP_Mesh,A,t_i);
  }
  
  exit(EXIT_SUCCESS);
  
}


