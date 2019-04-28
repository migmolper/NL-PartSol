#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "InOutFun/InOutFun.h"
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

  /* Initialize parser to read files */
  ParserDictionary Dict = InitParserDictionary();

  /* Declaration of the structure that contain the data of the mesh */
  ElementProperties ElemProp;
  MeshProperties MeshProp;

  /* Read mesh data */
  ReadGidMesh(argv[1],&ElemProp,&MeshProp,Dict);

  printf("%i ; %i \n",MeshProp.Nnodes,ElemProp.Nnodes);
  puts("paso");
  exit(0);

  for(int i = 0 ; i<MeshProp.Nnodes ; i++){
    printf("Nodes of the element (%i) = \n",i);
    for(int j = 0 ; j<ElemProp.Nnodes ; j++){
      printf(" %i ",MeshProp.Connectivity[i][j]);
    }
    printf("\n");
  }

  


  /* Read the mesh */

  /* Physical parameters */
  double PoissonRatio = 0.1;
  double YoungModulus = 0.1;

  /* Define Gauss point mesh */
  GaussPoint * GP_e = AllocateGaussPoints(4);
  /* Matrix with the intial coordinates */
  Matrix GP_gp0 = MatAlloc(4,2);
  /* Matrix with the coordiantes of the gauss point */
  Matrix RefCoords = MatAlloc(1,2);

  GP_gp0.nM[0][0] = - 1/pow(3,0.5);
  GP_gp0.nM[0][1] = - 1/pow(3,0.5);
  
  GP_gp0.nM[1][0] = + 1/pow(3,0.5);
  GP_gp0.nM[1][1] = - 1/pow(3,0.5);

  GP_gp0.nM[2][0] = + 1/pow(3,0.5);
  GP_gp0.nM[2][1] = + 1/pow(3,0.5);

  GP_gp0.nM[3][0] = - 1/pow(3,0.5);
  GP_gp0.nM[3][1] = + 1/pow(3,0.5);

  for(int i = 0 ; i<4 ; i++){

    /* Pass by reference the coordinates of the GP */
    RefCoords.nV = GP_gp0.nM[i];
    
    /* Initialize the Gauss points mesh */
    GP_e[i] = Initialize_GP(i,RefCoords,PoissonRatio,YoungModulus);
  }

  /* Define one single element */
  Element ElementMesh;

  /* Type of the element */
  char TypeElement[4] = "Q4";
  
  /* Coordinates of the nodes of the element */
  Matrix GlobalNodsCoords_e = MatAlloc(4,2);
  GlobalNodsCoords_e.nM[0][0] = -1;
  GlobalNodsCoords_e.nM[0][1] = -1;
  
  GlobalNodsCoords_e.nM[1][0] = 1;
  GlobalNodsCoords_e.nM[1][1] = -1;

  GlobalNodsCoords_e.nM[2][0] = 1;
  GlobalNodsCoords_e.nM[2][1] = 1;

  GlobalNodsCoords_e.nM[3][0] = -1;
  GlobalNodsCoords_e.nM[3][1] = 1;

  /* Conectivity mesh of the element */
  int * GlobalNodsId_e = (int *)Allocate_Array(4,sizeof(int));
  GlobalNodsId_e[0] = 0;
  GlobalNodsId_e[1] = 1;
  GlobalNodsId_e[2] = 2;
  GlobalNodsId_e[3] = 3;

  
  /* Initialize the element mesh */
  ElementMesh = Initialize_Element(GP_e,TypeElement,
				   GlobalNodsCoords_e,
				   GlobalNodsId_e);
  
  for(int i = 0 ; i<4 ; i++){

    Get_RefDeformation_Gradient(&GP_e[i],&ElementMesh);

  }

  Matrix K;

  K = Get_Stiffness_Matrix(&ElementMesh);

  printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n",
  	 K.nM[0][0],K.nM[0][1],K.nM[0][2],K.nM[0][3],K.nM[0][4],K.nM[0][5],K.nM[0][6],K.nM[0][7]);
  printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n",
  	 K.nM[1][0],K.nM[1][1],K.nM[1][2],K.nM[1][3],K.nM[1][4],K.nM[1][5],K.nM[1][6],K.nM[1][7]);
  printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n",
  	 K.nM[2][0],K.nM[2][1],K.nM[2][2],K.nM[2][3],K.nM[2][4],K.nM[2][5],K.nM[2][6],K.nM[2][7]);
  printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n",
  	 K.nM[3][0],K.nM[3][1],K.nM[3][2],K.nM[3][3],K.nM[3][4],K.nM[3][5],K.nM[3][6],K.nM[3][7]);
  printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n",
  	 K.nM[4][0],K.nM[4][1],K.nM[4][2],K.nM[4][3],K.nM[4][4],K.nM[4][5],K.nM[4][6],K.nM[4][7]);
  printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n",
  	 K.nM[5][0],K.nM[5][1],K.nM[5][2],K.nM[5][3],K.nM[5][4],K.nM[5][5],K.nM[5][6],K.nM[5][7]);
  printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n",
  	 K.nM[6][0],K.nM[6][1],K.nM[6][2],K.nM[6][3],K.nM[6][4],K.nM[6][5],K.nM[6][6],K.nM[6][7]);
  printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n",
  	 K.nM[7][0],K.nM[7][1],K.nM[7][2],K.nM[7][3],K.nM[7][4],K.nM[7][5],K.nM[7][6],K.nM[7][7]);

  
  exit(EXIT_SUCCESS);
}
