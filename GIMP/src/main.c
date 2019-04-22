#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "ToolsLib/ElementTools.h"
#include "ToolsLib/GaussPointsTools.h"
#include "ToolsLib/Solvers.h"


int main(void){

  /* Physical parameters */
  double PoissonRatio = 0.1;
  double YoungModulus = 0.1;

  /* Define Gauss point mesh */
  GaussPoint * GP_e = AllocateGaussPoints(4);
  /* Matrix with the intial coordinates */
  Matrix GP_gp0 = MatAlloc(2,4);
  /* Matrix with the coordiantes of the gauss point */
  Matrix RefCoords = MatAlloc(1,2);

  GP_gp0.nM[0][0] = 0.5 - pow(3,0.5);
  GP_gp0.nM[1][0] = 0.5 - pow(3,0.5);
  GP_gp0.nM[0][1] = 0.5 + pow(3,0.5);
  GP_gp0.nM[1][1] = 0.5 - pow(3,0.5);
  GP_gp0.nM[0][2] = 0.5 + pow(3,0.5);
  GP_gp0.nM[1][2] = 0.5 + pow(3,0.5);
  GP_gp0.nM[0][3] = 0.5 - pow(3,0.5);
  GP_gp0.nM[1][3] = 0.5 + pow(3,0.5);

  for(int i = 0 ; i<4 ; i++){
    
    RefCoords.nV[0] = GP_gp0.nM[0][i];
    RefCoords.nV[1] = GP_gp0.nM[1][i];
    
    /* Initialize the Gauss points mesh */
    Initialize_GP(&GP_e[i],i,RefCoords,PoissonRatio,YoungModulus);    
  }

  /* Define one single element */
  Element Elem;

  /* Type of the element */
  char TypeElement[4] = "Q4";
  
  /* Coordinates of the nodes of the element */
  Matrix GlobalNodsCoords_e = MatAlloc(4,2);
  GlobalNodsCoords_e.nM[0][0] = 0;
  GlobalNodsCoords_e.nM[0][1] = 0;    
  GlobalNodsCoords_e.nM[1][0] = 1;
  GlobalNodsCoords_e.nM[1][1] = 0;
  GlobalNodsCoords_e.nM[2][0] = 1;
  GlobalNodsCoords_e.nM[2][1] = 1;
  GlobalNodsCoords_e.nM[3][0] = 0;
  GlobalNodsCoords_e.nM[3][1] = 1;

  /* Conectivity mesh of the element */
  int * GlobalNodsId_e = (int *)Allocate_Array(4,sizeof(int));
  GlobalNodsId_e[0] = 0;
  GlobalNodsId_e[1] = 1;
  GlobalNodsId_e[2] = 2;
  GlobalNodsId_e[3] = 3;
      
  /* Initialize the element mesh */
  Initialize_Element(&Elem,GP_e,TypeElement,
		     GlobalNodsCoords_e,
		     GlobalNodsId_e);

  for(int i = 0 ; i<4 ; i++){

    Get_RefDeformation_Gradient(&GP_e[i],&Elem);

  }
 
  Matrix K;

  K = Get_Stiffness_Matrix(&Elem);

  /* printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n", */
  /* 	 K.nM[0][0],K.nM[0][1],K.nM[0][2],K.nM[0][3],K.nM[0][4],K.nM[0][5],K.nM[0][6],K.nM[0][7]); */
  /* printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n", */
  /* 	 K.nM[1][0],K.nM[1][1],K.nM[1][2],K.nM[1][3],K.nM[1][4],K.nM[1][5],K.nM[1][6],K.nM[1][7]); */
  /* printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n", */
  /* 	 K.nM[2][0],K.nM[2][1],K.nM[2][2],K.nM[2][3],K.nM[2][4],K.nM[2][5],K.nM[2][6],K.nM[2][7]); */
  /* printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n", */
  /* 	 K.nM[3][0],K.nM[3][1],K.nM[3][2],K.nM[3][3],K.nM[3][4],K.nM[3][5],K.nM[3][6],K.nM[3][7]); */
  /* printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n", */
  /* 	 K.nM[4][0],K.nM[4][1],K.nM[4][2],K.nM[4][3],K.nM[4][4],K.nM[4][5],K.nM[4][6],K.nM[4][7]); */
  /* printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n", */
  /* 	 K.nM[5][0],K.nM[5][1],K.nM[5][2],K.nM[5][3],K.nM[5][4],K.nM[5][5],K.nM[5][6],K.nM[5][7]); */
  /* printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n", */
  /* 	 K.nM[6][0],K.nM[6][1],K.nM[6][2],K.nM[6][3],K.nM[6][4],K.nM[6][5],K.nM[6][6],K.nM[6][7]); */
  /* printf("[%f ; %f ; %f ; %f ; %f ; %f ; %f ; %f ] \n", */
  /* 	 K.nM[7][0],K.nM[7][1],K.nM[7][2],K.nM[7][3],K.nM[7][4],K.nM[7][5],K.nM[7][6],K.nM[7][7]); */

  printf("Resolver sistema de ecuaciones \n");

  Matrix A = MatAlloc(2,2);
  A.nM[0][0] = 4;
  A.nM[0][1] = 1;
  A.nM[1][0] = 1;
  A.nM[1][1] = 3;
  Matrix b = MatAlloc(2,1);
  b.nV[0] = 1;
  b.nV[1] = 2;
  Matrix x0 = MatAllocZ(2,1);

  Matrix x = Conjugate_Gradient_Method(A,b,x0);

  printf("Sol Ax = b : %f ; %f \n",x.nV[0],x.nV[1]);
  
  return 0;
}
