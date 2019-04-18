#include <stdio.h>
#include <stdlib.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "ToolsLib/ElementTools.h"


int main(void){

  /* Physical parameters */
  double PoissonRatio = 0.1;
  double YoungModulus = 0.1;

  /* Define one single Gauss point */
  GaussPoint GP_e;
  /* Define one single element */
  Element Elem;

  /* Matrix with the coordiantes of the gauss point */
  Matrix RefCoords = MatAlloc(1,2);
  RefCoords.nV[0] = 0.5;
  RefCoords.nV[1] = 0.5;

  /* Type of the element */
  char TypeElement[4] = "Q4";
  
  /* Initialize the Gauss points mesh */
  Initialize_GP(&GP_e,RefCoords,PoissonRatio,YoungModulus);    

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
  Initialize_Element(&Elem,&GP_e,TypeElement,
		     GlobalNodsCoords_e,
		     GlobalNodsId_e);

  Get_RefDeformation_Gradient(&GP_e,&Elem);
 
  Matrix K;

  K = Get_Stiffness_Matrix(&Elem);

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
   
  return 0;
}
