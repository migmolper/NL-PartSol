#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../Matlib/Matlib.h"

/*******************************************************/

void UpdateGaussPointStress(GaussPoint MPM_Mesh){

  /* 0º Variable declaration  */
  Matrix StrainTensor_GP;
  Matrix StressTensor_GP;

  /* 1º Switch the dimensions of the aulixiar strain tensor */
  switch(NumberDimensions){
  case 1:
    StrainTensor_GP = MatAssign(1,1,NAN,NULL,NULL);
    break;
  case 2:
    StrainTensor_GP = MatAssign(3,1,NAN,NULL,NULL);
    break;
  case 3:
    StrainTensor_GP = MatAssign(6,1,NAN,NULL,NULL);
    break;
  default :
    printf("%s : %s \n","Error in UpdateGaussPointStress()",
	   "Wrong number of dimensions !!! ");
    exit(0);
  }

  /* 2º Iterate over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    /* 3º Store in an auxiliar variable the strain tensor in the GP */
    StrainTensor_GP.nV = MPM_Mesh.Phi.Strain.nM[i];
    /* 4º Get the new stress tensor (2D Linear elastic) */
    StressTensor_GP =
      MPM_Mesh.D.LE2D(StrainTensor_GP,PoissonModulus,ElasticModulus);
    /* 5º Update the stress tensor with the new-one */
    for(int j = 0 ; j<StrainTensor_GP.N_rows ; j++){
      MPM_Mesh.Phi.Stress.nM[i][j] = StressTensor_GP.nV[j];
    }
    /* 6º Free memory */
    FreeMat(StressTensor_GP);
  }
  
}

/*******************************************************/
