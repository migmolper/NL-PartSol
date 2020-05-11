#include "nl-partsol.h"

Tensor SolidRigid(Tensor Stress){

  int Ndim = NumberDimensions;
  
  /*Check in the input its is ok */
  if (Stress.Order == 2){
    for(int i = 0 ; i<Ndim ; i++){
      for(int j = 0 ; j<Ndim ; j++){
	Stress.N[i][j] = 0;
      }
    }
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in SolidRigid()",
	    "The input should be 2nd tensor and a 2nd tensor");
    exit(EXIT_FAILURE);
  }
  return Stress;  
}

