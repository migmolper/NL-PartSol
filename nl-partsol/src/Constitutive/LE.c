#include "grams.h"


Tensor LinearElastic(Tensor Strain, Tensor Stress, Material Mat){

  int Ndim = NumberDimensions;
  
  /*Check in the input its is ok */
  if ((Strain.Order == 2) && (Stress.Order == 2)){
    /* Define material and other properties */
    double mu = Mat.mu; 
    double E = Mat.E;
    double lambda = mu*E/((1-mu*2)*(1+mu));
    double G = E/(2*(1+mu));
    double traceStrain = get_I1_Of(Strain);
    Tensor I = get_I();

    for(int i = 0 ; i<Ndim ; i++){
      /* Add first component */
      for(int j = 0 ; j<Ndim ; j++){
	Stress.N[i][j] = 2*G*Strain.N[i][j];
      }
      /* Add diagonal component */
      Stress.N[i][i] += lambda*traceStrain*I.N[i][i];
    }

    free_Tensor(I);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in LinearElastic()",
	    "The input should be 2nd tensor and a 2nd tensor");
    exit(EXIT_FAILURE);
  }
  return Stress;  
}

