#include "nl-partsol.h"

double compute_InternalEnergy(Tensor Strain, Tensor Stress){

  /* Internal energy for the Gauss-Point */
  double W = 0; 
  /*Check in the input its is ok */
  if ((Strain.Order == 2) && (Stress.Order == 2)){    
    /* Calcule the internal work */
    W = 0.5*get_dotProduct_Of(Strain, Stress);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in compute_InternalEnergy()",
	    "The input should be 2nd tensor and a 2nd tensor");
    exit(EXIT_FAILURE);
  }
  return W;
}

