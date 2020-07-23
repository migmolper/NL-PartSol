#include "nl-partsol.h"

double internal_energy__Particles__(Tensor Strain, Tensor Stress){

  /* Internal energy for the Gauss-Point */
  double W = 0; 
  /*Check in the input its is ok */
  if ((Strain.Order == 2) && (Stress.Order == 2)){    
    /* Calcule the internal work */
    W = 0.5*inner_product__TensorLib__(Strain, Stress);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in internal_energy__Particles__()",
	    "The input should be 2nd tensor and a 2nd tensor");
    exit(EXIT_FAILURE);
  }
  return W;
}

