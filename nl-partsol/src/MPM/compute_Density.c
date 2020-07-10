#include "nl-partsol.h"

/*******************************************************/

double update_Density(double rho_n, double TimeStep,
		      Tensor Rate_Strain)
{
  /* Define output variable */
  double rho_n1;  
  /* Check if the input is ok*/
  if ((Rate_Strain.Order == 2)){
    /* Update the density */
    rho_n1 = rho_n/(1 + TimeStep*I1__TensorLib__(Rate_Strain));
  }
  else{
    fprintf(stderr,"%s : %s %s !!! \n",
	    "Error in update_Density()",
	    "The input should be",
	    "a 2nd order tensor and a scalar");
    exit(EXIT_FAILURE);    
  }
  return rho_n1;  
}

/*******************************************************/ 
