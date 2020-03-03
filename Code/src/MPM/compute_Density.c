#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

/*******************************************************/

double update_Density(double rho_n, double TimeStep,
		      Tensor Rate_Strain)
{
  /* Define output variable */
  double rho_n1;  
  /* Check if the input is ok*/
  if ((Rate_Strain.Order == 2)){
    /* Update the density */
    rho_n1 = rho_n/(1 + TimeStep*get_I1_Of(Rate_Strain));
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

/* double UpdateGaussPointDensity(double rho_n_GP, */
/* 			       double Delta_TraceStrain_GP) */
/* /\*! */
/*  * \brief Brief description of UpdateGaussPointDensity. */
/*  *        Update the density field of the Gauss Point .  */
/*  * */
/*  *  The parameters for this functions are  : */
/*  *  @param rho_n_GP : Density of the previous step. */
/*  *  @param Delta_TraceStrain_GP : Increment of the trace of the strain tensor. */
/*  * */
/*  *\/ */
/* { */
/*   double rho_n1_GP; /\* Density for the next step *\/ */

/*   /\* Update the density *\/ */
/*   rho_n1_GP = (double)rho_n_GP/(1 + Delta_TraceStrain_GP); */

/*   return rho_n1_GP;   */
/* } */

/*******************************************************/
