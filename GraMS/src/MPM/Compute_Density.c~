#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/TypeDefinitions.h"
#include "../GRAMS/GlobalVariables.h"

/*******************************************************/


double UpdateGaussPointDensity(double rho_n,
			       double Incr_TraceStrain){

  /* 1º Density for the next step */
  double rho_n1;

  /* 2º Update the density */
  rho_n1 = (double)rho_n/(1 + Incr_TraceStrain);

  /* 3º Return density updated */
  return rho_n1;  
}

/*******************************************************/
