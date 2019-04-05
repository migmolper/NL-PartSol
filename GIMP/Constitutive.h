#include <stdio.h>
#include <stdlib.h>
#include "Utils.h"

double ** LinearElastic(double PoissonRatio,double YoungModulus){

  double LameFirstParam = /* Lambda */
    PoissonRatio*YoungModulus/((1-PoissonRatio*2)*(1+PoissonRatio));
  double LameSecondParam = /* G */
    YoungModulus/(2*(1+PoissonRatio));

  double ** D = (double **)Allocate_Matrix(3,3,sizeof(double));

  D[0][0] = LameFirstParam + 2*LameSecondParam; 
  D[0][1] = LameFirstParam; 
  D[0][2] = 0; 
  D[1][0] = LameFirstParam; 
  D[1][1] = LameFirstParam + 2*LameSecondParam; 
  D[1][2] = 0; 
  D[2][0] = 0; 
  D[2][1] = 0;
  D[2][2] = 2*LameSecondParam;
  
  return D;
}
