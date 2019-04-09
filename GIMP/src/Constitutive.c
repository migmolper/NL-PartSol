#include <stdio.h>
#include <stdlib.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"


Matrix LinearElastic(double PoissonRatio,double YoungModulus){

  Matrix D;
  D.N_cols = 3;
  D.N_rows = 3;
  
  double LameFirstParam = /* Lambda */
    PoissonRatio*YoungModulus/((1-PoissonRatio*2)*(1+PoissonRatio));
  double LameSecondParam = /* G */
    YoungModulus/(2*(1+PoissonRatio));

  D.n = (double **)Allocate_Matrix(D.N_rows,D.N_cols,sizeof(double));

  D.n[0][0] = LameFirstParam + 2*LameSecondParam; 
  D.n[0][1] = LameFirstParam; 
  D.n[0][2] = 0; 
  D.n[1][0] = LameFirstParam; 
  D.n[1][1] = LameFirstParam + 2*LameSecondParam; 
  D.n[1][2] = 0; 
  D.n[2][0] = 0; 
  D.n[2][1] = 0;
  D.n[2][2] = 2*LameSecondParam;
  
  return D;
}
