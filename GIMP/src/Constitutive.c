#include <stdio.h>
#include <stdlib.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"


Matrix LinearElastic(double PoissonRatio,double YoungModulus){

  Matrix D = MatAlloc(3,3);
  
  double LameFirstParam = /* Lambda */
    PoissonRatio*YoungModulus/((1-PoissonRatio*2)*(1+PoissonRatio));
  double LameSecondParam = /* G */
    YoungModulus/(2*(1+PoissonRatio));

  D.nM[0][0] = LameFirstParam + 2*LameSecondParam;
  D.nM[0][1] = LameFirstParam;
  D.nM[0][2] = 0;
  D.nM[1][0] = LameFirstParam;
  D.nM[1][1] = LameFirstParam + 2*LameSecondParam;
  D.nM[1][2] = 0;
  D.nM[2][0] = 0;
  D.nM[2][1] = 0;
  D.nM[2][2] = 2*LameSecondParam;
  
  return D;
}
