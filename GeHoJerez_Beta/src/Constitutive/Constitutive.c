#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"

Matrix LinearElastic1D(double ElasticModulus_Mat){

  Matrix D  = MatAlloc(1,1);

  D.n = ElasticModulus_Mat;

  return D;

}

Matrix LinearElastic2D(double PoissonRatio_Mat,double ElasticModulus_Mat){

  Matrix D = MatAlloc(3,3);
  
  double LameFirstParam = /* Lambda */
    PoissonRatio_Mat*ElasticModulus_Mat/((1-PoissonRatio_Mat*2)*(1+PoissonRatio_Mat));
  double LameSecondParam = /* G */
    ElasticModulus_Mat/(2*(1+PoissonRatio_Mat));

  D.nM[0][0] = LameFirstParam + 2*LameSecondParam;
  D.nM[0][1] = LameFirstParam;
  D.nM[0][2] = 0;
  D.nM[1][0] = LameFirstParam;
  D.nM[1][1] = LameFirstParam + 2*LameSecondParam;
  D.nM[1][2] = 0;
  D.nM[2][0] = 0;
  D.nM[2][1] = 0;
  D.nM[2][2] = 2*LameSecondParam;

  
  /* Initialize name of material */
  strcpy(D.Info,"Elastic");
  
  return D;
  
}
