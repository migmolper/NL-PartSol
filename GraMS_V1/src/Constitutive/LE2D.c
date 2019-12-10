#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../Matlib/Matlib.h"

Matrix LinearElastic2D(Matrix Strain,
		       double PoissonRatio_Mat,
		       double ElasticModulus_Mat){

  Matrix Stress = MatAlloc(1,3);
  
  double LameFirstParam = /* Lambda */
    PoissonRatio_Mat*ElasticModulus_Mat/((1-PoissonRatio_Mat*2)*(1+PoissonRatio_Mat));
  double LameSecondParam = /* G */
    ElasticModulus_Mat/(2*(1+PoissonRatio_Mat));

  Stress.nV[0] = (LameFirstParam + 2*LameSecondParam)*Strain.nV[0] +
    LameFirstParam*Strain.nV[1];
  Stress.nV[1] = (LameFirstParam + 2*LameSecondParam)*Strain.nV[1] +
    LameFirstParam*Strain.nV[0];
  Stress.nV[2] = 2*LameSecondParam*Strain.nV[2];
    
  return Stress; 
}
