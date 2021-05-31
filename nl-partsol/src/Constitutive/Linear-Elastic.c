#include "nl-partsol.h"

/*************************************************************/

Tensor LinearElastic(Tensor Stress, Tensor Strain, Material Mat)
{

  int Ndim = NumberDimensions;
  double nu = Mat.nu; 
  double E = Mat.E;
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double traceStrain = I1__TensorLib__(Strain);

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      Stress.N[i][j] = 2*G*Strain.N[i][j] + (K - 2*G/3.0)*(i==j)*traceStrain;
    }
  }

  return Stress;  
}

/**************************************************************/

