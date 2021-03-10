#include "nl-partsol.h"

/**************************************************************/

Tensor compute_log_stress_tensor__Hencky__(Tensor T, Tensor E, Material MatProp_p)
{
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /*
    Compute auxiliar tensors
  */
  Tensor I = Identity__TensorLib__();
  
  /* Material parameters */
  double ElasticModulus = MatProp_p.E;
  double nu = MatProp_p.nu;
  double lambda = nu*ElasticModulus/((1-nu*2)*(1+nu));
  double G = ElasticModulus/(2*(1+nu));
  double trE = I1__TensorLib__(E);

  for(int i = 0 ; i < Ndim ; i++)
  {
  	for(int j = 0 ; j < Ndim ; j++)
	{
	  T.N[i][j] = lambda*trE*I.N[i][j] + 2*G*E.N[i][j];
	}
  }
  
  /*
    Free
  */
  free__TensorLib__(I);

  return T;
}

/**************************************************************/