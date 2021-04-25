#include "nl-partsol.h"

/**************************************************************/

Tensor compute_1PK_Stress_Tensor_Newtonian_Fluid(
  Tensor P, 
  Tensor F,
  Tensor d,
  double J, 
  Material MatProp_p)
{
  /* Number of dimensions */
  int Ndim = NumberDimensions;
  
  /* Material parameters */
  double ElasticModulus = MatProp_p.E;
  double nu = MatProp_p.nu;
  double G = ElasticModulus/(2*(1+nu));
  double lambda = nu*ElasticModulus/((1-nu*2)*(1+nu));
  double J2 = J*J;
  
  /*
    Auxiliar tensors
  */
  Tensor FT = transpose__TensorLib__(F);
  Tensor FmT = Inverse__TensorLib__(FT);

  for(int i = 0 ; i < Ndim ; i++)
  {
    for(int j = 0 ; j < Ndim ; j++)
    {
      P.N[i][j] = lambda*0.5*(J2 - 1)*FmT.N[i][j] + G*(F.N[i][j] - FmT.N[i][j]);
    }
  }
  
  /*
    Free tensors 
  */
  free__TensorLib__(FT);
  free__TensorLib__(FmT);

  return P;
}

/**************************************************************/