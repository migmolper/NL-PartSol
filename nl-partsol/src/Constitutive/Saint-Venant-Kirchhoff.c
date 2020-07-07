#include "nl-partsol.h"

/**************************************************************/

double energy_Saint_Venant_Kirchhoff(Tensor E, Material MatProp_p)
{
  
  /* Material parameters */
  double trE = get_I1_Of(E);
  double EE = get_innerProduct_Of(E,E);
  double ElasticModulus = MatProp_p.E;
  double mu = MatProp_p.mu;
  double lambda = mu*ElasticModulus/((1-mu*2)*(1+mu));
  double G = ElasticModulus/(2*(1+mu));

  return 0.5*lambda*DSQR(trE) + G*EE;
}

/**************************************************************/

Tensor grad_energy_Saint_Venant_Kirchhoff(Tensor grad_e, Tensor E, Material MatProp_p)
{
  /* Number of dimensions */
  int Ndim = NumberDimensions;
  
  /* Material parameters */

  double ElasticModulus = MatProp_p.E;
  double mu = MatProp_p.mu;
  double lambda = mu*ElasticModulus/((1-mu*2)*(1+mu));
  double G = ElasticModulus/(2*(1+mu));
  double trE = get_I1_Of(E);

  Tensor I = get_I();

  for(int i = 0 ; i < Ndim ; i++)
    {
      for(int j = 0 ; j < Ndim ; j++)
	{
	  grad_e.N[i][j] = 0.5*lambda*trE*I.N[i][j] + G*E.N[i][j];
	}
    }
  
  /* Free identity */
  free_Tensor(I);

  return grad_e;
}

/**************************************************************/

