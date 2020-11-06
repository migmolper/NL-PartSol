#include "nl-partsol.h"

/**************************************************************/

double energy_Saint_Venant_Kirchhoff(Tensor C, Material MatProp_p)
{
  /* Material parameters */
  double ElasticModulus = MatProp_p.E;
  double mu = MatProp_p.mu;
  double lambda = mu*ElasticModulus/((1-mu*2)*(1+mu));
  double G = ElasticModulus/(2*(1+mu));

  Tensor E = strain_Green_Lagrange__Particles__(C);
  double I1_E = I1__TensorLib__(E);

  double W = 0.5*lambda*I1_E*I1_E + G*inner_product__TensorLib__(E, E);

  free__TensorLib__(E);

  return W;
}

/**************************************************************/

Tensor grad_energy_Saint_Venant_Kirchhoff(Tensor grad_e, Tensor C, Material MatProp_p)
{
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /*
    Compute auxiliar tensors
  */
  Tensor E = strain_Green_Lagrange__Particles__(C);
  Tensor I = Identity__TensorLib__();
  
  /* Material parameters */
  double ElasticModulus = MatProp_p.E;
  double mu = MatProp_p.mu;
  double lambda = mu*ElasticModulus/((1-mu*2)*(1+mu));
  double G = ElasticModulus/(2*(1+mu));
  double trE = I1__TensorLib__(E);

  for(int i = 0 ; i < Ndim ; i++)
    {
      for(int j = 0 ; j < Ndim ; j++)
	{
	  grad_e.N[i][j] = lambda*trE*I.N[i][j] + G*E.N[i][j];
	}
    }
  
  /*
    Free
  */
  free__TensorLib__(I);
  free__TensorLib__(E);

  return grad_e;
}

/**************************************************************/

Tensor compute_stiffness_density_Saint_Venant_Kirchhoff(Tensor v, Tensor w, Material MatProp)
{

  /*
    Number of dimensions
  */
  int Ndim = NumberDimensions;
    
  /*
    Material parameters 
  */
  double ElasticModulus = MatProp.E;
  double mu = MatProp.mu;
  double lambda = mu*ElasticModulus/((1-mu*2)*(1+mu));
  double G = ElasticModulus/(2*(1+mu));

  /*
    Auxiliar variables
   */  
  double v_dot_w;
  Tensor v_o_w;
  Tensor C_mat = alloc__TensorLib__(2);
  
  v_o_w  = dyadic_Product__TensorLib__(v,w);
  v_dot_w = inner_product__TensorLib__(v,w);

  for(int A = 0 ; A<Ndim ; A++)
    {
      for(int B = 0 ; B<Ndim ; B++)
	{
	  C_mat.N[A][B] +=
	    lambda*v_o_w.N[A][B] + G*v_o_w.N[B][A] + (A == B ? 1 : 0)*G*v_dot_w;
	}
    }

  /*
    Free memory
   */
  free__TensorLib__(v_o_w);


  return C_mat;
}   

/**************************************************************/
