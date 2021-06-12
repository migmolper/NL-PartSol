#include "nl-partsol.h"

/**************************************************************/

double energy_Saint_Venant_Kirchhoff(Tensor C, Material MatProp_p)
{
  /* Material parameters */
  double ElasticModulus = MatProp_p.E;
  double nu  = MatProp_p.nu ;
  double lambda = nu *ElasticModulus/((1-nu *2)*(1+nu ));
  double G = ElasticModulus/(2*(1+nu ));

  Tensor E = strain_Green_Lagrange__Particles__(C);
  double I1_E = I1__TensorLib__(E);

  double W = 0.5*lambda*I1_E*I1_E + G*inner_product__TensorLib__(E, E);

  free__TensorLib__(E);

  return W;
}

/**************************************************************/

State_Parameters compute_1PK_Stress_Tensor_Saint_Venant_Kirchhoff(
  State_Parameters Intput_SP,
  Material MatProp_p)
{
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /* 
    Output state parameter
  */
  State_Parameters Output_SP;

  /* Get information from the state parameter */
  Tensor P = memory_to_tensor__TensorLib__(Intput_SP.Stress,2);
  const Tensor F = memory_to_tensor__TensorLib__(Intput_SP.F_n1_p,2);
  
  /* Material parameters */
  const double ElasticModulus = MatProp_p.E;
  const double nu = MatProp_p.nu;
  const double lambda = nu*ElasticModulus/((1-nu*2)*(1+nu));
  const double G = ElasticModulus/(2*(1+nu));
  
  /*
    Auxiliar tensors and variables
  */
  Tensor Ft = transpose__TensorLib__(F);
  Tensor Ft_x_F = matrix_product__TensorLib__(Ft,F);
  Tensor F_x_Ft_x_F = matrix_product__TensorLib__(F,Ft_x_F);
  double tr__Ft_x_F = I1__TensorLib__(Ft_x_F);

  for(int i = 0 ; i < Ndim ; i++)
  {
    for(int j = 0 ; j < Ndim ; j++)
    {
      P.N[i][j] = lambda*0.5*(tr__Ft_x_F - Ndim)*F.N[i][j] + G*(F_x_Ft_x_F.N[i][j] - F.N[i][j]);
    }
  }
  
  /*
    Free tensors 
  */
  free__TensorLib__(Ft);
  free__TensorLib__(Ft_x_F);
  free__TensorLib__(F_x_Ft_x_F);

  return Output_SP;
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
  double nu = MatProp_p.nu;
  double lambda = nu*ElasticModulus/((1-nu*2)*(1+nu));
  double G = ElasticModulus/(2*(1+nu));
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
  double nu = MatProp.nu;
  double lambda = nu*ElasticModulus/((1-nu*2)*(1+nu));
  double G = ElasticModulus/(2*(1+nu));

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
