#include "nl-partsol.h"

/**************************************************************/

Tensor grad_energy_Neo_Hookean(Tensor grad_e, Tensor C, double J, Material MatProp_p)
{
  /* Number of dimensions */
  int Ndim = NumberDimensions;
  
  /* Material parameters */
  double ElasticModulus = MatProp_p.E;
  double mu = MatProp_p.mu;
  double G = ElasticModulus/(2*(1+mu));
  double lambda = mu*ElasticModulus/((1-mu*2)*(1+mu));
  double J2 = J*J;
  
  /*
    Auxiliar tensors
  */
  Tensor I    = Identity__TensorLib__();
  Tensor C_m1 = Inverse__TensorLib__(C);

  for(int i = 0 ; i < Ndim ; i++)
    {
      for(int j = 0 ; j < Ndim ; j++)
	{
	  grad_e.N[i][j] =
	    lambda*0.5*(J2 - 1)*C_m1.N[i][j] +
	    G*(I.N[i][j] - C_m1.N[i][j]);
	}
    }
  
  /*
    Free tensors 
  */
  free__TensorLib__(I);
  free__TensorLib__(C_m1);

  return grad_e;
}

/**************************************************************/

Tensor compute_stiffness_density_Neo_Hookean(Tensor v, Tensor w, Tensor C,
					     double J, Material MatProp)
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
  double G = ElasticModulus/(2*(1+mu));
  double lambda = mu*ElasticModulus/((1-mu*2)*(1+mu));
  double J2 = J*J;
  
  /*
    Stifness density tensor
  */
  Tensor C_mat = alloc__TensorLib__(2);

  /*
    Auxiliar variables
  */  
  Tensor Cm1;
  Tensor Cm1T;
  Tensor Cm1_dot_v;
  Tensor Cm1_dot_w;
  Tensor Cm1T_dot_v;
  Tensor Cm1_dot_v_o_Cm1_dot_w;
  Tensor Cm1_dot_w_o_Cm1T_dot_v;
  double Cm1_dot_w_dot_v;
  
  Cm1                    = Inverse__TensorLib__(C);
  Cm1T                   = transpose__TensorLib__(Cm1);
  Cm1_dot_w              = vector_linear_mapping__TensorLib__(Cm1,w);
  Cm1_dot_v              = vector_linear_mapping__TensorLib__(Cm1,v);
  Cm1T_dot_v             = vector_linear_mapping__TensorLib__(Cm1T,v);
  Cm1_dot_w_dot_v        = inner_product__TensorLib__(Cm1_dot_w,v);
  Cm1_dot_v_o_Cm1_dot_w  = dyadic_Product__TensorLib__(Cm1_dot_v,Cm1_dot_w);    
  Cm1_dot_w_o_Cm1T_dot_v = dyadic_Product__TensorLib__(Cm1_dot_w,Cm1T_dot_v);  
  
  for(int A = 0 ; A<Ndim ; A++)
    {
      for(int B = 0 ; B<Ndim ; B++)
	{
	  C_mat.N[A][B] +=
	    lambda*J2*Cm1_dot_v_o_Cm1_dot_w.N[A][B] -
	    (0.5*(J2 - 1) - G)*Cm1_dot_w_dot_v*Cm1.N[A][B] -
	    (0.5*lambda*(J2 - 1) - G)*Cm1_dot_w_o_Cm1T_dot_v.N[A][B];
	}
    }

  /*
    Free memory
   */
  free__TensorLib__(Cm1);
  free__TensorLib__(Cm1T);
  free__TensorLib__(Cm1_dot_w);
  free__TensorLib__(Cm1_dot_v);
  free__TensorLib__(Cm1T_dot_v);
  free__TensorLib__(Cm1_dot_v_o_Cm1_dot_w);
  free__TensorLib__(Cm1_dot_w_o_Cm1T_dot_v);
  

  return C_mat;
}   

/**************************************************************/
