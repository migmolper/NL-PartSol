#include "nl-partsol.h"

/**************************************************************/

double energy_Neo_Hookean_Wriggers(Tensor C, double J, Material MatProp_p)
{
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /* Material parameters */
  double ElasticModulus = MatProp_p.E;
  double nu  = MatProp_p.nu ;
  double G = ElasticModulus/(2*(1+nu ));
  double lambda = nu *ElasticModulus/((1-nu *2)*(1+nu ));
  double I1_C = I1__TensorLib__(C);
  double f_J = 0.25*lambda*(J*J - 1) - 0.5*lambda*log(J) - G*log(J);

  double W =  f_J + 0.5*G*(I1_C - Ndim);

  return W;
}

/**************************************************************/

Tensor grad_energy_Neo_Hookean_Wriggers(Tensor grad_e, Tensor C,
					double J, Material MatProp_p)
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
  Tensor I    = Identity__TensorLib__();
  Tensor C_m1 = Inverse__TensorLib__(C);

  for(int i = 0 ; i < Ndim ; i++)
    {
      for(int j = 0 ; j < Ndim ; j++)
	{
	  grad_e.N[i][j] = lambda*0.5*(J2 - 1)*C_m1.N[i][j] + G*(I.N[i][j] - C_m1.N[i][j]);
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

Tensor compute_stiffness_density_Neo_Hookean_Wriggers(Tensor v, Tensor w,
						      Tensor C, double J,
						      Material MatProp)
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
  double G = ElasticModulus/(2*(1+nu));
  double lambda = nu*ElasticModulus/((1-nu*2)*(1+nu));
  double J2 = J*J;
  double alpha = lambda*J2;
  double beta = 0.5*lambda*(J2 - 1) - G;
  
  /*
    Stifness density tensor
  */
  Tensor C_mat = alloc__TensorLib__(2);

  /*
    Auxiliar variables
  */  
  Tensor Cm1;
  Tensor Cm1_dot_v;
  Tensor Cm1_dot_w;
  Tensor Cm1_dot_v_o_Cm1_dot_w;
  double Cm1_dot_w_dot_v;
  
  Cm1                     = Inverse__TensorLib__(C);
  Cm1_dot_v               = vector_linear_mapping__TensorLib__(Cm1,v);
  Cm1_dot_w               = vector_linear_mapping__TensorLib__(Cm1,w);
  Cm1_dot_w_dot_v         = inner_product__TensorLib__(Cm1_dot_w,v);
  Cm1_dot_v_o_Cm1_dot_w   = dyadic_Product__TensorLib__(Cm1_dot_v,Cm1_dot_w);    

  
  for(int A = 0 ; A<Ndim ; A++)
    {
      for(int B = 0 ; B<Ndim ; B++)
      {
	     C_mat.N[A][B] += alpha*Cm1_dot_v_o_Cm1_dot_w.N[A][B] -
                  	    beta*Cm1_dot_w_dot_v*Cm1.N[A][B] -
                  	    beta*Cm1_dot_v_o_Cm1_dot_w.N[B][A];
      }
    }

  /*
    Free memory
   */
  free__TensorLib__(Cm1);
  free__TensorLib__(Cm1_dot_w);
  free__TensorLib__(Cm1_dot_v);
  free__TensorLib__(Cm1_dot_v_o_Cm1_dot_w);
  

  return C_mat;
}   

/**************************************************************/

Matrix compute_D_matrix_Neo_Hookean_Wriggers(Tensor F, double J, Material MatProp)
{

  /*
    Material parameters 
  */
  double ElasticModulus = MatProp.E;
  double nu = MatProp.nu;
  double G = ElasticModulus/(2*(1+nu));
  double lambda = nu*ElasticModulus/((1-nu*2)*(1+nu));
  double J2 = J*J;

  Tensor C = right_Cauchy_Green__Particles__(F);
  Tensor Cm1 = Inverse__TensorLib__(C);


  Matrix D = allocZ__MatrixLib__(3,3);
  
  D.nM[0][0] = lambda*J2*Cm1.N[0][0]*Cm1.N[0][0]-(lambda*(J2-1)-2*G)*Cm1.N[0][0]*Cm1.N[0][0];
  D.nM[0][1] = lambda*J2*Cm1.N[0][0]*Cm1.N[1][1]-(lambda*(J2-1)-2*G)*Cm1.N[0][1]*Cm1.N[1][0];
  D.nM[1][0] = D.nM[0][1];
  D.nM[0][2] = lambda*J2*Cm1.N[0][0]*Cm1.N[0][1]-(lambda*(J2-1)-2*G)*Cm1.N[0][0]*Cm1.N[1][0];
  D.nM[2][0] = D.nM[0][2];
  D.nM[1][1] = lambda*J2*Cm1.N[1][1]*Cm1.N[1][1]-(lambda*(J2-1)-2*G)*Cm1.N[1][1]*Cm1.N[1][1];
  D.nM[1][2] = lambda*J2*Cm1.N[1][1]*Cm1.N[0][1]-(lambda*(J2-1)-2*G)*Cm1.N[1][0]*Cm1.N[1][1];
  D.nM[2][1] = D.nM[1][2];
  D.nM[2][2] = lambda*J2*Cm1.N[0][1]*Cm1.N[0][1]-0.5*(lambda*(J2-1)-2*G)*(Cm1.N[0][0]*Cm1.N[1][1]+Cm1.N[0][1]*Cm1.N[0][1]);

  free__TensorLib__(C);
  free__TensorLib__(Cm1);

  return D;
}


/**************************************************************/
