#include "nl-partsol.h"

/**************************************************************/

State_Parameters compute_1PK_Stress_Tensor_Newtonian_Fluid( 
  State_Parameters Intput_SP,
  Material MatProp_p)
{
  /*
    Number of dimensions
  */
  int Ndim = NumberDimensions;

  /* 
    Output state parameter
  */
  State_Parameters Output_SP;

  /*
    Take information from input state parameters
  */
  Tensor P = memory_to_tensor__TensorLib__(Intput_SP.Stress,2);
  Tensor F = memory_to_tensor__TensorLib__(Intput_SP.F_n1_p,2);
  Tensor dFdt = memory_to_tensor__TensorLib__(Intput_SP.dFdt,2);
  double J = Intput_SP.J;

  /*
    Material parameters
  */
  double p0 = MatProp_p.ReferencePressure;
  double mu = MatProp_p.Viscosity;
  double n = MatProp_p.n_Macdonald_model;
  double K = MatProp_p.Compressibility;

  /*
    Auxiliar tensors
  */
  Tensor Fm1 = Inverse__TensorLib__(F);
  Tensor FmT = transpose__TensorLib__(Fm1);
  Tensor dFdt__x__Fm1 = matrix_product__TensorLib__(dFdt,Fm1);
  Tensor d = symmetrise__TensorLib__(dFdt__x__Fm1);
  Tensor d__x__FmT = matrix_product__TensorLib__(d, FmT);

  /*
    Auxiliar parameters
  */
  double tr_d = I1__TensorLib__(d);


  for(int i = 0 ; i < Ndim ; i++)
  {
    for(int j = 0 ; j < Ndim ; j++)
    {
      P.N[i][j] = 
      - J*(p0 + (K/n)*(pow(J,-n) - 1))*FmT.N[i][j]
      + 2.0*J*mu*d__x__FmT.N[i][j]
      - (2.0/3.0)*J*mu*tr_d*FmT.N[i][j];
    }
  }

  #if NumberDimensions == 2
    Intput_SP.Stress[4] = - (J*p0 + J*(K/n)*(pow(J,-n) - 1));
  #endif
  /*
    Free tensors 
  */
  free__TensorLib__(Fm1);
  free__TensorLib__(FmT);
  free__TensorLib__(dFdt__x__Fm1);
  free__TensorLib__(d);
  free__TensorLib__(d__x__FmT);

  return Output_SP;
}

/**************************************************************/

Tensor compute_stiffness_density_Newtonian_Fluid(
  Tensor GRAD_I, 
  Tensor GRAD_J,
  Tensor F,
  Tensor dFdt,
  double J, 
  double alpha4,
  Material MatProp_p)
{

  /*
    Number of dimensions
  */
  int Ndim = NumberDimensions;
    
  /* Material parameters */
  double p0 = MatProp_p.ReferencePressure;
  double mu = MatProp_p.Viscosity;
  double n = MatProp_p.n_Macdonald_model;
  double K = MatProp_p.Compressibility;
  
  /*
    Stifness density tensor
  */
  Tensor A = alloc__TensorLib__(2);

  /*
    Auxiliar variables
  */    
  Tensor Fm1 = Inverse__TensorLib__(F);
  Tensor FmT = transpose__TensorLib__(Fm1);

  Tensor dFdt_Fm1 = matrix_product__TensorLib__(dFdt,Fm1);
  Tensor d = symmetrise__TensorLib__(dFdt_Fm1);

  Tensor FmTGRAD_I = vector_linear_mapping__TensorLib__(FmT,GRAD_I);
  Tensor FmTGRAD_J = vector_linear_mapping__TensorLib__(FmT,GRAD_J);

  Tensor Fm1GRAD_o_FmTGRAD_IJ = dyadic_Product__TensorLib__(FmTGRAD_I,FmTGRAD_J);
  Tensor Fm1GRAD_o_FmTGRAD_JI = dyadic_Product__TensorLib__(FmTGRAD_J,FmTGRAD_I);

  Tensor d_Fm1GRAD_o_FmTGRAD_IJ = matrix_product__TensorLib__(d,Fm1GRAD_o_FmTGRAD_IJ);
  Tensor d_Fm1GRAD_o_FmTGRAD_JI = matrix_product__TensorLib__(d,Fm1GRAD_o_FmTGRAD_JI);

  Tensor Fm1GRAD_o_FmTGRAD_IJ_dFdt_Fm1 = matrix_product__TensorLib__(Fm1GRAD_o_FmTGRAD_IJ,dFdt_Fm1);
  Tensor Fm1GRAD_o_FmTGRAD_JI_dFdt_Fm1 = matrix_product__TensorLib__(Fm1GRAD_o_FmTGRAD_JI,dFdt_Fm1);

  double GRAD_I_dot_GRAD_J = inner_product__TensorLib__(GRAD_I,GRAD_J);

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      A.N[i][j] = 
      + (- J*p0 - J*(K/n)*(pow(J,-n) - 1) + K*pow(J,1-n) - (2.0/3.0)*alpha4*J*mu)*Fm1GRAD_o_FmTGRAD_IJ.N[i][j]
      + 2*J*mu*d_Fm1GRAD_o_FmTGRAD_IJ.N[i][j]
      + (J*p0 + J*(K/n)*(pow(J,-n)-1) + alpha4*J*mu)*Fm1GRAD_o_FmTGRAD_JI.N[i][j]
      - 2*J*mu*d_Fm1GRAD_o_FmTGRAD_JI.N[i][j]
      + alpha4*J*mu*(i==j)*GRAD_I_dot_GRAD_J
      - J*mu*GRAD_I_dot_GRAD_J*dFdt_Fm1.N[i][j]
      - J*mu*Fm1GRAD_o_FmTGRAD_JI_dFdt_Fm1.N[i][j]
      + (2.0/3.0)*J*mu*Fm1GRAD_o_FmTGRAD_IJ_dFdt_Fm1.N[i][j];
    }
  }

  /*
    Free memory
   */
  free__TensorLib__(Fm1);
  free__TensorLib__(FmT);
  free__TensorLib__(dFdt_Fm1);
  free__TensorLib__(d);
  free__TensorLib__(FmTGRAD_I);
  free__TensorLib__(FmTGRAD_J);
  free__TensorLib__(Fm1GRAD_o_FmTGRAD_IJ);
  free__TensorLib__(Fm1GRAD_o_FmTGRAD_JI);
  free__TensorLib__(d_Fm1GRAD_o_FmTGRAD_IJ);
  free__TensorLib__(d_Fm1GRAD_o_FmTGRAD_JI);
  free__TensorLib__(Fm1GRAD_o_FmTGRAD_IJ_dFdt_Fm1);
  free__TensorLib__(Fm1GRAD_o_FmTGRAD_JI_dFdt_Fm1);

  return A;
}   

/**************************************************************/
