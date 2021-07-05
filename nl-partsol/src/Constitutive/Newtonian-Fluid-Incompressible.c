#include "nl-partsol.h"

/**************************************************************/

State_Parameters compute_1PK_Stress_Tensor_Newtonian_Fluid_Incompressible( 
  State_Parameters Input_SP,
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
  Tensor P = memory_to_tensor__TensorLib__(Input_SP.Stress,2);
  Tensor F = memory_to_tensor__TensorLib__(Input_SP.F_n1_p,2);
  Tensor dFdt = memory_to_tensor__TensorLib__(Input_SP.dFdt,2);
  double Pressure = Input_SP.Pressure;
  double J = Input_SP.J;
  
  /*
    Material parameters
  */
  double mu = MatProp_p.Viscosity;

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
      - Pressure*FmT.N[i][j]
      + 2*J*mu*d__x__FmT.N[i][j]
      - (2.0/((double)Ndim))*J*mu*tr_d*FmT.N[i][j];
    }
  }
  
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

Tensor compute_stiffness_density_Newtonian_Fluid_Incompressible(
  Tensor GRAD_pA, 
  Tensor GRAD_pB,
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
  double mu = MatProp_p.Viscosity;
  
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

  Tensor FmTGRAD_pA = vector_linear_mapping__TensorLib__(FmT,GRAD_pA);
  Tensor FmTGRAD_pB = vector_linear_mapping__TensorLib__(FmT,GRAD_pB);

  Tensor Fm1GRAD_o_FmTGRAD_pAB = dyadic_Product__TensorLib__(FmTGRAD_pA,FmTGRAD_pB);
  Tensor Fm1GRAD_o_FmTGRAD_pBA = dyadic_Product__TensorLib__(FmTGRAD_pB,FmTGRAD_pA);

  Tensor d_Fm1GRAD_o_FmTGRAD_pAB = matrix_product__TensorLib__(d,Fm1GRAD_o_FmTGRAD_pAB);
  Tensor d_Fm1GRAD_o_FmTGRAD_pBA = matrix_product__TensorLib__(d,Fm1GRAD_o_FmTGRAD_pBA);

  Tensor Fm1GRAD_o_FmTGRAD_pAB_dFdt_Fm1 = matrix_product__TensorLib__(Fm1GRAD_o_FmTGRAD_pAB,dFdt_Fm1);
  Tensor Fm1GRAD_o_FmTGRAD_pBA_dFdt_Fm1 = matrix_product__TensorLib__(Fm1GRAD_o_FmTGRAD_pBA,dFdt_Fm1);

  double GRAD_A_dot_GRAD_B = inner_product__TensorLib__(GRAD_pA,GRAD_pB);

  for(int i = 0 ; i<Ndim ; i++)
  {

    for(int j = 0 ; j<Ndim ; j++)
    {
      A.N[i][j] += 
      - ((2.0/((double)Ndim))*alpha4*J*mu)*Fm1GRAD_o_FmTGRAD_pAB.N[i][j]
      + 2*J*mu*d_Fm1GRAD_o_FmTGRAD_pAB.N[i][j]
      + (alpha4*J*mu)*Fm1GRAD_o_FmTGRAD_pBA.N[i][j]
      - 2*J*mu*d_Fm1GRAD_o_FmTGRAD_pBA.N[i][j]
      + alpha4*J*mu*(i==j)*GRAD_A_dot_GRAD_B
      - J*mu*GRAD_A_dot_GRAD_B*dFdt_Fm1.N[i][j]
      - J*mu*Fm1GRAD_o_FmTGRAD_pBA_dFdt_Fm1.N[i][j]
      + (2.0/((double)Ndim))*J*mu*Fm1GRAD_o_FmTGRAD_pAB_dFdt_Fm1.N[i][j];
    }
  }

  /*
    Free memory
   */
  free__TensorLib__(Fm1);
  free__TensorLib__(FmT);
  free__TensorLib__(dFdt_Fm1);
  free__TensorLib__(d);
  free__TensorLib__(FmTGRAD_pA);
  free__TensorLib__(FmTGRAD_pB);
  free__TensorLib__(Fm1GRAD_o_FmTGRAD_pAB);
  free__TensorLib__(Fm1GRAD_o_FmTGRAD_pBA);
  free__TensorLib__(d_Fm1GRAD_o_FmTGRAD_pAB);
  free__TensorLib__(d_Fm1GRAD_o_FmTGRAD_pBA);
  free__TensorLib__(Fm1GRAD_o_FmTGRAD_pAB_dFdt_Fm1);
  free__TensorLib__(Fm1GRAD_o_FmTGRAD_pBA_dFdt_Fm1);

  return A;
}   

/**************************************************************/
