/**
 * @file Mooney-Rivlin.c
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Constitutive/Hyperelastic/Mooney-Rivlin.h"
#include "Globals.h"

/**************************************************************/

double energy_Mooney_Rivlin(Tensor C, double J, Material MatProp_p) {
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /* Material parameters */
  //  double mu_1 = MatProp_p.mu_Ogden[0];
  //  double mu_2 = MatProp_p.mu_Ogden[1];
  double nu = MatProp_p.nu;
  double E = MatProp_p.E;
  double K = E / (3 * (1 - 2 * nu));

//  double I1_C = I1__TensorLib__(C);
//  double I2_C = I2__TensorLib__(C);
//  double I3_C = J * J;

  double f_J = 0.5 * K * log(J) * log(J);

  double W = 0.0; // 0.5*mu_1*(I1_C - 3.0) - 0.5*mu_2*(I2_C/I3_C - 3.0) + f_J;

  return W;
}

/**************************************************************/

State_Parameters
compute_1PK_Stress_Tensor_Mooney_Rivlin(State_Parameters Intput_SP,
                                        Material MatProp_p) {
  /* Number of dimensions */
  int Ndim = NumberDimensions;

  /*
    Output state parameter
  */
  State_Parameters Output_SP;

  /* Get information from the state parameter */
  Tensor F = memory_to_tensor__TensorLib__(Intput_SP.D_phi_n1, 2);
  Tensor P = memory_to_tensor__TensorLib__(Intput_SP.Stress, 2);
  double J = Intput_SP.J;

  /*
    Compute Auxiliar tensors
  */
  Tensor Fm1 = Inverse__TensorLib__(F);
  Tensor C = right_Cauchy_Green__Particles__(F);
  Tensor FC = matrix_product_old__TensorLib__(F, C);

  /* Material parameters */
  //  double mu_1 = MatProp_p.mu_Ogden[0];
  //  double mu_2 = MatProp_p.mu_Ogden[1];
  double nu = MatProp_p.nu;
  double E = MatProp_p.E;
  double K = E / (3 * (1 - 2 * nu));

//  double I1_C = I1__TensorLib__(C);
//  double I2_C = I2__TensorLib__(C);

//  for (int i = 0; i < Ndim; i++) {
//    for (int j = 0; j < Ndim; j++) {
      //      P.N[i][j] = (mu_1 - mu_2*I1_C)*F.N[i][j] - mu_2*FC.N[i][j] +
      //      (mu_2*I2_C + K*log(J))*Fm1.N[j][i];
//    }
//  }

  /*
    Free tensors
  */
  free__TensorLib__(Fm1);
  free__TensorLib__(C);
  free__TensorLib__(FC);

  return Output_SP;
}

/**************************************************************/

Tensor compute_stiffness_density_Mooney_Rivlin(Tensor GRAD_I, Tensor GRAD_J,
                                               Tensor F, double J,
                                               Material MatProp) {

  /*
    Number of dimensions
  */
  int Ndim = NumberDimensions;

  /*
    Stifness density tensor
  */
  Tensor A = alloc__TensorLib__(2);

  /*
    Compute Auxiliar tensors
  */
  Tensor Fm1 = Inverse__TensorLib__(F);
  Tensor C = right_Cauchy_Green__Particles__(F);
  //  Tensor b = left_Cauchy_Green__Particles__(F);
  Tensor FC = matrix_product_old__TensorLib__(F, C);

  /* Material parameters */
  //  double mu_1 = MatProp.mu_Ogden[0];
  //  double mu_2 = MatProp.mu_Ogden[1];
  double nu = MatProp.nu;
  double E = MatProp.E;
  double K = E / (3 * (1 - 2 * nu));

//  double I1_C = I1__TensorLib__(C);
//  double I2_C = I2__TensorLib__(C);

  /*
    Auxiliar variables
  */
  Tensor GRAD_o_GRAD = dyadic_Product__TensorLib__(GRAD_I, GRAD_J);
  double GRAD_dot_GRAD_IJ = inner_product__TensorLib__(GRAD_I, GRAD_J);

  Tensor FGRAD_I = vector_linear_mapping__TensorLib__(F, GRAD_I);
  Tensor FGRAD_J = vector_linear_mapping__TensorLib__(F, GRAD_J);
  Tensor FGRAD_o_FGRAD = dyadic_Product__TensorLib__(FGRAD_I, FGRAD_J);
  double FGRAD_dot_FGRAD_IJ = inner_product__TensorLib__(FGRAD_I, FGRAD_J);

  Tensor FmT = Inverse__TensorLib__(Fm1);
  Tensor FmTGRAD_I = vector_linear_mapping__TensorLib__(FmT, GRAD_I);
  Tensor FmTGRAD_J = vector_linear_mapping__TensorLib__(FmT, GRAD_J);
  Tensor Fm1GRAD_o_FmTGRAD = dyadic_Product__TensorLib__(FmTGRAD_I, FmTGRAD_J);

  //  for(int i = 0 ; i<Ndim ; i++)
  //  {
  //    for(int j = 0 ; j<Ndim ; j++)
  //    {
  //      A.N[i][j] = mu_1*GRAD_o_GRAD.N[i][j] - 2*mu_2*FGRAD_o_FGRADN[i][j] -
  //      mu_2*I1_C*GRAD_o_GRAD.N[i][j] - mu_2*FGRAD_dot_FGRAD_IJ*(i==j) -
  //      mu_2*FGRAD_o_FGRAD.N[j][i] + mu_2*GRAD_dot_GRAD_IJ*b.N[i][j] -
  //    }
  //  }

  /*
    Free memory
   */
  free__TensorLib__(Fm1);
  free__TensorLib__(FmT);
  free__TensorLib__(FmTGRAD_I);
  free__TensorLib__(FmTGRAD_J);
  free__TensorLib__(Fm1GRAD_o_FmTGRAD);

  return A;
}

/**************************************************************/
