#include "Constitutive/Fluid/Bingham-Fluid.h"

static double compute_Bingham_viscosity(Tensor, Material);

static double compute_p_Tait_Murnaghan(double, Material);

/**************************************************************/

State_Parameters
compute_1PK_Stress_Tensor_Bingham_Fluid(State_Parameters Intput_SP,
                                        Material MatProp_p) {
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
  Tensor P = memory_to_tensor__TensorLib__(Intput_SP.Stress, 2);
  Tensor F = memory_to_tensor__TensorLib__(Intput_SP.D_phi_n1, 2);
  Tensor dFdt = memory_to_tensor__TensorLib__(Intput_SP.dFdt, 2);
  double J = Intput_SP.J;

  /*
    Tensor variables
  */
  Tensor Fm1 = Inverse__TensorLib__(F);
  Tensor FmT = transpose__TensorLib__(Fm1);
  Tensor dFdt__x__Fm1 = matrix_product_old__TensorLib__(dFdt, Fm1);
  Tensor d;// = symmetrise__TensorLib__(dFdt__x__Fm1);
  Tensor d__x__FmT = matrix_product_old__TensorLib__(d, FmT);
  double tr_d = I1__TensorLib__(d);

  /*
    Material parameters
  */
  double mu = compute_Bingham_viscosity(d, MatProp_p);
  double p = compute_p_Tait_Murnaghan(J, MatProp_p);

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      P.N[i][j] = -J * p * FmT.N[i][j] + 2 * J * mu * d__x__FmT.N[i][j] -
                  (2.0 / ((double)Ndim)) * J * mu * tr_d * FmT.N[i][j];
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

static double compute_p_Tait_Murnaghan(double J, Material MatProp_p) {
  double p0 = MatProp_p.ReferencePressure;
  double n = MatProp_p.n_Macdonald_model;
  double K = MatProp_p.Compressibility;

  return p0 + (K / n) * (pow(J, -n) - 1);
}

/**************************************************************/

static double compute_Bingham_viscosity(Tensor d, Material MatProp_p) {
  double mu_0 = MatProp_p.Viscosity;
  double taub_yield = MatProp_p.kappa_0;
  double m = MatProp_p.fluidity_param;
  double EVPS = sqrt(2 * inner_product__TensorLib__(d, d));

  return mu_0 + (taub_yield / EVPS) * (1 - exp(-m * EVPS));
}

/**************************************************************/

Tensor compute_stiffness_density_Bingham_Fluid(Tensor GRAD_I, Tensor GRAD_J,
                                               Tensor F, Tensor dFdt, double J,
                                               double alpha4,
                                               Material MatProp_p) {

  /*
    Number of dimensions
  */
  int Ndim = NumberDimensions;

  /*
    Stifness density tensor
  */
  Tensor A = alloc__TensorLib__(2);

  /*
    Tensor variables
  */
  Tensor Fm1 = Inverse__TensorLib__(F);
  Tensor FmT = transpose__TensorLib__(Fm1);

  Tensor dFdt_Fm1 = matrix_product_old__TensorLib__(dFdt, Fm1);
  Tensor d;// = symmetrise__TensorLib__(dFdt_Fm1);

  Tensor FmTGRAD_I = vector_linear_mapping__TensorLib__(FmT, GRAD_I);
  Tensor FmTGRAD_J = vector_linear_mapping__TensorLib__(FmT, GRAD_J);

  Tensor Fm1GRAD_o_FmTGRAD_IJ =
      dyadic_Product__TensorLib__(FmTGRAD_I, FmTGRAD_J);
  Tensor Fm1GRAD_o_FmTGRAD_JI =
      dyadic_Product__TensorLib__(FmTGRAD_J, FmTGRAD_I);

  Tensor d_Fm1GRAD_o_FmTGRAD_IJ =
      matrix_product_old__TensorLib__(d, Fm1GRAD_o_FmTGRAD_IJ);
  Tensor d_Fm1GRAD_o_FmTGRAD_JI =
      matrix_product_old__TensorLib__(d, Fm1GRAD_o_FmTGRAD_JI);

  Tensor Fm1GRAD_o_FmTGRAD_IJ_dFdt_Fm1 =
      matrix_product_old__TensorLib__(Fm1GRAD_o_FmTGRAD_IJ, dFdt_Fm1);
  Tensor Fm1GRAD_o_FmTGRAD_JI_dFdt_Fm1 =
      matrix_product_old__TensorLib__(Fm1GRAD_o_FmTGRAD_JI, dFdt_Fm1);

  double GRAD_I_dot_GRAD_J = inner_product__TensorLib__(GRAD_I, GRAD_J);

  /*
    Material parameters
  */
  double p0 = MatProp_p.ReferencePressure;
  double mu = compute_Bingham_viscosity(d, MatProp_p);
  double n = MatProp_p.n_Macdonald_model;
  double K = MatProp_p.Compressibility;

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      A.N[i][j] =
          +(-J * p0 - J * (K / n) * (pow(J, -n) - 1) + K * pow(J, 1 - n) -
            (2.0 / ((double)Ndim)) * alpha4 * J * mu) *
              Fm1GRAD_o_FmTGRAD_IJ.N[i][j] +
          2 * J * mu * d_Fm1GRAD_o_FmTGRAD_IJ.N[i][j] +
          (J * p0 + J * (K / n) * (pow(J, -n) - 1) + alpha4 * J * mu) *
              Fm1GRAD_o_FmTGRAD_JI.N[i][j] -
          2 * J * mu * d_Fm1GRAD_o_FmTGRAD_JI.N[i][j] +
          alpha4 * J * mu * (i == j) * GRAD_I_dot_GRAD_J -
          J * mu * GRAD_I_dot_GRAD_J * dFdt_Fm1.N[i][j] -
          J * mu * Fm1GRAD_o_FmTGRAD_JI_dFdt_Fm1.N[i][j] +
          (2.0 / ((double)Ndim)) * J * mu *
              Fm1GRAD_o_FmTGRAD_IJ_dFdt_Fm1.N[i][j];
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
