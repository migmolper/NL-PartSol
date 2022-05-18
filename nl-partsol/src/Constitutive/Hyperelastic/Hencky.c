

#include "Constitutive/Hyperelastic/Hencky.h"

/**************************************************************/

double energy_Hencky__Constitutive__() {}

/**************************************************************/

int compute_Kirchhoff_Stress_Hencky(State_Parameters IO_State,
                                    Material MatProp) {
  int STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  double b[4];
  double eigvec_b[4] = {0.0, 0.0, 0.0, 0.0};
  double eigval_b[3] = {0.0, 0.0, 1.0};
  double T_aux[4] = {0.0, 0.0, 0.0, 0.0};
#else
  double b[9];
  double eigvec_b[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double eigval_b[3] = {0.0, 0.0, 0.0};
  double T_aux[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif
  double E_hencky[3] = {0.0, 0.0, 0.0};
  double tr_E_hencky;
  const double *D_phi_n1 = IO_State.D_phi_n1;
  double *T = IO_State.Stress;
  double T_A;

  // Material parameters

  double Lambda =
      MatProp.E * MatProp.nu / ((1.0 - 2.0 * MatProp.nu) * (1.0 + MatProp.nu));
  double G = MatProp.E / (2.0 * (1.0 + MatProp.nu));

  left_Cauchy_Green__Particles__(b, D_phi_n1);

  STATUS = sym_eigen_analysis__TensorLib__(eigval_b, eigvec_b, b);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,
            "" RED "Error in sym_eigen_analysis__TensorLib__()" RESET "\n");
    return EXIT_FAILURE;
  }

  E_hencky[0] = 0.5 * log(eigval_b[0]);
  E_hencky[1] = 0.5 * log(eigval_b[1]);
  E_hencky[2] = 0.5 * log(eigval_b[2]);

  tr_E_hencky = E_hencky[0] + E_hencky[1] + E_hencky[2];

  for (unsigned A = 0; A < Ndim; A++) {

    T_A = Lambda * tr_E_hencky + 2 * G * E_hencky[A];

    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        T_aux[i * Ndim + j] +=
            T_A * eigvec_b[A * Ndim + i] * eigvec_b[A * Ndim + j];
      }
    }
  }

#if NumberDimensions == 2
  T[0] = T_aux[0];
  T[1] = T_aux[1];
  T[2] = T_aux[2];
  T[3] = T_aux[3];
  T[4] = Lambda * tr_E_hencky + 2 * G * E_hencky[2];
#else
  T[0] = T_aux[0];
  T[1] = T_aux[1];
  T[2] = T_aux[2];
  T[3] = T_aux[3];
  T[4] = T_aux[4];
  T[5] = T_aux[5];
  T[6] = T_aux[6];
  T[7] = T_aux[7];
  T[8] = T_aux[8];
#endif

  return EXIT_SUCCESS;
}

/**************************************************************/
