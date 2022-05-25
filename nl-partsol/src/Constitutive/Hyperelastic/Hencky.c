/**
 * @file Hencky.c
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Constitutive/Hyperelastic/Hencky.h"
#include "Globals.h"

/**************************************************************/

/**
 * @brief
 *
 * @param[out] AA Elastic matrix
 * @param[in] Lame Lamé parameter
 * @param[in] G Shear modulus
 */
static void __elastic_tangent(double *AA, double Lame, double G);
/**************************************************************/

/**
 * @brief
 *
 * @param[out] T_ppal Elastic stress tensor
 * @param[in] E_hencky Henky strain
 * @param[in] AA Elastic matrix
 */
static void __elastic_stress_ppal(double *T_ppal, const double *E_hencky,
                                  const double *AA);
/**************************************************************/

/**
 * @brief
 *
 * @param T_xyz Stress tensor in the xyz coordinates
 * @param T_ppal Stress tensor in the ppal coordinates
 * @param eigvec_b Eigenvectors of b
 */
static void __elastic_stress_xyz(double *T_xyz, const double *T_ppal,
                                 const double *eigvec_b);

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
  double T_ppal[3] = {0.0, 0.0, 0.0};
  const double *D_phi_n1 = IO_State.D_phi_n1;

  // Material parameters
  double E = MatProp.E;
  double nu = MatProp.nu;
  double Lame = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  double G = E / (2.0 * (1.0 + nu));
  double AA[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  //! Compute the elastic matrix and the compliance for the tangent
  __elastic_tangent(AA, Lame, G);

  //! Left cauchy-green tensor
  left_Cauchy_Green__Particles__(b, D_phi_n1);

  //! Compute Hencky strain in the ppal directions
  STATUS = sym_eigen_analysis__TensorLib__(eigval_b, eigvec_b, b);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,
            "" RED "Error in sym_eigen_analysis__TensorLib__()" RESET "\n");
    return EXIT_FAILURE;
  }

  E_hencky[0] = 0.5 * log(eigval_b[0]);
  E_hencky[1] = 0.5 * log(eigval_b[1]);
  E_hencky[2] = 0.5 * log(eigval_b[2]);

  //! Compute elastic stress tensor in the ppal directions
  __elastic_stress_ppal(T_ppal, E_hencky, AA);

  //! Rotate stress tensor ppal --> xyz
  __elastic_stress_xyz(IO_State.Stress, T_ppal, eigvec_b);

  //! Tangent matrix (ppal)
#if NumberDimensions == 2
  IO_State.C_ep[0] = AA[0];
  IO_State.C_ep[1] = AA[1];
  IO_State.C_ep[2] = AA[3];
  IO_State.C_ep[3] = AA[4];
#else
  IO_State.C_ep[0] = AA[0];
  IO_State.C_ep[1] = AA[1];
  IO_State.C_ep[2] = AA[2];
  IO_State.C_ep[3] = AA[3];
  IO_State.C_ep[4] = AA[4];
  IO_State.C_ep[5] = AA[5];
  IO_State.C_ep[6] = AA[6];
  IO_State.C_ep[7] = AA[7];
  IO_State.C_ep[8] = AA[8];
#endif

  //! Compute deformation energy
  *(IO_State.W) = 0.5 * (T_ppal[0] * E_hencky[0] + T_ppal[1] * E_hencky[1] +
                         T_ppal[2] * E_hencky[2]);

  return EXIT_SUCCESS;
}

/**************************************************************/

static void __elastic_tangent(double *AA, double Lame, double G) {

  AA[0] = Lame + 2 * G;
  AA[1] = Lame;
  AA[2] = Lame;

  AA[3] = Lame;
  AA[4] = Lame + 2 * G;
  AA[5] = Lame;

  AA[6] = Lame;
  AA[7] = Lame;
  AA[8] = Lame + 2 * G;
}

/**************************************************************/

static void __elastic_stress_ppal(double *T_ppal, const double *E_hencky,
                                  const double *AA) {

  T_ppal[0] = AA[0] * E_hencky[0] + AA[1] * E_hencky[1] + AA[2] * E_hencky[2];
  T_ppal[1] = AA[3] * E_hencky[0] + AA[4] * E_hencky[1] + AA[5] * E_hencky[2];
  T_ppal[2] = AA[6] * E_hencky[0] + AA[7] * E_hencky[1] + AA[8] * E_hencky[2];
}

/**************************************************************/

static void __elastic_stress_xyz(double *T_xyz, const double *T_ppal,
                                 const double *eigvec_b) {

  unsigned Ndim = NumberDimensions;

#if NumberDimensions == 2
  double T_aux[4] = {0.0, 0.0, 0.0, 0.0};
#else
  double T_aux[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#endif

  double T_A;

  for (unsigned A = 0; A < Ndim; A++) {

    T_A = T_ppal[A];

    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        T_aux[i * Ndim + j] +=
            T_A * eigvec_b[A * Ndim + i] * eigvec_b[A * Ndim + j];
      }
    }
  }

#if NumberDimensions == 2
  T_xyz[0] = T_aux[0];
  T_xyz[1] = T_aux[1];
  T_xyz[2] = T_aux[2];
  T_xyz[3] = T_aux[3];
  T_xyz[4] = T_ppal[2];
#else
  T_xyz[0] = T_aux[0];
  T_xyz[1] = T_aux[1];
  T_xyz[2] = T_aux[2];
  T_xyz[3] = T_aux[3];
  T_xyz[4] = T_aux[4];
  T_xyz[5] = T_aux[5];
  T_xyz[6] = T_aux[6];
  T_xyz[7] = T_aux[7];
  T_xyz[8] = T_aux[8];
#endif
}

/**************************************************************/
