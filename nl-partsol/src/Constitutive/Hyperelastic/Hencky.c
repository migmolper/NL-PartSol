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

int compute_Kirchhoff_Stress_Hencky__Constitutive__(State_Parameters IO_State,
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
  double AA[9] = {Lame + 2 * G, Lame, Lame, Lame,        Lame + 2 * G,
                  Lame,         Lame, Lame, Lame + 2 * G};

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

  //! Compute deformation energy
  *(IO_State.W) = 0.5 * (T_ppal[0] * E_hencky[0] + T_ppal[1] * E_hencky[1] +
                         T_ppal[2] * E_hencky[2]);

  return EXIT_SUCCESS;
}

/**************************************************************/

int compute_stiffness_density_Hencky__Constitutive__(
    double *Stiffness_density, const double *dN_alpha_n1,
    const double *dN_beta_n1, const State_Parameters IO_State,
    Material MatProp) {

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;

  // Material parameters
  double E = MatProp.E;
  double nu = MatProp.nu;
  double Lame = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  double G = E / (2.0 * (1.0 + nu));

#if NumberDimensions == 2
  double b[4];
  double eigvec_b[4] = {0.0, 0.0, 0.0, 0.0};
  double eigval_b[2] = {0.0, 0.0};
  double eigvec_T[4] = {0.0, 0.0, 0.0, 0.0};
  double eigval_T[2] = {0.0, 0.0};
  double u__o__v[2][2];
  double AA[4] = {Lame + 2 * G, Lame, Lame, Lame + 2 * G};
  Stiffness_density[0] = 0.0;
  Stiffness_density[1] = 0.0;
  Stiffness_density[2] = 0.0;
  Stiffness_density[3] = 0.0;
#else
  double b[9];
  double eigvec_b[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double eigval_b[3] = {0.0, 0.0, 0.0};
  double eigvec_T[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double eigval_T[3] = {0.0, 0.0, 0.0};
  double T_aux[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double u__o__v[3][3];
  double AA[9] = {Lame + 2 * G, Lame, Lame, Lame,        Lame + 2 * G,
                  Lame,         Lame, Lame, Lame + 2 * G};
  Stiffness_density[0] = 0.0;
  Stiffness_density[1] = 0.0;
  Stiffness_density[2] = 0.0;
  Stiffness_density[3] = 0.0;
  Stiffness_density[4] = 0.0;
  Stiffness_density[5] = 0.0;
  Stiffness_density[6] = 0.0;
  Stiffness_density[7] = 0.0;
  Stiffness_density[8] = 0.0;
#endif
  const double *D_phi_n1 = IO_State.D_phi_n1;

  //! Left cauchy-green tensor
  left_Cauchy_Green__Particles__(b, D_phi_n1);

  //! Compute Hencky strain in the ppal directions
  STATUS = sym_eigen_analysis__TensorLib__(eigval_b, eigvec_b, b);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,
            "" RED "Error in sym_eigen_analysis__TensorLib__()" RESET "\n");
    return EXIT_FAILURE;
  }

  STATUS = sym_eigen_analysis__TensorLib__(eigval_T, eigvec_T, IO_State.Stress);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,
            "" RED "Error in sym_eigen_analysis__TensorLib__()" RESET "\n");
    return EXIT_FAILURE;
  }


  // Do the diadic product of gradient directions
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      u__o__v[i][j] = dN_beta_n1[i] * dN_alpha_n1[j];
    }
  }


  for (unsigned A = 0; A < Ndim; A++) {

    double u_A = 0.0;
    double v_A = 0.0;
    for (unsigned i = 0; i < Ndim; i++) 
    {
      u_A += dN_alpha_n1[i]*eigvec_b[A + i * Ndim];
      v_A += dN_beta_n1[i]*eigvec_b[A + i * Ndim];
    }

    for (unsigned B = 0; B < Ndim; B++) {

      double u_B = 0.0;
      double v_B = 0.0;
      for (unsigned i = 0; i < Ndim; i++) 
      {
        u_B += dN_alpha_n1[i]*eigvec_b[B + i * Ndim];        
        v_B += dN_beta_n1[i]*eigvec_b[B + i * Ndim];
      }

      double C_ep_AB = AA[A * Ndim + B];
      double v_A__dot__u_B = u_B*v_A;
      double u_A__dot__v_B = u_A*v_B;    
      double u_B__dot__v_B = u_B*v_B;
      
      for (unsigned i = 0; i < Ndim; i++) {
        for (unsigned j = 0; j < Ndim; j++) {
          Stiffness_density[i * Ndim + j] +=
              C_ep_AB * u_A__dot__v_B * eigvec_b[A + i * Ndim] * eigvec_b[B + j * Ndim];

          if (A != B) {
            if (fabs(eigval_b[B] - eigval_b[A]) > 1E-14) {
              Stiffness_density[i * Ndim + j] +=
                  0.5 *
                  ((eigval_T[B] - eigval_T[A]) /
                   (eigval_b[B] - eigval_b[A])) *
                  (eigval_b[B] * u_B__dot__v_B*(eigvec_b[A + i * Ndim] * eigvec_b[A + j * Ndim]) +
                   eigval_b[A] * v_A__dot__u_B*(eigvec_b[A + i * Ndim] * eigvec_b[B + j * Ndim]));
            } 
          }          
        }
      }
    }
  }

  // Assemble the geometrical contribution to the tanget matrix
  for (unsigned i = 0; i < Ndim; i++) {
    for (unsigned j = 0; j < Ndim; j++) {
      for (unsigned k = 0; k < Ndim; k++) {
        Stiffness_density[i * Ndim + j] +=
            -IO_State.Stress[i * Ndim + k] * u__o__v[k][j];
      }
    }
  }

  return EXIT_SUCCESS;
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
            T_A * eigvec_b[A + i * Ndim] * eigvec_b[A + j * Ndim];
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
