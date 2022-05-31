/**
 * @file EigenSoftening.c
 * @author Miguel Molinos (@migmolper)
 * @brief   Pedro Navas, Rena C. Yu, Bo Li & Gonzalo Ruiz.
 * Modeling the dynamic fracture in concrete:
 * an eigensoftening meshfree approach.
 * International Journal of Impact Engineering. 113 (2018) 9-20
 *  NOTE : Here notation is the same as in the paper.
 * Inputs :
 * -> Ji_k0 : Matrix with the value of the damage parameter.
 * -> Mass : Matrix with the mass of the GP.
 * -> StrainF : Value of the strain field at the failure init.
 * -> Beps : Table with the list of neighbours per GP.
 * -> Neps : Number of neighbours per GP
 * -> Num_GP : Number of GP of the mesh.
 * @version 0.1
 * @date 2022-05-30
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Constitutive/Fracture/EigenSoftening.h"
#include "Globals.h"

/*******************************************************/

int Eigensoftening__Constitutive__(unsigned p, const double *Damage_n,
                                    double *Damage_n1, const double *Strain_p,
                                    const double *StrainF_n_p, double *StrainF_n1_p,
                                    const double *Mass, const double *Stress,
                                    Material MatProp_p, const ChainPtr Beps_p) {
  int STATUS = EXIT_SUCCESS;

  //! Define auxiliar variable
  double ft_p, wcrit_p, heps_p, Teps_p, Damage_n1_aux_p;
  double m_p;
  double m_q;
  double sum_m_x_T_ppal, sum_m;
  unsigned q;

  //! For the current particle get first principal stress
  double eigval_stress_p[3] = {0.0, 0.0, 0.0};
  double eigval_stress_q[3];
  double eigval_strain_p[3] = {0.0, 0.0, 0.0};
#if NumberDimensions == 2
  double eigvec_stress_p[4] = {0.0, 0.0, 0.0, 0.0};
  double eigvec_stress_q[4];
  double eigvec_strain_p[4] = {0.0, 0.0, 0.0, 0.0};
  const double * Stress_p = &(Stress[p * 5]);
#else
  double eigvec_stress_p[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double eigvec_stress_q[9];
  double eigvec_strain_p[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double * Stress_p = &(Stress[p * 9]);
#endif

  STATUS = sym_eigen_analysis__TensorLib__(eigval_stress_p, eigvec_stress_p,
                                           Stress_p);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr,
            "" RED "Error in sym_eigen_analysis__TensorLib__()" RESET "\n");
    return EXIT_FAILURE;
  }

#if NumberDimensions == 2
  eigval_stress_p[2] = Stress_p[4];
#endif

  //! Material parameters
  ft_p = MatProp_p.ft;     // Tensile strengt of the material
  heps_p = MatProp_p.heps; // Bandwidth of the cohesive fracture (Bazant)
  wcrit_p = MatProp_p.wcrit;     // critical opening displacement

  if ((Damage_n[p] == 0.0) && (eigval_stress_p[0] > 0.0)) {

    //! For the current particle get the mass */
    m_p = Mass[p];

    //! Add contribution of p to the individual sum
    sum_m = m_p;
    sum_m_x_T_ppal = m_p * eigval_stress_p[0];

    /*! Iterate in the neighbours of the particle */
    ChainPtr idx_Beps_p = Beps_p;

    while (idx_Beps_p != NULL) {

      q = idx_Beps_p->Idx;

      //! Mass of the particle q
      m_q = Mass[q];

      //! Add contribution of p to the individual sum
      sum_m += m_q;
      if (Damage_n[q] < 1.0) {

#if NumberDimensions == 2
        const double * Stress_q = &(Stress[q * 5]);
#else
        const double * Stress_q = &(Stress[q * 9]);
#endif

        STATUS = sym_eigen_analysis__TensorLib__(eigval_stress_q,
                                                 eigvec_stress_q, Stress_q);
        if (STATUS == EXIT_FAILURE) {
          fprintf(stderr,
                  "" RED "Error in sym_eigen_analysis__TensorLib__()" RESET
                  "\n");
          return EXIT_FAILURE;
        }

        sum_m_x_T_ppal = m_q * eigval_stress_q[0];
      }

      //! Next particle
      idx_Beps_p = idx_Beps_p->next;
    }

    /* Get the equivalent critical stress */
    Teps_p = sum_m_x_T_ppal / sum_m;

    //! Check fracture criterium
    if (Teps_p > ft_p) {

      STATUS = sym_eigen_analysis__TensorLib__(eigval_strain_p, eigvec_strain_p,
                                               Strain_p);
      if (STATUS == EXIT_FAILURE) {
        fprintf(stderr,
                "" RED "Error in sym_eigen_analysis__TensorLib__()" RESET "\n");
        return EXIT_FAILURE;
      }

      //! Strain during fracture
      *(StrainF_n1_p) = eigval_strain_p[0];
    }

  } else if ((Damage_n[p] != 1.0) && (*(StrainF_n_p) > 0)) {

    STATUS = sym_eigen_analysis__TensorLib__(eigval_strain_p, eigvec_strain_p,
                                             Strain_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in sym_eigen_analysis__TensorLib__()" RESET "\n");
      return EXIT_FAILURE;
    }

    //! Compute particle damage
    Damage_n1_aux_p = (eigval_strain_p[0] - *(StrainF_n_p)) * heps_p / wcrit_p;
    Damage_n1[p] = DMIN(1.0, DMAX(Damage_n1_aux_p, Damage_n[p]));
  }

  return EXIT_SUCCESS;
}