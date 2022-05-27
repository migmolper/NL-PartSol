/**
 * @file EigenErosion.c
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Constitutive/Fracture/EigenErosion.h"
#include "Globals.h"

/**************************************************************/

int Eigenerosion__Constitutive__(unsigned p, const double *Damage_n,
                                 double *Damage_n1, const double *W,
                                 const double *J_n1, const double *Vol_0,
                                 const double *Stress_p, Material MatProp_p,
                                 const ChainPtr Beps_p, double DeltaX) {

  int STATUS = EXIT_SUCCESS;

  //! Define auxiliar variable
  double Ceps_p, G_p, Gf_p;
  double V0_p, Jacobian_p, V_p, W_p;
  double V0_q, Jacobian_q, V_q, W_q;
  double sum_V_x_W, sum_V;
  unsigned q;

  //! For the current particle get first principal stress
  double eigval_stress_p[3] = {0.0, 0.0, 0.0};
#if NumberDimensions == 2
  double eigvec_stress_p[4] = {0.0, 0.0, 0.0, 0.0};
#else
  double eigvec_stress_p[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
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

  if ((Damage_n[p] < 1.0) && (eigval_stress_p[0] > 0.0)) {

    //! Material parameters
    Ceps_p = MatProp_p.Ceps; // Normalizing constant
    Gf_p = MatProp_p.Gf;     // Griffit fracture criterium

    //! Compute current volume of the particle p
    Jacobian_p = J_n1[p];
    V0_p = Vol_0[p];
    V_p = V0_p * Jacobian_p;

    //! Internal work of particle p
    W_p = W[p];

    //! Add contribution of p to the individual sum
    sum_V = V_p;
    sum_V_x_W = V_p * W_p;

    /*! Iterate in the neighbours of the particle */
    ChainPtr idx_Beps_p = Beps_p;

    while (idx_Beps_p != NULL) {

      q = idx_Beps_p->Idx;

      //! Compute current volume of the particle p
      Jacobian_q = J_n1[q];
      V0_q = Vol_0[q];
      V_q = V0_q * Jacobian_q;

      //! Internal work of particle q
      W_q = W[q];

      //! Add contribution of p to the individual sum
      sum_V += V_q;
      if (Damage_n[q] < 1.0)
        sum_V_x_W += V_q * W_q;

      //! Next particle
      idx_Beps_p = idx_Beps_p->next;
    }

    //!  Compute energy-release rate for the particle
    G_p = (Ceps_p * DeltaX / sum_V) * sum_V_x_W;

    //! Fracture criterium
    if (G_p > Gf_p) {
      Damage_n1[p] = 1.0;
    }
  }

  return EXIT_SUCCESS;
}

/**************************************************************/