
#include "Formulations/Displacements-WaterPressure/U-pw-Analisys.h"

/*********************************************************************/

Fields allocate_upw_vars__Fields__(int NumParticles) {
  int Ndim = NumberDimensions;
  Fields Phi;

  /*!
    Global coordinates
  */
  Phi.x_GC = allocZ__MatrixLib__(NumParticles, Ndim);

  /*!
    Natural coordinates (Vectorial)
  */
  Phi.x_EC = allocZ__MatrixLib__(NumParticles, Ndim);

  /*!
    Displacement and increment of displacement fields (Vectorial)
  */
  Phi.dis = allocZ__MatrixLib__(NumParticles, Ndim);
  Phi.D_dis = allocZ__MatrixLib__(NumParticles, Ndim);
 
  /*!
    Velocity field (Vectorial)
  */
  Phi.vel = allocZ__MatrixLib__(NumParticles, Ndim);
  
  /*!
    Acceleration field (Vectorial)
  */
  Phi.acc = allocZ__MatrixLib__(NumParticles, Ndim);

  /*!
    Strain field (Tensor)
  */
  Phi.Strain = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);

  /*!
    Deformation gradient field (Tensor) + Initialise it with the indentity
  */
  Phi.F_n = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
  Phi.F_n1 = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
  Phi.DF = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);

  for (int p = 0; p < NumParticles; p++) {
    for (int i = 0; i < Ndim; i++) {
      Phi.F_n.nM[p][i + i * Ndim] = 1.0;
      Phi.F_n1.nM[p][i + i * Ndim] = 1.0;
      Phi.DF.nM[p][i + i * Ndim] = 1.0;
    }
  }

  Phi.dt_F_n = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
  Phi.dt_F_n1 = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
  Phi.dt_DF = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);

  /*!
    Allocate and initialise the Jacobian of the deformation gradient
    ant its rate
  */
  Phi.J_n = allocZ__MatrixLib__(NumParticles, 1);
  Phi.J_n1 = allocZ__MatrixLib__(NumParticles, 1);
  Phi.dJ_dt = allocZ__MatrixLib__(NumParticles, 1);

  for (int p = 0; p < NumParticles; p++) {
    Phi.J_n.nV[p] = 1.0;
    Phi.J_n1.nV[p] = 1.0;
  }

  /*!
    Elastic left Cauchy-Green tensor
  */
#if NumberDimensions == 2
  Phi.b_e_n = allocZ__MatrixLib__(NumParticles, 5);
  Phi.b_e_n1 = allocZ__MatrixLib__(NumParticles, 5);
#else
  Phi.b_e_n = allocZ__MatrixLib__(NumParticles, 9);
  Phi.b_e_n1 = allocZ__MatrixLib__(NumParticles, 9);
#endif

  for (int p = 0; p < NumParticles; p++) {
#if NumberDimensions == 2
  Phi.b_e_n.nM[p][0] = 0.0;
  Phi.b_e_n.nM[p][3] = 0.0;
  Phi.b_e_n.nM[p][4] = 0.0;
#else
  Phi.b_e_n.nM[p][0] = 0.0;
  Phi.b_e_n.nM[p][4] = 0.0;
  Phi.b_e_n.nM[p][9] = 0.0;
#endif
  }


  /*!
   * F-bar variables
   * */
  Phi.Fbar = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);

  for (int p = 0; p < NumParticles; p++) {
    for (int i = 0; i < Ndim; i++) {
      Phi.Fbar.nM[p][i + i * Ndim] = 1.0;
    }
  }

  Phi.Jbar = allocZ__MatrixLib__(NumParticles, 1);

  for (int p = 0; p < NumParticles; p++) {
    Phi.Jbar.nV[p] = 1.0;
  }

  /*!
    Strain_If field (Scalar)
  */
  Phi.Strain_If = allocZ__MatrixLib__(NumParticles, 1);

  /*!
    Stress field (Tensor)
  */
  Phi.Stress =
      allocZ__MatrixLib__(NumParticles, Ndim * Ndim + (Ndim == 2 ? 1 : 0));
  /*!
   * Pore water preassure and its rate
   */
  Phi.Pw = allocZ__MatrixLib__(NumParticles, 1);
  Phi.Pw_0 = allocZ__MatrixLib__(NumParticles, 1);
  Phi.Pw_n1 = allocZ__MatrixLib__(NumParticles, 1);
  Phi.D_Pw = allocZ__MatrixLib__(NumParticles, 1);
  Phi.d_Pw_dt_n = allocZ__MatrixLib__(NumParticles, 1);
  Phi.d_Pw_dt_n1 = allocZ__MatrixLib__(NumParticles, 1);
  Phi.d2_Pw_dt2 = allocZ__MatrixLib__(NumParticles, 1);

  /*!
    Deformation Energy (Scalar)
  */
  Phi.W = allocZ__MatrixLib__(NumParticles, 1);

  /*!
    Mass
  */
  Phi.mass = allocZ__MatrixLib__(NumParticles, 1);

  /*!
    Density of the mixture
  */
  Phi.rho = allocZ__MatrixLib__(NumParticles, 1);

  /*!
    Density of each phase (intrinsic)
  */
  Phi.rho_s = allocZ__MatrixLib__(NumParticles, 1);
  Phi.rho_f = allocZ__MatrixLib__(NumParticles, 1);

  /*!
    Inital volume
  */
  Phi.Vol_0 = allocZ__MatrixLib__(NumParticles, 1);

  /*!
    Relative volume fraction for each phase
  */
  Phi.phi_s = allocZ__MatrixLib__(NumParticles, 1);
  Phi.phi_f = allocZ__MatrixLib__(NumParticles, 1);


  Phi.Chi = (double *)calloc(NumParticles,sizeof(double));

  Phi.EPS_n = (double *)calloc(NumParticles,sizeof(double));

  Phi.EPS_n1 = (double *)calloc(NumParticles,sizeof(double));

  /*!
   * Back stress for kinematic hardening (plasticity)
   */
  Phi.Back_stress = allocZ__MatrixLib__(NumParticles, 3);
  return Phi;
}

/*********************************************************************/

void free_upw_vars__Fields__(Fields Phi) {
  free__MatrixLib__(Phi.rho);
  free__MatrixLib__(Phi.rho_s);
  free__MatrixLib__(Phi.rho_f);
  free__MatrixLib__(Phi.Vol_0);
  free__MatrixLib__(Phi.mass);
  free__MatrixLib__(Phi.x_GC);
  free__MatrixLib__(Phi.x_EC);
  free__MatrixLib__(Phi.dis);
  free__MatrixLib__(Phi.D_dis);
  free__MatrixLib__(Phi.vel);
  free__MatrixLib__(Phi.acc);
  free__MatrixLib__(Phi.Stress);
  free__MatrixLib__(Phi.Pw);
  free__MatrixLib__(Phi.Pw_0);
  free__MatrixLib__(Phi.Pw_n1);
  free__MatrixLib__(Phi.D_Pw);
  free__MatrixLib__(Phi.d_Pw_dt_n);
  free__MatrixLib__(Phi.d_Pw_dt_n1);
  free__MatrixLib__(Phi.d2_Pw_dt2);
  free__MatrixLib__(Phi.Strain);
  free__MatrixLib__(Phi.Strain_If);
  free__MatrixLib__(Phi.b_e_n);
  free__MatrixLib__(Phi.b_e_n1);
  free__MatrixLib__(Phi.F_n);
  free__MatrixLib__(Phi.F_n1);
  free__MatrixLib__(Phi.DF);
  free__MatrixLib__(Phi.dt_F_n);
  free__MatrixLib__(Phi.dt_F_n1);
  free__MatrixLib__(Phi.dt_DF);
  free__MatrixLib__(Phi.J_n);
  free__MatrixLib__(Phi.J_n1);
  free__MatrixLib__(Phi.dJ_dt);
  free__MatrixLib__(Phi.Fbar);
  free__MatrixLib__(Phi.Jbar);
  free__MatrixLib__(Phi.W);
  free__MatrixLib__(Phi.phi_s);
  free__MatrixLib__(Phi.phi_f);
  free(Phi.Chi);
  free(Phi.EPS_n);
  free(Phi.EPS_n1);  
  free__MatrixLib__(Phi.Back_stress);
}

/*********************************************************************/
