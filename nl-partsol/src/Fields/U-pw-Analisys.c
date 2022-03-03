#include <string.h>
#include "nl-partsol.h"

/*********************************************************************/

Fields allocate_upw_vars__Fields__(int NumParticles) {
  int Ndim = NumberDimensions;
  Fields Phi;

  /*!
    Global coordinates
  */
  Phi.x_GC = allocZ__MatrixLib__(NumParticles, Ndim);
  strcpy(Phi.x_GC.Info, "Global Coordinates");

  /*!
    Natural coordinates (Vectorial)
  */
  Phi.x_EC = allocZ__MatrixLib__(NumParticles, Ndim);
  strcpy(Phi.x_EC.Info, "Element Coordinates GP");

  /*!
    Displacement and increment of displacement fields (Vectorial)
  */
  Phi.dis = allocZ__MatrixLib__(NumParticles, Ndim);
  strcpy(Phi.dis.Info, "Displacement field GP");
  Phi.D_dis = allocZ__MatrixLib__(NumParticles, Ndim);
  strcpy(Phi.D_dis.Info, "Increment of displacement field GP");

  /*!
    Velocity field (Vectorial)
  */
  Phi.vel = allocZ__MatrixLib__(NumParticles, Ndim);
  strcpy(Phi.vel.Info, "Velocity field GP");

  /*!
    Acceleration field (Vectorial)
  */
  Phi.acc = allocZ__MatrixLib__(NumParticles, Ndim);
  strcpy(Phi.acc.Info, "Acceleration field GP");

  /*!
    Strain field (Tensor)
  */
  Phi.Strain = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
  strcpy(Phi.Strain.Info, "Strain field GP");

  /*!
    Deformation gradient field (Tensor) + Initialise it with the indentity
  */
  Phi.F_n = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
  strcpy(Phi.F_n.Info, "Deformation gradient at t = n");
  Phi.F_n1 = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
  strcpy(Phi.F_n1.Info, "Deformation gradient at t = n + 1");
  Phi.DF = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
  strcpy(Phi.DF.Info, "Increment of deformation gradient");

  for (int p = 0; p < NumParticles; p++) {
    for (int i = 0; i < Ndim; i++) {
      Phi.F_n.nM[p][i + i * Ndim] = 1.0;
      Phi.F_n1.nM[p][i + i * Ndim] = 1.0;
      Phi.DF.nM[p][i + i * Ndim] = 1.0;
    }
  }

  Phi.dt_F_n = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
  strcpy(Phi.dt_F_n.Info, "Rate of deformation gradient at t = n");
  Phi.dt_F_n1 = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
  strcpy(Phi.dt_F_n1.Info, "Rate of deformation gradient at t = n + 1");
  Phi.dt_DF = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);
  strcpy(Phi.dt_DF.Info, "Rate of increment deformation gradient");

  /*!
    Allocate and initialise the Jacobian of the deformation gradient
    ant its rate
  */
  Phi.J_n = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.J_n.Info, "Jacobian of the particle");
  Phi.J_n1 = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.J_n1.Info, "Jacobian of the deformation gradient at t = n + 1");
  Phi.dJ_dt = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.dJ_dt.Info,
         "Rate of the Jacobian of the deformation gradient at t = n + 1");

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
  strcpy(Phi.Fbar.Info, "Fbar deformation gradient");

  for (int p = 0; p < NumParticles; p++) {
    for (int i = 0; i < Ndim; i++) {
      Phi.Fbar.nM[p][i + i * Ndim] = 1.0;
    }
  }

  Phi.Jbar = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.Jbar.Info, "Jacobian of the particle");

  for (int p = 0; p < NumParticles; p++) {
    Phi.Jbar.nV[p] = 1.0;
  }

  /*!
    Strain_If field (Scalar)
  */
  Phi.Strain_If = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.Strain_If.Info, "Strain in fracture GP");

  /*!
    Stress field (Tensor)
  */
  Phi.Stress =
      allocZ__MatrixLib__(NumParticles, Ndim * Ndim + (Ndim == 2 ? 1 : 0));
  strcpy(Phi.Stress.Info, "Stress field GP");

  /*!
   * Pore water preassure and its rate
   */
  Phi.Pw = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.Pw.Info, "Pore water preassure");
  Phi.Pw_0 = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.Pw_0.Info, "Initial pore water preassure");
  Phi.Pw_n1 = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.Pw_n1.Info, "Pore water preassure at t = n + 1");
  Phi.D_Pw = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.D_Pw.Info, "Increment of pore water preassure");
  Phi.d_Pw_dt_n = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.d_Pw_dt_n.Info,
         "First time derivative of pore water preassure at t = n");
  Phi.d_Pw_dt_n1 = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.d_Pw_dt_n1.Info,
         "First time derivative of pore water preassure at t = n + 1");
  Phi.d2_Pw_dt2 = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.d2_Pw_dt2.Info, "Second time derivative of pore water preassure");

  /*!
    Deformation Energy (Scalar)
  */
  Phi.W = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.W.Info, "Deformation Energy GP");

  /*!
    Mass
  */
  Phi.mass = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.mass.Info, "Mass GP");

  /*!
    Density of the mixture
  */
  Phi.rho = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.rho.Info, "Density GP");

  /*!
    Density of each phase (intrinsic)
  */
  Phi.rho_s = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.rho_s.Info, "Intrinsic soil density");
  Phi.rho_f = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.rho_f.Info, "Intrinsic water density");

  /*!
    Inital volume
  */
  Phi.Vol_0 = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.Vol_0.Info, "Inital volume GP");

  /*!
    Relative volume fraction for each phase
  */
  Phi.phi_s = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.phi_s.Info, "Volume fraction soil");
  Phi.phi_f = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.phi_s.Info, "Volume fraction water");

  /*!
    Damage parameter (Fracture)
  */
  Phi.chi = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.chi.Info, "Damage parameter GP");

  /*!
    Equivalent plastic strain (Plasticity)
  */
  Phi.Equiv_Plast_Str = allocZ__MatrixLib__(NumParticles, 1);
  strcpy(Phi.Equiv_Plast_Str.Info, "EPS GP");

  /*!
   * Back stress for kinematic hardening (plasticity)
   */
  Phi.Back_stress = allocZ__MatrixLib__(NumParticles, 3);
  strcpy(Phi.Back_stress.Info, "Back stress GP");

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
  free__MatrixLib__(Phi.chi);
  free__MatrixLib__(Phi.Equiv_Plast_Str);
  free__MatrixLib__(Phi.Back_stress);
}

/*********************************************************************/
