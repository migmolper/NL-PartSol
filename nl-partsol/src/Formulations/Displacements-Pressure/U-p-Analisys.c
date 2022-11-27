#include "Formulations/Displacements-Pressure/U-p-Analisys.h"

/*********************************************************************/

Fields allocate_Up_vars__Fields__(int NumParticles) {
  int Ndim = NumberDimensions;
  Fields Phi;

  Phi.x_GC = allocZ__MatrixLib__(NumParticles, Ndim);

  Phi.x_EC = allocZ__MatrixLib__(NumParticles, Ndim);

  Phi.dis = allocZ__MatrixLib__(NumParticles, Ndim);

  Phi.D_dis = allocZ__MatrixLib__(NumParticles, Ndim);

  Phi.vel = allocZ__MatrixLib__(NumParticles, Ndim);
  
  Phi.acc = allocZ__MatrixLib__(NumParticles, Ndim);

  Phi.Strain = allocZ__MatrixLib__(NumParticles, Ndim * Ndim);

#if NumberDimensions == 2
  Phi.F_n = allocZ__MatrixLib__(NumParticles, 5);
  Phi.F_n1 = allocZ__MatrixLib__(NumParticles, 5);
  Phi.DF = allocZ__MatrixLib__(NumParticles, 5);
#else
  Phi.F_n = allocZ__MatrixLib__(NumParticles, 9);
  Phi.F_n1 = allocZ__MatrixLib__(NumParticles, 9);
  Phi.DF = allocZ__MatrixLib__(NumParticles, 9);
#endif

  for (int p = 0; p < NumParticles; p++) {
#if NumberDimensions == 2
  Phi.F_n.nM[p][0] = Phi.F_n1.nM[p][0] = Phi.DF.nM[p][0] = 1.0;
  Phi.F_n.nM[p][3] = Phi.F_n1.nM[p][3] = Phi.DF.nM[p][3] = 1.0;
  Phi.F_n.nM[p][4] = Phi.F_n1.nM[p][4] = Phi.DF.nM[p][4] = 1.0;
#else
  Phi.F_n.nM[p][0] = Phi.F_n1.nM[p][0] = Phi.DF.nM[p][0] = 1.0;
  Phi.F_n.nM[p][4] = Phi.F_n1.nM[p][4] = Phi.DF.nM[p][4] = 1.0;
  Phi.F_n.nM[p][9] = Phi.F_n1.nM[p][9] = Phi.DF.nM[p][9] = 1.0;
#endif
  }

#if NumberDimensions == 2
  Phi.dt_F_n = allocZ__MatrixLib__(NumParticles, 5);
  Phi.dt_F_n1 = allocZ__MatrixLib__(NumParticles, 5);
  Phi.dt_DF = allocZ__MatrixLib__(NumParticles, 5);
#else
  Phi.dt_F_n = allocZ__MatrixLib__(NumParticles, 9);
  Phi.dt_F_n1 = allocZ__MatrixLib__(NumParticles, 9);
  Phi.dt_DF = allocZ__MatrixLib__(NumParticles, 9);
#endif
  
#if NumberDimensions == 2
  Phi.b_e_n = allocZ__MatrixLib__(NumParticles, 5);
  Phi.b_e_n1 = allocZ__MatrixLib__(NumParticles, 5);
#else
  Phi.b_e_n = allocZ__MatrixLib__(NumParticles, 9);
  Phi.b_e_n1 = allocZ__MatrixLib__(NumParticles, 9);
#endif

  for (int p = 0; p < NumParticles; p++) {
#if NumberDimensions == 2
  Phi.b_e_n.nM[p][0] = 1.0;
  Phi.b_e_n.nM[p][3] = 1.0;
  Phi.b_e_n.nM[p][4] = 1.0;
#else
  Phi.b_e_n.nM[p][0] = 1.0;
  Phi.b_e_n.nM[p][4] = 1.0;
  Phi.b_e_n.nM[p][9] = 1.0;
#endif
  }

  Phi.J_n = (double *)calloc(NumParticles,sizeof(double));
  
  Phi.J_n1 = (double *)calloc(NumParticles,sizeof(double));
  
  for (int p = 0; p < NumParticles; p++) {
    Phi.J_n[p] = 1.0;
    Phi.J_n1[p] = 1.0;
  }

#if NumberDimensions == 2
  Phi.Fbar = allocZ__MatrixLib__(NumParticles, 5);
#else
  Phi.Fbar = allocZ__MatrixLib__(NumParticles, 9);
#endif

  for (int p = 0; p < NumParticles; p++) {
#if NumberDimensions == 2
  Phi.Fbar.nM[p][0] = 0.0;
  Phi.Fbar.nM[p][3] = 0.0;
  Phi.Fbar.nM[p][4] = 0.0;
#else
  Phi.Fbar.nM[p][0] = 0.0;
  Phi.Fbar.nM[p][4] = 0.0;
  Phi.Fbar.nM[p][9] = 0.0;
#endif
  }

  Phi.Jbar = allocZ__MatrixLib__(NumParticles, 1);

  for (int p = 0; p < NumParticles; p++) {
    Phi.Jbar.nV[p] = 1.0;
  }

#if NumberDimensions == 2
  Phi.Stress = allocZ__MatrixLib__(NumParticles, 5);
#else
  Phi.Stress = allocZ__MatrixLib__(NumParticles, 9);
#endif

  Phi.lambda_pressure_n = allocZ__MatrixLib__(NumParticles, 1);
  
  Phi.lambda_pressure_n1 = allocZ__MatrixLib__(NumParticles, 1);
  
  Phi.W = (double *)calloc(NumParticles,sizeof(double));
  
  Phi.mass = allocZ__MatrixLib__(NumParticles, 1);

  Phi.rho = allocZ__MatrixLib__(NumParticles, 1);
  
  Phi.Vol_0 = allocZ__MatrixLib__(NumParticles, 1);
  
if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
  
  Phi.Damage_n = (double *)calloc(NumParticles,sizeof(double));

  Phi.Damage_n1 = (double *)calloc(NumParticles,sizeof(double));
}

if (Driver_EigenSoftening == true) {
  Phi.Strain_f_n = (double *)calloc(NumParticles,sizeof(double));

  Phi.Strain_f_n1 = (double *)calloc(NumParticles,sizeof(double));
}


  Phi.EPS_n = (double *)calloc(NumParticles,sizeof(double));

  Phi.EPS_n1 = (double *)calloc(NumParticles,sizeof(double));

  Phi.Kappa_n = (double *)calloc(NumParticles,sizeof(double));

  Phi.Kappa_n1 = (double *)calloc(NumParticles,sizeof(double));
  
  Phi.Back_stress = allocZ__MatrixLib__(NumParticles, 3);
  
  return Phi;
}

/*********************************************************************/

void free_Up_vars__Fields__(Fields Phi) {
  free__MatrixLib__(Phi.rho);
  free__MatrixLib__(Phi.mass);
  free__MatrixLib__(Phi.x_GC);
  free__MatrixLib__(Phi.x_EC);
  free__MatrixLib__(Phi.dis);
  free__MatrixLib__(Phi.D_dis);
  free__MatrixLib__(Phi.vel);
  free__MatrixLib__(Phi.acc);
  free__MatrixLib__(Phi.Stress);
  free__MatrixLib__(Phi.lambda_pressure_n);
  free__MatrixLib__(Phi.lambda_pressure_n1);
  free__MatrixLib__(Phi.Strain);
  free__MatrixLib__(Phi.b_e_n);
  free__MatrixLib__(Phi.b_e_n1);
  free__MatrixLib__(Phi.F_n);
  free__MatrixLib__(Phi.F_n1);
  free(Phi.J_n);
  free(Phi.J_n1);
  free__MatrixLib__(Phi.dt_F_n);
  free__MatrixLib__(Phi.dt_F_n1);
  free__MatrixLib__(Phi.dt_DF);
  free__MatrixLib__(Phi.Fbar);
  free__MatrixLib__(Phi.Jbar);
  free__MatrixLib__(Phi.Jbar);
  free__MatrixLib__(Phi.DF);
  free(Phi.W);
  free__MatrixLib__(Phi.Vol_0);
if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
    free(Phi.Damage_n);
  free(Phi.Damage_n1);
}
if (Driver_EigenSoftening == true) {
  free(Phi.Strain_f_n);
  free(Phi.Strain_f_n1);
}
  free(Phi.EPS_n);
  free(Phi.EPS_n1);  
  free(Phi.Kappa_n);
  free(Phi.Kappa_n1);
  free__MatrixLib__(Phi.Back_stress);
}

/*********************************************************************/