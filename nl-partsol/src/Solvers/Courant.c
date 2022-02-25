#include "nl-partsol.h"

/*
  Call global variables
*/
Mixture *Soil_Water_Mixtures;
int Number_Soil_Water_Mixtures;

/*********************************************************************/

double U_DeltaT__SolversLib__(Particle MPM_Mesh, double h,
                              Time_Int_Params Parameters_Solver) {

  double DeltaT;
  double CEL_MAX = 0;
  double C[3] = {0, 0, 0};
  int Ndim = NumberDimensions;
  int Nmat = MPM_Mesh.NumberMaterials;
  bool DynamicTimeStep = false;

  /*
    Read paramter solver
  */
  double CFL = Parameters_Solver.CFL;
  double CEL_MAT = Parameters_Solver.Cel;

  /*
    Consider the velocity of the particles for the courant. In
    some cases, for instance Fr>1 is important.
   */
  if (DynamicTimeStep) {
    /*
      Get the maximum wave speed in any direction
    */
    for (int i = 0; i < MPM_Mesh.NumGP; i++) {
      for (int j = 0; j < Ndim; j++) {
        C[j] = DMAX(C[j], CEL_MAT + fabs(MPM_Mesh.Phi.vel.nM[i][j]));
      }
    }

    /*
      Get the absolute maximum value of the celerity
    */
    for (int j = 0; j < Ndim; j++) {
      CEL_MAX = DMAX(CEL_MAX, C[j]);
    }

  } else {
    CEL_MAX = CEL_MAT;
  }

  /*
    Get the minimum value of the time step
  */
  DeltaT = CFL * h / CEL_MAX;

  /*
    Return new time step
  */
  return DeltaT;
}

/*********************************************************************/

double DeltaT_Coussy__SolversLib__(Particle MPM_Mesh, double h, double xi,
                                   double CFL) {

  bool DynamicTimeStep = false;
  double DeltaT;
  double E;
  double rho_s;
  double rho_f;
  double phi_s;
  double phi_f;
  double rho_mixture;
  double CEL_MAT = 0;
  double CEL_MAX = 0;
  double C[3] = {0, 0, 0};
  int Ndim = NumberDimensions;
  int Nmixt = Number_Soil_Water_Mixtures;
  int Np = MPM_Mesh.NumGP;
  int Mixture_idx;
  int Material_Soil_idx;

  /*
    Get the maximum material celerity
  */
  for (int p = 0; p < Np; p++) {
    /*
      Load the elastic modulus os the soil phase
    */
    Mixture_idx = MPM_Mesh.MixtIdx[p];
    Material_Soil_idx = Soil_Water_Mixtures[Mixture_idx].Soil_Idx;
    E = MPM_Mesh.Mat[Material_Soil_idx].E;

    /*
      Get the relative density of each phase
    */
    rho_s = MPM_Mesh.Phi.rho_s.nV[p];
    rho_f = MPM_Mesh.Phi.rho_f.nV[p];
    phi_s = MPM_Mesh.Phi.phi_s.nV[p];
    phi_f = MPM_Mesh.Phi.phi_f.nV[p];

    /*
      Compute the mixture density
    */
    rho_mixture = rho_s * phi_s + xi * rho_f * phi_f;

    /*
      Compute the celerity of the material for the particle p
    */
    CEL_MAT = DMAX(sqrt(E / rho_mixture), CEL_MAT);
  }

  /*
    Consider the velocity of the particles for the courant. In
    some cases, for instance Fr>1 is important.
   */
  if (DynamicTimeStep) {
    /*
      Get the maximum wave speed in any direction
    */
    for (int i = 0; i < MPM_Mesh.NumGP; i++) {
      for (int j = 0; j < Ndim; j++) {
        C[j] = DMAX(C[j], CEL_MAT + fabs(MPM_Mesh.Phi.vel.nM[i][j]));
      }
    }

    /*
      Get the absolute maximum value of the celerity
    */
    for (int j = 0; j < Ndim; j++) {
      CEL_MAX = DMAX(CEL_MAX, C[j]);
    }
  } else {
    CEL_MAX = CEL_MAT;
  }

  /*
    Get the minimum value of the time step
  */
  DeltaT = CFL * h / CEL_MAX;

  /*
    Return new time step
  */
  return DeltaT;
}

/*********************************************************************/
