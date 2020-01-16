#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*******************************************************/

void UpdateGaussPointStress(GaussPoint MPM_Mesh){

  /* 1º Variable declaration  */
  Matrix Strain_k1;
  Matrix Stress_k0;
  Matrix Stress_k1;
  double mu;
  double E;
  int N_GP = MPM_Mesh.NumGP;
  int Mat_GP;

  /* 2º Iterate over the Gauss-Points */
  for(int i = 0 ; i<N_GP ; i++){
    if(MPM_Mesh.Phi.ji.nV[i] != 1.0){

    Mat_GP = MPM_Mesh.MatIdx[i];
    
    /* 3º Asign materials of the GP */
    mu = MPM_Mesh.Mat[Mat_GP].mu; /* Elastic modulus */
    E =  MPM_Mesh.Mat[Mat_GP].E; /* Poisson ratio */
    
    /* 4º Use pointers to memory manage */
    Strain_k1.nV = MPM_Mesh.Phi.Strain.nM[i];
    Stress_k0.nV = MPM_Mesh.Phi.Stress.nM[i];
    Stress_k1.nV = MPM_Mesh.Phi.Stress.nM[i];
    
    /* 5º Get the new stress tensor (2D Linear elastic) */
    Stress_k1 =
      MPM_Mesh.Mat[Mat_GP].D.LE(Strain_k1,Stress_k0,mu,E);

    /* 6º Get the deformation energy */
    MPM_Mesh.Phi.W.nV[i] =
      W_LinearElastic(Strain_k1,Stress_k1,MPM_Mesh.Phi.ji.nV[i]);
    }
    else{
      for(int j = 0 ; j<MPM_Mesh.Phi.Stress.N_cols ; j++){
	MPM_Mesh.Phi.Stress.nM[i][j] = 0.0;
      }
    }
  }
}

/*******************************************************/
