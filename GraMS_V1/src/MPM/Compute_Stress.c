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

  /* 2º Iterate over the Gauss-Points */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    /* 3º Use pointers to memory manage */
    Strain_k1.nV = MPM_Mesh.Phi.Strain.nM[i];
    Stress_k0.nV = MPM_Mesh.Phi.Stress.nM[i];
    Stress_k1.nV = MPM_Mesh.Phi.Stress.nM[i];
    /* 4º Get the new stress tensor (2D Linear elastic) */
    Stress_k1 =
      MPM_Mesh.D.LE(Strain_k1,Stress_k0,PoissonModulus,ElasticModulus);
    /* 5º Get the deformation energy */
    MPM_Mesh.Phi.W.nV[i] = W_LinearElastic(Strain_k1,Stress_k1);
  }
  
}

/*******************************************************/
