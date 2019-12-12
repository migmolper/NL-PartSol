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
  double Damage;
  int N_GP = MPM_Mesh.NumGP;

  /* 2º Iterate over the Gauss-Points */
  for(int i = 0 ; i<N_GP ; i++){
    
    /* 3º Asign materials of the GP */
    mu = MPM_Mesh.Mat.mu.nV[i]; /* Elastic modulus */
    E =  MPM_Mesh.Mat.E.nV[i]; /* Poisson ratio */
    
    /* 4º Use pointers to memory manage */
    Strain_k1.nV = MPM_Mesh.Phi.Strain.nM[i];
    Stress_k0.nV = MPM_Mesh.Phi.Stress.nM[i];
    Stress_k1.nV = MPM_Mesh.Phi.Stress.nM[i];
    
    /* 5º Get the new stress tensor (2D Linear elastic) */
    Stress_k1 =
      MPM_Mesh.D.LE(Strain_k1,Stress_k0,mu,E);
    
    /* 6º Get the deformation energy */
    Damage = MPM_Mesh.Mat.ji.nV[i];
    MPM_Mesh.Phi.W.nV[i] = W_LinearElastic(Strain_k1,Stress_k1,Damage);
    
  }

  /* /\* 7º Calcule damage parameter *\/ */
  /* MPM_Mesh.Mat.ji = EigenerosionAlgorithm(MPM_Mesh.Mat.ji, MPM_Mesh.Mat.mass, */
  /* 					  MPM_Mesh.Phi.W, MPM_Mesh.Mat.Ceps, */
  /* 					  double G_F, MPM_Mesh.ListNeighbours); */
  
}

/*******************************************************/
