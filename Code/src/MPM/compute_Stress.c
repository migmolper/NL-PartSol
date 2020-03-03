#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"


/*************************************************************/

Tensor compute_Stress(Tensor Strain, Tensor Stress, Material Mat)
{
  /* Variable definition  */
  Tensor Strain_n1;
    
  /* Select the constitutive model */
  if(strcmp(Mat.Type,"LE") == 0){
    Stress = LinearElastic(Strain,Stress,Mat);
  }
  else{
    exit(0);
  }
  
  /* Return the stress tensor */
  return Stress;
}

/* /\*************************************************************\/ */

/* void UpdateGaussPointStress(GaussPoint MPM_Mesh) */
/* /\*! */
/*  * \brief Brief description of UpdateGaussPointStress. */
/*  *        Update the stress state of the body evaluating the */
/*  *        strains in each Gauss-Point to get the stress state. */
/*  * */
/*  *  The parameters for this functions are  : */
/*  *  @param MPM_Mesh : Mesh with the material points. */
/*  * */
/*  *\/ */
/* { */
/*   /\* Variable definition  *\/ */
/*   Matrix Strain_n1; /\* Strain tensor of the GP in the next step *\/ */
/*   Matrix Stress_n; /\* Stress tensor of the GP in the previous step *\/ */
/*   Matrix Stress_n1; /\* Stress tensor of the GP in the next step *\/ */
/*   Material Material_GP; /\* Material properties of the GP *\/ */
/*   int N_GP = MPM_Mesh.NumGP; */
/*   int Mat_GP; */

/*   /\* 2º Iterate over the Gauss-Points *\/ */
/*   for(int i = 0 ; i<N_GP ; i++){ */
/*     if(MPM_Mesh.Phi.ji.nV[i] != 1.0){ */
      
/*       /\* 3º Asign materials of the GP *\/ */
/*       Mat_GP = MPM_Mesh.MatIdx[i]; */
/*       Material_GP = MPM_Mesh.Mat[Mat_GP];  */
      
/*       /\* 4º Use pointers to memory manage *\/ */
/*       Strain_n1.nV = MPM_Mesh.Phi.Strain.nM[i]; */
/*       Stress_n.nV = MPM_Mesh.Phi.Stress.nM[i]; */
/*       Stress_n1.nV = MPM_Mesh.Phi.Stress.nM[i]; */
    
/*       /\* 5º Get the new stress tensor *\/ */
/*       if(strcmp(Material_GP.Type,"LE") == 0){ */
/* 	Stress_n1 = LinearElastic(Strain_n1,Stress_n,Material_GP); */
/*       } */
/*       else{ */
/* 	exit(0); */
/*       } */
      
/*       /\* 6º Get the deformation energy *\/ */
/*       MPM_Mesh.Phi.W.nV[i] = */
/* 	W_LinearElastic(Strain_n1,Stress_n1,MPM_Mesh.Phi.ji.nV[i]); */
/*     } */
/*     else{ */
/*       for(int j = 0 ; j<MPM_Mesh.Phi.Stress.N_cols ; j++){ */
/* 	MPM_Mesh.Phi.Stress.nM[i][j] = 0.0; */
/*       } */
/*     } */
/*   } */
/* } */

/* /\*******************************************************\/ */

