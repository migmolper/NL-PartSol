#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

Matrix EigenerosionAlgorithm(Matrix ji, Matrix W, Matrix Mass,
			     int * MatIdx, Material * MatPro,
			     ChainPtr * Beps, double DeltaX)
/*
  A.Pandolfi & M.Ortiz.
  An eigenerosion approach to brittle fracture. 
  International Journal for Numerical Methods in Enginnering.
  92:694-714, 2012.
  NOTE : Here the notation is the same as in the paper.
  Inputs :
  -> Ji_k0 : Matrix with the value of the damage parameter.
  -> Mass : Matrix with the mass of the GP.
  -> W : Incremental free-energy density per unit mass. 
  -> Ceps : Matrix with the normalizing parameter.
  -> Gf : Failure value for the energy-release rate.
  -> Beps : Table with the list of neighbours per GP.
  -> Neps : Number of neighbours per GP
  -> Num_GP : Number of GP of the mesh.
*/
{
  /* Define auxiliar variable */
  int Num_GP = W.N_rows*W.N_cols;
  
  double Ceps_p, m_p, sum_p, G_p, Gf_p;
  double m_q, W_q;
  int * Beps_p;
  int Neps_p;
  int q;
  int Mat_p;
  
  for(int p = 0 ; p < Num_GP ; p++){
    /* Calcule damage if the GP is not broken */
    if(ji.nV[p] < 1){

      /* Kind of material */
      Mat_p = MatIdx[p];
      /* Normalizing constant */
      Ceps_p = MatPro[Mat_p].Ceps;
      /* Normalizing constant */
      Gf_p = MatPro[Mat_p].Gf;

      /* Include the main GP in the calculus */
      m_p = Mass.nV[p];
      sum_p = Mass.nV[p]*W.nV[p];

      /* neighbours */
      Neps_p = LenghtChain(Beps[p]);
      Beps_p = ChainToArray(Beps[p],Neps_p);

      /* Loop over the neighbours */
      for(int j = 0; j < Neps_p ; j++){
	/* Index */
	q = Beps_p[j];
	m_q = Mass.nV[q];
	W_q = W.nV[q];
	/* Get m_p */
	m_p += m_q;
	/* Get sum_p */
	sum_p += m_q*W_q;
      }

      /* Free memory */
      free(Beps_p);

      /* Compute energy-release rate for the GP */
      G_p = (Ceps_p*DeltaX/m_p)*sum_p;

      /* Fracture criterium */
      if(G_p > Gf_p){
	ji.nV[p] = 1.0;
      }
      
    }    
  }
  
  return ji;
}

/*******************************************************/

Matrix ComputeDamage(Matrix ji, Matrix W, Matrix Mass,
		     int * MatIdx, Material * MatProp,
		     ChainPtr * Beps, double DeltaX){

  /* Choose the damage model */
  Matrix Damage_n1 = EigenerosionAlgorithm(ji, W, Mass,
					   MatIdx, MatProp,
					   Beps, DeltaX);

  return Damage_n1;
}

/*******************************************************/
