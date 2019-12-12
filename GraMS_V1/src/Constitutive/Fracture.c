#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

Matrix EigenerosionAlgorithm(Matrix Ji_k0, Matrix Mass,
			     Matrix W, Matrix Ceps,
			     double G_F, ChainPtr * Beps)
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
  -> G_F : Failure value for the energy-release rate.
  -> Beps : Table with the list of neighbours per GP.
  -> Neps : Number of neighbours per GP
  -> Num_GP : Number of GP of the mesh.
*/
{
  /* Define auxiliar variable */
  int Num_GP = W.N_rows*W.N_cols;
  Matrix Ji_k1 = MatAssign(Num_GP,1,NAN,Ji_k0.nV,NULL);
  double Ceps_p, m_p, sum_p, G_p;
  double m_q, W_q;
  int * Beps_p;
  int Neps_p;
  int q;
  
  for(int p = 0 ; p < Num_GP ; p++){
    /* Calcule damage if the GP is not broken */
    if(Ji_k0.nV[p] < 1){

      /* Include the main GP in the calculus */
      m_p = Mass.nV[p];
      sum_p = Mass.nV[p]*W.nV[p];

      /* neighbours */
      Neps_p = Neps[p];
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

      /* Normalizing constant */
      Ceps_p = Ceps.nV[p];

      /* Compute energy-release rate for the GP */
      G_p = (Ceps_p/m_p)*sum_p;

      /* Fracture criterium */
      if(G_p > G_F){
	Ji_k1.nV[p] = 1.0;
      }
      
    }    
  }
  
  return Ji_k1;
}
