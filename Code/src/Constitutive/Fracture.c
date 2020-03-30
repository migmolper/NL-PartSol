#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

void EigenerosionAlgorithm(Matrix ji, Matrix W,
			   Matrix Mass, Matrix Stress,
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
    if((ji.nV[p] < 1) && /* Non broken GP */
       (Stress.nM[p][0]>0) && (Stress.nM[p][1]>0)){ /* Only traction */

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
  
}

/*******************************************************/

void EigensofteningAlgorithm(Matrix ji, Matrix RateStrain,
			     Matrix Mass, Matrix Stress,
			     int * MatIdx, Material * MatPro,
			     ChainPtr * Beps, double DeltaT)
/*
  Pedro Navas, Rena C. Yu, Bo Li & Gonzalo Ruiz.
  Modeling the dynamic fracture in concrete: 
  an eigensoftening meshfree approach.
  International Journal of Impact Engineering.
  113 (2018) 9-20
  NOTE : Here notation is the same as in the paper.
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
  int Num_GP = RateStrain.N_rows;
  
  
  /* Material properties of the eigensoftening algorithm */
  double ft_p, Wc_p, heps_p;

  /* Invariants for stress and rate of strain */
  double Stress_I, Stress_II, RateStrain_I, RateStrain_II;
    
  double m_p, Stress_1p, sum_p, Seps_p, RateStrain_1p;
  double m_q, Stress_1q;
  int * Beps_p;
  int Neps_p;
  int q;
  int Mat_p;
  
  for(int p = 0 ; p < Num_GP ; p++){

    /* Kind of material */
    Mat_p = MatIdx[p];
    /* Get the tensile strengt of the material */
    ft_p = MatPro[Mat_p].ft;
    /* Get the bandwidth of the cohesive fracture (Bazant) */
    heps_p = MatPro[Mat_p].heps;
    /* Get the critical opening displacement */
    Wc_p = MatPro[Mat_p].Wc;

    /* Only for intact particles */
    if(ji.nV[p] == 0){
          
      /* For the current particle get the mass and first principal stress */
      m_p = Mass.nV[p];
      Stress_I = Stress.nM[p][0]+Stress.nM[p][1];
      Stress_II = Stress.nM[p][2]*Stress.nM[p][2];
      Stress_1p = 0.5*(- Stress_I + sqrt(Stress_I*Stress_I +
					 4*Stress_II*Stress_II));

      /* Add the first term to the sumation */
      sum_p = m_p*Stress_1p;

      /* Get the neighbours */
      Neps_p = LenghtChain(Beps[p]);
      Beps_p = ChainToArray(Beps[p],Neps_p);    

      /* Loop over the neighbours */
      for(int j = 0; j < Neps_p ; j++){

	/* Get the indedx of each particle q close to p */
	q = Beps_p[j];

	/* For the current particle get the mass and first principal stress */
	m_q = Mass.nV[q];      
	Stress_I = Stress.nM[q][0]+Stress.nM[q][1];
	Stress_II = Stress.nM[q][2]*Stress.nM[q][2];
	Stress_1q = 0.5*(- Stress_I + sqrt(Stress_I*Stress_I +
					   4*Stress_II*Stress_II));
      
	/* Get sum_p */
	sum_p += m_q*Stress_1q;
	m_p += m_q;
      }
    
      /* Free memory */
      free(Beps_p);

      /* Get the equivalent critical stress */
      Seps_p = sum_p/m_p;
    }
    
    /* Compute the damage parameter if the particle is damaged */
    if((ji.nV[p] > 0.0) || (Seps_p > ft_p )){ 

      /* Get the principal rate of strain */
      RateStrain_I = RateStrain.nM[p][0]+RateStrain.nM[p][1];
      RateStrain_II = RateStrain.nM[p][2]*RateStrain.nM[p][2];
      RateStrain_1p = 0.5*(- RateStrain_I + sqrt(RateStrain_I*RateStrain_I +
						 RateStrain_II*RateStrain_II));
            
      /* Fracture criterium */
      ji.nV[p] = RateStrain_1p*DeltaT*heps_p/Wc_p;
      
    }    
  }
 
}

/*******************************************************/

