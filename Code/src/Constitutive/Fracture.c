#include "grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

void EigenerosionAlgorithm(int p, Matrix ji, Matrix W,  Matrix Mass,
			   Matrix Rho, Matrix Stress,
			   Material MatPro, ChainPtr * Beps, double DeltaX)
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
  Tensor Stress_p, EV_Stress_p; /* Stress tensor */
  double Ceps_p, m_p, rho_p, V_p, sum_p, G_p, Gf_p;
  double m_q, rho_q, V_q, W_q;
  int * Beps_p;
  int NumBeps_p;
  int q;
  
  /* For the current particle get first principal stress */	
  Stress_p = memory_to_Tensor(Stress.nM[p], 2);
  EV_Stress_p = get_Eigenvalues_Of(Stress_p);

  /* Non broken GP Only traction */ 
  if((ji.nV[p] < 1) && (EV_Stress_p.n[0]>0)){ 

    /* Normalizing constant */
    Ceps_p = MatPro.Ceps;
    /* Normalizing constant */
    Gf_p = MatPro.Gf;

    /* Include the main GP in the calculus */
    m_p = Mass.nV[p];
    rho_p = Rho.nV[p];
    V_p = m_p/rho_p;
      
    sum_p = V_p*W.nV[p];

    /* Get a pointer with the list of neighbours */
    NumBeps_p = get_Lenght_Set(Beps[p]);
    Beps_p = Set_to_Pointer(Beps[p],NumBeps_p);

    /* Loop over the neighbours */
    for(int j = 0; j < NumBeps_p ; j++){
      /* Index */
      q = Beps_p[j];

      /* Get volume of particle q */
      m_q = Mass.nV[q];
      rho_q = Rho.nV[q];
      V_q = m_q/rho_q;

      /* Get m_p */
      V_p += V_q;

      if(ji.nV[q] < 1){
	/* Get internal work of particle q */
	W_q = W.nV[q]; 	
	/* Add to sum_p */
	sum_p += V_q*W_q;
      }
	
    }

    /* Free memory */
    free(Beps_p);

    /* Compute energy-release rate for the particle */
    G_p = (Ceps_p*DeltaX/V_p)*sum_p;

    /* Fracture criterium */
    if(G_p > Gf_p){
      ji.nV[p] = 1.0;
    }
      
  }

  /* Free eigenvalues */
  free_Tensor(EV_Stress_p);
    
  
}

/*******************************************************/

void EigensofteningAlgorithm(int p, Matrix ji, Matrix Strain,
			     Matrix StrainF, Matrix Mass,
			     Matrix Stress, Material MatPro,
			     ChainPtr * Beps)
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
  -> StrainF : Value of the strain field at the failure init. 
  -> Beps : Table with the list of neighbours per GP.
  -> Neps : Number of neighbours per GP
  -> Num_GP : Number of GP of the mesh.
*/
{
  /* Define auxiliar variable */

  Tensor Stress_p, Stress_q, EV_Stress_p, EV_Stress_q; /* Stress tensor */
  Tensor Strain_p, EV_Strain_p; /* Stress tensor */
  
  /* Material properties of the eigensoftening algorithm */
  double ft_p, Wc_p, heps_p;
    
  double m_p, sum_p, Seps_p, ji_p;
  double m_q;
  int * Beps_p;
  int NumBeps_p;
  int q;

  /* Get the tensile strengt of the material */
  ft_p = MatPro.ft;
  /* Get the bandwidth of the cohesive fracture (Bazant) */
  heps_p = MatPro.heps;
  /* Get the critical opening displacement */
  Wc_p = MatPro.Wc;
    
  /* Only for intact particles */
  if((ji.nV[p] == 0.0) && (StrainF.nV[p] == 0.0)){

    /* Get the number of neighbours */
    NumBeps_p = get_Lenght_Set(Beps[p]);
      
    if(NumBeps_p > 0){

      /* Get the neighbours */      
      Beps_p = Set_to_Pointer(Beps[p],NumBeps_p);    
      
      /* For the current particle get the mass */
      m_p = Mass.nV[p];

      /* For the current particle get first principal stress */	
      Stress_p = memory_to_Tensor(Stress.nM[p], 2);
      EV_Stress_p = get_Eigenvalues_Of(Stress_p);
      
      /* Add the first term to the sumation */
      sum_p = m_p*EV_Stress_p.n[0];

      /* Free eigenvalues */
      free_Tensor(EV_Stress_p);

      /* Loop over the neighbours */
      for(int j = 0; j < NumBeps_p ; j++){

	/* Get the indedx of each particle q close to p */
	q = Beps_p[j];
  
	if(ji.nV[q] < 1.0){
	  /* For the current particle get the mass and first principal stress */
	  m_q = Mass.nV[q];

	  Stress_q = memory_to_Tensor(Stress.nM[q], 2);
	  EV_Stress_q = get_Eigenvalues_Of(Stress_q);
	    
	  /* Get sum_p */
	  sum_p += m_q*EV_Stress_q.n[0];

	  /* Free eigenvalues */
	  free_Tensor(EV_Stress_q);	    
	}
	/* Add mass contribution */
	m_p += m_q;
      }
    
      /* Free memory */
      free(Beps_p);

      /* Get the equivalent critical stress */
      Seps_p = sum_p/m_p;

      /* Store the principal strain when crack start */
      if(Seps_p>ft_p){

	Strain_p = memory_to_Tensor(Strain.nM[p], 2);
	EV_Strain_p = get_Eigenvalues_Of(Strain_p);
	  	  
	/* Strain during fracture */
	StrainF.nV[p] = EV_Strain_p.n[0];

	/* Free eigenvalues */
	free_Tensor(EV_Strain_p);	    	  
      }
	
    }
      
  }    
  /* Compute the damage parameter if the particle is damaged */
  else if((ji.nV[p] != 1.0) && (StrainF.nV[p] > 0 )){

    Strain_p = memory_to_Tensor(Strain.nM[p], 2);
    EV_Strain_p = get_Eigenvalues_Of(Strain_p);
	  	  
    /* Fracture criterium */
    ji_p = (EV_Strain_p.n[0]-StrainF.nV[p])*heps_p/Wc_p;
    ji.nV[p] = MINVAL(1,MAXVAL(ji_p,ji.nV[p]));       

    /* Free eigenvalues */
    free_Tensor(EV_Strain_p);	    	  
      
  }
 
}

/*******************************************************/

