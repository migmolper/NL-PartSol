#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../GRAMS/grams.h"

Matrix LinearElastic(Matrix Strain_n1,
		     Matrix Stress_n,
		     Material Material_GP){

  double mu = Material_GP.mu; /* Poisson ratio */
  double E = Material_GP.E; /* Elastic modulus */ 
  Matrix Stress_n1;  
  double LameFirstParam = /* Lambda */
    mu*E/((1-mu*2)*(1+mu));
  double LameSecondParam = /* G */
    E/(2*(1+mu));

  /* Assign memory position */
  Stress_n1.nV = Stress_n.nV;

  switch(NumberDimensions){
  case 1 :
    puts("Not implemented");
    break;
  case 2 :
    Stress_n1.nV[0] = (LameFirstParam + 2*LameSecondParam)*Strain_n1.nV[0] +
      LameFirstParam*Strain_n1.nV[1];
    Stress_n1.nV[1] = (LameFirstParam + 2*LameSecondParam)*Strain_n1.nV[1] +
      LameFirstParam*Strain_n1.nV[0];
    Stress_n1.nV[2] = 2*LameSecondParam*Strain_n1.nV[2];
    break;
  case 3 :
    puts("Not implemented");
    break;    
  default :
    printf("%s : %s \n",
	   "Error in LinearElastic",
	   "Incorrect number of dimensions");
    exit(0);
  }
    
  return Stress_n1; 
}

double W_LinearElastic(Matrix Strain,
		       Matrix Stress,
		       double Damage){

  double W = 0; /* Internal energy for the GP */
    
  /* Calcule the internal work */
  switch(NumberDimensions){
  case 1 :
    if(Damage < 1){
      W = 0.5*Strain.n*Stress.n;
    }
    return W;
  case 2 :
    if(Damage < 1){
      for(int i = 0 ; i<3 ; i++){
	W += Strain.nV[i]*Stress.nV[i];
      }
      W *= 0.5;
    }
    return W;
  case 3 :
    if(Damage < 1){
      for(int i = 0 ; i<6 ; i++){
	W += Strain.nV[i]*Stress.nV[i];
      }
      W *= 0.5;
    }
    return W; 
  default :
    printf("%s : %s \n",
	   "Error in W_LinearElastic",
	   "Incorrect number of dimensions");
    exit(0);
  }
}
