#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grams.h"


Tensor LinearElastic(Tensor Strain, Tensor Stress, Material Mat){

  /*Check in the input its is ok */
  if ((Strain.Order == 2) && (Stress.Order == 2)){
    /* Define material and other properties */
    double mu = Mat.mu; 
    double E = Mat.E;
    double lambda = mu*E/((1-mu*2)*(1+mu));
    double G = E/(2*(1+mu));
    double traceStrain = get_I1_Of(Strain);
    Tensor I = get_I();

    for(int i = 0 ; i<3 ; i++){
      for(int j = 0 ; j<3 ; j++){
	Stress.N[i][j] = lambda*traceStrain*I.N[i][j] + 2*G*Strain.N[i][j];
      }
    }

    free_Tensor(I);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in LinearElastic()",
	    "The input should be 2nd tensor and a 2nd tensor");
    exit(EXIT_FAILURE);
  }
  return Stress;  
}

/* Matrix LinearElastic(Matrix Strain_n1, */
/* 		     Matrix Stress_n, */
/* 		     Material Material_GP){ */

/*   double mu = Material_GP.mu; /\* Poisson ratio *\/ */
/*   double E = Material_GP.E; /\* Elastic modulus *\/  */
/*   Matrix Stress_n1;   */
/*   double LameFirstParam = /\* Lambda *\/ */
/*     mu*E/((1-mu*2)*(1+mu)); */
/*   double LameSecondParam = /\* G *\/ */
/*     E/(2*(1+mu)); */

/*   /\* Assign memory position *\/ */
/*   Stress_n1.nV = Stress_n.nV; */

/*   switch(NumberDimensions){ */
/*   case 1 : */
/*     puts("Not implemented"); */
/*     break; */
/*   case 2 : */
/*     Stress_n1.nV[0] = (LameFirstParam + 2*LameSecondParam)*Strain_n1.nV[0] + */
/*       LameFirstParam*Strain_n1.nV[1]; */
/*     Stress_n1.nV[1] = (LameFirstParam + 2*LameSecondParam)*Strain_n1.nV[1] + */
/*       LameFirstParam*Strain_n1.nV[0]; */
/*     Stress_n1.nV[2] = 2*LameSecondParam*Strain_n1.nV[2]; */
/*     break; */
/*   case 3 : */
/*     puts("Not implemented"); */
/*     break;     */
/*   default : */
/*     printf("%s : %s \n", */
/* 	   "Error in LinearElastic", */
/* 	   "Incorrect number of dimensions"); */
/*     exit(0); */
/*   } */
    
/*   return Stress_n1;  */
/* } */

