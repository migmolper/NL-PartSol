#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

double compute_InternalEnergy(Tensor Strain, Tensor Stress){

  /* Internal energy for the Gauss-Point */
  double W = 0; 
  /*Check in the input its is ok */
  if ((Strain.Order == 2) && (Stress.Order == 2)){    
    /* Calcule the internal work */
    W = 0.5*get_innerProduct_Of(Strain, Stress);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in compute_InternalEnergy()",
	    "The input should be 2nd tensor and a 2nd tensor");
    exit(EXIT_FAILURE);
  }
  return W;
}

/* double W_LinearElastic(Matrix Strain, */
/* 		       Matrix Stress, */
/* 		       double Damage){ */

/*   double W = 0; /\* Internal energy for the GP *\/ */
    
/*   /\* Calcule the internal work *\/ */
/*   switch(NumberDimensions){ */
/*   case 1 : */
/*     if(Damage < 1){ */
/*       W = 0.5*Strain.n*Stress.n; */
/*     } */
/*     return W; */
/*   case 2 : */
/*     if(Damage < 1){ */
/*       for(int i = 0 ; i<3 ; i++){ */
/* 	W += Strain.nV[i]*Stress.nV[i]; */
/*       } */
/*       W *= 0.5; */
/*     } */
/*     return W; */
/*   case 3 : */
/*     if(Damage < 1){ */
/*       for(int i = 0 ; i<6 ; i++){ */
/* 	W += Strain.nV[i]*Stress.nV[i]; */
/*       } */
/*       W *= 0.5; */
/*     } */
/*     return W;  */
/*   default : */
/*     printf("%s : %s \n", */
/* 	   "Error in W_LinearElastic", */
/* 	   "Incorrect number of dimensions"); */
/*     exit(0); */
/*   } */
/* } */

