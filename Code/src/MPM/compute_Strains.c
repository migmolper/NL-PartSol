#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

/*************************************************************/

Tensor compute_RateOfStrain(Matrix Velocity, Matrix Gradient)
{
  int Ndim = NumberDimensions;
  Tensor Rate_Strain = alloc_Tensor(2);
  Tensor Velocity_I;
  Tensor Gradient_I;
  Tensor VoG_I;

  int NodesElem = Gradient.N_rows;

  /* Compute rate of strain */
  for(int I = 0 ; I<NodesElem ; I++){
    /* Assign from matrix to tensor */
    Velocity_I = memory_to_Tensor(Velocity.nM[I], 1);
    Gradient_I = memory_to_Tensor(Gradient.nM[I], 1);
   
    /* Compute the dyadic product of the nodal velocity and the
       gradient of the shape functions */
    VoG_I = get_dyadicProduct_Of(Velocity_I, Gradient_I);
    
    /* Ad the nodal contribution to the train tensor */
    for(int i = 0 ; i<Ndim ; i++){
      for(int j = 0 ; j<Ndim ; j++){
	Rate_Strain.N[i][j] +=
	  0.5*(VoG_I.N[i][j] + VoG_I.N[j][i]);
      }
    }
    /* Free memory */
    free_Tensor(VoG_I);
  }
  
  return Rate_Strain;
}

/*******************************************************/

Tensor update_Strain(Tensor Strain, Tensor Rate_Strain, double TimeStep)
{
  int Ndim = NumberDimensions;
  /* Check in the input its is ok */
  if ((Strain.Order == 2) && (Rate_Strain.Order == 2)){
    /* Update strain tensor with the rate of strain tensor */
    for(int i = 0 ; i<Ndim ; i++){
      for(int j = 0 ; j<Ndim ; j++){
	Strain.N[i][j] += TimeStep*Rate_Strain.N[i][j];
      }
    }
  }
  else{
    fprintf(stderr,"%s : %s %s !!! \n",
	    "Error in update_Strain()",
	    "The input should be",
	    "two tensors of 2nd order and a scalar");
    exit(EXIT_FAILURE);    
  }

  return Strain;
}

/*******************************************************/

