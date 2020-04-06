#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"


/*************************************************************/

Tensor compute_Stress(Tensor Strain, Tensor Stress, Material Mat)
{   
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

/*************************************************************/
