#include "nl-partsol.h"

/*************************************************************/

Tensor compute_Stress(Tensor Strain, Tensor Stress, Material Mat)
{   
  /* Select the constitutive model */
  if(strcmp(Mat.Type,"SR") == 0){
    Stress = SolidRigid(Stress);
  }  
  else if(strcmp(Mat.Type,"LE") == 0){
    Stress = LinearElastic(Strain,Stress,Mat);
  }
  else{
    exit(0);
  }
  
  /* Return the stress tensor */
  return Stress;
}

/*************************************************************/
