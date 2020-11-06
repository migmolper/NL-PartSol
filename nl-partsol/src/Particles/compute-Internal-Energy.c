#include "nl-partsol.h"

double internal_energy__Particles__(Tensor Strain, Tensor Stress){

  /* Internal energy for the Gauss-Point */
  double W = 0; 
  /*Check in the input its is ok */
  if ((Strain.Order == 2) && (Stress.Order == 2)){    
    /* Calcule the internal work */
    W = 0.5*inner_product__TensorLib__(Strain, Stress);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in internal_energy__Particles__()",
	    "The input should be 2nd tensor and a 2nd tensor");
    exit(EXIT_FAILURE);
  }
  return W;
}

/*********************************************************************/

double finite_strains_internal_energy__Particles__(Tensor F_p, Tensor S_p)
{
  /* Internal energy for the Gauss-Point */
  double W = 0; 
  
  /*Check in the input its is ok */
  if ((F_p.Order == 2) && (S_p.Order == 2))
    {
      /* Strain measure */
      Tensor C_p = right_Cauchy_Green__Particles__(F_p);

      /* Calcule the internal work */
      W = 0.5*inner_product__TensorLib__(S_p, C_p);

      /*
	Free auxiliar variables
      */
      free__TensorLib__(C_p);
    }
  else
    {
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in finite_strains_internal_energy__Particles__()",
	      "The input should be 2nd tensor and a 2nd tensor");
      exit(EXIT_FAILURE);
    }
  
  return W;
}

