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

double finite_strains_internal_energy__Particles__(Tensor F_p, Material MatProp_p, double V0_p)
{
  /* Internal energy for the Gauss-Point */
  double W = 0; 

  double J = I3__TensorLib__(F_p);

  /* Strain measure */
  Tensor C_p = right_Cauchy_Green__Particles__(F_p);

  if(strcmp(MatProp_p.Type,"Saint-Venant-Kirchhoff") == 0)
    {
      W = V0_p*energy_Saint_Venant_Kirchhoff(C_p, MatProp_p);
    }
  else if(strcmp(MatProp_p.Type,"Neo-Hookean-Wriggers") == 0)
    {
      W = V0_p*energy_Neo_Hookean_Wriggers(C_p, J, MatProp_p);
    }     
  else
    {
      fprintf(stderr,"%s : %s %s %s \n",
        "Error in finite_strains_internal_energy__Particles__()",
        "The material",MatProp_p.Type,"has not been yet implemnented");
      exit(EXIT_FAILURE);
    }
    /*
	     Free auxiliar variables
    */
    free__TensorLib__(C_p);
  
  return W;
}

/*********************************************************************/

