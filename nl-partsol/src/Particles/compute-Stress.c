#include "nl-partsol.h"

/*************************************************************/

Tensor explicit_integration_stress__Particles__(Tensor Strain,
						Tensor Stress,
						Material Mat)
{   
  /*
    Select the constitutive model 
  */
  if(strcmp(Mat.Type,"SR") == 0){
    Stress = SolidRigid(Stress);
  }  
  else if(strcmp(Mat.Type,"LE") == 0){
    Stress = LinearElastic(Strain,Stress,Mat);
  }
  else{
    exit(EXIT_FAILURE);
  }
  
  /* Return the stress tensor */
  return Stress;
}

/**************************************************************/

Tensor forward_integration_Stress__Particles__(Tensor S_p,Tensor F_n1_p,Material MatProp_p)
{


  /*
    Compute the right Cauchy Green tensor
   */
  Tensor C_n1_p = right_Cauchy_Green__Particles__(F_n1_p);

  /*
    Compute the Stress tensor
   */
  if(strcmp(MatProp_p.Type,"Saint-Venant-Kirchhoff") == 0)
    {
      S_p = grad_energy_Saint_Venant_Kirchhoff(S_p, C_n1_p, MatProp_p);
    }
  else if(strcmp(MatProp_p.Type,"Neo-Hookean-Wriggers") == 0)
    {
      double J_n1_p = I3__TensorLib__(F_n1_p);
      S_p = grad_energy_Neo_Hookean_Wriggers(S_p, C_n1_p, J_n1_p, MatProp_p);
    }
  else
    {
      fprintf(stderr,"%s : %s %s %s \n",
	      "Error in forward_integration_Stress__Particles__()",
	      "The material",MatProp_p.Type,"has not been yet implemnented");
      exit(EXIT_FAILURE);
    }

  /*
    Free auxiliar variables
   */
  free__TensorLib__(C_n1_p);

  return S_p; 
  
}

/**************************************************************/

Tensor configurational_midpoint_integration_Stress__Particles__(Tensor S_p,
								Tensor F_n1_p,
								Tensor F_n_p,
								Material MatProp_p)
{

  /*
    Compute the midpoint deformation gradient 
   */
  Tensor F_n12_p = Convex_combination__TensorLib__(F_n1_p,F_n_p,0.5);

  /*
    Compute the right Cauchy Green tensor
   */
  Tensor C_n12_p = right_Cauchy_Green__Particles__(F_n12_p);

  /*
    Compute the Stress tensor
   */
  if(strcmp(MatProp_p.Type,"Saint-Venant-Kirchhoff") == 0)
    {
      S_p = grad_energy_Saint_Venant_Kirchhoff(S_p, C_n12_p, MatProp_p);
    }
  else if(strcmp(MatProp_p.Type,"Neo-Hookean-Wriggers") == 0)
    {
      double J_n12_p = I3__TensorLib__(F_n12_p);
      S_p = grad_energy_Neo_Hookean_Wriggers(S_p, C_n12_p, J_n12_p, MatProp_p);
    }
  else
    {
      fprintf(stderr,"%s : %s %s %s \n",
	      "Error in configurational_midpoint_integration_Stress__Particles__()",
	      "The material",MatProp_p.Type,"has not been yet implemnented");
      exit(EXIT_FAILURE);
    }

  /*
    Free auxiliar variables
   */
  free__TensorLib__(F_n12_p);
  free__TensorLib__(C_n12_p);

  return S_p;  
  
}

/**************************************************************/

Tensor average_strain_integration_Stress__Particles__(Tensor S_p,
						      Tensor F_n1_p,
						      Tensor F_n_p,
						      Material MatProp_p)
{
  
  /*
    Compute the right Cauchy-Green tensor in diferent time steps
  */
  Tensor C_n_p  = right_Cauchy_Green__Particles__(F_n_p);
  Tensor C_n1_p = right_Cauchy_Green__Particles__(F_n1_p);
  Tensor C_n12_p = Convex_combination__TensorLib__(C_n1_p,C_n_p,0.5);
  
  /*
    Compute the Stress tensor
  */
  if(strcmp(MatProp_p.Type,"Saint-Venant-Kirchhoff") == 0)
    {
      S_p = grad_energy_Saint_Venant_Kirchhoff(S_p, C_n12_p, MatProp_p);
    }
  else if(strcmp(MatProp_p.Type,"Neo-Hookean-Wriggers") == 0)
    {
      double J_n12_p = 0.5*(I3__TensorLib__(F_n_p) + I3__TensorLib__(F_n1_p));
      S_p = grad_energy_Neo_Hookean_Wriggers(S_p, C_n12_p, J_n12_p, MatProp_p);
    }
  else
    {
      fprintf(stderr,"%s : %s %s %s \n",
	      "Error in average_strain_integration_Stress__Particles__()",
	      "The material",MatProp_p.Type,"has not been yet implemnented");
      exit(EXIT_FAILURE);
    }
  
  /*
    Free auxiliar variables
   */
  free__TensorLib__(C_n_p);
  free__TensorLib__(C_n1_p);
  free__TensorLib__(C_n12_p);

  return S_p;
}

/**************************************************************/

Tensor average_itegration_Stress__Particles__(Tensor S_p,
					      Tensor F_n1_p,
					      Tensor F_n_p,
					      Material MatProp_p)
{

  int Ndim = NumberDimensions;
  Tensor C_n_p;
  Tensor C_n1_p;
  Tensor S_n_p;
  Tensor S_n1_p;
  Tensor S_n12_p;
  
  /*
    Compute the right Cauchy-Green tensor
  */
  C_n_p  = right_Cauchy_Green__Particles__(F_n_p);
  C_n1_p = right_Cauchy_Green__Particles__(F_n1_p);

  /*
    Compute the Stress tensor
  */
  S_n_p = alloc__TensorLib__(2);
  S_n1_p = alloc__TensorLib__(2);
  
  if(strcmp(MatProp_p.Type,"Saint-Venant-Kirchhoff") == 0)
    {
      S_n_p = grad_energy_Saint_Venant_Kirchhoff(S_n_p, C_n_p, MatProp_p);
      S_n1_p = grad_energy_Saint_Venant_Kirchhoff(S_n1_p, C_n1_p, MatProp_p);
      S_n12_p = Convex_combination__TensorLib__(S_n1_p,S_n_p,0.5);
    }
  else if(strcmp(MatProp_p.Type,"Neo-Hookean-Wriggers") == 0)
    {
      double J_n_p = I3__TensorLib__(F_n_p);
      double J_n1_p = I3__TensorLib__(F_n1_p);
      S_n_p  = grad_energy_Neo_Hookean_Wriggers(S_n_p, C_n_p, J_n_p, MatProp_p);
      S_n1_p = grad_energy_Neo_Hookean_Wriggers(S_n1_p, C_n1_p, J_n1_p, MatProp_p);
      S_n_p = Convex_combination__TensorLib__(S_n1_p,S_n_p,0.5);
    }
  else
    {
      fprintf(stderr,"%s : %s %s %s \n",
	      "Error in average_strain_integration_Stress__Particles__()",
	      "The material",MatProp_p.Type,"has not been yet implemnented");
      exit(EXIT_FAILURE);
    }
  
  /*
    Compute the Stress tensor
  */

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
	{
	  S_p.N[i][j] = S_n12_p.N[i][j];	  
	}
    }
  
  /*
    Free auxiliar variables
  */
  free__TensorLib__(C_n_p);
  free__TensorLib__(C_n1_p);
  free__TensorLib__(S_n_p);
  free__TensorLib__(S_n1_p);
  free__TensorLib__(S_n12_p);

  return S_p;  
}

/**************************************************************/


