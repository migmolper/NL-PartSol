#include "nl-partsol.h"

/**************************************************************/

Tensor Itegration_Stress_Configurational_Midpoint(Tensor S_p,
						  Tensor F_n1_p, Tensor F_n_p,
						  Material MatProp_p)
{
  /* Auxiliar variables */
  Tensor F_n12_p;
  Tensor C_n12_p;
  Tensor E_n12_p;

  /*
    Compute the midpoint deformation gradient 
   */
  F_n12_p = compute_midpoint_Tensor(F_n1_p,F_n_p,0.5);

  /*
    Compute the right Cauchy Green tensor
   */
  C_n12_p = compute_RightCauchyGreen(F_n_p);

  /*
    Compute the lagrangian strain tensor
   */
  E_n12_p =  compute_LagrangianStrain(C_n12_p);

  /*
    Compute the Stress tensor
   */
  S_p = grad_energy_Saint_Venant_Kirchhoff(S_p, E_n12_p, MatProp_p);

  /*
    Free auxiliar variables
   */
  free_Tensor(F_n12_p);
  free_Tensor(C_n12_p);
  free_Tensor(E_n12_p);

  return S_p;  
  
}

/**************************************************************/

Tensor Itegration_Stress_Average_Strain(Tensor S_p,
					Tensor F_n1_p, Tensor F_n_p,
					Material MatProp_p)
{

  /* Auxiliar variables */
  Tensor C_n_p;
  Tensor C_n1_p;
  Tensor C_n12_p;
  Tensor E_n12_p;
  
  /*
    Compute the right Cauchy-Green tensor
  */
  C_n_p  = compute_RightCauchyGreen(F_n_p);
  C_n1_p = compute_RightCauchyGreen(F_n1_p);
  C_n12_p = compute_midpoint_Tensor(C_n1_p,C_n_p,0.5);
  
  /*
    Compute the lagrangian strain tensor
   */
  E_n12_p =  compute_LagrangianStrain(C_n12_p);
  
  /*
    Compute the Stress tensor
   */
  S_p = grad_energy_Saint_Venant_Kirchhoff(S_p, E_n12_p, MatProp_p);
  
  /*
    Free auxiliar variables
   */
  free_Tensor(C_n_p);
  free_Tensor(C_n1_p);
  free_Tensor(C_n12_p);
  free_Tensor(E_n12_p);

  return S_p;
}

/**************************************************************/

Tensor Itegration_Stress_Average(Tensor S_p,
				 Tensor F_n1_p, Tensor F_n_p,
				 Material MatProp_p)
{

  int Ndim = NumberDimensions;
  Tensor C_n_p;
  Tensor C_n1_p;
  Tensor E_n_p;
  Tensor E_n1_p;
  Tensor S_n_p;
  Tensor S_n1_p;
  Tensor S_n12_p;
  
  /*
    Compute the right Cauchy-Green tensor
  */
  C_n_p  = compute_RightCauchyGreen(F_n_p);
  C_n1_p = compute_RightCauchyGreen(F_n1_p);

  /*
    Compute the lagrangian strain tensor
  */
  E_n_p =  compute_LagrangianStrain(C_n_p);
  E_n1_p =  compute_LagrangianStrain(C_n1_p);
    
  /*
    Compute the Stress tensor
  */
  S_n_p = alloc_Tensor(2);
  S_n1_p = alloc_Tensor(2);
  S_n_p = grad_energy_Saint_Venant_Kirchhoff(S_n_p, E_n_p, MatProp_p);
  S_n1_p = grad_energy_Saint_Venant_Kirchhoff(S_n1_p, E_n1_p, MatProp_p);
  S_n12_p = compute_midpoint_Tensor(S_n1_p,S_n_p,0.5);

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
  free_Tensor(C_n_p);
  free_Tensor(C_n1_p);
  free_Tensor(E_n_p);
  free_Tensor(E_n1_p);
  free_Tensor(S_n_p);
  free_Tensor(S_n1_p);
  free_Tensor(S_n12_p);

  return S_p;  
}

/**************************************************************/

/* Tensor Itegration_Stress_Simo(); */


/* double compute_beta_Simo(Tensor C_n, Tensor C_n1) */
/* { */
/*   /\* Number of dimensions *\/ */
/*   int Ndim = NumberDimensions; */
  
/*   /\* Initialice beta to 0.5 *\/ */
/*   double beta = 0.5; */
/*   double dg_dbeta, ddg_ddbeta; */
/*   Tensor C_nbeta; */

    
/*   dg_dbeta = */

  
/* } */

/**************************************************************/

/* double compute_dg_dbeta(Tensor C_n, Tensor C_nbeta, Tensor C_n1, */
/* 			double (* energy)(Tensor, Material), */
/* 			Tensor (* grad_energy)(Tensor, Material)) */
/* { */
/* } */

/**************************************************************/
