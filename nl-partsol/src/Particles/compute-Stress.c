#include "nl-partsol.h"

/*************************************************************/

Tensor explicit_integration_stress__Particles__(
  int p,
  Particle MPM_Mesh,
  Material MatProp)
{   
  Tensor Stress = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2); 
  Tensor Strain = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Strain.nM[p], 2);

  Plastic_status Input_Plastic_Parameters;
  Plastic_status Output_Plastic_Parameters;

  /*
    Select the constitutive model 
  */
  if(strcmp(MatProp.Type,"SR") == 0)
  {
    Stress = SolidRigid(Stress);
  }  
  else if(strcmp(MatProp.Type,"LE") == 0)
  {
    Stress = LinearElastic(Strain,Stress,MatProp);
  }
  else if(strcmp(MatProp.Type,"Von-Mises") == 0)
  {

    Input_Plastic_Parameters.EPS = MPM_Mesh.Phi.EPS.nV[p];

    Output_Plastic_Parameters = infinitesimal_strains_plasticity_Von_Mises(Stress, Strain, Input_Plastic_Parameters, MatProp);

    MPM_Mesh.Phi.EPS.nV[p] = Output_Plastic_Parameters.EPS;
  }
  else
  {
    exit(EXIT_FAILURE);
  }
  
  /* Return the stress tensor */
  return Stress;
}

/**************************************************************/

Tensor forward_integration_Stress__Particles__(
  int p,
  Particle MPM_Mesh,
  Material MatProp_p)
{
  
  Tensor P_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2); 

  // Variables for the constitutive model
  double J_p;
  Tensor F_n1_p;
  Tensor dFdt_n1_p;
  Tensor F_plastic_p;
  Plastic_status Input_Plastic_Parameters;
  Plastic_status Output_Plastic_Parameters;

  if(strcmp(MatProp_p.Type,"Saint-Venant-Kirchhoff") == 0)
  {
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    P_p = compute_1PK_Stress_Tensor_Saint_Venant_Kirchhoff(P_p, F_n1_p, MatProp_p);
  }
  else if(strcmp(MatProp_p.Type,"Neo-Hookean-Wriggers") == 0)
  {
    J_p = MPM_Mesh.Phi.J.nV[p];
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    P_p = compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(P_p, F_n1_p, J_p, MatProp_p);
  }
  else if(strcmp(MatProp_p.Type,"Newtonian-Fluid-Compressible") == 0)
  {
    J_p = MPM_Mesh.Phi.J.nV[p];
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    dFdt_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.dt_F_n1.nM[p],2);
    P_p = compute_1PK_Stress_Tensor_Newtonian_Fluid(P_p,F_n1_p,dFdt_n1_p,J_p,MatProp_p);
  }
  else if(strcmp(MatProp_p.Type,"Von-Mises") == 0)
  {
    J_p = MPM_Mesh.Phi.J.nV[p];
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    F_plastic_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_plastic.nM[p],2);
    Input_Plastic_Parameters.EPS = MPM_Mesh.Phi.EPS.nV[p];

    Output_Plastic_Parameters = finite_strains_plasticity_Von_Mises(P_p,F_plastic_p,F_n1_p,Input_Plastic_Parameters,MatProp_p,J_p);

    MPM_Mesh.Phi.EPS.nV[p] = Output_Plastic_Parameters.EPS;
  }
  else if((strcmp(MatProp_p.Type,"Drucker-Prager-Plane-Strain") == 0) || (strcmp(MatProp_p.Type,"Drucker-Prager-Outer-Cone") == 0))
  {
    J_p = MPM_Mesh.Phi.J.nV[p];
    F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
    F_plastic_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_plastic.nM[p],2);
    Input_Plastic_Parameters.Cohesion = MPM_Mesh.Phi.cohesion.nV[p];
    Input_Plastic_Parameters.EPS = MPM_Mesh.Phi.EPS.nV[p];

    Output_Plastic_Parameters = finite_strains_plasticity_Drucker_Prager_Sanavia(P_p,F_plastic_p,F_n1_p,Input_Plastic_Parameters,MatProp_p,J_p);

    MPM_Mesh.Phi.cohesion.nV[p] = Output_Plastic_Parameters.Cohesion;
    MPM_Mesh.Phi.EPS.nV[p] = Output_Plastic_Parameters.EPS;

  }
  else
  {
    fprintf(stderr,"%s : %s %s %s \n","Error in forward_integration_Stress__Particles__()",
    "The material",MatProp_p.Type,"has not been yet implemnented");
    exit(EXIT_FAILURE);
  }

  return P_p;
}

/**************************************************************/

Tensor configurational_midpoint_integration_Stress__Particles__(
  Tensor S_p,
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
      S_p = compute_2PK_Stress_Tensor_Neo_Hookean_Wriggers(S_p, C_n12_p, J_n12_p, MatProp_p);
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
      S_p = compute_2PK_Stress_Tensor_Neo_Hookean_Wriggers(S_p, C_n12_p, J_n12_p, MatProp_p);
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
      S_n_p  = compute_2PK_Stress_Tensor_Neo_Hookean_Wriggers(S_n_p, C_n_p, J_n_p, MatProp_p);
      S_n1_p = compute_2PK_Stress_Tensor_Neo_Hookean_Wriggers(S_n1_p, C_n1_p, J_n1_p, MatProp_p);
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

void compute_Piola_transformation__Particles__(Tensor S_p, Tensor sigma_k1, Tensor F_total, double J)
{
  
  int Ndim = NumberDimensions;

  contravariant_pull_back_tensor__TensorLib__(S_p, sigma_k1, F_total);
  
  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        S_p.N[i][j] = J*S_p.N[i][j];
      }
    }
}

/**************************************************************/
