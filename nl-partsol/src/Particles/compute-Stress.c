#include "nl-partsol.h"

/*************************************************************/

Tensor explicit_integration_stress__Particles__(
  int p,
  Particle MPM_Mesh,
  Material MatProp)
{   
  int Ndim = NumberDimensions;

  Tensor Stress = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p],2); 
  Tensor Strain = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Strain.nM[p], 2);

  State_Parameters Input_SP;
  State_Parameters Output_SP;

  /*
    Select the constitutive model 
  */
  if(strcmp(MatProp.Type,"SR") == 0)
  {
    Stress = SolidRigid(Stress);
  }  
  else if(strcmp(MatProp.Type,"LE") == 0)
  {
    Input_SP.Stress = MPM_Mesh.Phi.Stress.nM[p];
    Input_SP.Strain = MPM_Mesh.Phi.Strain.nM[p];

    Output_SP = compute_kirchhoff_isotropic_linear_elasticity(Input_SP,MatProp);
  }
  else if(strcmp(MatProp.Type,"Von-Mises") == 0)
  {
    Output_SP = compute_kirchhoff_isotropic_linear_elasticity(Input_SP,MatProp);

    Input_SP.Stress = MPM_Mesh.Phi.Stress.nM[p];
    Input_SP.EPS = MPM_Mesh.Phi.EPS.nV[p];
    Input_SP.Back_stress = MPM_Mesh.Phi.Back_stress.nM[p];

    if(strcmp(MatProp.Plastic_Solver,"Backward-Euler") == 0)
    {
      Output_SP = Von_Mises_backward_euler(Input_SP,MatProp);
    }
    else if(strcmp(MatProp.Plastic_Solver,"Forward-Euler") == 0)
    {
      Output_SP = Von_Mises_forward_euler(Input_SP,MatProp);
    }
    else
    {
      fprintf(stderr,"%s : %s %s %s \n","Error in stress_integration__Particles__()",
    "The solver",MatProp.Type,"has not been yet implemented");
      exit(EXIT_FAILURE);
    }

    MPM_Mesh.Phi.EPS.nV[p] = Output_SP.EPS;

    free(Output_SP.Increment_E_plastic);
  }
  else
  {
    exit(EXIT_FAILURE);
  }
  
  /* Return the stress tensor */
  return Stress;
}

/**************************************************************/

void Stress_integration__Particles__(
  int p,
  Particle MPM_Mesh,
  Mesh FEM_Mesh,
  Material MatProp_p)
{
  int Ndim = NumberDimensions;

  State_Parameters Input_SP;
  State_Parameters Output_SP;

  if(strcmp(MatProp_p.Type,"Saint-Venant-Kirchhoff") == 0)
  {
    Input_SP.Stress = MPM_Mesh.Phi.Stress.nM[p];
    
    if(MatProp_p.Locking_Control_Fbar)
    {
      Input_SP.F_n1_p = MPM_Mesh.Phi.Fbar.nM[p];
    }
    else
    {
      Input_SP.F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
    }

    Output_SP = compute_1PK_Stress_Tensor_Saint_Venant_Kirchhoff(Input_SP, MatProp_p);

  }
  else if(strcmp(MatProp_p.Type,"Neo-Hookean-Wriggers") == 0)
  {
    Input_SP.Stress = MPM_Mesh.Phi.Stress.nM[p];
        
    if(MatProp_p.Locking_Control_Fbar)
    {
      Input_SP.F_n1_p = MPM_Mesh.Phi.Fbar.nM[p];
      Input_SP.J = MPM_Mesh.Phi.Jbar.nV[p];
    }
    else
    {
      Input_SP.F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
      Input_SP.J = MPM_Mesh.Phi.J_n1.nV[p];
    }

    Output_SP = compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(Input_SP, MatProp_p);

  }
  else if(strcmp(MatProp_p.Type,"Newtonian-Fluid-Compressible") == 0)
  {
    Input_SP.Stress = MPM_Mesh.Phi.Stress.nM[p];
    Input_SP.dFdt = MPM_Mesh.Phi.dt_F_n1.nM[p];
    
    if(MatProp_p.Locking_Control_Fbar)
    {
//      Input_SP.J = MPM_Mesh.Phi.Jbar.nV[p];
      Input_SP.F_n1_p = MPM_Mesh.Phi.Fbar.nM[p];
    }
    else
    {
      Input_SP.J = MPM_Mesh.Phi.J_n1.nV[p];
      Input_SP.F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
    }

    Output_SP = compute_1PK_Stress_Tensor_Newtonian_Fluid(Input_SP,MatProp_p);

  }
  else if(strcmp(MatProp_p.Type,"Newtonian-Fluid-Incompressible") == 0)
  {
    Input_SP.Stress = MPM_Mesh.Phi.Stress.nM[p];
    Input_SP.dFdt = MPM_Mesh.Phi.dt_F_n1.nM[p];
    Input_SP.J = MPM_Mesh.Phi.J_n1.nV[p];
    Input_SP.F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
    Input_SP.Pressure = MPM_Mesh.Phi.lambda_pressure_n1.nV[p];

    Output_SP = compute_1PK_Stress_Tensor_Newtonian_Fluid_Incompressible(Input_SP,MatProp_p);

  }
  else if(strcmp(MatProp_p.Type,"Von-Mises") == 0)
  {

    Input_SP.Stress = MPM_Mesh.Phi.Stress.nM[p];
    Input_SP.F_m1_plastic_p = MPM_Mesh.Phi.F_m1_plastic.nM[p];
    Input_SP.EPS = MPM_Mesh.Phi.EPS.nV[p];
    Input_SP.Back_stress = MPM_Mesh.Phi.Back_stress.nM[p];

    if(MatProp_p.Locking_Control_Fbar)
    {
      Input_SP.F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
      Input_SP.Fbar = MPM_Mesh.Phi.Fbar.nM[p];
    }
    else
    {
      Input_SP.F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
    }
  
    if(strcmp(MatProp_p.Plastic_Solver,"Backward-Euler") == 0)
    {
      Output_SP = finite_strain_plasticity(Input_SP,MatProp_p,Von_Mises_backward_euler);
    }
    else if(strcmp(MatProp_p.Plastic_Solver,"Forward-Euler") == 0)
    {
      Output_SP = finite_strain_plasticity(Input_SP,MatProp_p,Von_Mises_forward_euler);
    }
    else
    {
      fprintf(stderr,"%s : %s %s %s \n","Error in stress_integration__Particles__()",
    "The solver",MatProp_p.Type,"has not been yet implemented");
      exit(EXIT_FAILURE);
    }

    MPM_Mesh.Phi.EPS.nV[p] = Output_SP.EPS;
  }
    else if(strcmp(MatProp_p.Type,"Granular") == 0)
    {
      Input_SP.Stress = MPM_Mesh.Phi.Stress.nM[p];
      Input_SP.F_m1_plastic_p = MPM_Mesh.Phi.F_m1_plastic.nM[p];
      Input_SP.Kappa = MPM_Mesh.Phi.Kappa_hardening.nV[p];
      Input_SP.EPS = MPM_Mesh.Phi.EPS.nV[p];

 //     if(MatProp_p.Locking_Control_Fbar)
 //     {
 //       Input_SP.F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
 //       Input_SP.Fbar = MPM_Mesh.Phi.Fbar.nM[p];
 //     }
 //     else
 //     {
        Input_SP.F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
 //     }
  
      Output_SP = finite_strain_plasticity(Input_SP,MatProp_p,Frictional_Monolithic);

      MPM_Mesh.Phi.Kappa_hardening.nV[p] = Output_SP.Kappa;
      MPM_Mesh.Phi.EPS.nV[p] = Output_SP.EPS;
      
    }

  else
  {
    fprintf(stderr,"%s : %s %s %s \n","Error in forward_integration_Stress__Particles__()",
    "The material",MatProp_p.Type,"has not been yet implemnented");
    exit(EXIT_FAILURE);
  }

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

Tensor average_strain_integration_Stress__Particles__(
  Tensor S_p,
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
