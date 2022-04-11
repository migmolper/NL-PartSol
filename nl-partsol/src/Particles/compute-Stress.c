#include <string.h>
#include "nl-partsol.h"

/*************************************************************/

Tensor explicit_integration_stress__Particles__(int p, Particle MPM_Mesh,
                                                Material MatProp) {
  int Ndim = NumberDimensions;

  Tensor Stress = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[p], 2);
  Tensor Strain = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Strain.nM[p], 2);

  State_Parameters Input_SP;
  State_Parameters Output_SP;

  /*
    Select the constitutive model
  */
  if (strcmp(MatProp.Type, "SR") == 0) {
    Stress = SolidRigid(Stress);
  } else if (strcmp(MatProp.Type, "LE") == 0) {
    Input_SP.Stress = MPM_Mesh.Phi.Stress.nM[p];
    Input_SP.Strain = MPM_Mesh.Phi.Strain.nM[p];

    Output_SP =
        compute_kirchhoff_isotropic_linear_elasticity(Input_SP, MatProp);
  } else {
    exit(EXIT_FAILURE);
  }

  /* Return the stress tensor */
  return Stress;
}

/**************************************************************/

int Stress_integration__Particles__(int p, Particle MPM_Mesh, Mesh FEM_Mesh,
                                     Material MatProp_p) {
  int Ndim = NumberDimensions;
  int STATUS = EXIT_SUCCESS;

  State_Parameters Input_SP;
  State_Parameters Output_SP;
  State_Parameters IO_State;

  if (strcmp(MatProp_p.Type, "Saint-Venant-Kirchhoff") == 0) {
    Input_SP.Particle_Idx = p;
    Input_SP.Stress = MPM_Mesh.Phi.Stress.nM[p];

    if (MatProp_p.Locking_Control_Fbar) {
      Input_SP.D_phi_n1 = MPM_Mesh.Phi.Fbar.nM[p];
      Input_SP.J = MPM_Mesh.Phi.Jbar.nV[p];
    } else {
      Input_SP.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
      Input_SP.J = MPM_Mesh.Phi.J_n1.nV[p];
    }

    Output_SP =
        compute_1PK_Stress_Tensor_Saint_Venant_Kirchhoff(Input_SP, MatProp_p);

  } else if (strcmp(MatProp_p.Type, "Neo-Hookean-Wriggers") == 0) {
    IO_State.Particle_Idx = p;
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];

    if (MatProp_p.Locking_Control_Fbar) {
      IO_State.D_phi_n1 = MPM_Mesh.Phi.Fbar.nM[p];
      IO_State.J = MPM_Mesh.Phi.Jbar.nV[p];
    } else {
      IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
      IO_State.J = MPM_Mesh.Phi.J_n1.nV[p];
    }

    STATUS = compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(IO_State, MatProp_p);
    if(STATUS == EXIT_FAILURE){
      fprintf(stderr, ""RED"Error in compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(,)"RESET" \n");
      return EXIT_FAILURE;
    }

  } else if (strcmp(MatProp_p.Type, "Newtonian-Fluid-Compressible") == 0) {
    IO_State.Particle_Idx = p;
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];
    IO_State.dFdt = MPM_Mesh.Phi.dt_F_n1.nM[p];

    if (MatProp_p.Locking_Control_Fbar) {
      IO_State.D_phi_n1 = MPM_Mesh.Phi.Fbar.nM[p];
      IO_State.J = MPM_Mesh.Phi.Jbar.nV[p];
    } else {
      IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
      IO_State.J = MPM_Mesh.Phi.J_n1.nV[p];
    }

    STATUS = compute_1PK_Stress_Tensor_Newtonian_Fluid(IO_State, MatProp_p);
    if(STATUS == EXIT_FAILURE){
      fprintf(stderr, ""RED"Error in compute_1PK_Stress_Tensor_Newtonian_Fluid(,)"RESET" \n");
      return EXIT_FAILURE;
    }

  } else if (strcmp(MatProp_p.Type, "Newtonian-Fluid-Incompressible") == 0) {
    Input_SP.Particle_Idx = p;
    Input_SP.Stress = MPM_Mesh.Phi.Stress.nM[p];
    Input_SP.dFdt = MPM_Mesh.Phi.dt_F_n1.nM[p];
    Input_SP.J = MPM_Mesh.Phi.J_n1.nV[p];
    Input_SP.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    Input_SP.Pressure = MPM_Mesh.Phi.lambda_pressure_n1.nV[p];

    Output_SP = compute_1PK_Stress_Tensor_Newtonian_Fluid_Incompressible(
        Input_SP, MatProp_p);

  } else if (strcmp(MatProp_p.Type, "Von-Mises") == 0) {

    // Asign variables to the solver
    IO_State.Particle_Idx = p;
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];
    IO_State.Back_stress = MPM_Mesh.Phi.Back_stress.nM[p];
    IO_State.b_e = MPM_Mesh.Phi.b_e_n1.nM[p];
    IO_State.EPS = &MPM_Mesh.Phi.EPS_n1[p];
    IO_State.d_phi = MPM_Mesh.Phi.DF.nM[p];
    IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    IO_State.compute_C_ep = true;
    IO_State.C_ep = MPM_Mesh.Phi.C_ep.nM[p];

#if NumberDimensions == 2
  for (unsigned i = 0 ; i<5 ; i++) IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i]; 
#else
  for (unsigned i = 0 ; i<9 ; i++) IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i]; 
#endif
  *(IO_State.EPS) = MPM_Mesh.Phi.EPS_n[p];
    
      STATUS = compute_1PK_Von_Mises(IO_State, MatProp_p);
  if(STATUS == EXIT_FAILURE){
    fprintf(stderr, ""RED"Error in compute_1PK_Von_Mises(,)"RESET" \n");
    return EXIT_FAILURE;
  }


  } 
  
  else if (strcmp(MatProp_p.Type, "Drucker-Prager") == 0) {

    // Asign variables to the solver
    IO_State.Particle_Idx = p;
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];
    IO_State.b_e = MPM_Mesh.Phi.b_e_n1.nM[p];
    IO_State.EPS = &MPM_Mesh.Phi.EPS_n1[p];
    IO_State.Kappa = &MPM_Mesh.Phi.Kappa_n1[p];
    IO_State.d_phi = MPM_Mesh.Phi.DF.nM[p];
    IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    IO_State.Failure = &(MPM_Mesh.Phi.Status_particle[p]);
    IO_State.compute_C_ep = true;
    IO_State.C_ep = MPM_Mesh.Phi.C_ep.nM[p];

#if NumberDimensions == 2
  for (unsigned i = 0 ; i<5 ; i++) IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i]; 
#else
  for (unsigned i = 0 ; i<9 ; i++) IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i]; 
#endif
  *(IO_State.Kappa) = MPM_Mesh.Phi.Kappa_n[p];
  *(IO_State.EPS) = MPM_Mesh.Phi.EPS_n[p];

  STATUS = compute_1PK_Drucker_Prager(IO_State, MatProp_p);
  if(STATUS == EXIT_FAILURE){
    fprintf(stderr, ""RED"Error in compute_1PK_Drucker_Prager(,)"RESET" \n");
    return EXIT_FAILURE;
  }

  }  
  else if (strcmp(MatProp_p.Type, "Matsuoka-Nakai") == 0) {

    // Asign variables to the solver
    IO_State.Particle_Idx = p;
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];
    IO_State.b_e = MPM_Mesh.Phi.b_e_n1.nM[p];
    IO_State.EPS = &MPM_Mesh.Phi.EPS_n1[p];
    IO_State.Kappa = &MPM_Mesh.Phi.Kappa_n1[p];
    IO_State.d_phi = MPM_Mesh.Phi.DF.nM[p];
    IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    IO_State.Failure = &(MPM_Mesh.Phi.Status_particle[p]);
    IO_State.compute_C_ep = true;
    IO_State.C_ep = MPM_Mesh.Phi.C_ep.nM[p];

#if NumberDimensions == 2
  for (unsigned i = 0 ; i<5 ; i++) IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i]; 
#else
  for (unsigned i = 0 ; i<9 ; i++) IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i]; 
#endif
    *(IO_State.Kappa) = MPM_Mesh.Phi.Kappa_n[p];
    *(IO_State.EPS) = MPM_Mesh.Phi.EPS_n[p];

  STATUS = compute_1PK_Matsuoka_Nakai(IO_State, MatProp_p);
  if(STATUS == EXIT_FAILURE){
    fprintf(stderr, ""RED"Error in compute_1PK_Matsuoka_Nakai(,)"RESET" \n");
    return EXIT_FAILURE;
  }

  }
  else if (strcmp(MatProp_p.Type, "Lade-Duncan") == 0) {

    // Asign variables to the solver
    IO_State.Particle_Idx = p;
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];
    IO_State.b_e = MPM_Mesh.Phi.b_e_n1.nM[p];
    IO_State.EPS = &MPM_Mesh.Phi.EPS_n1[p];
    IO_State.Kappa = &MPM_Mesh.Phi.Kappa_n1[p];
    IO_State.d_phi = MPM_Mesh.Phi.DF.nM[p];
    IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    IO_State.Failure = &(MPM_Mesh.Phi.Status_particle[p]);
    IO_State.compute_C_ep = true;
    IO_State.C_ep = MPM_Mesh.Phi.C_ep.nM[p];

#if NumberDimensions == 2
  for (unsigned i = 0 ; i<5 ; i++) IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i]; 
#else
  for (unsigned i = 0 ; i<9 ; i++) IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i]; 
#endif
    *(IO_State.Kappa) = MPM_Mesh.Phi.Kappa_n[p];
    *(IO_State.EPS) = MPM_Mesh.Phi.EPS_n[p];

  STATUS = compute_1PK_Lade_Duncan(IO_State, MatProp_p);
  if(STATUS == EXIT_FAILURE){
    fprintf(stderr, ""RED"Error in compute_1PK_Lade_Duncan(,)"RESET" \n");
    return EXIT_FAILURE;
  }

  }
  else {
    fprintf(stderr, "%s : %s %s %s \n",
            "Error in forward_integration_Stress__Particles__()",
            "The material", MatProp_p.Type, "has not been yet implemnented");
    exit(EXIT_FAILURE);
  }

  return STATUS;
}

/**************************************************************/

Tensor tangent_matrix__Particles__(Tensor GRADIENT_pA,
 Tensor GRADIENT_pB, 
 Tensor F_n1_p,
 Tensor dFdt_n1_p,
 double J_p,
 double alpha_4,
 Material MatProp_p)
{

Tensor Stiffness_density_p;

if (strcmp(MatProp_p.Type, "Newtonian-Fluid-Incompressible") == 0) {
  Stiffness_density_p = compute_stiffness_density_Newtonian_Fluid_Incompressible(GRADIENT_pA, GRADIENT_pB, F_n1_p, dFdt_n1_p, J_p, alpha_4,MatProp_p);
}
else {
  fprintf(stderr, "%s : %s %s %s \n",
   "Error in assemble_Nodal_Tangent_Stiffness()", "The material",
   MatProp_p.Type, "has not been yet implemnented");
   exit(EXIT_FAILURE);
}

return Stiffness_density_p;

}


/**************************************************************/

Tensor average_strain_integration_Stress__Particles__(Tensor S_p, Tensor F_n1_p,
                                                      Tensor F_n_p,
                                                      Material MatProp_p) {

  /*
    Compute the right Cauchy-Green tensor in diferent time step
  */
  Tensor C_n_p = right_Cauchy_Green__Particles__(F_n_p);
  Tensor C_n1_p = right_Cauchy_Green__Particles__(F_n1_p);
  Tensor C_n12_p = Convex_combination__TensorLib__(C_n1_p, C_n_p, 0.5);

  /*
    Compute the Stress tensor
  */
  if (strcmp(MatProp_p.Type, "Saint-Venant-Kirchhoff") == 0) {
    S_p = grad_energy_Saint_Venant_Kirchhoff(S_p, C_n12_p, MatProp_p);
  } else if (strcmp(MatProp_p.Type, "Neo-Hookean-Wriggers") == 0) {
    double J_n12_p = 0.5 * (I3__TensorLib__(F_n_p) + I3__TensorLib__(F_n1_p));
    S_p = compute_2PK_Stress_Tensor_Neo_Hookean_Wriggers(S_p, C_n12_p, J_n12_p,
                                                         MatProp_p);
  } else {
    fprintf(stderr, "%s : %s %s %s \n",
            "Error in average_strain_integration_Stress__Particles__()",
            "The material", MatProp_p.Type, "has not been yet implemnented");
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
