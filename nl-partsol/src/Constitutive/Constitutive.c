
/**
 * @file Constitutive.c
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Constitutive/Constitutive.h"
#include "Globals.h"

/**************************************************************/

int Stress_integration__Constitutive__(int p, Particle MPM_Mesh,
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

  } else if (strcmp(MatProp_p.Type, "Hencky") == 0) {

    IO_State.Particle_Idx = p;
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];
    IO_State.W = &(MPM_Mesh.Phi.W[p]);
    IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    IO_State.J = MPM_Mesh.Phi.J_n1.nV[p];
    
    STATUS = compute_Kirchhoff_Stress_Hencky__Constitutive__(IO_State, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in compute_Kirchhoff_Stress_Hencky__Constitutive__(,)" RESET
              " \n");
      return EXIT_FAILURE;
    }
    
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

    STATUS = compute_Kirchhoff_Stress_Neo_Hookean(IO_State, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in compute_Kirchhoff_Stress_Neo_Hookean(,)" RESET
              " \n");
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

    //    IO_State.Pressure = MPM_Mesh.Phi.lambda_pressure_n1.nV[p];

    STATUS =
        compute_Kirchhoff_Stress_Tensor_Newtonian_Fluid(IO_State, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(
          stderr,
          "" RED
          "Error in compute_Kirchhoff_Stress_Tensor_Newtonian_Fluid(,)" RESET
          " \n");
      return EXIT_FAILURE;
    }

  } else if (strcmp(MatProp_p.Type, "Von-Mises") == 0) {

    // Asign variables to the solver
    IO_State.Particle_Idx = p;
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];
    IO_State.W = &(MPM_Mesh.Phi.W[p]);
    IO_State.Back_stress = MPM_Mesh.Phi.Back_stress.nM[p];
    IO_State.b_e = MPM_Mesh.Phi.b_e_n1.nM[p];
    IO_State.EPS = &MPM_Mesh.Phi.EPS_n1[p];
    IO_State.d_phi = MPM_Mesh.Phi.DF.nM[p];
    IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    IO_State.compute_C_ep = true;
    IO_State.C_ep = MPM_Mesh.Phi.C_ep.nM[p];

#if NumberDimensions == 2
    for (unsigned i = 0; i < 5; i++)
      IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i];
#else
    for (unsigned i = 0; i < 9; i++)
      IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i];
#endif
    *(IO_State.EPS) = MPM_Mesh.Phi.EPS_n[p];

    STATUS = compute_1PK_Von_Mises(IO_State, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in compute_1PK_Von_Mises(,)" RESET " \n");
      return EXIT_FAILURE;
    }

  }

  else if (strcmp(MatProp_p.Type, "Drucker-Prager") == 0) {

    // Asign variables to the solver
    IO_State.Particle_Idx = p;
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];
    IO_State.W = &(MPM_Mesh.Phi.W[p]);
    IO_State.b_e = MPM_Mesh.Phi.b_e_n1.nM[p];
    IO_State.EPS = &MPM_Mesh.Phi.EPS_n1[p];
    IO_State.Kappa = &MPM_Mesh.Phi.Kappa_n1[p];
    IO_State.d_phi = MPM_Mesh.Phi.DF.nM[p];
    IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    IO_State.Failure = &(MPM_Mesh.Phi.Status_particle[p]);
    IO_State.compute_C_ep = true;
    IO_State.C_ep = MPM_Mesh.Phi.C_ep.nM[p];

#if NumberDimensions == 2
    for (unsigned i = 0; i < 5; i++)
      IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i];
#else
    for (unsigned i = 0; i < 9; i++)
      IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i];
#endif
    *(IO_State.Kappa) = MPM_Mesh.Phi.Kappa_n[p];
    *(IO_State.EPS) = MPM_Mesh.Phi.EPS_n[p];

    STATUS = compute_Kirchhoff_Stress_Drucker_Prager__Constitutive__(IO_State, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in compute_Kirchhoff_Stress_Drucker_Prager__Constitutive__(,)" RESET " \n");
      return EXIT_FAILURE;
    }

  } else if (strcmp(MatProp_p.Type, "Matsuoka-Nakai") == 0) {

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
    for (unsigned i = 0; i < 5; i++)
      IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i];
#else
    for (unsigned i = 0; i < 9; i++)
      IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i];
#endif
    *(IO_State.Kappa) = MPM_Mesh.Phi.Kappa_n[p];
    *(IO_State.EPS) = MPM_Mesh.Phi.EPS_n[p];

    STATUS = compute_1PK_Matsuoka_Nakai(IO_State, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in compute_1PK_Matsuoka_Nakai(,)" RESET " \n");
      return EXIT_FAILURE;
    }

  } else if (strcmp(MatProp_p.Type, "Lade-Duncan") == 0) {

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
    for (unsigned i = 0; i < 5; i++)
      IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i];
#else
    for (unsigned i = 0; i < 9; i++)
      IO_State.b_e[i] = MPM_Mesh.Phi.b_e_n.nM[p][i];
#endif
    *(IO_State.Kappa) = MPM_Mesh.Phi.Kappa_n[p];
    *(IO_State.EPS) = MPM_Mesh.Phi.EPS_n[p];

    STATUS = compute_1PK_Lade_Duncan(IO_State, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr, "" RED "Error in compute_1PK_Lade_Duncan(,)" RESET " \n");
      return EXIT_FAILURE;
    }

  } else {
    fprintf(stderr, "%s : %s %s %s \n",
            "Error in forward_integration_Stress__Particles__()",
            "The material", MatProp_p.Type, "has not been yet implemnented");
    exit(EXIT_FAILURE);
  }

  return STATUS;
}

/**************************************************************/

int stiffness_density__Constitutive__(int p, double *Stiffness_density,
                                      const double *dN_alpha_n1,
                                      const double *dN_beta_n1,
                                      const double *dN_alpha_n,
                                      const double *dN_beta_n, double alpha_4,
                                      Particle MPM_Mesh, Material MatProp_p) {

  State_Parameters IO_State;
  int STATUS = EXIT_SUCCESS;

  // Compute local stiffness density
  if (strcmp(MatProp_p.Type, "Neo-Hookean-Wriggers") == 0) {
    IO_State.D_phi_n = MPM_Mesh.Phi.F_n.nM[p];
    IO_State.J = MPM_Mesh.Phi.J_n1.nV[p];
    STATUS = compute_stiffness_density_Neo_Hookean(
        Stiffness_density, dN_alpha_n1, dN_beta_n1, dN_alpha_n, dN_beta_n,
        IO_State, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in compute_stiffness_density_Neo_Hookean" RESET
              "\n");
      return EXIT_FAILURE;
    }
  } else if (strcmp(MatProp_p.Type, "Hencky") == 0) {
    IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];
    
    STATUS = compute_stiffness_density_Hencky__Constitutive__(
        Stiffness_density, dN_alpha_n1, dN_beta_n1, IO_State, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in compute_stiffness_density_Hencky__Constitutive__" RESET
              "\n");
      return EXIT_FAILURE;
    }
    
 //   printf(" %e, %e \n %e, %e\n",Stiffness_density[0],Stiffness_density[1],Stiffness_density[2],Stiffness_density[3]);


  } else if (strcmp(MatProp_p.Type, "Newtonian-Fluid-Compressible") == 0) {
    IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    IO_State.D_phi_n = MPM_Mesh.Phi.F_n.nM[p];
    IO_State.dFdt = MPM_Mesh.Phi.dt_F_n1.nM[p];
    IO_State.J = MPM_Mesh.Phi.J_n1.nV[p];
    IO_State.alpha_4 = alpha_4;
    STATUS = compute_stiffness_density_Newtonian_Fluid(
        Stiffness_density, dN_alpha_n1, dN_beta_n1, dN_alpha_n, dN_beta_n,
        IO_State, MatProp_p);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in compute_stiffness_density_Newtonian_Fluid" RESET
              "\n");
      return EXIT_FAILURE;
    }
  } else if (strcmp(MatProp_p.Type, "Drucker-Prager") == 0) {
    IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    IO_State.D_phi_n = MPM_Mesh.Phi.F_n.nM[p];
    IO_State.b_e = MPM_Mesh.Phi.b_e_n1.nM[p];
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];
    IO_State.C_ep = MPM_Mesh.Phi.C_ep.nM[p];
    STATUS = compute_1PK_elastoplastic_tangent_matrix(
        Stiffness_density, dN_alpha_n1, dN_beta_n1, IO_State);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in compute_1PK_elastoplastic_tangent_matrix" RESET
              "\n");
      return EXIT_FAILURE;
    }

  } else if (strcmp(MatProp_p.Type, "Matsuoka-Nakai") == 0) {
    IO_State.D_phi_n1 = MPM_Mesh.Phi.F_n1.nM[p];
    IO_State.D_phi_n = MPM_Mesh.Phi.F_n.nM[p];
    IO_State.b_e = MPM_Mesh.Phi.b_e_n1.nM[p];
    IO_State.Stress = MPM_Mesh.Phi.Stress.nM[p];
    IO_State.C_ep = MPM_Mesh.Phi.C_ep.nM[p];
    STATUS = compute_1PK_elastoplastic_tangent_matrix(
        Stiffness_density, dN_alpha_n1, dN_beta_n1, IO_State);
    if (STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in compute_1PK_elastoplastic_tangent_matrix" RESET
              "\n");
      return EXIT_FAILURE;
    }
  } else {
    fprintf(stderr,
            "" RED "The material %s has not been yet implemnented" RESET "\n",
            MatProp_p.Type);
    return EXIT_FAILURE;
  }

  return STATUS;
}
