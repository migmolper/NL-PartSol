/**
 * @file U-pw-Functions.c
 * @author Miguel Molinos (@migmolper)
 * @brief Functions to compute the residual and the jacobian
 * for the equations of balance of momentum of the mixture
 * and the darcy equation
 * @version 0.1
 * @date 2022-05-15
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Formulations/Displacements-WaterPressure/U-pw-Functions.h"

/**************************************************************/

void local_compatibility_conditions__upw__(const Nodal_Field D_upw,
                                           Mask ActiveNodes, Particle MPM_Mesh,
                                           Mesh FEM_Mesh, int *STATUS) {

  /*
    Auxiliar variables
  */
  *STATUS = EXIT_SUCCESS;
  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned Nnodes_mask = ActiveNodes.Nactivenodes;
  unsigned NumNodes_p;
  unsigned Ndof_U_p;
  unsigned Ndof_pw_p;
  int Mixture_idx;
  int Material_Soil_idx;
  int Material_Water_idx;

#ifdef USE_PETSC
  const double *dU;
  VecGetArrayRead(D_U.value, &dU);
  const double *dU_dt;
  VecGetArrayRead(D_U.d_value_dt, &dU_dt);
  const double *dPw;
  VecGetArrayRead(D_U.value, &dPw);
  const double *dPw_dt;
  VecGetArrayRead(D_U.d_value_dt, &dPw_dt);
#else
  const double *dU = D_U.value;
  const double *dU_dt = D_U.d_value_dt;
  const double *dPw = D_U.value;
  const double *dPw_dt = D_U.d_value_dt;
#endif

  double Jacobian;      /* Jacobian of the deformation gradient at t = n + 1 */
  double rate_Jacobian; /* Rate of the Jacobian of the deformation gradient at t
                        = n + 1 */
  double K_f;           /* Compressibility (fluid) */
  double rho_f_0;       /* Initial density of the fluid */
  double phi_s_0;       /* Initial volume fraction (solid) */
  double phi_f_0;       /* Initial volume fraction (fluid) */
  Element Nodes_p;      /* Element for each particle */
  Material MatProp_Soil_p; /* Variable with the material properties of the solid
                              phase */
  Material MatProp_Water_p; /* Variable with the material properties of the
                               fluid phase */
  Matrix Nodal_D_Displacement_p;
  Matrix Nodal_D_Velocity_p;
  Matrix Nodal_D_theta_p;
  Matrix Nodal_D_theta_dt;
  double *F_n_p;  /* Deformation gradient of the soil skeleton (t = n) */
  double *F_n1_p; /* Deformation gradient of the soil skeleton (t = n + 1) */
  double *DF_p; /* Increment of the deformation gradient of the soil skeleton */
  double *dFdt_n_p;  /* Rate of the deformation gradient of the soil skeleton (t
                      =  n) */
  double *dFdt_n1_p; /* Rate of the deformation gradient of the soil skeleton (t
                      = n + 1) */
  double *dt_DF_p;   /* Rate of the increment of the deformation gradient of the
                      soil skeleton */

  for (int p = 0; p < Np; p++) {

    //!  Define nodal connectivity for each particle
    //!  and compute gradient of the shape function
    NumNodes_p = MPM_Mesh.NumberNodes[p];
    Element Nodes_p =
        nodal_set__Particles__(p, MPM_Mesh.ListNodes[p], NumNodes_p);
    Matrix shapefunction_n_p =
        compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    Matrix d_shapefunction_n_p =
        compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    Ndof_U_p = NumNodes_p * Ndim;
    Ndof_pw_p = NumNodes_p;

    //  Get the nodal increment of displacement using the mask
    double *dU_Ap = (double *)calloc(Ndof_U_p, __SIZEOF_DOUBLE__);
    double *dU_dt_Ap = (double *)calloc(Ndof_U_p, __SIZEOF_DOUBLE__);
    double *dPw_Ap = (double *)calloc(Ndof_pw_p, __SIZEOF_DOUBLE__);
    double *dPw_dt_Ap = (double *)calloc(Ndof_pw_p, __SIZEOF_DOUBLE__);

    if ((dU_Ap == NULL) || (dU_dt_Ap == NULL) || (dPw_Ap == NULL) ||
        (dPw_dt_Ap == NULL)) {
      fprintf(stderr, "" RED "Error in calloc(): Out of memory" RESET " \n");
      *STATUS = EXIT_FAILURE;
    }

    // Get nodal field (displacement-pressure) which interacts with the particle
    get_set_vectorial_field__MeshTools__(dU_Ap, dU, Nodes_p, ActiveNodes);
    get_set_vectorial_field__MeshTools__(dU_dt_Ap, dU_dt, Nodes_p, ActiveNodes);

    get_set_scalar_field__MeshTools__(dPw_Ap, dPw, Nodes_p, ActiveNodes);
    get_set_scalar_field__MeshTools__(dPw_dt_Ap, dPw_dt, Nodes_p, ActiveNodes);

    //  Evaluate the shape function and its gradient in the coordinates of the particle
    Matrix d_shapefunction_n_p = compute_dN__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);
    Matrix shapefunction_n_p = compute_N__MeshTools__(Nodes_p, MPM_Mesh, FEM_Mesh);

    //  Take the values of the deformation gradient from the previous step
    const double *F_n_p = MPM_Mesh.Phi.F_n.nM[p];
    double *F_n1_p = MPM_Mesh.Phi.F_n1.nM[p];
    double *DF_p = MPM_Mesh.Phi.DF.nM[p];
    const double *dFdt_n_p = MPM_Mesh.Phi.dt_F_n.nM[p];
    double *dFdt_n1_p = MPM_Mesh.Phi.dt_F_n1.nM[p];
    double *dt_DF_p = MPM_Mesh.Phi.dt_DF.nM[p];

    //  Update the deformation gradient in t = n + 1 with the information
    //  from t = n and the increment of deformation gradient.
    update_increment_Deformation_Gradient__Particles__(
        DF_p, dU_Ap, d_shapefunction_n_p.nV, NumNodes_p);

    update_rate_increment_Deformation_Gradient__Particles__(
        dt_DF_p, dU_dt_Ap, d_shapefunction_n_p.nV, NumNodes_p);

    update_Deformation_Gradient_n1__Particles__(F_n1_p, F_n_p, DF_p);

    update_rate_Deformation_Gradient_n1__Particles__(dFdt_n1_p, dt_DF_p, F_n_p,
                                                     DF_p, dFdt_n_p);

    //  Compute the Jacobian of the deformation gradient and its rate
    Jacobian = I3__TensorLib__(F_n1_p);
    if (MPM_Mesh.Phi.J_n1.nV[p] <= 0.0) {
      fprintf(stderr, "" RED "Negative jacobian in particle %i" RESET " \n", p);
      *STATUS = EXIT_FAILURE;
    }

    *STATUS = compute_Jacobian_Rate__Particles__(&rate_Jacobian, Jacobian,
                                                 F_n1_p, dFdt_n1_p);
    if (*STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in compute_Jacobian_Rate__Particles__()" RESET
              " \n");
    }

    // Update soil skeleton jacobian and its rate
    MPM_Mesh.Phi.J_n1.nV[p] = Jacobian;
    MPM_Mesh.Phi.dJ_dt.nV[p] = rate_Jacobian;

    // Take the values of the pore water pressure from previous steps and update
    double Pw_n_p = MPM_Mesh.Phi.Pw.nV[p];
    double Pw_dt_n_p = MPM_Mesh.Phi.d_Pw_dt_n.nV[p];
    double * Pw_n1_p = &MPM_Mesh.Phi.Pw_n1.nV[p];
    double * Pw_dt_n1_p = &MPM_Mesh.Phi.d_Pw_dt_n1.nV[p];

    *Pw_n1_p = Pw_n_p;
    *Pw_dt_n1_p = Pw_dt_n_p;
    for (unsigned A = 0; A < NumNodes_p; A++) {
      *Pw_n1_p += dPw_Ap[A] * shapefunction_n_p.nV[A];
      *Pw_dt_n1_p += dPw_dt_Ap[A] * shapefunction_n_p.nV[A];
    }

    /*
      Free memory
    */
    free(dU_Ap);
    free(dU_dt_Ap);
    free(dPw_Ap);
    free(dPw_dt_Ap);
    free__MatrixLib__(d_shapefunction_n_p);
    free__MatrixLib__(shapefunction_n_p);
    free(Nodes_p.Connectivity);
  }
}

/**************************************************************/

void constitutive_update__upw__(Particle MPM_Mesh, Mesh FEM_Mesh,
                                  int *STATUS) {
  *STATUS = EXIT_SUCCESS;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned MixtIdx_p;
  unsigned MatIndx_Soil_p;
  unsigned MatIndx_Water_p;
  unsigned p;

#pragma omp for private(p, MixtIdx_p, MatIndx_Soil_p, MatIndx_Water_p)
  for (p = 0; p < Np; p++) {

    /*
      Load material properties for each phase
    */
    MixtIdx_p = MPM_Mesh.MixtIdx[p];
    MatIndx_Soil_p = Soil_Water_Mixtures[MixtIdx_p].Soil_Idx;
    MatIndx_Water_p = Soil_Water_Mixtures[MixtIdx_p].Water_Idx;
    Material MatProp_Soil_p = MPM_Mesh.Mat[MatIndx_Soil_p];
    Material MatProp_Water_p = MPM_Mesh.Mat[MatIndx_Water_p];

    //  Update the Kirchhoff stress tensor with an apropiate
    //  integration rule.
    *STATUS = Stress_integration__Constitutive__(p, MPM_Mesh, MatProp_Soil_p);
    if (*STATUS == EXIT_FAILURE) {
      fprintf(stderr,
              "" RED "Error in Stress_integration__Constitutive__(,)" RESET
              " \n");
    }

    // Get the Jacobian
    double Jacobian = MPM_Mesh.Phi.J_n1.nV[p];

    //  Read intrinsic material properties (fluid)
    double K_f = MatProp_Water_p.Compressibility;
    double rho_f_0 = MatProp_Water_p.rho;

    //  Read reference volume fraction for each phase
    double phi_f_0 = Soil_Water_Mixtures[MixtIdx_p].phi_f_0;
    double phi_s_0 = Soil_Water_Mixtures[MixtIdx_p].phi_s_0;

    // Get the initial pressure
    double Pw_0 = MPM_Mesh.Phi.Pw_0.nV[p];

    // From the Kirchhoff pressure compure the cauchy pore water pressure
    double Pw_n1 = MPM_Mesh.Phi.Pw_n1.nV[p] / Jacobian;

    // Update the intrinsic fluid density
    MPM_Mesh.Phi.rho_f.nV[p] = rho_f_0 * exp((Pw_n1 - Pw_0) / K_f);

    // Update the volume fraction of the solid/fluid phase
    MPM_Mesh.Phi.phi_s.nV[p] = phi_s_0 / Jacobian;
    MPM_Mesh.Phi.phi_f.nV[p] = 1.0 - (1.0 - phi_f_0) / Jacobian;

    // Update density of the mixture
    MPM_Mesh.Phi.rho.nV[p] =
        MPM_Mesh.Phi.rho_s.nV[p] * MPM_Mesh.Phi.phi_s.nV[p] +
        MPM_Mesh.Phi.rho_f.nV[p] * MPM_Mesh.Phi.phi_f.nV[p];
  }
}

/**************************************************************/

void compute_gradient_kirchoff_pw__upw__(
  double * gradient_theta_n1_p,
  const double * dPw,
  const double * Pw_n,
  const double * d_shapefunction_n1_p,
  unsigned Nnodes_p){


  unsigned Ndim = NumberDimensions;
  double Pw_n1;

  for (unsigned A = 0; A < Nnodes_p; A++) {

    //  Get nodal value (A) of the kirchhoff pore water pressure at n+1
    Pw_n1 = Pw_n[A] + dPw[A];

    //  Compute nodal contribution
    for (unsigned i = 0; i < Ndim; i++) {
      gradient_theta_n1_p[i] += Pw_n1 * d_shapefunction_n1_p[A*Ndim + i];
    }
  }

}

/**************************************************************/

void compute_total_acceleration__upw__(double *b_n1_p, const double *a_n_p,
                                         const double *da_I,
                                         const double *shapefunction_n_p,
                                         const ChainPtr ListNodes_p,
                                         Mask ActiveNodes) {

  unsigned Ndim = NumberDimensions;

  //  Compute mixture/fluid dynamics in the next time step
  for (unsigned i = 0; i < Ndim; i++) {
    b_n1_p[i] = a_n_p[i] - gravity_field.Value[i].Fx[TimeStep];
  }

  ChainPtr Idx = ListNodes_p;
  unsigned A = 0;
  while (Idx != NULL) {
    // Get the shape function evaluation in node A
    // and the masked index of the node A
    double shapefunction_n_pA = shapefunction_n_p[A];
    int Ap = Idx->Idx;
    int Mask_node_A = ActiveNodes.Nodes2Mask[Ap];

    for (unsigned i = 0; i < Ndim; i++) {
      b_n1_p[i] += shapefunction_n_pA * da_I[Mask_node_A * Ndim + i];
    }

    Idx = Idx->next;
    A++;
  }
}


/**************************************************************/

void compute_L0__upw__(double *L0, const double *dN_alpha_u_n1,
                       const double *kichhoff_stress, double V0) {

#if NumberDimensions == 2
  L0[0] = (kichhoff_stress[0] * dN_alpha_u_n1[0] +
           kichhoff_stress[1] * dN_alpha_u_n1[1]) *
          V0;
  L0[1] = (kichhoff_stress[2] * dN_alpha_u_n1[0] +
           kichhoff_stress[3] * dN_alpha_u_n1[1]) *
          V0;
#else
  L0[0] = (kichhoff_stress[0] * dN_alpha_u_n1[0] +
           kichhoff_stress[1] * dN_alpha_u_n1[1] +
           kichhoff_stress[2] * dN_alpha_u_n1[2]) *
          V0;
  L0[1] = (kichhoff_stress[3] * dN_alpha_u_n1[0] +
           kichhoff_stress[4] * dN_alpha_u_n1[1] +
           kichhoff_stress[5] * dN_alpha_u_n1[2]) *
          V0;
  L0[2] = (kichhoff_stress[6] * dN_alpha_u_n1[0] +
           kichhoff_stress[7] * dN_alpha_u_n1[1] +
           kichhoff_stress[8] * dN_alpha_u_n1[2]) *
          V0;
#endif
}

/**************************************************************/

void compute_L1__upw__(double *L1, const double *dN_alpha_u_n1,
                       double kichhoff_pressure, double V0) {

#if NumberDimensions == 2
  L1[0] = -kichhoff_pressure * dN_alpha_u_n1[0] * V0;
  L1[1] = -kichhoff_pressure * dN_alpha_u_n1[1] * V0;
#else
  L1[0] = -kichhoff_pressure * dN_alpha_u_n1[0] * V0;
  L1[1] = -kichhoff_pressure * dN_alpha_u_n1[1] * V0;
  L1[2] = -kichhoff_pressure * dN_alpha_u_n1[2] * V0;
#endif
}

/**************************************************************/

void compute_DL1_Du__upw__(double *DL1_Du, const double *dN_alpha_u_n1,
                           const double *dN_beta_u_n1, double kichhoff_pressure,
                           double V0) {

#if NumberDimensions == 2
  DL1_Du[0] = kichhoff_pressure * dN_beta_u_n1[0] * dN_alpha_u_n1[0] * V0;
  DL1_Du[1] = kichhoff_pressure * dN_beta_u_n1[0] * dN_alpha_u_n1[1] * V0;

  DL1_Du[2] = kichhoff_pressure * dN_beta_u_n1[1] * dN_alpha_u_n1[0] * V0;
  DL1_Du[3] = kichhoff_pressure * dN_beta_u_n1[1] * dN_alpha_u_n1[1] * V0;
#else
  DL1_Du[0] = kichhoff_pressure * dN_beta_u_n1[0] * dN_alpha_u_n1[0] * V0;
  DL1_Du[1] = kichhoff_pressure * dN_beta_u_n1[0] * dN_alpha_u_n1[1] * V0;
  DL1_Du[2] = kichhoff_pressure * dN_beta_u_n1[0] * dN_alpha_u_n1[2] * V0;

  DL1_Du[3] = kichhoff_pressure * dN_beta_u_n1[1] * dN_alpha_u_n1[0] * V0;
  DL1_Du[4] = kichhoff_pressure * dN_beta_u_n1[1] * dN_alpha_u_n1[1] * V0;
  DL1_Du[5] = kichhoff_pressure * dN_beta_u_n1[1] * dN_alpha_u_n1[2] * V0;

  DL1_Du[6] = kichhoff_pressure * dN_beta_u_n1[2] * dN_alpha_u_n1[0] * V0;
  DL1_Du[7] = kichhoff_pressure * dN_beta_u_n1[2] * dN_alpha_u_n1[1] * V0;
  DL1_Du[8] = kichhoff_pressure * dN_beta_u_n1[2] * dN_alpha_u_n1[2] * V0;
#endif
}

/**************************************************************/

void compute_DL1_Dpw__upw__(double *DL1_Dpw, const double *dN_alpha_u_n1,
                            double N_beta_pw_n1, double V0) {

#if NumberDimensions == 2
  DL1_Dpw[0] = -dN_alpha_u_n1[0] * N_beta_pw_n1 * V0;
  DL1_Dpw[1] = -dN_alpha_u_n1[1] * N_beta_pw_n1 * V0;
#else
  DL1_Dpw[0] = -dN_alpha_u_n1[0] * N_beta_pw_n1 * V0;
  DL1_Dpw[1] = -dN_alpha_u_n1[1] * N_beta_pw_n1 * V0;
  DL1_Dpw[2] = -dN_alpha_u_n1[2] * N_beta_pw_n1 * V0;
#endif
}

/**************************************************************/

void compute_L2__upw__(double *L2, double N_alpha_u_n1, const double *b,
                       double m) {
#if NumberDimensions == 2
  L2[0] = -N_alpha_u_n1 * m * b[0];
  L2[1] = -N_alpha_u_n1 * m * b[1];
#else
  L2[0] = -N_alpha_u_n1 * m * b[0];
  L2[1] = -N_alpha_u_n1 * m * b[1];
  L2[2] = -N_alpha_u_n1 * m * b[2];
#endif
}

/**************************************************************/

void compute_DL2_Du__upw__(double *DL2_Du, double N_alpha_u_n1,
                           double N_beta_u_n1, const double *dN_beta_u_n1,
                           const double *b, double kichhoff_pressure,
                           double phi_f, double intrinsic_rho_f, double kappa_f,
                           double Jacobian, double m, double alpha_1,
                           double V0) {

  double c0 =
      intrinsic_rho_f * (Jacobian - kichhoff_pressure * phi_f / kappa_f);

#if NumberDimensions == 2
  DL2_Du[0] = -c0 * N_alpha_u_n1 * b[0] * dN_beta_u_n1[0] * V0 -
              alpha_1 * N_alpha_u_n1 * N_beta_u_n1 * m;
  DL2_Du[1] = -c0 * N_alpha_u_n1 * b[0] * dN_beta_u_n1[1] * V0;

  DL2_Du[2] = -c0 * N_alpha_u_n1 * b[1] * dN_beta_u_n1[0] * V0;
  DL2_Du[3] = -c0 * N_alpha_u_n1 * b[1] * dN_beta_u_n1[1] * V0 -
              alpha_1 * N_alpha_u_n1 * N_beta_u_n1 * m;
#else
  DL2_Du[0] = -c0 * N_alpha_u_n1 * b[0] * dN_beta_u_n1[0] * V0 -
              alpha_1 * N_alpha_u_n1 * N_beta_u_n1 * m;
  DL2_Du[1] = -c0 * N_alpha_u_n1 * b[0] * dN_beta_u_n1[1] * V0;
  DL2_Du[2] = -c0 * N_alpha_u_n1 * b[0] * dN_beta_u_n1[2] * V0;

  DL2_Du[3] = -c0 * N_alpha_u_n1 * b[1] * dN_beta_u_n1[0] * V0;
  DL2_Du[4] = -c0 * N_alpha_u_n1 * b[1] * dN_beta_u_n1[1] * V0 -
              alpha_1 * N_alpha_u_n1 * N_beta_u_n1 * m;
  DL2_Du[5] = -c0 * N_alpha_u_n1 * b[1] * dN_beta_u_n1[2] * V0;

  DL2_Du[6] = -c0 * N_alpha_u_n1 * b[2] * dN_beta_u_n1[0] * V0;
  DL2_Du[7] = -c0 * N_alpha_u_n1 * b[2] * dN_beta_u_n1[1] * V0;
  DL2_Du[8] = -c0 * N_alpha_u_n1 * b[2] * dN_beta_u_n1[2] * V0 -
              alpha_1 * N_alpha_u_n1 * N_beta_u_n1 * m;
#endif
}

/**************************************************************/

void compute_DL2_Dpw__upw__(double *DL2_Dpw, double N_alpha_u_n1,
                            double N_beta_pw_n1, const double *b, double phi_f,
                            double intrinsic_rho_f, double kappa_f, double V0) {

  double c0 = phi_f * intrinsic_rho_f / kappa_f;

#if NumberDimensions == 2
  DL2_Dpw[0] = -c0 * b[0] * N_alpha_u_n1 * N_beta_pw_n1 * V0;
  DL2_Dpw[1] = -c0 * b[1] * N_alpha_u_n1 * N_beta_pw_n1 * V0;
#else
  DL2_Dpw[0] = -c0 * b[0] * N_alpha_u_n1 * N_beta_pw_n1 * V0;
  DL2_Dpw[1] = -c0 * b[1] * N_alpha_u_n1 * N_beta_pw_n1 * V0;
  DL2_Dpw[2] = -c0 * b[2] * N_alpha_u_n1 * N_beta_pw_n1 * V0;
#endif
}

/**************************************************************/

void compute_G0__upw__(double *G0, double N_alpha_pw_n1, double intrinsic_rho_f,
                       double relative_rho_f_p, double rate_kichhoff_pressure,
                       double rate_Jacobian, double kappa_f, double V0) {
  *G0 = N_alpha_pw_n1 *
        (relative_rho_f_p * rate_kichhoff_pressure / kappa_f +
         intrinsic_rho_f * rate_Jacobian) *
        V0;
}

/**************************************************************/

void compute_DG0_Du__upw__(double *DG0_Du, double N_alpha_pw_n1,
                           const double *dN_beta_u_n1, const double *grad_v,
                           double kichhoff_pressure,
                           double rate_kichhoff_pressure, double Jacobian,
                           double rate_Jacobian, double phi_s, double phi_f,
                           double intrinsic_rho_f, double kappa_f,
                           double alpha_4, double V0) {

  double c0 = (rate_kichhoff_pressure * intrinsic_rho_f / kappa_f) *
              (phi_s - kichhoff_pressure * phi_f / (Jacobian * kappa_f));
  double c1 = rate_Jacobian * intrinsic_rho_f * kichhoff_pressure /
              (Jacobian * kappa_f);
  double c2 = intrinsic_rho_f * (rate_Jacobian + alpha_4 * Jacobian);
  double c3 = intrinsic_rho_f * Jacobian;

#if NumberDimensions == 2
  DG0_Du[0] =
      N_alpha_pw_n1 *
      ((c0 - c1 + c2) * dN_beta_u_n1[0] -
       c3 * (grad_v[0] * dN_beta_u_n1[0] + grad_v[1] * dN_beta_u_n1[1])) *
      V0;

  DG0_Du[1] =
      N_alpha_pw_n1 *
      ((c0 - c1 + c2) * dN_beta_u_n1[1] -
       c3 * (grad_v[2] * dN_beta_u_n1[0] + grad_v[3] * dN_beta_u_n1[1])) *
      V0;

#else
  DG0_Du[0] = N_alpha_pw_n1 *
              ((c0 - c1 + c2) * dN_beta_u_n1[0] -
               c3 * (grad_v[0] * dN_beta_u_n1[0] + grad_v[1] * dN_beta_u_n1[1] +
                     grad_v[2] * dN_beta_u_n1[2])) *
              V0;
  DG0_Du[1] = N_alpha_pw_n1 *
              ((c0 - c1 + c2) * dN_beta_u_n1[1] -
               c3 * (grad_v[3] * dN_beta_u_n1[0] + grad_v[4] * dN_beta_u_n1[1] +
                     grad_v[5] * dN_beta_u_n1[2])) *
              V0;
  DG0_Du[2] = N_alpha_pw_n1 *
              ((c0 - c1 + c2) * dN_beta_u_n1[2] -
               c3 * (grad_v[6] * dN_beta_u_n1[0] + grad_v[7] * dN_beta_u_n1[1] +
                     grad_v[8] * dN_beta_u_n1[2])) *
              V0;
#endif
}

/**************************************************************/

void compute_DG0_Dpw__upw__(double *DG0_Dpw, double N_alpha_pw_n1,
                            double N_beta_pw_n1, double rate_kichhoff_pressure,
                            double Jacobian, double rate_Jacobian, double phi_f,
                            double intrinsic_rho_f, double kappa_f, double alpha_4, double V0) {
  
  double relative_rho_f = intrinsic_rho_f*phi_f;

  double c0 = (rate_kichhoff_pressure * intrinsic_rho_f * phi_f) /
              (kappa_f * kappa_f * Jacobian);

  double c1 = relative_rho_f * alpha_4 / kappa_f;

  double c2 = rate_Jacobian * intrinsic_rho_f / (Jacobian * kappa_f);

  *DG0_Dpw = N_alpha_pw_n1 * (c0 + c1 + c2) * N_beta_pw_n1 * V0;
}

/**************************************************************/

void compute_G1__upw__(double *G1, const double *dN_alpha_pw_n1,
                       const double *K, const double *grad_kichhoff_pressure,
                       double g, double V0) {
#if NumberDimensions == 2
  *G1 = -(1.0 / g) *
        (dN_alpha_pw_n1[0] * (K[0] * grad_kichhoff_pressure[0] +
                              K[1] * grad_kichhoff_pressure[1]) +
         dN_alpha_pw_n1[1] * (K[2] * grad_kichhoff_pressure[0] +
                              K[3] * grad_kichhoff_pressure[1])) *
        V0;
#else
  *G1 = -(1.0 / g) *
        (dN_alpha_pw_n1[0] * (K[0] * grad_kichhoff_pressure[0] +
                              K[1] * grad_kichhoff_pressure[1] +
                              K[2] * grad_kichhoff_pressure[2]) +
         dN_alpha_pw_n1[1] * (K[3] * grad_kichhoff_pressure[0] +
                              K[4] * grad_kichhoff_pressure[1] +
                              K[5] * grad_kichhoff_pressure[2]) +
         dN_alpha_pw_n1[2] * (K[6] * grad_kichhoff_pressure[0] +
                              K[7] * grad_kichhoff_pressure[1] +
                              K[8] * grad_kichhoff_pressure[2])) *
        V0;
#endif
}

/**************************************************************/

void compute_DG1_Du__upw__(double *DG1_Du, const double *dN_alpha_pw_n1,
                           const double *dN_beta_u_n1,
                           const double *grad_kichhoff_pressure,
                           const double *K, double g, double V0) {
#if NumberDimensions == 2
  DG1_Du[0] = (1.0 / g) * ();

  dN_alpha_pw_n1[0] * dN_beta_u_n1[0] * +dN_alpha_pw_n1[0] * dN_beta_u_n1[1] *

      DG1_Du[1] = (1.0 / g) * ();
#else
  DG1_Du[0] = (1.0 / g) *;
  DG1_Du[1] = (1.0 / g) *;
  DG1_Du[2] = (1.0 / g) *;
#endif
}

/**************************************************************/

void compute_G2__upw__(double *G2, const double *dN_alpha_pw_n1,
                       const double *K, const double *b, double Jacobian,
                       double intrinsic_rho_f, double g, double V0) {
#if NumberDimensions == 2
  *G2 = -(Jacobian * intrinsic_rho_f / g) *
        (dN_alpha_pw_n1[0] * (K[0] * b[0] + K[1] * b[1]) +
         dN_alpha_pw_n1[1] * (K[2] * b[0] + K[3] * b[1])) *
        V0;
#else
  *G2 = -(Jacobian * intrinsic_rho_f / g) *
        (dN_alpha_pw_n1[0] * (K[0] * b[0] + K[1] * b[1] + K[2] * b[2]) +
         dN_alpha_pw_n1[1] * (K[3] * b[0] + K[4] * b[1] + K[5] * b[2]) +
         dN_alpha_pw_n1[2] * (K[6] * b[0] + K[7] * b[1] + K[8] * b[2])) *
        V0;
#endif
}

/**************************************************************/