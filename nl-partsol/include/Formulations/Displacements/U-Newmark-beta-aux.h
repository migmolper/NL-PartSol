// Global libs
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"

typedef struct {

#ifdef USE_PETSC
  Vec value;
  Vec d_value_dt;
  Vec d2_value_dt2;
#else
  double *value;
  double *d_value_dt;
  double *d2_value_dt2;
#endif

} Nodal_Field;

typedef struct {
  double alpha_1;
  double alpha_2;
  double alpha_3;
  double alpha_4;
  double alpha_5;
  double alpha_6;
  double epsilon;
  double DeltaTimeStep;
} Newmark_parameters;

static double __compute_deltat(Particle MPM_Mesh /**< */, double h /**< */,
                               Time_Int_Params Parameters_Solver /**< */);
/**************************************************************/

/*!
  \brief Finite strains Newmark-beta

  \param[in] beta: First Newmark-beta parameter
  \param[in] gamma: Second Newmark-beta parameter
  \param[in] DeltaTimeStep: Timestep increment
  \param[in] epsilon: Preconditioner parameter for the mass matrix
*/
static Newmark_parameters __compute_Newmark_parameters(double beta,
                                                       double gamma,
                                                       double DeltaTimeStep,
                                                       double epsilon);
/**************************************************************/

/*!
  \brief Compute the contribution of particle p to the lumped \n
  mass matrix correspoding to node A

  \param[out] Local_Mass_Matrix_p Local lumped Mass matrix
  \param[in] Na_p Shape function evaluation at node A
  \param[in] m_p Particle mass
*/
static void __compute_local_mass_matrix(double *Local_Mass_Matrix_p,
                                        double Na_p, double m_p);
/**************************************************************/

/*!
  \brief Returns a mask with the position of the dofs \n
  with contributions to the lumped mass matrix for node A \n
  and particle p.

  \param[out] Mask_active_dofs_A Global dofs positions for node A.
  \param[in] Mask_node_A Index of the node A with mask.
*/
static void __get_assembling_locations_lumped_mass(int *Mask_active_dofs_A,
                                                   int Mask_node_A);
/**************************************************************/

/*!
  \brief This function returns the lumped mass matrix \n
  in order to reduce storage, this matrix is presented as \n
  a vector.

  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[out] STATUS Returns failure or success
  \return The lumped mass matrix for the active system
*/
#ifdef USE_PETSC
static Vec __compute_nodal_lumped_mass(Particle MPM_Mesh, Mesh FEM_Mesh,
                                       Mask ActiveNodes, int *STATUS);
#else
static double *__compute_nodal_lumped_mass(Particle MPM_Mesh, Mesh FEM_Mesh,
                                           Mask ActiveNodes, int *STATUS);
#endif
/**************************************************************/

static void __local_projection_displacement(double *U_N_m_IP,
                                            const double *dis_p,
                                            double m__x__ShapeFunction_pA);
/**************************************************************/

static void __local_projection_velocity(double *V_N_m_IP, const double *vel_p,
                                        double m__x__ShapeFunction_pA);
/**************************************************************/

static void __local_projection_acceleration(double *A_N_m_IP,
                                            const double *acc_p,
                                            double m__x__ShapeFunction_pA);
/**************************************************************/

static void __get_assembling_locations_nodal_kinetics(int *Mask_active_dofs_A,
                                                      int Mask_node_A,
                                                      Mask ActiveDOFs);
/**************************************************************/

/*!
  \brief Project the kinematic information of the particles towards the nodes \n
  of the mesh to get the nodal fields at t = n
  \param[in] Lumped_Mass Lumped mass matrix used for the projection
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] Params Time integration parameters
  \param[out] STATUS Returns failure or success
  \return Nodal kinetics information from the step n
*/
#ifdef USE_PETSC
static Nodal_Field __get_nodal_field_tn(Vec Lumped_Mass, Particle MPM_Mesh,
                                        Mesh FEM_Mesh, Mask ActiveNodes,
                                        Mask ActiveDOFs,
                                        Newmark_parameters Params, int *STATUS);
#else
static Nodal_Field __get_nodal_field_tn(double *Lumped_Mass, Particle MPM_Mesh,
                                        Mesh FEM_Mesh, Mask ActiveNodes,
                                        Mask ActiveDOFs,
                                        Newmark_parameters Params, int *STATUS);
#endif
/**************************************************************/

/*!

*/
static Nodal_Field __initialise_nodal_increments(Nodal_Field U_n, Mesh FEM_Mesh,
                                                 Mask ActiveNodes,
                                                 Newmark_parameters Params,
                                                 int *STATUS);
/**************************************************************/

/*!
  \brief Update the local deformation of the particles \n
  ensuring the local compatibility conditions
  \param[in] D_U Increment of nodal kinetics
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[out] STATUS Returns failure or success
*/
static void __local_compatibility_conditions(const Nodal_Field D_U,
                                             Mask ActiveNodes,
                                             Particle MPM_Mesh, Mesh FEM_Mesh,
                                             int *STATUS);
/**************************************************************/

/*!
  \brief Update the stress tensor of the particle
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[out] STATUS Returns failure or success
*/
static void __constitutive_update(Particle MPM_Mesh, Mesh FEM_Mesh,
                                  int *STATUS);
/**************************************************************/

/**
 * @brief Function used to compute the equilibrium residual
 * 
 * @param [in] U_n Nodal kinetics information from the step n 
 * @param [in] D_U Increment of nodal kinetics
 * @param [in] Lumped_Mass Effective mass matrix (vector)
 * @param [in] ActiveNodes List of nodes which takes place in the computation
 * @param [in] ActiveDOFs List of dofs which takes place in the computation
 * @param [in] MPM_Mesh Information of the particles
 * @param [in] FEM_Mesh Information of the background nodes
 * @param [in] Params Time integration parameters
 * @param [in] TimeStep Current time step for time-dependent loads
 * @param [in] Is_compute_Residual The function computes the residual
 * @param [in] Is_compute_Reactions The function computes the reaction
 * @param [out] STATUS Returns failure or success
 * @return Vec/double* 
 */
#ifdef USE_PETSC
static Vec __assemble_residual(Nodal_Field U_n, Nodal_Field D_U,
                               Vec Lumped_Mass, Mask ActiveNodes,
                               Mask ActiveDOFs, Particle MPM_Mesh,
                               Mesh FEM_Mesh, Newmark_parameters Params,
                               unsigned TimeStep,
                               bool Is_compute_Residual,
                               bool Is_compute_Reactions, int *STATUS);
#else
static double *__assemble_residual(Nodal_Field U_n, Nodal_Field D_U,
                                   double *Lumped_Mass, Mask ActiveNodes,
                                   Mask ActiveDOFs, Particle MPM_Mesh,
                                   Mesh FEM_Mesh, Newmark_parameters Params,
                                   unsigned TimeStep,
                                   bool Is_compute_Residual,
                                   bool Is_compute_Reactions, int *STATUS);
#endif
/**************************************************************/

/*!
  \brief Returns a mask with the position of the residual values.

  \param[out] Mask_active_dofs_A Global dofs positions for node A.
  \param[in] Mask_node_A Index of the node A with mask.
  \param[in] ActiveDOFs List of dofs which takes place in the computation
*/
static void __get_assembling_locations_residual(int *Mask_active_dofs_A,
                                                int Mask_node_A,
                                                Mask ActiveDOFs);
/**************************************************************/

/*!
  \brief Returns a mask with the position of the reactions values.

  \param[out] Mask_active_dofs_A Global dofs positions for node A.
  \param[in] Mask_node_A Index of the node A with mask.
*/
static void __get_assembling_locations_reactions(int *Mask_active_dofs_A,
                                                 int Mask_node_A);
/**************************************************************/

/*!
  \brief Do the assembly process in the global residual/reaction vector. \n
  Replaces a PETSc function with the same name in case PETSc is not defined.

  \param[in,out] Residual Global residual/reaction vector
  \param[in] Local_Value Local contribution
  \param[in] Mask_active_dofs_A Global dofs positions for node A
*/
#ifndef USE_PETSC
static void VecSetValues(double *Residual, const double *Local_Value,
                         const int *Mask_active_dofs_A);
#endif
/**************************************************************/

/*!
  \brief Function used to compute the contribution of the \n
  internal forces to the residual

  \param[in,out] Residual Residual/reaction vector
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] Is_compute_Residual The function computes the residual
  \param[in] Is_compute_Reactions The function computes the reaction
  \param[out] STATUS Returns failure or success
*/
#ifdef USE_PETSC
static void __Nodal_Internal_Forces(Vec Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions, int *STATUS);
#else
static void __Nodal_Internal_Forces(double *Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions, int *STATUS);
#endif
/**************************************************************/

static void __internal_force_density(double *InternalForcesDensity_Ap,
                                     const double *kirchhoff_p,
                                     const double *gradient_n1_pA, double V0_p);
/**************************************************************/

/*!
  \brief Function used to compute the contribution of the \n
  contact forces to the residual

  \param[in,out] Residual Residual vector
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] Is_compute_Residual The function computes the residual
  \param[in] Is_compute_Reactions The function computes the reaction
*/
#ifdef USE_PETSC
static void __Nodal_Traction_Forces(Vec Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions);
#else
static void __Nodal_Traction_Forces(double *Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions);
#endif
/**************************************************************/

static void __local_traction_force(double *LocalTractionForce_Ap,
                                   const double *T, double N_n1_pA,
                                   double A0_p);
/**************************************************************/


/**
 * @brief Function used to compute the contribution of the \n 
 * body (distance) forces to the residual
 * 
 * @param[in,out] Residual Residual vector
 * @param[in] ActiveNodes List of nodes which takes place in the computation
 * @param[in] ActiveDOFs List of dofs which takes place in the computation
 * @param[in] MPM_Mesh Information of the particles
 * @param[in] FEM_Mesh Information of the background nodes 
 * @param TimeStep 
 * @param[in] Is_compute_Residual The function computes the residual
 * @param[in] Is_compute_Reactions The function computes the reaction
 */
#ifdef USE_PETSC
static int __Nodal_Body_Forces(Vec Residual, Mask ActiveNodes, Mask ActiveDOFs,
                                Particle MPM_Mesh, Mesh FEM_Mesh,
                                unsigned TimeStep,
                                bool Is_compute_Residual,
                                bool Is_compute_Reactions);
#else
static int __Nodal_Body_Forces(double *Residual, Mask ActiveNodes,
                                Mask ActiveDOFs, Particle MPM_Mesh,
                                Mesh FEM_Mesh, 
                                unsigned TimeStep,
                                bool Is_compute_Residual,
                                bool Is_compute_Reactions);
#endif
/**************************************************************/

static void __local_body_force(double *LocalBodyForce_Ap, const double *b,
                               double N_n1_pA, double m_p);

/**************************************************************/

/*!
  \brief Function used to compute the contribution of the \n
  inertial forces to the residual

  \param[in,out] Residual Residual vector
  \param[in] Mass Mass matrix
  \param[in] U_n Nodal kinetics information from the step n
  \param[in] D_U Increment of nodal kinetics
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] Params Time integration parameters
*/
#ifdef USE_PETSC
static void __Nodal_Inertial_Forces(Vec Residual, Vec Mass, Nodal_Field U_n,
                                    Nodal_Field D_U, Mask ActiveNodes,
                                    Mask ActiveDOFs, Newmark_parameters Params);
#else
static void __Nodal_Inertial_Forces(double *Residual, double *Mass,
                                    Nodal_Field U_n, Nodal_Field D_U,
                                    Mask ActiveNodes, Mask ActiveDOFs,
                                    Newmark_parameters Params);
#endif
/**************************************************************/

/*!
  \brief Compute the 2-norm of the residual
  \param[in] Residual
  \param[in] Total_dof Size of the vector
  \return The norm of the residual
*/
#ifdef USE_PETSC
static double __error_residual(const Vec Residual, unsigned Total_dof);
#else
static double __error_residual(const double *Residual, unsigned Total_dof);
#endif
/**************************************************************/

/*!

*/
#ifdef USE_PETSC
static Mat __preallocation_tangent_matrix(Mask ActiveNodes, Mask ActiveDOFs,
                                          Particle MPM_Mesh, int *STATUS);
#else
static double *__preallocation_tangent_matrix(Mask ActiveNodes, Mask ActiveDOFs,
                                              Particle MPM_Mesh, int *STATUS);
#endif
/**************************************************************/

static void compute_local_intertia(double *Inertia_density_p /**< */,
                                   double Na_p /**< */, double Nb_p /**< */,
                                   double m_p /**< */, double alpha_1 /**< */,
                                   double epsilon /**< */, unsigned A /**< */,
                                   unsigned B /**< */);
/**************************************************************/

/*!
  \brief Returns a mask with the position of the stifness in the global tangent
  \n matrix for an eeficient assembly process.

  \param[out] Mask_active_dofs_A Global dofs positions for node A.
  \param[in] Mask_node_A Index of the node A with mask.
  \param[out] Mask_active_dofs_B Global dofs positions for node B
  \param[in] Mask_node_B Index of the node B with mask
  \param[in] ActiveDOFs List of dofs which takes place in the computation
*/
static void __get_tangent_matrix_assembling_locations(int *Mask_active_dofs_A,
                                                      int Mask_node_A,
                                                      int *Mask_active_dofs_B,
                                                      int Mask_node_B,
                                                      Mask ActiveDOFs);
/**************************************************************/

/*!
  \brief Computes the local tangent stiffness as the addition of the \n
  local tangent stiffness density and the intertia density matrix

  \param[out] Tangent_Stiffness_p Local tangent stiffness
  \param[in] Stiffness_density_p Local stiffness density
  \param[in] Inertia_density_p Inertia density matrix, a.k.a mass matrix
  \param[in] V0_p Particle volume
*/
static void local_tangent_stiffness(double *Tangent_Stiffness_p,
                                    const double *Stiffness_density_p,
                                    const double *Inertia_density_p,
                                    double V0_p);
/**************************************************************/

/*!
  \brief Do the assembly process in the global tangent matrix. \n
  Replaces a PETSc function with the same name in case PETSc is not defined.

  \param[in,out] Tangent_Stiffness Global tangent stiffness
  \param[in] Tangent_Stiffness_p Local tangent stiffness
  \param[in] Mask_active_dofs_A Global dofs positions for node A
  \param[in] Mask_active_dofs_B Global dofs positions for node B
  \param[in] Nactivedofs Number of active dofs in the Global tangent matrix
*/
#ifndef USE_PETSC
static void MatSetValues(double *Tangent_Stiffness,
                         const double *Tangent_Stiffness_p,
                         const int *Mask_active_dofs_A,
                         const int *Mask_active_dofs_B, int Nactivedofs);
#endif
/**************************************************************/

/*!
  \param[in,out] Tangent_Stiffness The tangent matrix for the problem
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] Params Time integration parameters
  \param[in] Iter Current iteration of the solver
  \param[out] STATUS Returns failure or success
*/
#ifdef USE_PETSC
static void __assemble_tangent_stiffness(Mat Tangent_Stiffness,
                                         Mask ActiveNodes, Mask ActiveDOFs,
                                         Particle MPM_Mesh, Mesh FEM_Mesh,
                                         Newmark_parameters Params,
                                         unsigned Iter, int *STATUS);
#else
static void __assemble_tangent_stiffness(double *Tangent_Stiffness,
                                         Mask ActiveNodes, Mask ActiveDOFs,
                                         Particle MPM_Mesh, Mesh FEM_Mesh,
                                         Newmark_parameters Params,
                                         unsigned Iter, int *STATUS);
#endif
/**************************************************************/

/*!
  \brief This function takes the nodal increments comming from \n
  the linear solver and computes the updated nodal increments \n
  for the nodal kinetics
  \param[in] Residual Increments during the k iteration
  \param[in,out] D_U Nodal incremented (displacement,velocity,acceleration)
  \param[in] U_n Previously converged nodal kinetics
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] Params Time integration parameters
  \param[in] Ntotaldofs Number of dofs
*/
#ifdef USE_PETSC
static void __update_Nodal_Increments(const Vec Residual, Nodal_Field D_U,
                                      Nodal_Field U_n, Mask ActiveDOFs,
                                      Newmark_parameters Params,
                                      unsigned Ntotaldofs);
#else
static void __update_Nodal_Increments(const double *Residual, Nodal_Field D_U,
                                      Nodal_Field U_n, Mask ActiveDOFs,
                                      Newmark_parameters Params,
                                      unsigned Ntotaldofs);
#endif
/**************************************************************/

static void __update_Particles(Nodal_Field D_U /**< */,
                               Particle MPM_Mesh /**< */, Mesh FEM_Mesh /**< */,
                               Mask ActiveNodes /**< */);