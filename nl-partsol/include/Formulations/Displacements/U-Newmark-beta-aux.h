// Global libs
#include "Macros.h"
#include "Matlib.h"
#include "Particles.h"
#include "Types.h"

typedef struct {

  Vec value;
  Vec d_value_dt;
  Vec d2_value_dt2;

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
static Vec __compute_nodal_lumped_mass(Particle MPM_Mesh, Mesh FEM_Mesh,
                                       Mask ActiveNodes, int *STATUS);
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

static void __local_projection_increment_displacement(
    double *DU_N_m_IP, const double *vel_p, const double *acc_p,
    double m__x__ShapeFunction_pA, double DeltaTimeStep);

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
static Nodal_Field __get_nodal_field_tn(Vec Lumped_Mass, Particle MPM_Mesh,
                                        Mesh FEM_Mesh, Mask ActiveNodes,
                                        Mask ActiveDOFs,
                                        Newmark_parameters Params, int *STATUS);
/**************************************************************/

/**
 * @brief
 *
 * @param Lumped_Mass
 * @param MPM_Mesh
 * @param FEM_Mesh
 * @param U_n
 * @param ActiveNodes
 * @param ActiveDOFs
 * @param Params
 * @param STATUS
 * @return Nodal_Field
 */
static Nodal_Field __initialise_nodal_increments(
    Vec Lumped_Mass, Particle MPM_Mesh, Mesh FEM_Mesh, Nodal_Field U_n,
    Mask ActiveNodes, Mask ActiveDOFs, Newmark_parameters Params, int *STATUS);
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
static int __constitutive_update(Particle MPM_Mesh, Mesh FEM_Mesh);
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
 * @param [in] Is_compute_Residual The function computes the residual
 * @param [in] Is_compute_Reactions The function computes the reaction
 * @param [out] STATUS Returns failure or success
 * @return Vec/double*
 */
static Vec __assemble_residual(Nodal_Field U_n, Nodal_Field D_U,
                               Vec Lumped_Mass, Mask ActiveNodes,
                               Mask ActiveDOFs, Particle MPM_Mesh,
                               Mesh FEM_Mesh, Newmark_parameters Params,
                               bool Is_compute_Residual,
                               bool Is_compute_Reactions, int *STATUS);
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
static void __Nodal_Internal_Forces(Vec Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions, int *STATUS);
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
static void __Nodal_Traction_Forces(Vec Residual, Mask ActiveNodes,
                                    Mask ActiveDOFs, Particle MPM_Mesh,
                                    Mesh FEM_Mesh, bool Is_compute_Residual,
                                    bool Is_compute_Reactions);
/**************************************************************/

static void __local_traction_force(double *LocalTractionForce_Ap,
                                   const double *T, double N_n1_pA,
                                   double A0_p);
/**************************************************************/

/**
 * @brief Function used to compute the contribution of the \n
 * inertial forces to the residual (a - b)
 *
 * @param[in,out] Residual Residual vector
 * @param[in] Mass Mass matrix
 * @param[in] U_n Nodal kinetics information from the step n
 * @param[in] D_U Increment of nodal kinetics
 * @param[in] ActiveNodes List of nodes which takes place in the computation
 * @param[in] ActiveDOFs List of dofs which takes place in the computation
 * @param[in] Params Time integration parameters
 * @param[in] Is_compute_Residual The function computes the residual
 * @param[in] Is_compute_Reactions The function computes the reaction
 */
static void __Nodal_Inertial_Forces(Vec Residual, Vec Mass, Nodal_Field U_n,
                                    Nodal_Field D_U, Mask ActiveNodes,
                                    Mask ActiveDOFs, Newmark_parameters Params,
                                    bool Is_compute_Residual,
                                    bool Is_compute_Reactions);
/**************************************************************/

/*!
  \brief Compute the 2-norm of the residual
  \param[in] Residual
  \return The norm of the residual
*/
static double __error_residual(const Vec Residual);
/**************************************************************/

/*!

*/
static Mat __preallocation_tangent_matrix(Mask ActiveNodes, Mask ActiveDOFs,
                                          Particle MPM_Mesh, int *STATUS);
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
  \param[in,out] Tangent_Stiffness The tangent matrix for the problem
  \param[in] ActiveNodes List of nodes which takes place in the computation
  \param[in] ActiveDOFs List of dofs which takes place in the computation
  \param[in] MPM_Mesh Information of the particles
  \param[in] FEM_Mesh Information of the background nodes
  \param[in] Params Time integration parameters
  \param[in] Iter Current iteration of the solver
  \param[out] STATUS Returns failure or success
*/
static void __assemble_tangent_stiffness(Mat Tangent_Stiffness,
                                         Mask ActiveNodes, Mask ActiveDOFs,
                                         Particle MPM_Mesh, Mesh FEM_Mesh,
                                         Newmark_parameters Params,
                                         unsigned Iter, int *STATUS);
/**************************************************************/

/**
 * @brief 
 * 
 * @param D_U Nodal incremented (displacement,velocity,acceleration)
 * @param U_n Previously converged nodal kinetics
 * @param ActiveDOFs List of dofs which takes place in the computation
 * @param Params Time integration parameters
 * @param Ntotaldofs Number of dofs
 */
static void __trial_Nodal_Increments(Nodal_Field D_U,
                                      Nodal_Field U_n, Mask ActiveDOFs,
                                      Newmark_parameters Params,
                                      unsigned Ntotaldofs);

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
static void __update_Nodal_Increments(const Vec Residual, Nodal_Field D_U,
                                      Nodal_Field U_n, Mask ActiveDOFs,
                                      Newmark_parameters Params,
                                      unsigned Ntotaldofs);
/**************************************************************/

static void __update_Particles(Nodal_Field D_U /**< */,
                               Particle MPM_Mesh /**< */, Mesh FEM_Mesh /**< */,
                               Mask ActiveNodes /**< */);