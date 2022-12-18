#ifndef _TYPES_H_
#define _TYPES_H_

#include <stdbool.h>
#include "Macros.h"

/** @brief \file Types.h
 *  \brief File with the definitions of the user-defined variables.
 */

/*******************************************************/

/** @brief \struct Matrix
 * This structure is devoted to store in matricial
 * the allocated memory
 */
typedef struct {

  /** @brief
   * Number of rows
   */
  int N_rows;

  /** @brief
   * Number of columns
   */
  int N_cols;

  /** @brief
   * Value Is an scalar
   */
  double n;

  /** @brief
   * Pointer for a vector
   */
  double *nV;

  /** @brief
   * Table of pointers for a matrix
   */
  double **nM;

  /** @brief
   * Aditional information
   */
  char Info[100];

} Matrix;

/*******************************************************/

/** @brief \struct Chain
 * This structure is devoted to define
 * each component of a set
 */
typedef struct Chain {

  /** @brief
   * Index of a node in the set
   */
  int Idx;

  /** @brief
   * Pointer to the next node in the set
   */
  struct Chain *next;

} Chain;

/** @brief Pointer to the header of the set
 */
typedef Chain *ChainPtr;

/*******************************************************/

/** @brief Table
 * This structure is devoted to store a table with
 * values of the type double.
 */
typedef struct {

  /** @brief
   * Number of rows
   */
  int N_rows;

  /** @brief
   * Number of columns
   */
  int N_cols;

  /** @brief
   * Scalar
   */
  int n;

  /** @brief
   * 1D list
   */
  int *nV;

  /** @brief
   * 2D list
   */
  int **nM;

  /** @brief
   * Aditional information
   */
  char Info[100];

} Table;

/*******************************************************/

/** @brief \struct Curve
 *  Curve definition
 */
typedef struct {

  /** @brief
   * Number of items in the curve */
  int Num;

  /** @brief
   * Values for each time */
  double *Fx;

  /** @brief
   * Aditional information */
  char Info[100];

} Curve;

/*******************************************************/

/** @brief \struct Tensor
 * This structure is devoted
 * to help the code to deal with high level functions
 * for the B-free approach
 */
typedef struct {

  /** @brief
   * Order of the tensor
   */
  int Order;

  /** @brief
   * First order tensor
   */
  double *n;

  /** @brief
   * Second order tensor
   */
  double *N[NumberDimensions];

  /** @brief
   * Aditional information
   */
  char Info[100];

} Tensor;

/*******************************************************/

/** @brief \struct EigenTensor
  Structure for the output of a function wich compute eigenvalues
    and eigenvectors
*/
typedef struct {

  Tensor Value;

  Tensor Vector;

} EigenTensor;

/*******************************************************/

/** @brief \struct Fields
 * Variable devoted to store in the memory the physicl
 * information of each particle
 */
typedef struct {

  /**
   * @brief Reference volume
   *
   */
  Matrix Vol_0;

  /**
   * @brief Reference area
   *
   */
  Matrix Area_0;

  Matrix rho; /** @brief Density field */

  Matrix rho_s, rho_f; /** @brief Intrinsic density field (solid/fluid) */

  Matrix phi_s, phi_f; /** @brief Volume fraction (solid/fluid) */

  Matrix mass; /** @brief Mass field */

  Matrix x_GC; /** @brief Position in global coordinates */

  Matrix x_EC; /** @brief Position in element coordinates */

  Matrix dis, D_dis; /** @brief Displacement field */

  Matrix vel; /** @brief Velocity field */

  Matrix acc; /** @brief Acceleration fields */

  Matrix Stress; /** @brief Stress field */

  /** @brief
   * Lagrange multiplier for incompressible formulations
   * */
  Matrix lambda_pressure_n;
  Matrix lambda_pressure_n1;

  /** @brief
   * Pore water pressure, initial state, and the increment
   */
  Matrix Pw;
  Matrix d_Pw;
  Matrix Pw_0;
  Matrix Pw_n1;
  Matrix D_Pw;

  /** @brief
   * Rates of the pore water pressure
   */
  Matrix d_Pw_dt_n;
  Matrix d_Pw_dt_n1;

  Matrix d2_Pw_dt2;

  Matrix Strain; /** @brief Strain field */

  Matrix F_n, F_n1; /** @brief Total deformation gradient */

  Matrix DF; /** @brief Incremental deformation gradient */

  Matrix dt_F_n, dt_F_n1, dt_DF; /** @brief Rate deformation gradient */

  double *J_n, *J_n1, *dJ_dt; /** @brief Jacobian of the deformation gradient */

  Matrix Fbar; /** @brief Locking-free deformation gradient (F-bar) */

  Matrix Jbar; /** @brief Jacobian of the F-bar */

  double *Temperature_n, *Temperature_n1; /** @brief Particle temperature */

  double *W; /** @brief Deformation Energy */

  double *Damage_n, *Damage_n1; /** @brief Damage parameter (Fracture) */

  double *Strain_f_n, *Strain_f_n1; /** @brief Strain during crack (Fracture) */

  double *EPS_n, *EPS_n1; /** @brief Equivalent plastic strain. */

  double *Kappa_n,
      *Kappa_n1; /** @brief Isotropic hardening (stress like) variable. */

  Matrix Back_stress; /** @brief Kinematic hardening variable. */

  Matrix b_e_n, b_e_n1; /** @brief Elastic left deformation gradient */

  Matrix C_ep; /** @brief Elastoplastic tangent matrix */

  bool *Status_particle; /** @brief Check if the particle is consider failed or
                            not */

#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  Matrix PU; /** @brief Partition of unity property. */

#endif
#endif

} Fields;

/*******************************************************/

/** @brief \struct Load
 *   An important aclaration about this code, a force is a load,
 * or even a boundary condition is defined as a load from the
 * point of view of the code (as a structure). The idea is
 * a force and a boundary condition both of them has a direction,
 * a value, a list of nodes/GPs where it is applied,... So it's
 * a good idea to use this structure for both porpouses.
 */

typedef struct {

  /** @brief
   * Number of nodes/GP with this load
   */
  int NumNodes;

  /** @brief
   * Number of dimensions of the load
   */
  int Dim;

  /** @brief
   * Direction of the load {0,0} {1,0} {0,1} {1,1}
   */
  int *Dir;

  /** @brief
   * List of nodes with this load
   */
  int *Nodes;

  /** @brief
   * Curve for each dimension with the evolution
   * value with the time
   */
  Curve *Value;

  // On/Off variable
  bool STATUS;

  /** @brief
   * Some information about this load
   */
  char Info[100];

} Load;

/*******************************************************/

/** @brief \struct Boundaries
 * Boundary conditions definition
 */
typedef struct {

  /** @brief
   * Number of boundaries of the domain
   */
  int NumBounds;

  /** @brief
   * Table with all the boundaries and its values
   */
  Load *BCC_i;

} Boundaries;

/*******************************************************/

/**
 * @brief Properties of a general material model
 *
 */
typedef struct {

  int Id; /** @brief Index of the material */

  char Type[100]; /** @brief Name of the material */

  double Cel; /** @brief Material celerity */

  double rho; /** @brief Initial density (mixture/fluid/solid) */

  double E; /** @brief Elastic modulus */

  double nu; /** @brief Poisson ratio */

  double Compressibility; /** @brief Compressibility */

  /** @brief
   * Fluid parameters
   */
  double ReferencePressure;
  double Viscosity;
  double n_Macdonald_model;

  double Ceps; /** @brief Normalizing constant (Eigenerosion/Eigensoftening) */

  double Gf; /** @brief Failure energy (Eigenerosion) */

  double ft; /** @brief Tensile strengt of the material (Eigensoftening) */

  double heps; /** @brief Bandwidth of the cohesive fracture (Eigensoftening) */

  double wcrit; /** @brief Critical opening displacement (Eigensoftening) */

  char Plastic_Solver[100]; /** @brief Integration algorithm for plasticity */

  /** @brief
   * General plastic parameters
   */
  double kappa_0;
  double Hardening_modulus;
  double Plastic_Strain_0;
  double atmospheric_pressure;
  double J2_degradated;
  double Cohesion;

  /** @brief
   * Frictional material (Borja et al. 2003)
   * */
  char Yield_Function_Frictional[100];
  double m_Frictional;
  double c0_Frictional;
  double phi_Frictional;
  double psi_Frictional;

  /** @brief
   * Hardening Hughes
   */
  bool Hardening_Hughes;
  double Parameter_Hardening_Hughes;

  /** @brief
   * Hardening Cervera
   * */
  bool Hardening_Cervera;

  /** @brief
   * Hardening Ortiz
   * */
  bool Hardening_Ortiz;
  double Exponent_Hardening_Ortiz;

  /** @brief
   * Hardening Voce
   * */
  bool Hardening_Voce;
  double K_0_Hardening_Voce;
  double K_inf_Hardening_Voce;
  double delta_Hardening_Voce;
  double theta_Hardening_Voce;

  /** @brief
   * Hardening Borja et al. 2003
   * */
  bool Hardening_Borja;
  double a_Hardening_Borja[3];
  double alpha_Hardening_Borja;

  /** @brief
   * Viscoplasticity parameters
   * */
  bool Viscous_regularization;
  double fluidity_param;

  /** @brief
   * Activate auxiliar techniques
   * */
  bool Locking_Control_Fbar;
  double alpha_Fbar;

} Material;

/*******************************************************/

/** @brief \struct Mixture
 * This structure is devoted to store information
 * for a general kind of mixtures
 */
typedef struct {

  /** @brief
   * Index of the mixture
   */
  int Id;

  /** @brief
   * Name of the mixture
   */
  char Type[100];

  /*
    Index for the constitutive description of each phase
    for the soil-water mixture
  */
  int Soil_Idx;
  int Water_Idx;

  /** @brief
   * Permeability of the soil skleleton
   */
  Tensor Permeability;

  /** @brief
   * Initial volume fractions (Soil/Water)
   */
  double phi_s_0;
  double phi_f_0;

} Mixture;

/*******************************************************/

/** @brief
 * \struct State_Parameters
 */
typedef struct {
  /** @brief
   * Particle identifier
   * */
  int Particle_Idx;

  /** @brief
   * Stress/strain parameters
   * */
  double *Stress;
  double *Strain;
  double Pressure;
  double *W; /** @brief Deformation energy */

  double *D_phi_n1; /** @brief Total deformation gradient (t = n + 1) */
  double *D_phi_n;  /** @brief Total deformation gradient (t = n) */
  double *d_phi;    /** @brief Incremental deformation gradient */
  double *b_e;      /** @brief Elastic left Cauchy-Green tensor */
  double *Fbar;
  double J; /**<z Jocobian */

  double *dFdt;         /** @brief Rate of the deformation gradient */
  double DeltaTimeStep; /** @brief Increment of the time step */
  double alpha_4;       /** @brief Time integation paramter (Newmark-beta) */

  double *Back_stress; /** @brief State varible for kinematic hardening*/
  double *Kappa;       /** @brief  Hardening Parameter */
  double *EPS;         /** @brief  Equivalent plastic strain */
  double *C_ep;        /** @brief  Elastoplastic tangent matrix */
  bool compute_C_ep;

  double Cohesion;
  double Yield_stress;

  // Failed material point
  bool *Failure;

} State_Parameters;

/*******************************************************/

/** @brief \struct Particle
 * This structure is devoted to store all the information
 * of a list of particles
 */
typedef struct {

  /** @brief
   * Number of particles
   * */
  int NumGP;

  /** @brief
   * Index with the closest node to each particle
   * */
  int *I0;

  /** @brief
   * Index of the element
   * */
  int *Element_p;

  /** @brief
   * Tributary nodes variables
   * */
  int *NumberNodes;
  ChainPtr *ListNodes;

  /** @brief
   * Set of particles close to each particle
   * */
  ChainPtr *Beps;

  /** @brief
   * Store the values of each field in the current time step
   * */
  Fields Phi;

  /** @brief
   * Values from the previous step
   * */
  Fields Phi_n0;

  /** @brief
   * Material variables
   * */
  int NumberMaterials;
  int *MatIdx;
  int *MixtIdx;
  Material *Mat;

  /** @brief
   * Neumann boundary conditions
   * */
  int NumNeumannBC;
  Load *F;

  /** @brief
   * Structure to store Neumann boundary conditions will replace NumNeumannBC
   * and F;
   * */
  Boundaries Neumann_Contours;

  /** @brief
   * uGIMP shape function parameter
   * */
  Matrix lp;

  /** @brief
   * LME shape function parameters:
   * */
  Matrix lambda; // Lagrange multiplier
  Matrix Beta;   // Thermalization or regularization parameter (scalar/matrix)
  Matrix Cut_off_Ellipsoid;

  /** @brief
   * Function to compute the stress state of the particle
   * */
  State_Parameters (*constitutive)(State_Parameters, Material);

} Particle;

/*******************************************************/

/** @brief \struct Mesh
 * This structure is devoted to store all the information
 * of a list of nodes
 */
typedef struct {

  /** @brief
   * Number of nodes in the mesh
   */
  int NumNodesMesh;

  /** @brief
   * Number of elements in the mesh
   */
  int NumElemMesh;

  /** @brief
   * Table with the coordinates of the nodes of the mesh
   */
  Matrix Coordinates;

  /** @brief
   * Number of nodes in a element
   */
  int *NumNodesElem;

  /** @brief
   * List of nodes for each element (Connectivity)
   */
  ChainPtr *Connectivity;

  /** @brief
   * Number of elements close to each node
   */
  int *NumNeighbour;

  /** @brief
   * List of elements close to each node
   */
  ChainPtr *NodeNeighbour;

  /** @brief
   * Number of nodes close to a node
   */
  int *SizeNodalLocality_0;
  int *SizeNodalLocality;

  /** @brief
   * List of nodes close to a node
   */
  ChainPtr *NodalLocality_0;
  ChainPtr *NodalLocality;

  /** @brief
   * Defines if a node is activated or not
   */
  bool *ActiveNode;

  /** @brief
   * Defines if a node belongs to a boundary or not
   */
  bool *BoundaryNode;

  /** @brief
   * List of boundaries of the domain
   */
  Boundaries Bounds;

  /** @brief
   * Number of dimensions of the element (OLD)
   */
  int Dimension;

  /** @brief
   * Minimum distance between nodes
   */
  double DeltaX;

  /** @brief
   * Nodal spacing
   */
  double *h_avg;

  /** @brief
   * Name of the element (OLD)
   */
  char TypeElem[20];

  /** @brief
   * Function with the interpolation technique
   */
  Matrix (*N_ref)(Matrix);

  /** @brief
   * Function with the gradient of the interpolation technique
   */
  Matrix (*dNdX_ref)(Matrix);

  /** @brief
   * Function with the gradient of the interpolation technique
   */
  Matrix (*dNdX)(Matrix, Matrix);

  /*
   * Function of compute natural coordinate
   */
  //  void (* X_to_Xi)(Matrix, Matrix, Matrix);

  /* * @brief
   * Function to compute the volume of an element
   */
  double (*volume_Element)(Matrix);

  /** @brief
   * Function to check if a point is inside or outside of a elemnt
   */
  bool (*In_Out_Element)(Matrix, Matrix);

  /** @brief
   * List of particles adjacent to a node
   */
  int *Num_Particles_Node;
  ChainPtr *List_Particles_Node;

  /** @brief
   * Variables and function for F-bar calculation
   * */
  bool Locking_Control_Fbar;
  int Num_Patch_Mesh;
  int *Idx_Patch;
  double *Vol_Patch_n;
  double *Vol_Patch_n1;

} Mesh;

/*******************************************************/

/** @brief \struct Mask
 *  Structure with the current "element" of the particle
 */
typedef struct {

  int Nactivenodes;

  int *Nodes2Mask;

} Mask;

/*******************************************************/

/** @brief \struct Element
 *  Structure with the current "element" of the particle
 */
typedef struct {

  /** @brief
   * Index of the particle
   */
  int i_GP;

  /** @brief
   * Number of nodes close to the particle
   */
  int NumberNodes;

  /** @brief
   * List of tributary nodes for the particle
   */
  int *Connectivity;

} Element;

/*******************************************************/

/** @brief \struct Time_Int_Params
 * Parameters of the time integration scheme
 */
typedef struct {

  /** @brief
   * Courant number
   */
  double CFL;

  /** @brief
   * Material celerity
   */
  double Cel;

  /** @brief
   * Initial step
   */
  int InitialTimeStep;

  /** @brief
   * Número de pasos de tiempo
   */
  int NumTimeStep;

  /** @brief
   * Simulation time
   */
  double FinalTime;

  /** @brief
   * Parameter to generate a mass matrix
   */
  double epsilon_Mass_Matrix;

  /** @brief
   * Conserving Energy-Momentum parameters
   */
  double TOL_Conserving_Energy_Momentum;

  /** @brief
   * Generalized alpha parameters
   */
  double rb_Generalized_alpha;
  double TOL_Generalized_alpha;

  /** @brief
   * Newmark parameters
   */
  double beta_Newmark_beta;
  double gamma_Newmark_beta;
  double TOL_Newmark_beta;
  bool Use_explicit_trial;

  /** @brief
   * Maximum number of interations
   */
  int MaxIter;

  /** @brief
   * Time integration scheme
   */
  char TimeIntegrationScheme[100];

} Time_Int_Params;

/*******************************************************/

/** @brief \struct Event
  Structure with output control
 */
typedef struct {

  /** @brief
    Number of files of the output
  */
  int NumFiles;

  /** @brief
    Index for the particle or node in the case of single
    node/particle analysis
  */
  int Idx;

  /** @brief
    List of nodes or particles to get the information from
    a path of nodes/particles.
  */
  int *Idx_Path;
  int Lenght_Path;

  /** @brief
    Physical time to start the event. Default = 0.0
   */
  double start;

  /** @brief
    Numerical step to start the event. Default = 0
   */
  int i_start;

  /** @brief
    Physical time step
   */
  double step;

  /** @brief
    Numerical time step
   */
  int i_step;

  /** @brief
    Physical time to finish the event
   */
  double end;

  /** @brief
    Numerical step to finish the event
   */
  int i_end;

  /** @brief
    Name of the output directory
   */
  char Directory[MAXC];

  /** @brief
    Output selector CSV
  */
  bool Out_csv_nodes_path_Velocity;
  bool Out_csv_nodes_path_Acceleration;
  bool Out_csv_nodes_path_D_Displacement;
  bool Out_csv_nodes_path_Forces;
  bool Out_csv_nodes_path_Reactions;
  bool Out_csv_nodes_path_Residual;

  bool Out_csv_particles_path_Damage;
  bool Out_csv_particles_path_Velocity;
  bool Out_csv_particles_path_Acceleration;
  bool Out_csv_particles_path_Displacement;
  bool Out_csv_particles_path_Stress;
  bool Out_csv_particles_path_Strain;
  bool Out_csv_particles_path_Deformation_gradient;

  bool Out_csv_Gauss_Point_evolution_Stress;
  bool Out_csv_Gauss_Point_evolution_Strain;
  bool Out_csv_Gauss_Point_evolution_Deformation_gradient;
  bool Out_csv_Gauss_Point_evolution_Plastic_Deformation_gradient;
  bool Out_csv_Gauss_Point_evolution_EPS;
  bool Out_csv_Gauss_Point_evolution_Cohesion;

  /** @brief
    Output selector Vtk
  */
  bool Out_vtk_nodes;
  bool Out_vtk_global_coordinates;
  bool Out_vtk_element_coordinates;
  bool Out_vtk_mass;
  bool Out_vtk_density;
  bool Out_vtk_damage;
  bool Out_vtk_nodal_idx;
  bool Out_vtk_material_idx;
  bool Out_vtk_velocity;
  bool Out_vtk_acceleration;
  bool Out_vtk_displacement;
  bool Out_vtk_stress;
  bool Out_vtk_eigenvalues_stress;
  bool Out_vtk_volumetric_stress;
  bool Out_vtk_strain;
  bool Out_vtk_eigenvalues_strain;
  bool Out_vtk_deformation_gradient;
  bool Out_vtk_energy;
  bool Out_vtk_Von_Mises;

} Event;

/*******************************************************/

#endif
