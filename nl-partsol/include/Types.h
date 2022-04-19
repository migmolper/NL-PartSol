#ifndef _TYPES_H_
#define _TYPES_H_


/*! \file Types.h
 *  \brief File with the definitions of the user-defined variables.    
 */

/*******************************************************/

/*! \struct Matrix
 * This structure is devoted to store in matricial
 * the allocated memory
 */
typedef struct {
  
  /*!
   * Number of rows 
   */
  int N_rows;

  /*! 
   * Number of columns 
   */
  int N_cols;

  /*! 
   * Value Is an scalar 
   */
  double n;
  
  /*!
   * Pointer for a vector 
   */
  double * nV;

  /*!
   * Table of pointers for a matrix
   */
  double ** nM;

  /*!
   * Aditional information 
   */
  char Info [100]; 

} Matrix;

/*******************************************************/

/*! \struct Chain
 * This structure is devoted to define
 * each component of a set
 */
typedef struct Chain {
  
  /*! 
   * Index of a node in the set 
   */
  int Idx;

  /*!
   * Pointer to the next node in the set
   */
  struct Chain * next;

  /*!
    
  */
  
} Chain; 

/*! \struct ChainPtr
 * Pointer to the header of the set
 */
typedef Chain * ChainPtr;

/*******************************************************/

/*! \struct Table
 * This structure is devoted to store a table with
 * values of the type double.
 */
typedef struct{

  /*!
   * Number of rows 
   */
  int N_rows;

  /*!
   * Number of columns 
   */
  int N_cols;

  /*!
   * Scalar
   */
  int n;

  /*!
   * 1D list 
   */
  int * nV;

  /*!
   * 2D list 
   */
  int ** nM;

  /*!
   * Aditional information 
   */
  char Info [100];
  
} Table;

/*******************************************************/

/*! \struct Curve
 *  Curve definition
 */
typedef struct{

  /*! Number of items in the curve */
  int Num; 
  /*! Values for each time */
  double * Fx; 
  /*! Aditional information */
  char Info [100];
  
} Curve;

/*******************************************************/

/*! \struct Tensor
 * This structure is devoted
 * to help the code to deal with high level functions
 * for the B-free approach
 */
typedef struct{
  
  /*!
   * Order of the tensor 
   */
  int Order;

  /*! 
   * First order tensor 
   */
  double *n;

  /*! 
   * Second order tensor 
   */
  double *N[NumberDimensions];

  /*! 
   * Aditional information 
   */
  char Info [100];
  
} Tensor;

/*******************************************************/

/*! \struct EigenTensor
  Structure for the output of a function wich compute eigenvalues 
    and eigenvectors
*/
typedef struct{

  Tensor Value;

  Tensor Vector;

} EigenTensor;

/*******************************************************/


/*! \struct Fields
 * Variable devoted to store in the memory the physicl 
 * information of each particle
 */
typedef struct {


  Matrix Vol_0; /**< Reference volume */

  Matrix Area_0; /**< Reference area */

  Matrix rho; /**< Density field */

  Matrix rho_s, rho_f; /**< Intrinsic density field (solid/fluid) */

  Matrix phi_s, phi_f; /**< Volume fraction (solid/fluid) */
  
  Matrix mass; /**< Mass field */
  
  Matrix x_GC; /**< Position in global coordinates */
  
  Matrix x_EC; /**< Position in element coordinates */
  
  Matrix dis, D_dis; /**< Displacement field */
  
  Matrix vel; /**< Velocity field */
  
  Matrix acc; /**< Acceleration fields */
  
  Matrix Stress; /**< Stress field */

  /*!
  * Lagrange multiplier for incompressible formulations
  * */
  Matrix lambda_pressure_n;
  Matrix lambda_pressure_n1;

  /*!
  * Pore water pressure, initial state, and the increment
  */
  Matrix Pw;
  Matrix d_Pw;
  Matrix Pw_0;
  Matrix Pw_n1;
  Matrix D_Pw;

  /*!
  * Rates of the pore water pressure
  */
  Matrix d_Pw_dt_n;
  Matrix d_Pw_dt_n1;

  Matrix d2_Pw_dt2;
  
  Matrix Strain; /**< Strain field */

  Matrix F_n, F_n1; /**< Total deformation gradient */

  Matrix DF; /**< Incremental deformation gradient */

  Matrix dt_F_n, dt_F_n1, dt_DF; /**< Rate deformation gradient */

  Matrix J_n, J_n1, dJ_dt; /**< Jacobian of the deformation gradient */

  Matrix Fbar;  /**< Locking-free deformation gradient (F-bar) */
  
  Matrix Jbar; /**< Jacobian of the F-bar */
  
  Matrix Strain_If; /**< Strain during crack */
  
  Matrix W; /**< Deformation Energy */

  double * Chi; /**< Damage parameter (Fracture) */

  double * EPS_n, * EPS_n1; /**< Equivalent plastic strain. */

  double * Kappa_n, * Kappa_n1; /**< Isotropic hardening (stress like) variable. */

  Matrix Back_stress;  /**< Kinematic hardening variable. */

  Matrix b_e_n, b_e_n1; /**< Elastic left deformation gradient */

  Matrix C_ep; /**< Elastoplastic tangent matrix */

  bool * Status_particle; /**< Check if the particle is consider failed or not */


#ifdef DEBUG_MODE
#if DEBUG_MODE + 0

  Matrix PU; /**< Partition of unity property. */

  #endif
#endif

} Fields;

/*******************************************************/

/*! \struct Load
 *   An important aclaration about this code, a force is a load,
 * or even a boundary condition is defined as a load from the
 * point of view of the code (as a structure). The idea is
 * a force and a boundary condition both of them has a direction, 
 * a value, a list of nodes/GPs where it is applied,... So it's
 * a good idea to use this structure for both porpouses.
 */

typedef struct {

  /*!
   * Number of nodes/GP with this load 
   */
  int NumNodes;
  
  /*! 
   * Number of dimensions of the load 
   */
  int Dim;

  /*!
   * Direction of the load {0,0} {1,0} {0,1} {1,1} 
   */
  int * Dir;

  /*!
   * List of nodes with this load 
   */
  int * Nodes;

  /*! 
   * Curve for each dimension with the evolution
   * value with the time 
   */
  Curve * Value;

  /*!
   * Some information about this load 
   */
  char Info [100];

} Load;

/*******************************************************/

/*! \struct Boundaries
 * Boundary conditions definition 
 */
typedef struct {

  /*!
   * Number of boundaries of the domain 
   */
  int NumBounds;
  
  /*!
   * Table with all the boundaries and its values 
   */
  Load * BCC_i;
  
} Boundaries;

/*******************************************************/

/*! \struct Material
 * Properties of a general material model 
 */
typedef struct {

  /*!
   * Index of the material 
   */
  int Id;

  /*!
   * Name of the material
   */
  char Type [100];
  
  /*!
   * Material celerity 
   */
  double Cel;
  
  /*!
   * Initial density (mixture/fluid/solid)
   */
  double rho;
  
  /*!
   * Elastic modulus 
   */
  double E;
  
  /*!
   * Poisson ratio 
   */
  double nu;

  /*!
  * Compressibility
  */
  double Compressibility;

  /*!
  * Fluid parameters
  */
  double ReferencePressure;
  double Viscosity;
  double n_Macdonald_model;

  /*!
   * Activate eigenerosion-fracture modulus (Eigenerosion/Eigensoftening)
   */
  bool Eigenerosion;
  bool Eigensoftening;
  double Ceps; /*! Normalizing constant (Eigenerosion/Eigensoftening) */
  double Gf; /*! Failure energy (Eigenerosion) */
  double ft; /*! Tensile strengt of the material (Eigensoftening) */
  double heps; /*! Bandwidth of the cohesive fracture (Eigensoftening) */
  double Wc; /*! Critical opening displacement (Eigensoftening) */

  /*!
   * Integration algorithm for plasticity
  */
  char Plastic_Solver [100];

  /*!
  * General plastic parameters
  */
  double kappa_0;
  double Hardening_modulus;
  double Plastic_Strain_0;
  double atmospheric_pressure;
  double J2_degradated;
  double Cohesion;

  /*!
   * Frictional material (Borja et al. 2003)
   * */
  char Yield_Function_Frictional [100];
  double m_Frictional;
  double c0_Frictional;
  double phi_Frictional;
  double psi_Frictional;

  /*!
  * Hardening Hughes
  */
  bool Hardening_Hughes;
  double Parameter_Hardening_Hughes;

  /*!
   * Hardening Cervera
   * */
  bool Hardening_Cervera;

  /*!
   * Hardening Ortiz
   * */
  bool Hardening_Ortiz;
  double Exponent_Hardening_Ortiz;

  /*!
   * Hardening Voce
   * */
  bool Hardening_Voce;
  double K_0_Hardening_Voce;
  double K_inf_Hardening_Voce;
  double delta_Hardening_Voce;
  double theta_Hardening_Voce;

  /*!
   * Hardening Borja et al. 2003
   * */
  bool Hardening_Borja;
  double a_Hardening_Borja[3];
  double alpha_Hardening_Borja;

  /*!
   * Viscoplasticity parameters
   * */
  bool Viscous_regularization;
  double fluidity_param;

  /*!
   * Activate auxiliar techniques
   * */
  bool Locking_Control_Fbar;
  double alpha_Fbar;

} Material;

/*******************************************************/

/*! \struct Mixture
 * This structure is devoted to store information 
 * for a general kind of mixtures
 */
typedef struct {

  /*!
   * Index of the mixture 
   */
  int Id;

  /*!
   * Name of the mixture
   */
  char Type [100];

  /* 
    Index for the constitutive description of each phase
    for the soil-water mixture
  */
  int Soil_Idx;
  int Water_Idx;

  /*!
  * Permeability of the soil skleleton
  */
  Tensor Permeability;

  /*!
  * Initial volume fractions (Soil/Water)
  */
  double phi_s_0;
  double phi_f_0;

} Mixture;

/*******************************************************/

/*! 
 * \struct State_Parameters
 */
typedef struct
{
  /*!
   * Particle identifier
   * */
  int Particle_Idx;
  
  /*!
   * Stress/strain parameters
   * */
  double * Stress;
  double * Strain;
  double Pressure;
  double pw;

  double * D_phi_n1;  /**< Total deformation gradient (t = n + 1) */
  double * D_phi_n; /**< Total deformation gradient (t = n) */
  double * d_phi; /**< Incremental deformation gradient */
  double * b_e; /**< Elastic left Cauchy-Green tensor */
  double * Fbar;
  double J; /**<z Jocobian */

  double * dFdt; /**< Rate of the deformation gradient */
  double DeltaTimeStep; /**< Increment of the time step */
  double alpha_4; /**< Time integation paramter (Newmark-beta) */

  double * Back_stress; /**< State varible for kinematic hardening*/
  double * Kappa; /**<  Hardening Parameter */
  double * EPS; /**<  Equivalent plastic strain */
  double * C_ep; /**<  Elastoplastic tangent matrix */
  bool compute_C_ep;

  double Cohesion; 
  double Yield_stress;
  
  // Failed material point
  bool *Failure;

} State_Parameters;

/*******************************************************/

/*! \struct Particle
 * This structure is devoted to store all the information 
 * of a list of particles
 */
typedef struct {

  /*!
   * Number of particles
   * */
  int NumGP;
  
  /*!
   * Index with the closest node to each particle 
   * */
  int * I0;

  /*!
   * Index of the element
   * */
  int * Element_p;

  /*!
   * Tributary nodes variables
   * */
  int * NumberNodes;
  ChainPtr * ListNodes;

  /*! 
   * Set of particles close to each particle 
   * */
  ChainPtr * Beps;

  /*! 
   * Store the values of each field in the current time step
   * */
  Fields Phi;
  
  /*!
   * Values from the previous step 
   * */
  Fields Phi_n0; 
  
  /*!
   * Material variables
   * */
  int NumberMaterials;
  int * MatIdx;
  int * MixtIdx;
  Material * Mat;

  /*!
   * Neumann boundary conditions 
   * */
  int NumNeumannBC;
  Load * F;

  /*!
   * Body forces 
   * */
  int NumberBodyForces;
  Load * B;


  /*!
   * Structure to store Neumann boundary conditions will replace NumNeumannBC and F;
   * */
  Boundaries Neumann_Contours;

  /*!
   * Current vector of distance accelerations
   * */
  Tensor b;


  /*!
   * uGIMP shape function parameter
   * */
  Matrix lp;

  /*!
   * LME shape function parameters:
   * */
  Matrix lambda; // Lagrange multiplier
  Matrix Beta; // Thermalization or regularization parameter (scalar/matrix)
  Matrix Cut_off_Ellipsoid;

  /*!
   * Function to compute the stress state of the particle
   * */
  State_Parameters (* constitutive)(State_Parameters,Material);


} Particle;


/*******************************************************/

/*! \struct Mesh
 * This structure is devoted to store all the information 
 * of a list of nodes
 */
typedef struct {

  /*!
   * Number of nodes in the mesh 
   */
  int NumNodesMesh;
  
  /*!
   * Number of elements in the mesh 
   */
  int NumElemMesh;
  
  /*!
   * Table with the coordinates of the nodes of the mesh 
   */
  Matrix Coordinates;
  
  /*!
   * Number of nodes in a element
   */
  int * NumNodesElem;
    
  /*!
   * List of nodes for each element (Connectivity) 
   */
  ChainPtr * Connectivity;
  
  /*!
   * Number of elements close to each node 
   */
  int * NumNeighbour;
  
  /*!
   * List of elements close to each node 
   */
  ChainPtr * NodeNeighbour;

  /*!
   * Number of nodes close to a node 
   */
  int * SizeNodalLocality_0;
  int * SizeNodalLocality;

  /*!
   * List of nodes close to a node 
   */  
  ChainPtr * NodalLocality_0;
  ChainPtr * NodalLocality;

  /*!
  * Defines if a node is activated or not
  */
  bool * ActiveNode;

  /*!
  * Defines if a node belongs to a boundary or not
  */
  bool * BoundaryNode;

  /*!
   * List of boundaries of the domain
   */
  Boundaries Bounds;
  
  /*! 
   * Number of dimensions of the element (OLD)
   */
  int Dimension;

  /*!
   * Minimum distance between nodes 
   */
  double DeltaX;

  /*! 
  * Nodal spacing
  */
  double * h_avg;

  /*! 
   * Name of the element (OLD)
   */
  char TypeElem [20];

  /*!
   * Function with the interpolation technique
   */
  Matrix (* N_ref)(Matrix );

  /*!
   * Function with the gradient of the interpolation technique
   */
  Matrix (* dNdX_ref)(Matrix );

  /*!
   * Function with the gradient of the interpolation technique
   */
  Matrix (* dNdX)(Matrix, Matrix);

  /*
   * Function of compute natural coordinate
  */
//  void (* X_to_Xi)(Matrix, Matrix, Matrix);

  /* !
  * Function to compute the volume of an element
  */
  double (* volume_Element)(Matrix);

  /*! 
  * Function to check if a point is inside or outside of a elemnt
  */
  bool (* In_Out_Element)(Matrix, Matrix);

  /*!
   * List of particles adjacent to a node 
   */
  int * Num_Particles_Node;
  ChainPtr * List_Particles_Node;

  /*!
   * Variables and function for F-bar calculation
   * */
//  int * Num_Particles_Element;
//  ChainPtr * List_Particles_Element;
  bool Locking_Control_Fbar;
  int Num_Patch_Mesh;
  int * Idx_Patch;
  double * Vol_Patch_n;
  double * Vol_Patch_n1;

} Mesh;

/*******************************************************/

/*! \struct Mask
 *  Structure with the current "element" of the particle
 */
typedef struct {

  int Nactivenodes;

  int * Nodes2Mask;
  
} Mask;

/*******************************************************/

/*! \struct Element
 *  Structure with the current "element" of the particle
 */
typedef struct {

  /*!
   * Index of the particle
   */ 
  int i_GP;
  
  /*! 
   * Number of nodes close to the particle 
   */
  int NumberNodes;
  
  /*!
   * List of tributary nodes for the particle 
   */ 
  int * Connectivity;
  
} Element;

/*******************************************************/

/*! \struct Time_Int_Params
 * Parameters of the time integration scheme
 */
typedef struct {

  /*!
  * Courant number
  */ 
  double CFL;

  /*!
   * Material celerity 
   */
  double Cel;

  /*!
  * Initial step
  */
  int InitialTimeStep;

  /*!
  * Número de pasos de tiempo
  */
  int NumTimeStep;

  /*!
  * Simulation time
  */
  double FinalTime;

  /*!
  * Parameter to generate a mass matrix
  */
  double epsilon_Mass_Matrix;

  /*!
  * Conserving Energy-Momentum parameters 
  */
  double TOL_Conserving_Energy_Momentum;

  /*!
   * Generalized alpha parameters 
   */
  double rb_Generalized_alpha;
  double TOL_Generalized_alpha;

  /*!
  * Newmark parameters
  */
  double beta_Newmark_beta;   
  double gamma_Newmark_beta;
  double TOL_Newmark_beta;


  /*!
  * Maximum number of interations
  */
  int MaxIter;

  /*!
  * Time integration scheme
  */
  char TimeIntegrationScheme [100]; 

} Time_Int_Params;

/*******************************************************/

/*! \struct Event
  Structure with output control
 */
typedef struct {

  /*!
    Number of files of the output
  */
  int NumFiles;

  /*!
    Index for the particle or node in the case of single
    node/particle analysis
  */
  int Idx;

  /*!
    List of nodes or particles to get the information from
    a path of nodes/particles.
  */
  int * Idx_Path;
  int Lenght_Path;

  /*!
    Physical time to start the event. Default = 0.0
   */
  double start;

  /*!
    Numerical step to start the event. Default = 0
   */
  int i_start;

  /*!
    Physical time step
   */
  double step;

  /*!
    Numerical time step
   */
  int i_step;

  /*!
    Physical time to finish the event    
   */
  double end;

  /*!
    Numerical step to finish the event
   */  
  int i_end;

  /*!
    Name of the output directory
   */
  char Directory[MAXC];


  /*!
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

  /*!
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
