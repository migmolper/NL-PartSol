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
  int I;

  /*!
   * Pointer to the next node in the set
   */
  struct Chain * next;
  
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
  double *N[3];

  /*! 
   * Aditional information 
   */
  char Info [100];
  
} Tensor;

/*******************************************************/


/*! \struct Fields
 * Variable devoted to store in the memory the physicl 
 * information of each particle
 */
typedef struct {

  /*!
   * Density field 
   */
  Matrix rho;
  
  /*!
   * Mass field 
   */
  Matrix mass;
  
  /*!
   * Position in global coordinates 
   */
  Matrix x_GC;
  
  /*!
   * Position in element coordiantes 
   */
  Matrix x_EC;
  
  /*!
   * Displacement field 
   */
  Matrix dis;
  
  /*! 
   * Velocity field 
   */
  Matrix vel;
  
  /*!
   * Acceleration field 
   */
  Matrix acc;
  
  /*!
   * Stress field
   */
  Matrix Stress;
  
  /*!
   * Strain field 
   */
  Matrix Strain;

  /*!
   * Deformation gradient at t = n and t = n + 1
   */
  Matrix F_n;
  Matrix F_n1;

  /*!
   * Plastic deformation gradient
   */
  Matrix F_plastic;
  
  /*!
   * Strain during crack 
   */
  Matrix Strain_If;
  
  /*!
   * Deformation Energy 
   */
  Matrix W;
  
  /*! 
   * Damage parameter (Fracture) 
   */
  Matrix chi;

  /*!
   * Cohesion of the particle (plasticity) 
   */
  Matrix cohesion;

  /*!
   * Equivalent plastic strain of the particle (plasticity) 
   */
  Matrix EPS;

  
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
 * Properties of a material model 
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
   * Thickness of the Material 
   */
  double thickness;
  
  /*!
   * Material celerity 
   */
  double Cel;
  
  /*!
   * Initial density 
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
   * General plastic parameters
   */
  double yield_stress;
  double cohesion_reference;
  double friction_angle;
  double dilatancy_angle;
  double hardening_modulus;
  
  /*!
   * Parameters of the Drucker-Prager Sanavia
   */
  double alpha_F_Drucker_Prager;
  double alpha_Q_Drucker_Prager;
  double beta_Drucker_Prager;
  double hardening_exp_Drucker_Prager; /*! Hardening exponent */
  double E_plastic_reference_Drucker_Prager; /*! Reference plastic strain */
  
} Material;

/*******************************************************/

/*! \struct GaussPoint
 * This structure is devoted to store all the information 
 * of a list of particles
 */
typedef struct {

  /*!
   * Number of particles  
   */
  int NumGP;
  
  /*!
   * Index with the closest node to each particle 
   */
  int * I0;

  /*! Tributary nodes variables */
  int * NumberNodes;
  ChainPtr * ListNodes;

  /*! 
   * Set of particles close to each particle 
   */
  ChainPtr * Beps;

  /*! 
   * Store the values of each field in the current time step
   */
  Fields Phi;
  
  /*! 
   * Values from the previous step 
   */
  Fields Phi_n0; 
  
  /*! 
   * Number of materials 
   */
  int NumberMaterials;
  
  /*! 
   * Index of the material for each particle
   */
  int * MatIdx;
  
  /*!
   * Library of materials 
   */
  Material * Mat;

  /*!
   * Number of Neumann boundary conditions 
   */
  int NumNeumannBC;

  /*!
   * Load case of Neumann boundary conditions 
   */  
  Load * F;

  /*! 
   * Number of body forces 
   */
  int NumberBodyForces;

  /*!
   * Load case for the body forces
   */
  Load * B; 

  /*!
   * Size of the voxel for each particle. Variable for the uGIMP 
   * shape function
   */
  Matrix lp;

  /*!
   * Lagrange multiplier for the LME shape functions 
   */
  Matrix lambda;

  /*! 
   * Thermalization or regularization parameter for the LME shape functions
   */
  Matrix Beta;  

} GaussPoint;


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
  int * SizeNodalLocality;

  /*!
   * List of nodes close to a node 
   */  
  ChainPtr * NodalLocality;

  /*!
   * List with the number of particles close to a node
   */
  int * NumParticles;

  /*!
   * List of particles in a node 
   */
  ChainPtr * I_particles;

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
    
} Mesh;

/*******************************************************/

/*! \struct Mask
 *  Structure with the current "element" of the particle
 */
typedef struct{

  int Nactivenodes;

  int * Nodes2Mask;
  
} Mask;

/*******************************************************/

/*! \struct Element
 *  Structure with the current "element" of the particle
 */
typedef struct{

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
   * Generalized alpha parameter alpha 
   */
  double GA_alpha;

  /*!
   * Generalized alpha parameter alpha 
   */
  double GA_beta;

  /*!
   * Generalized alpha parameter gamma 
   */
  double GA_gamma;
  
} Time_Int_Params;

/*******************************************************/

/*! \struct Event
  Structure with output control
 */
typedef struct {

  /*!
    Physical time to start the event. Default = 0.0
   */
  double start;

  /*!
    Numerical step to start the event. Default = 0
   */
  int k_start;

  /*!
    Physical time step
   */
  double step;

  /*!
    Numerical time step
   */
  int k_step;

  /*!
    Physical time to finish the event    
   */
  double end;

  /*!
    Numerical step to finish the event
   */  
  int k_end;

  /*!
    Name of the output file
   */
  char * File;

} Event;

/*******************************************************/

#endif
