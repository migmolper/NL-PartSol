#ifndef _VARIABLES_H_
#define _VARIABLES_H_

/*******************************************************/

/* Matrix definition */
typedef struct{
  int N_rows; /* Number of rows */
  int N_cols; /* Number of columns */
  double n; /* Value if is an scalar */
  double * nV; /* Pointer for a vector */
  double ** nM; /* Table of pointers for a matrix */
  char Info [100]; /* Aditional information */
} Matrix;

/*******************************************************/

/* Chain of nodes */
typedef struct Chain { 
  int I; /* Index of the node */
  struct Chain * next;  /* Pointer to the next element */
} Chain; 

/* Pointer to a chain */
typedef Chain * ChainPtr;

/*******************************************************/

/* Table definition */
typedef struct{
  int N_rows; /* Number of rows */
  int N_cols; /* Number of columns */
  int n; /* Scalar*/
  int * nV; /* 1D list */
  int ** nM; /* 2D list */
  char Info [100]; /* Aditional information */
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

/* Tensor definition */
typedef struct{
  int Order; /* Order of the tensor */
  double *n; /* First order tensor */
  double *N[3]; /* Second order tensor */
  char Info [100]; /* Aditional information */
} Tensor;

/*******************************************************/


/*! \struct Fields
 * Physical fields 
*/
typedef struct {

  /*! Density field */
  Matrix rho;
  /*! Mass field */
  Matrix mass;
  /*! Position in global coordinates */
  Matrix x_GC;  
  /*! Position in element coordiantes */
  Matrix x_EC;
  /*! Displacement field */
  Matrix dis;
  /*! Velocity field */
  Matrix vel;
  /*! Acceleration field */
  Matrix acc;
  /*! Stress field */
  Matrix Stress;
  /*! Strain field */
  Matrix Strain;
  /*! Strain during crack */
  Matrix StrainF;
  /*! Deformation Energy */
  Matrix W;
  /*! Damage parameter (Fracture) */
  Matrix ji;
  
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

  /*! Number of nodes/GP with this load */
  int NumNodes;
  /*! Number of dimensions of the load */
  int Dim;
  /*! Direction of the load {0,0} {1,0} {0,1} {1,1} */
  int * Dir;
  /*! List of nodes with this load */
  int * Nodes;
  /*! Curve for each dimension with the evolution
     value with the time */
  Curve * Value;
  /*! Some information about this load */
  char Info [100];

} Load;

/*******************************************************/

/*! \struct Boundaries
 * Boundary conditions definition 
 */
typedef struct {

  /*! Number of boundaries of the domain */
  int NumBounds;
  /*! Table with all the boundaries and its values */
  Load * BCC_i;
  
} Boundaries;

/*******************************************************/

typedef struct {

  /* Name and id of the material */
  int Id;
  char Type [100];
  /* Material celerity */
  double Cel;
  /* Initial density */
  double rho;
  /* Elastic modulus */
  double E;
  /* Poisson ratio */
  double mu;
  /* Thickness of the Material */
  double thickness;
  /* Activate fracture modulus */
  bool Eigenerosion;
  bool Eigensoftening;
  /*Normalizing constant (Eigenerosion/Eigensoftening) */
  double Ceps;
  /* Failure energy (Eigenerosion)  */
  double Gf;
  /* Tensile strengt of the material (Eigensoftening) */
  double ft;
  /* Bandwidth of the cohesive fracture (Eigensoftening) */
  double heps;
  /* Critical opening displacement */
  double Wc;
  
} Material;

/*******************************************************/

/*! \struct GaussPoint
 * Gauss points definitions
 */
typedef struct {

  /*! Total number of the GP */
  int NumGP;
  /*! Identification of the element where it is */
  int * I0;

  /*! Tributary nodes variables */
  int * NumberNodes;
  ChainPtr * ListNodes;

  /*! GPs near to each GP */
  ChainPtr * Beps;

  /*! Values from the actual step */
  Fields Phi;
  /*! Values from the previous step */
  Fields Phi_n0; 
  
  /*! Number of materials */
  int NumberMaterials;
  /*! Index of the material for each GP */
  int * MatIdx;
  /*! Array wit the number of materials */
  Material * Mat;

  /*! Neumann boundary conditions */
  int NumNeumannBC;
  Load * F;

  /*! Body forces */
  int NumberBodyForces;
  Load * B; 

  /*! Shape functions variables */
  Matrix lp; /*! Voxel shape  */
  Matrix lambda; /*! Lagrange multiplier */
  Matrix Beta; /*! Norm parameter */

} GaussPoint;

/*******************************************************/

/* Mesh definition */
typedef struct {

  /*** GENERAL MESH PROPERTIES ***/
  /* Number of nodes in the mesh */
  int NumNodesMesh;
  /* Number of elements in the mesh */
  int NumElemMesh;
  /* Table with the coordinates of the nodes of the mesh */
  /* double ** Coordinates; */
  Matrix Coordinates;
  /* List of nodes for each element (Connectivity) */
  int * NumNodesElem;
  ChainPtr * Connectivity;
  /* Number of elements that shares a node and list 
     of elements that share a node */
  int * NumNeighbour;
  ChainPtr * NodeNeighbour;

  /* List of nodes close to a node */
  int * SizeNodalLocality;
  ChainPtr * NodalLocality;

  /* Number of particles close to this node */
  int * NumParticles;
  /* List of particles in a node */
  ChainPtr * I_particles;

  /*** BOUNDARIES ***/
  Boundaries Bounds;
  
  /*** ELEMENT PROPERTIES ***/
  /* Number of dimensions of the element */
  int Dimension;
  /* Size of the element */
  double DeltaX;
  /* Name of the element */
  char TypeElem [20];
  /* Shape function for the mesh, for interpolation porpouses */
  Matrix (* N_ref)(Matrix );
  /* Derivative shape function, for interpolation porpouses */
  Matrix (* dNdX_ref)(Matrix );
    
} Mesh;

/*******************************************************/

/*! \struct Element
 *  Individual mesh properties of a GPs
 */
typedef struct{

  /*! Index of the GP */ 
  int i_GP;
  /*! Number of nodes near to the GP */
  int NumberNodes;
  /*! Nodal connectivity of the GP */ 
  int * Connectivity;
  
} Element;

/*******************************************************/

/*! \struct Time_Int_Params
 */
typedef struct {

  /*! Generalized alpha parameters */
  double GA_alpha;
  double GA_beta;
  double GA_gamma;
  
} Time_Int_Params;

/*******************************************************/

#endif
