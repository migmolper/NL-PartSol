/*

  Library with all the structures adopted in the code to extend the
  basic functionalities of C with them : 

  -> Matrix : Usefull structure to deal with algebraic operations.
  -> Table : Table of integers.
  -> Chain and ChainPtr : Chain of nodes and pointer to a chain of nodes.
  -> Curve : Structure created for dealing with complex boundary
  conditions and loads.
  -> BoundaryConditions : Structure that store all the information
  of a boundary condition and allow to me and others programers to
  programe/imposse easily new ones.
  -> Load : Structure that the store all the information about 
  one kind of load over the GPs*.
  -> LoadCase : Structure created to join multiple kind of loads in
  one single structure to deal with complex load cases.
  -> Fields : Structure that store all the physical variables of
  the MPM** simulation related with the GPs*.
  -> GaussPoint : Structure created for store general properties of
  a set of GPs*, this allows to generate several groups of them, with
  diferent properties. Ie : water, soil, air.
  -> Mesh : Structure created for store general properties of a FEM*** 
  mesh, this allows to generate differents mesh for diferents sets of 
  GPs*.

  Note* : GPs -> Gauss Points
  Note** : MPM -> Material Point Method
  Note*** : FEM -> Finite Element Method

*/


/*******************************************************/

/* Matrix definition */
typedef struct{
  int N_rows; /* Number of rows */
  int N_cols; /* Number of columns */
  double n; /* Value if is an scalar */
  double * nV; /* Value if is a vector */
  double ** nM; /* Value if is a matrix */
  char Info [100]; /* Aditional information */
} Matrix;

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

/* Chain of nodes */
typedef struct Chain { 
  int I; /* Index of the node */
  struct Chain * next;  /* Pointer to the next element */
} Chain; 

/* Pointer to a chain */
typedef Chain * ChainPtr;

/*******************************************************/

/* Curve definition */
typedef struct{

  /* Number of items in the curve */
  int Num; 
  /* Values for each time */
  double * Fx; 
  /* Aditional information */
  char Info [100];
  
} Curve;

/*******************************************************/

/* Loads definition :
 An important aclaration about this code, a force is a load,
 or even a boundary condition is defined as a load from the
 point of view of the code (as a structure). The idea is
 a force and a boundary condition both of them has a direction, 
 a value, a list of nodes/GPs where it is applied,... So it's
 a good idea to use this structure for both porpouses.
*/
typedef struct {

  /* Number of nodes/GP with this load */
  int NumNodes;
  /* Number of dimensions of the load */
  int Dim;
  /* Direction of the load {0,0} {1,0} {0,1} {1,1} */
  int * Dir;
  /* List of nodes with this load */
  int * Nodes;
  /* Curve for each dimension with the evolution
     value with the time */
  Curve * Value;
  /* Some information about this load */
  char Info [100];

} Load;

/*******************************************************/

/* Boundary conditions definition */
typedef struct {

  /* Number of boundaries of the domain */
  int NumBounds;
  /* Table with all the boundaries and its values */
  Load * BCC_i;
  /* Some information about the boundaries */
  char Info [100];
  
} Boundaries;


/*******************************************************/

/* Load case : multiple loads */
typedef struct {

  /* Number of nodes/GP with this load */
  int NumLoads;
  /* List of loads */
  Load * Load_i;
  /* Some information about this load */
  char Info [100];

} LoadCase;

/*******************************************************/

/* Physical fields */
typedef struct {

  /* Position in global coordinates */
  Matrix x_GC;  
  /* Position in element coordiantes */
  Matrix x_EC;
  /* Voxel shape */
  Matrix lp;
  /* Density field */
  Matrix rho;
  /* Mass field */
  Matrix mass;
  /* Displacement field */
  Matrix dis;
  /* Velocity field */
  Matrix vel;
  /* Acceleration field */
  Matrix acc;
  /* Stress field */
  Matrix Stress;
  /* Strain field */
  Matrix Strain;
  
} Fields;

/*******************************************************/

/* Gauss points definitions */
typedef struct {

  /* Total number of the GP */
  int NumGP;
  /* Identification of the element where it is */
  int * Element_id;
  /* Tributary nodes variables */
  int * NumberNodes;
  ChainPtr * ListNodes;
  /* List of Fields */
  Fields Phi;
  /* Constitutive response */
  Matrix D;
  /* External forces */
  LoadCase F;
  /* Body forces */
  LoadCase B;

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
  /* Active node : 
     Boolean variable that set the node ative (1>) or not (0) */
  int * ActiveNode;
  /* Number of elements that shares a node and list 
     of elements that share a node */
  int * NumNeighbour;
  ChainPtr * NodeNeighbour;

  /*** BOUNDARIES ***/
  Boundaries Bounds;
  
  /*** ELEMENT PROPERTIES ***/
  /* Number of dimensions of the element */
  int Dimension;
  /* Size of the element (COURANT) */
  double DeltaX;
  /* Name of the element */
  char TypeElem [20];
  /* Shape function of the reference element evaluated in a GP */
  Matrix (* N_ref)(Matrix );
  /* Derivative shape function of the reference element evaluated in a GP */
  Matrix (* dNdX_ref)(Matrix );
    
} Mesh;

