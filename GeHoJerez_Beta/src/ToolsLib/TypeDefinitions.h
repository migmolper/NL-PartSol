

/* Math structures :
   As in matlab, arrays and matrix are the same,
   note that this definition in not very eficient,
   but give you more flexibility.
*/

/*******************************************************/

typedef struct{
  int N_rows; /* Number of rows */
  int N_cols; /* Number of columns */
  double n; /* Value if is an scalar */
  double * nV; /* Value if is a vector */
  double ** nM; /* Value if is a matrix */
  char Info [100]; /* Aditional information */
} Matrix;

/*******************************************************/

typedef struct{
  int Num; /* Number of items in the curve */
  int Dim; /* Number of dimensions of the curve */
  double ** Fx; /* Values for each time and dimensions */
  char Info [100]; /* Aditional information */
} Curve;

/*******************************************************/

/* Boundary conditions definition */
typedef struct {

  /* Number of nodes/GP with this BCC */
  int NumNodes;
  /* List of nodes with this BCC */
  int * Nodes;
  /* Curve with the value in the time for each direction */
  Curve Value;
  /* Some information about this BCC */
  char Info [100];
  
} BoundaryConditions;

/*******************************************************/

typedef struct {

  /* Number of nodes/GP with this load */
  int NumNodes;
  /* List of nodes with this load */
  int * Nodes;
  /* Curve with the evolution value with the time */
  Curve Value;
  /* Some information about this load */
  char Info [100];

} Load;

/*******************************************************/

typedef struct {

  /* Position in global coordinates */
  Matrix x_GC;  
  /* Position in element coordiantes */
  Matrix x_EC;
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
  /* List of Fields */
  Fields Phi;
  /* Constitutive response */
  Matrix D;
  /* External forces */
  Load F;
  /* Body forces */
  Load B;

} GaussPoint;


/*******************************************************/

/* Element type definition */
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
  int **  Connectivity;
  /* Active node : 
     Boolean variable that set the node ative (1>) or not (0) */
  int * ActiveNode;
  /* Number of elements that shares a node and list 
     of elements that share a node */
  int ** NodeNeighbour;

  /*** BOUNDARY NODES (Only square domains) ***/
  /* TOP */
  BoundaryConditions TOP;
  /* BOTTOM */
  BoundaryConditions BOTTOM;
  /* RIGHT */
  BoundaryConditions RIGHT;
  /* LEFT */
  BoundaryConditions LEFT;
  
  /*** INDIVIDUAL ELEMENT PROPERTIES ***/
  /* Number of dimensions of the element */
  int Dimension;
  /* Name of the element */
  char TypeElem [20];
  /* Number of nodes of the element */
  int NumNodesElem;
  /* Shape function of the reference element evaluated in a GP */
  Matrix (* N_ref)(Matrix );
  /* Derivative shape function of the reference element evaluated in a GP */
  Matrix (* dNdX_ref)(Matrix );
    
} Mesh;

