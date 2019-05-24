

/* Math structures :
   As in matlab, arrays and matrix are the same,
   note that this definition in not very eficient,
   but give you more flexibility.
*/

/*******************************************************/

typedef struct{
  int N_rows;
  int N_cols;
  double n;
  double * nV;
  double ** nM;
  char Info [100];
} Matrix;


/*******************************************************/

typedef struct{
  
  /* Integer identificator for the separator */
  int NumberSeparators;
  char ** sep;
  /* List of keyword por the parser engine */
  int NumberKeyWords;
  char ** KeyWords;
  
} ParserDictionary;

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
  /* Deformation gradient */
  Matrix F;
  
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
  
} GaussPoint;

/*******************************************************/

/* Element type definition */
typedef struct {

  /*** GENERAL MESH PROPERTIES ***/
  /* Number of nodes in the mesh */
  int NumNodesMesh;
  /* Number of elements in the mesh */
  int NumElemMesh;
  /* Number of boundary nodes */
  int NumNodesBound;
  /* Table with the coordinates of the nodes of the mesh */
  /* double ** Coordinates; */
  Matrix Coordinates;
  /* List of nodes for each element (Connectivity) */
  int **  Connectivity;
  /* List with the boundary nodes */
  int * NodesBound;
  /* Active element : 
     Boolean variable that set the element ative (1) or not (0) */
  int * ActiveElem;

  /*** INDIVIDUAL ELEMENT PROPERTIES ***/
  /* Number of dimensions of the element */
  int Dimension;
  /* Name of the element */
  char TypeElem [20];
  /* Number of nodes of the element */
  int NumNodesElem;
  /* Degrees of freedom for each node*/
  int NumberDOF;
  /* Shape function of the reference element evaluated in a GP */
  Matrix (* N_ref)(Matrix );
  /* Derivative shape function of the reference element evaluated in a GP */
  Matrix (* dNdX_ref)(Matrix );
    
} Mesh;

