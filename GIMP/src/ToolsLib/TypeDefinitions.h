
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
  char Info [80];
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

  /* Number of nodes in the mesh */
  int Nnodes;
  /* Number of elements in the mesh */
  int Nelem;
  /* Table with the coordinates of the nodes of the mesh */
  double ** Coordinates;
  /* List of nodes for each element (Connectivity) */
  int **  Connectivity;
  
} MeshProperties;

/*******************************************************/

typedef struct {

  /* Number of dimensions of the element */
  int Dimension;
  /* Name of the element */
  char Type [20];
  /* Number of nodes of the element */
  int Nnodes;
  
} ElementProperties;

/*******************************************************/

typedef struct {

  /* Position in global coordinates */
  Matrix x_GC;  
  /* Position in element coordiantes */
  Matrix x_EC;
  /* Density field */
  double rho;
  /* Mass field */
  double mass;
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

  /* Identification number of the element */
  int id;

  /* Number of nodes of the element */
  int NumberNodes;

  /* Degrees of freedom for each node*/
  int NumberDOF;

  /* Number of Gauss points inside of the element */
  int NumberGP;

  /* Global index of the nodes (connectivity) */
  int * N_id;
  
  /* Shape function of the reference element evaluated in a GP */
  Matrix (* N_ref)(Matrix );

  /* Derivative shape function of the reference element evaluated in a GP */
  Matrix (* dNdX_ref)(Matrix );
  
  /* Ful derivative matrix to get the deformation of the element */
  Matrix dNdx;

  /* Reference deformation tensor */
  Matrix F_ref;

  /* Deformation tensor */
  Matrix F;

  /* Eulerian Cauchy-Green tensor */
  Matrix B;
  
  /* Lagrangian Cauchy-Green tensor */
  Matrix C;
  
} Element;

/*******************************************************/
