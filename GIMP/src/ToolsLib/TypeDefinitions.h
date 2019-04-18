
/* Math structures :
   As in matlab, arrays and matrix are the same,
   note that this definition in not very eficient,
   but give you more flexibility.
*/

typedef struct{
  int N_rows;
  int N_cols;
  double n;
  double * nV;
  double ** nM;
} Matrix;

/* Gauss points definitions */
typedef struct {

  /* Identification number of the GP */
  int id;

  /* Identification of the element where it is */
  int Element_id;

  /* Name of the Material */
  char Material [20];

  /* Position of the GP with global coordiantes */
  Matrix x_GC;

  /* Position of the GP with element coordiantes */
  Matrix x_EC;

  /* Velocity field */
  Matrix v;
  
  /* Acceleration field */
  Matrix a;
  
  /* Stress field */
  Matrix Stress;
  
  /* Strain field */
  Matrix Strain;

  /* Constitutive response */
  Matrix D;
  
  /* Reference deformation gradient */
  Matrix F_ref;
  
  /* Deformation gradient */
  Matrix F;
  
} GaussPoint;


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

  /* Global coordinates of the element nodes */
  Matrix X_g;
  /*
    ^
    |_____ Future updates : remove this field and use only the connectivity mesh
  */
  
  /* Shape function of the reference element evaluated in a GP */
  Matrix (* N_ref)(Matrix );

  /* Derivative shape function of the reference element evaluated in a GP */
  Matrix (* dNdX_ref)(Matrix );

  /* List of GP inside of the element */
  GaussPoint * GP_e;

  /* It is not the Eulerian Cauchy-Green tensor, it is the operator matrix 
     to get the deformation of the mesh */
  Matrix B;
  
} Element;

