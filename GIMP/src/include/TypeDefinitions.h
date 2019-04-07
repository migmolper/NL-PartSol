

/* Math structures */
typedef struct{
  int N_elem;
  double * n;
} Array;

typedef struct{
  int N_rows;
  int N_cols;
  double ** n;
} Matrix;


/* Fields declaration */
typedef struct {
  int Size;
  double * n;
} Vector; /* Vectorial field */

typedef struct{
  int Size;
  double ** n;
} Tensor; /* Tensorial field */

/* Gauss points definitions */
typedef struct {
  
  int id; /* Identification number of the GP */
  int Element_id; /* Identification of the element where it is */
  char Material [20]; /* Id of the Material */
  Vector x_GC; /* Position of the GP with global coordiantes */
  Vector x_EC; /* Position of the GP with element coordiantes */
  Vector v; /* Velocity field */
  Vector a; /* Acceleration field */
  Tensor Stress; /* Stress field */
  Tensor Strain; /* Strain field */
  Tensor F_ref; /* Reference deformation gradient */
  Tensor F; /* Deformation gradient */
  Tensor C; /* Lagrangian Cauchy-Green tensor (right) */
  Tensor B; /* Eulerian Cauchy-Green tensor (left) */
  
} GaussPoint;


/* Element type definition */
typedef struct {
  
  int id; /* Identification number of the element */
  int NumberNodes; /* Number of nodes of the element */
  int NumberDOF; /* Degrees of freedom for each node*/
  int NumberGP; /* Number of Gauss points inside of the element */
  int * N_id; /* Global index of the nodes (connectivity) */
  double ** X_g; /* Global coordiantes of the element nodes */
  double * (* n)(Vector *); /* Shape function evaluated in a GP */
  double ** (* dn)(Vector *); /* Derivative shape function evaluated in a GP */
  GaussPoint * GP_e; /* List of GP inside of the element */
  Matrix B; /* It is not the Eulerian Cauchy-Green tensor, it is the operator matrix
	     to get the deformation of the mesh */
  
} Element;

