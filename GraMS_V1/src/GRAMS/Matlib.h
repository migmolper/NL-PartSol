#include <stdbool.h>
  
#ifndef TypeDefinitions
#define TypeDefinitions
#endif

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

/* Solvers library */
/* Matrix Nelder_Mead(); */
Matrix Newton_Rapson(Matrix(* Function)(Matrix, Matrix), Matrix,
		     Matrix(* Jacobian)(Matrix, Matrix), Matrix,
		     Matrix,Matrix);
Matrix Conjugate_Gradient_Method(Matrix,Matrix,Matrix);
Matrix Jacobi_Conjugate_Gradient_Method(Matrix,Matrix,Matrix);
Matrix One_Iteration_Lumped(Matrix, Matrix, Matrix);

/*******************************************************/

/* Matrix library */
void * Allocate_Array(int,int);
void * Allocate_ArrayZ(int,int);
void ** Allocate_Matrix(int,int,int);
void ** Allocate_MatrixZ(int,int,int);
Matrix MatAlloc(int,int);
Matrix MatAllocZ(int,int);
Matrix MatAssign(int,int,double,double *,double **);
void FreeMat(Matrix);
void PrintMatrix(Matrix, int, int);
double StatsDouMatrix(double *, int, char *);
double StatsIntMatrix(int *, int, char *);
Matrix CopyMat(Matrix);
double Norm_Mat(Matrix,int);
double Cond_Mat(Matrix);
double Get_Determinant(Matrix);
Matrix Get_Inverse(Matrix);
Matrix Transpose_Mat(Matrix);
Matrix Scalar_prod(Matrix,Matrix);
Matrix Vectorial_prod(Matrix a, Matrix b);
Matrix Tensorial_prod(Matrix,Matrix);
Matrix Incr_Mat(Matrix, Matrix);
Matrix Add_Mat(Matrix,Matrix);
Matrix Sub_Mat(Matrix,Matrix);
Matrix Get_Lumped_Matrix(Matrix);
double Area_Poligon(Matrix);
Matrix Centroid_Poligon(Matrix);
int InOut_Poligon(Matrix, Matrix);
double SignumFunct(double x);
Matrix SolvePolynomial(Matrix);

/*******************************************************/

/* Chain library */
ChainPtr ArrayToChain(int *, int);
int * ChainToArray(ChainPtr, int);
void FreeChain(ChainPtr);
int LenghtChain(ChainPtr);
bool IsPresentNode (ChainPtr, int);
void PushNodeTop (ChainPtr *, int);
void PopNode (ChainPtr *, int);
ChainPtr CopyChain(ChainPtr);
ChainPtr ChainUnion(ChainPtr, ChainPtr);
ChainPtr ChainIntersection(ChainPtr, ChainPtr);
void printList (ChainPtr);

/*******************************************************/

/* Library */

typedef struct { /* Matrix operators */

  /* Allocate matrix */
  Matrix (* Alloc)(int,int);
  /* Allocate matrix of zeros */
  Matrix (* AllocZ)(int,int);
  /* Asignation of properties */
  Matrix (* Assign)(int,int,double,double *,double **);
  /* Free matrix */
  void (* FreeMat)(Matrix);
  /* Print matrix */
  void (* Print)(Matrix, int, int);
  /* Copy matrix */
  Matrix (* Copy)(Matrix);
  /* Get the norm of a matrix */
  double (* Norm)(Matrix,int);
  /* Get the condition number */
  double (* Cond)(Matrix);
  /* Get the determinant of a matrix */
  double (* Det)(Matrix);
  /* Get the inverse of a matrix */
  Matrix (* Inv)(Matrix);
  /* Get the transpose of a matrix */
  Matrix (* Trans)(Matrix);
  /* Scalar product */
  Matrix (* Sprod)(Matrix,Matrix);
  /* Vectorial product */
  Matrix (* Vprod)(Matrix a, Matrix b);
  /* Tensorial product */
  Matrix (* Tprod)(Matrix,Matrix);
  /* Increment a matrix */
  Matrix (* Incr)(Matrix, Matrix);
  /* Add two matrix */
  Matrix (* Add)(Matrix,Matrix);
  /* Subtract two matrix */
  Matrix (* Sub)(Matrix,Matrix);
  /* Get the lumped matrix */
  Matrix (* Lumped)(Matrix);
  
} MatLib;

MatLib MatrixOperators(void);
