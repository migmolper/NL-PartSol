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

/* Tensor definition */
typedef struct{
  int Order; /* Order of the tensor */
  double *n; /* First order tensor */
  double *N[3]; /* Second order tensor */
  char Info [100]; /* Aditional information */
} Tensor;

/*******************************************************/

/* Math macros from numerical recipies */
static float sqr_arg;
#define SQR(a) ((sqr_arg=(a)) == 0.0 ? 0.0 : sqr_arg*sqr_arg)
static double dsqr_arg;
#define DSQR(a) ((dsqr_arg=(a)) == 0.0 ? 0.0 : dsqr_arg*dsqr_arg)
static double dmax_arg1, dmax_arg2;
#define DMAX(a,b) (dmax_arg1=(a),dmax_arg2=(b),(dmax_arg1) > (dmax_arg2) ? \
		   (dmax_arg1) : (dmax_arg2))
static double dmin_arg1, dmin_arg2;
#define DMIN(a,b) (dmin_arg1=(a),dmin_arg2=(b),(dmin_arg1) < (dmin_arg2) ? \
		   (dmin_arg1) : (dmin_arg2))
static float max_arg1, max_arg2;
#define FMAX(a,b) (max_arg1=(a),max_arg2=(b),(max_arg1) > (max_arg2) ? \
		   (max_arg1) : (max_arg2))
static float min_arg1, min_arg2;
#define FMIN(a,b) (min_arg1=(a),min_arg2=(b),(min_arg1) < (min_arg2) ? \
		   (min_arg1) : (min_arg2))
static long lmax_arg1, lmax_arg2;
#define LMAX(a,b) (lmax_arg1=(a),lmax_arg2=(b),(lmax_arg1) > (lmax_arg2) ? \
		   (lmax_arg1) : (lmax_arg2))
static long lmin_arg1, lmin_arg2;
#define LMIN(a,b) (lmin_arg1=(a),lmin_arg2=(b),(lmin_arg1) < (lmin_arg2) ? \
		   (lmin_arg1) : (lmin_arg2))
static int imax_arg1, imax_arg2;
#define IMAX(a,b) (imax_arg1=(a),imax_arg2=(b),(imax_arg1) > (imax_arg2) ?	\
		   (imax_arg1) : (imax_arg2))
static int imin_arg1, imin_arg2;
#define IMIN(a,b) (imin_arg1=(a),imin_arg2=(b),(imin_arg1) < (imin_arg2) ?	\
		   (imin_arg1) : (imin_arg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


/*******************************************************/

/* Solvers library */
/* Matrix Nelder_Mead(); */
Matrix Newton_Rapson(Matrix(* Function)(Matrix, Matrix), Matrix,
		     Matrix(* Jacobian)(Matrix, Matrix), Matrix,
		     Matrix,Matrix);
Matrix Solve_Linear_Sistem(Matrix, Matrix, Matrix);
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
double Norm_Mat(Matrix, int);
double Cond_Mat(Matrix, double);
double Get_Determinant(Matrix);
Matrix Get_Inverse(Matrix);
Matrix Transpose_Mat(Matrix);
Matrix Scalar_prod(Matrix, Matrix);
Matrix Vectorial_prod(Matrix, Matrix);
Matrix Tensorial_prod(Matrix,Matrix);
Matrix Incr_Mat(Matrix, Matrix);
Matrix Add_Mat(Matrix, Matrix);
Matrix Sub_Mat(Matrix, Matrix);
Matrix Get_Lumped_Matrix(Matrix);
Matrix Matrix_x_Scalar(Matrix, double);
double Area_Poligon(Matrix);
Matrix Centroid_Poligon(Matrix);
int InOut_Poligon(Matrix, Matrix);
double SignumFunct(double x);
Matrix SolvePolynomial(Matrix);
Matrix get_nurbs_distance(Matrix);
double Distance(Matrix, Matrix);
void get_SVD_Of(Matrix A, Matrix W, Matrix V);

/*******************************************************/

/* Chain library */
ChainPtr Pointer_to_Set(int *, int);
int * Set_to_Pointer(ChainPtr, int);
ChainPtr RangeChain(int, int);
void free_Set(ChainPtr);
bool is_in_Set(ChainPtr, int);
void push_to_Set(ChainPtr *, int);
void pop_from_Set(ChainPtr *, int);
ChainPtr CopyChain(ChainPtr);
int get_Lenght_Set(ChainPtr);
ChainPtr get_Union_Of(ChainPtr *, int);
ChainPtr get_Intersection_Of(ChainPtr, ChainPtr);
void print_Set(ChainPtr);
Matrix get_set_Coordinates(ChainPtr, Matrix, Matrix);
void order_Set(ChainPtr *, ChainPtr *, Matrix);

/*******************************************************/

/* Tensor libray */
Tensor alloc_Tensor(int Order);
Tensor memory_to_Tensor(double * A_mem, int Order);
void free_Tensor(Tensor A);
double get_I1_Of(Tensor A);
double get_I2_Of(Tensor A);
double get_I3_Of(Tensor A);
double get_J1_Of(Tensor A);
double get_J2_Of(Tensor A);
double get_J3_Of(Tensor A);
Tensor get_Eigenvalues_Of(Tensor);
double get_EuclideanNorm_Of(Tensor A);
Tensor get_I();
Tensor get_Inverse_Of(Tensor A);
Tensor get_Transpose_Of(Tensor A);
double get_innerProduct_Of(Tensor A, Tensor B);
Tensor get_vectorProduct_Of(Tensor a, Tensor b);
Tensor get_dyadicProduct_Of(Tensor a, Tensor b);
Tensor get_firstOrderContraction_Of(Tensor A, Tensor b);

/*******************************************************/
