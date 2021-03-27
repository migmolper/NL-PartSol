/*! \file Matlib.h
    \brief File with the prototype of math library
*/

#ifndef _MATLIB_H_
#define _MATLIB_H_


#define PI__MatrixLib__ 3.14159265358979323846

/*
  Math macros from numerical recipies 
*/
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
Matrix Solve_Linear_Sistem(Matrix, Matrix);
Matrix Conjugate_Gradient_Method(Matrix,Matrix,Matrix);
Matrix Jacobi_Conjugate_Gradient_Method(Matrix,Matrix,Matrix);
Matrix One_Iteration_Lumped(Matrix, Matrix, Matrix);

/*******************************************************/

void * Allocate_Array(int,int);
void * Allocate_ArrayZ(int,int);
void ** Allocate_Matrix(int,int,int);
void ** Allocate_MatrixZ(int,int,int);
double StatsDouMatrix(double *, int, char *);
double StatsIntMatrix(int *, int, char *);

/*******************************************************/

/*
  Matrix library 
*/
Matrix alloc__MatrixLib__(int,int);
Matrix allocZ__MatrixLib__(int,int);
Matrix memory_to_matrix__MatrixLib__(int,int,double *);
void   free__MatrixLib__(Matrix);
void   print__MatrixLib__(Matrix, int, int);
Matrix copy__MatrixLib__(Matrix);
double norm__MatrixLib__(Matrix, int);
double generalised_Euclidean_distance__MatrixLib__(Matrix, Matrix);
double conditioning__MatrixLib__(Matrix, double);
Matrix inverse__MatrixLib__(Matrix);
Matrix transpose__MatrixLib__(Matrix);
double I3__MatrixLib__(Matrix);
Matrix matrix_product__MatrixLib__(Matrix, Matrix);
double scalar_product__MatrixLib__(Matrix, Matrix);
Matrix vectorial_product__MatrixLib__(Matrix, Matrix);
Matrix dyadic_product__MatrixLib__(Matrix,Matrix);
Matrix increment__MatrixLib__(Matrix, Matrix);
Matrix addition__MatrixLib__(Matrix, Matrix);
Matrix substraction__MatrixLib__(Matrix, Matrix);
Matrix lumped__MatrixLib__(Matrix);
double area__MatrixLib__(Matrix);
Matrix centroid__MatrixLib__(Matrix);
int    inout__MatrixLib__(Matrix, Matrix);
Matrix solve_polynomial__MatrixLib__(Matrix);
Matrix nurbs_distance__MatrixLib__(Matrix);
double point_distance__MatrixLib__(Matrix, Matrix);
void   single_value_descomposition__MatrixLib__(Matrix,Matrix,Matrix);

/*******************************************************/

/*
  Set library 
*/
ChainPtr memory_to_set__SetLib__(int *, int);
int * set_to_memory__SetLib__(ChainPtr, int);
ChainPtr range__SetLib__(int, int);
void free__SetLib__(ChainPtr *);
ChainPtr * alloc_table__SetLib__(int);
void free_table__SetLib__(ChainPtr *, int);
bool inout__SetLib__(ChainPtr, int);
void push__SetLib__(ChainPtr *, int);
void pop__SetLib__(ChainPtr *, int);
ChainPtr copy__SetLib__(ChainPtr);
int lenght__SetLib__(ChainPtr);
ChainPtr union__SetLib__(ChainPtr *, int);
ChainPtr intersection__SetLib__(ChainPtr, ChainPtr);
void print__SetLib__(ChainPtr);
void order__SetLib__(ChainPtr *, ChainPtr *, Matrix);

/*******************************************************/

/* 
   Tensor libray 
*/
Tensor alloc__TensorLib__(int);
Tensor memory_to_tensor__TensorLib__(double *, int);
void   free__TensorLib__(Tensor);
double I1__TensorLib__(Tensor);
double I2__TensorLib__(Tensor);
double I3__TensorLib__(Tensor);
double J1__TensorLib__(Tensor);
double J2__TensorLib__(Tensor);
double J3__TensorLib__(Tensor);
Tensor Eigenvalues__TensorLib__(Tensor);
Tensor Eigenvectors__TensorLib__(Tensor,Tensor);
double EuclideanNorm__TensorLib__(Tensor);
double Generalised_norm__TensorLib__(Tensor, Tensor);
Tensor Identity__TensorLib__();
Tensor Inverse__TensorLib__(Tensor);
Tensor Solve_system__TensorLib__(Tensor, Tensor);
Tensor transpose__TensorLib__(Tensor);
Tensor subtraction__TensorLib__(Tensor, Tensor);
double inner_product__TensorLib__(Tensor, Tensor);
Tensor vector_product__TensorLib__(Tensor, Tensor);
Tensor dyadic_Product__TensorLib__(Tensor, Tensor);
Tensor vector_linear_mapping__TensorLib__(Tensor, Tensor);
Tensor matrix_product__TensorLib__(Tensor, Tensor);
Tensor Convex_combination__TensorLib__(Tensor, Tensor, double);
double volumetric_component__TensorLib__(Tensor);
Tensor deviatoric_component__TensorLib__(Tensor, double);
Tensor rotate__TensorLib__(Tensor, Tensor);
void   covariant_push_forward_tensor__TensorLib__(Tensor, Tensor, Tensor);
void   contravariant_push_forward_tensor__TensorLib__(Tensor, Tensor, Tensor);
void   covariant_pull_back_tensor__TensorLib__(Tensor, Tensor, Tensor);
void   contravariant_pull_back_tensor__TensorLib__(Tensor, Tensor, Tensor);
void   print__TensorLib__(Tensor);
/*******************************************************/


/* LAPACK interfase of the program */
Matrix solve_system_LAPACK(Matrix,Matrix);


#endif
