/*! \file Matlib.h
    \brief File with the prototype of math library
*/

#ifndef _MATLIB_H_
#define _MATLIB_H_

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
Matrix Identity__MatrixLib__(int);
Matrix memory_to_matrix__MatrixLib__(int,int,double *);
void   free__MatrixLib__(Matrix);
void   print__MatrixLib__(Matrix, int, int);
Matrix copy__MatrixLib__(Matrix);
double norm__MatrixLib__(Matrix, int);
double Euclidean_distance__MatrixLib__(Matrix);
double generalised_Euclidean_distance__MatrixLib__(Matrix, Matrix);
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
Matrix solve_polynomial__MatrixLib__(Matrix);
double rcond__MatrixLib__(Matrix);
Matrix solve__MatrixLib__(Matrix, Matrix);
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
ChainPtr create_circular_set__SetLib__(ChainPtr);
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
double I1__TensorLib__(const Tensor);
double I2__TensorLib__(const Tensor);
double I3__TensorLib__(const Tensor);
double J1__TensorLib__(const Tensor);
double J2__TensorLib__(const Tensor);
double J3__TensorLib__(const Tensor);

EigenTensor Eigen_analysis__TensorLib__(const Tensor);

double EuclideanNorm__TensorLib__(const Tensor);
double Generalised_norm__TensorLib__(const Tensor, const Tensor);
Tensor Identity__TensorLib__();
Tensor Inverse__TensorLib__(const Tensor);
Tensor Solve_system__TensorLib__(Tensor, Tensor);
Tensor transpose__TensorLib__(const Tensor);
Tensor addition__TensorLib__(Tensor, Tensor);
Tensor subtraction__TensorLib__(Tensor, Tensor);
double inner_product__TensorLib__(Tensor, Tensor);
Tensor vector_product__TensorLib__(Tensor, Tensor);
Tensor dyadic_Product__TensorLib__(Tensor, Tensor);
Tensor vector_linear_mapping__TensorLib__(Tensor, Tensor);
Tensor matrix_product__TensorLib__(Tensor, Tensor);
Tensor Convex_combination__TensorLib__(Tensor, Tensor, double);

void covariant_push_forward_tensor__TensorLib__(Tensor, Tensor, Tensor);
void contravariant_push_forward_tensor__TensorLib__(Tensor, Tensor, Tensor);
void covariant_pull_back_tensor__TensorLib__(Tensor, Tensor, Tensor);
void contravariant_pull_back_tensor__TensorLib__(Tensor, Tensor, Tensor);
void print__TensorLib__(Tensor);


/*******************************************************/

/*!
  \fn int compute_inverse__TensorLib__(double * A_m1, const double * A);

  \brief Compute the inverse of a tensor

  \param b : [out] Inverse of the tensor
  \param F : [in] Tensor
*/
int compute_inverse__TensorLib__(double * A_m1, const double * A);
/*******************************************************/

/*!
  \fn int compute_adjunt__TensorLib__(double * A_mT, const double * A);

  \brief Compute the adjunt of a tensor

  \param b : [out] Adjunt of the tensor
  \param F : [in] Tensor
*/
int compute_adjunt__TensorLib__(double * A_mT, const double * A);
/*******************************************************/

/*!
  \fn void symmetrise__TensorLib__(double * symA, const double * A);

  \brief Compute the symmetric conterpart of a tensor

  \param symA : [out] Symmetric part of a tensor
  \param A : [in] Tensor
*/
void symmetrise__TensorLib__(double * symA, const double * A);
/*******************************************************/

/* LAPACK interfase of the program */
Matrix solve_system_LAPACK(Matrix,Matrix);


#endif
