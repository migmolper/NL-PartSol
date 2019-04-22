#ifndef TypeDefinitions
#define TypeDefinitions
#endif

void * Allocate_Array(int,int);

double * Allocate_ArrayZ(int);

void ** Allocate_Matrix(int,int,int);

double ** Allocate_MatrixZ(int,int);

Matrix MatAlloc(int,int);

Matrix MatAllocZ(int,int);

Matrix CopyMat(Matrix);

double Get_Determinant(Matrix);

Matrix Get_Inverse(Matrix);

Matrix Transpose_Mat(Matrix);

Matrix Scalar_prod(Matrix,Matrix);

Matrix Tensorial_prod(Matrix,Matrix);

Matrix Add_Mat(Matrix,Matrix);

Matrix Sub_Mat(Matrix,Matrix);

Matrix Norm_Mat(Matrix,int);
