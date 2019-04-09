#ifndef TypeDefinitions
#define TypeDefinitions
#endif

void * Allocate_Array(int, int);

void ** Allocate_Matrix(int,int, int);

double Get_Determinant(Tensor);

Tensor Get_Inverse(Tensor);

Matrix Mat_Mul(Matrix, Matrix);

Matrix Mat_Sum(Matrix, Matrix);

Matrix Mat_Sub(Matrix, Matrix);

/*

diadic product



*/
