#ifndef TypeDefinitions
#define TypeDefinitions
#endif

void * Allocate_Array(int,int);

void * Allocate_ArrayZ(int,int);

void ** Allocate_Matrix(int,int,int);

void ** Allocate_MatrixZ(int,int,int);

Matrix MatAlloc(int,int);

Matrix MatAllocZ(int,int);

Matrix CopyMat(Matrix);

double Get_Determinant(Matrix);

Matrix Get_Inverse(Matrix);

Matrix Transpose_Mat(Matrix);

Matrix Scalar_prod(Matrix,Matrix);

Matrix Vectorial_prod(Matrix a, Matrix b);

Matrix Tensorial_prod(Matrix,Matrix);

Matrix Add_Mat(Matrix,Matrix);

Matrix Sub_Mat(Matrix,Matrix);

double Norm_Mat(Matrix,int);

int InOut_Poligon(Matrix, Matrix);

Matrix Newton_Rapson(Matrix(* Function)(Matrix, Matrix), Matrix,
		     Matrix(* Jacobian)(Matrix, Matrix), Matrix,
		     Matrix,Matrix);

Matrix Get_Lumped_Matrix(Matrix);
