
#ifndef TypeDefinitions
#define TypeDefinitions
#endif


/* SOLVERS */
Matrix Newton_Rapson(Matrix(* Function)(Matrix, Matrix), Matrix,
		     Matrix(* Jacobian)(Matrix, Matrix), Matrix,
		     Matrix,Matrix);
Matrix Conjugate_Gradient_Method(Matrix,Matrix,Matrix);
Matrix Jacobi_Conjugate_Gradient_Method(Matrix,Matrix,Matrix);
Matrix One_Iteration_Lumped(Matrix, Matrix, Matrix);

/* MATRIX UTILITIES */
void * Allocate_Array(int,int);
void * Allocate_ArrayZ(int,int);
void ** Allocate_Matrix(int,int,int);
void ** Allocate_MatrixZ(int,int,int);
Matrix MatAlloc(int,int);
Matrix MatAllocZ(int,int);
void FreeMat(Matrix);
void PrintMatrix(Matrix, int, int);
double StatsDouMatrix(double *, int, char *);
double StatsIntMatrix(int *, int, char *);
Matrix CopyMat(Matrix);
double Norm_Mat(Matrix,int);
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
double SignumFunct(double x);

/* MATH */
double Area_Poligon(Matrix);
Matrix Centroid_Poligon(Matrix);
int InOut_Poligon(Matrix, Matrix);

