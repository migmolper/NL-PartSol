#include "grams.h"

/*********************************************************************/

/* Function for arrays declaration */
void * Allocate_Array(int SizeArray, int SizeType)
/*
  Function for array declaration
  Inputs : Number of rows, number of columns and kind of element 
  (double, integer, ...)
  Outpus : Matrix
*/
{
  
  void * V;  
  V = (void *)malloc(SizeArray*SizeType);  
  if (V == NULL){
    puts("Error in the array declaration");
    exit(EXIT_FAILURE);
  }  
  return V;
  
}

/*********************************************************************/


/* Function for arrays declaration */
void * Allocate_ArrayZ(int SizeArray, int SizeType)
/*
  Function for array of zeros declaration (double)
  Inputs : Number of rows, number of columns and kind of element 
  (double, integer, ...)
  Outpus : Matrix
*/
{
  
  void * V;
  V = (void*)calloc(SizeArray, SizeType);
  if (V == NULL){puts("Error in the array declaration");
    exit(EXIT_FAILURE);}  
  return V;
  
}

/*********************************************************************/

void ** Allocate_Matrix(int NumberRows,int NumberColumns, int SizeType)
/*
  Function for matrix declaration
  Inputs : Number of rows, number of columns and kind of element 
  (double, integer, ...)
  Outpus : Matrix
 */
{

  void ** M;
  M = (void **)malloc((unsigned) NumberRows * sizeof(void *));

  if (M == NULL){
    puts("Error in matrix declaration");
    exit(EXIT_FAILURE);
  }
  
  for(int i = 0 ; i<NumberRows ; i++){
    M[i] = malloc((unsigned) NumberColumns*SizeType);
    
    if (M[i] == NULL){
      puts("Error in matrix declaration");
      exit(EXIT_FAILURE);
    }
    
  }
  return M;
  
}


/*********************************************************************/

void ** Allocate_MatrixZ(int NumberRows,int NumberColumns, int SizeType)
/*
  Function for matrix of zeros declaration (double)
  Inputs : Number of rows, number of columns.
  Outpus : Matrix
 */
{

  void ** M;
  M = (void **)calloc(NumberRows,sizeof(void *));
  if (M == NULL){puts("Error in matrix declaration");  exit(EXIT_FAILURE);}
  for(int i = 0 ; i<NumberRows ; i++){
    M[i] = (void *)calloc(NumberColumns,SizeType);
    if (M[i] == NULL){puts("Error in matrix declaration");  exit(EXIT_FAILURE);}
  }
  return M;
  
}


/*********************************************************************/

Matrix MatAlloc(int NumberRows,int NumberColumns){

  Matrix M;
  M.N_rows = NumberRows;
  M.N_cols = NumberColumns;

  if( (NumberRows == 1) || (NumberColumns == 1) ){ /* It is an array */
    M.nV = (double *)Allocate_Array(NumberRows*NumberColumns,sizeof(double));
    M.nM = NULL;
    M.n = NAN;
  }
  else if( (NumberRows != 1) && (NumberColumns != 1)  ){ /* It is a matrix */
    M.nM = (double **)Allocate_Matrix(NumberRows,NumberColumns,sizeof(double));
    M.nV = NULL;
    M.n = NAN;
  }
  
  return M;
}

/*********************************************************************/


Matrix MatAllocZ(int NumberRows,int NumberColumns)
/*
  Allocate the matrix structure with zeros 
*/
{

  Matrix M;
  M.N_rows = NumberRows;
  M.N_cols = NumberColumns;

  if( (NumberRows == 1) || (NumberColumns == 1) ){ /* It is an array */
    M.nV = (double *)Allocate_ArrayZ(NumberRows*NumberColumns,
				     sizeof(double));
    M.nM = NULL;
    M.n = NAN;
  }
  else if( (NumberRows != 1) && (NumberColumns != 1)  ){ /* It is a matrix */
    M.nM = (double **)Allocate_MatrixZ(NumberRows,
				      NumberColumns,
				      sizeof(double));
    M.nV = NULL;
    M.n = NAN;
  }
  
  return M;
}

/*********************************************************************/
Matrix MatAssign(int N_rows,int N_cols,
		 double n, double * nV,
		 double ** nM)
/*
  Define general parameters of a matrix for auxiliar
  matrix.
*/
{
  Matrix M;
  M.N_rows = N_rows;
  M.N_cols = N_cols;
  M.nV = nV;
  M.nM = nM;
  M.n = n;
  
  return M;
}

/*********************************************************************/

void FreeMat(Matrix Input){

  if( ((Input.N_cols == 1) && (Input.N_rows > 1)) ||
      ((Input.N_cols > 1) && (Input.N_rows == 1)) ){ /* It is a vector */
    if(Input.nV == NULL){
      printf("Error in FreeMat() : This vector is NULL");
       exit(EXIT_FAILURE);
    }
    free(Input.nV);
  }
  else if( (Input.N_cols > 1) && (Input.N_rows > 1) ){ /* It is a matrix */
    if(Input.nM == NULL){
      printf("Error in FreeMat() : This matrix is NULL !!! \n");
       exit(EXIT_FAILURE);
    }
    for(int i = 0 ; i<Input.N_rows ; i++){
      free(Input.nM[i]);
    }
    free(Input.nM);
  }
  else{
    printf("Error in FreeMat() : You are trying to free a scalar !!! \n");
     exit(EXIT_FAILURE);
  }
  
}


/*********************************************************************/

void PrintMatrix(Matrix In, int NumCols, int NumRows)
/*
  Print term by term the input matrix
*/
{

  int Nvalues;
  
  if( In.nV != NULL ){
    Nvalues = In.N_cols*In.N_rows;
    Nvalues = IMIN(NumRows,Nvalues);
    printf("%s : \n",In.Info);
    for(int i = 0 ; i<Nvalues ; i++){
      printf("%f \n",In.nV[i]);
    }
  }
  
  if( In.nM != NULL ){
    
    printf("%s : \n",In.Info);
    for(int i = 0 ; i<IMAX(In.N_cols,In.N_rows) ; i++){
      for(int j = 0 ; j<IMIN(In.N_cols,In.N_rows) ; j++){
	if(In.N_cols<In.N_rows)
	  printf(" %f ",In.nM[i][j]);
	else
	  printf(" %f ",In.nM[j][i]);
      }
      printf("\n");
    }
  }
}

/*********************************************************************/

Matrix CopyMat(Matrix In)
/* 
   Copy the input matrix in a auxiliar matrix, this is necessary because the operators 
   for the linear algebra are destructive for the input data
 */
{
  if( (In.N_rows==0) || (In.N_cols==0)){
    puts("Error in CopyMat() : Input matrix has 0 dimensions !");
     exit(EXIT_FAILURE);
  }
  Matrix Out = MatAlloc(In.N_rows,In.N_cols);

  if(In.nM != NULL){
    for(int i = 0; i<In.N_rows; i++){
      memcpy(&Out.nM[i], &In.nM[i], sizeof(In.nM[0]));
    }
  }
  if(In.nV != NULL){
    for(int i = 0; i<In.N_rows; i++){
      memcpy(&Out.nV[i], &In.nV[i], sizeof(In.nV[0]));
    }

  }
  if( (In.nM == NULL) &&(In.nV == NULL) ){
    puts("Error in CopyMat() : The input matrix is NULL !");
  }

  return Out;
}


/*********************************************************************/

double Get_Determinant(Matrix M_in)
/* 
   Function to get the determinant of a matrix (max 3x3) :
   Inputs : M_in Matrix (Tensor type)
   Outputs : Determinant
*/
{
  /* Check if we dont have a null matrix */
  if ( (M_in.nM == NULL) && (M_in.n != M_in.n) ){
    if(M_in.nV != NULL){
      puts("Error in Get_Determinant() : An array does not have determinant !");
       exit(EXIT_FAILURE);
    }
    puts("Error in Get_Determinant() : Input matrix = NULL !");    
     exit(EXIT_FAILURE);
  }
  
  /* Check if the matrix is square */
  if(M_in.N_cols != M_in.N_rows){
    puts(" Error in Get_Determinant() : Non square matrix !");
     exit(EXIT_FAILURE);
  }
   
  double Det_out;

  /* Get the determinant */
  switch(M_in.N_cols){

  case 1:
    Det_out = M_in.n;
    break;

  case 2:
    Det_out =
      M_in.nM[0][0]*M_in.nM[1][1] -
      M_in.nM[0][1]*M_in.nM[1][0];
    break;
    
  case 3:
    Det_out =
      M_in.nM[0][0]*M_in.nM[1][1]*M_in.nM[2][2] - /* + a11*a22*a33 */
      M_in.nM[0][0]*M_in.nM[1][2]*M_in.nM[2][1] + /* - a11*a23*a32 */
      M_in.nM[0][1]*M_in.nM[1][2]*M_in.nM[2][0] - /* + a12*a23*a31 */
      M_in.nM[0][1]*M_in.nM[1][0]*M_in.nM[2][2] + /* - a12*a33*a21 */
      M_in.nM[0][2]*M_in.nM[1][0]*M_in.nM[2][1] - /* + a13*a21*a32 */
      M_in.nM[0][2]*M_in.nM[1][1]*M_in.nM[2][0] ; /* - a13*a22*a31 */
    break;

  default :
    fprintf(stderr,"%s : %s \n",
	    "Error in Get_Determinant()",
	    "I am only able to get the determinat of a 3x3 matrix !");
     exit(EXIT_FAILURE); 
  }
  
  return Det_out;
}


/*********************************************************************/

Matrix Get_Inverse(Matrix M_in)
/* 
   Function to get the inverse of the matrix :
   Inputs : Matrix in, dimension os the matrix
   Outpus : Matrix out
*/
{ 
  /* Check if we dont have a null matrix */
  if (M_in.nM == NULL){
    if(M_in.nV != NULL){
      puts("Error in Get_Inverse : An array does not have inverse !");
       exit(EXIT_FAILURE);
    }
    else{
      puts("Error in Get_Inverse : Input matrix = NULL !");    
       exit(EXIT_FAILURE);
    }
  }

  /* Check if the Matrix is invertible */

  /* Is the matrix square ? */
  if(M_in.N_cols != M_in.N_rows){
    puts(" Error in Get_Inverse() : Non square matrix !");
     exit(EXIT_FAILURE);
  }
  
  /* Get the determinant of the matrix */
  double Det = Get_Determinant(M_in);
  if(Det == 0){
    puts("Error in Get_Inverse() : Determinant null !");
     exit(EXIT_FAILURE);
  }
  double Det_m1 = (double)1/Det;

  /* Allocate the output */
  Matrix M_out = MatAlloc(M_in.N_rows,M_in.N_cols);
 
  /* Do in a different fashion if it is 2D or 3D */
  if(M_in.N_cols == 2){
    M_out.nM[0][0] = Det_m1*M_in.nM[1][1];

    M_out.nM[0][1] = - Det_m1*M_in.nM[0][1];

    M_out.nM[1][0] = - Det_m1*M_in.nM[1][0];

    M_out.nM[1][1] = Det_m1*M_in.nM[0][0];
  }
  else if(M_in.N_cols == 3){
    M_out.nM[0][0] = Det_m1*(M_in.nM[1][1]*M_in.nM[2][2] -
			     M_in.nM[1][2]*M_in.nM[2][1]);
    
    M_out.nM[0][1] = -Det_m1*(M_in.nM[0][1]*M_in.nM[2][2] -
			      M_in.nM[0][2]*M_in.nM[2][1]);
    
    M_out.nM[0][2] = Det_m1*(M_in.nM[0][1]*M_in.nM[1][2] -
			     M_in.nM[0][2]*M_in.nM[1][1]);
    
    M_out.nM[1][0] = -Det_m1*(M_in.nM[1][0]*M_in.nM[2][2] -
			      M_in.nM[1][2]*M_in.nM[2][0]);
    
    M_out.nM[1][1] = Det_m1*(M_in.nM[0][0]*M_in.nM[2][2] -
			     M_in.nM[0][2]*M_in.nM[2][0]);
    
    M_out.nM[1][2] = -Det_m1*(M_in.nM[0][0]*M_in.nM[1][2] -
			      M_in.nM[0][2]*M_in.nM[1][0]);
    
    M_out.nM[2][0] = Det_m1*(M_in.nM[1][0]*M_in.nM[2][1] -
			     M_in.nM[1][1]*M_in.nM[2][0]);
    
    M_out.nM[2][1] = -Det_m1*(M_in.nM[0][0]*M_in.nM[2][1] -
			      M_in.nM[0][1]*M_in.nM[2][0]);
    
    M_out.nM[2][2] = Det_m1*(M_in.nM[0][0]*M_in.nM[1][1] -
			     M_in.nM[0][1]*M_in.nM[1][0]);
  }
  
  /* Return the inverse matrix */
  return M_out;
  
}

/*********************************************************************/

Matrix Transpose_Mat(Matrix M)
/* 
   It gives to you the transpose of a matrix or an array matrix.
 */
{

  /* Define and allocate the output matrix */
  Matrix M_T = MatAlloc(M.N_cols,M.N_rows);  

  /* Do it in a different fashion if it is a matrix o an array */
  if ( M.nM != NULL ){ /* It is a matrix */    
    for(int i = 0 ; i < M_T.N_rows ; i++){
      for(int j = 0 ; j < M_T.N_cols ; j++){
	  M_T.nM[i][j] = M.nM[j][i];
      }/* for i*/
    }/* for j */ 
  }
  if( M.nV != NULL){
    M_T.nV = M.nV;
  }
  if(M.n == M.n){
    puts("Warning : Your are transposing a scalar ");
  }
  
  return M_T;
}

/*********************************************************************/

Matrix Scalar_prod(Matrix A,Matrix B)
/*
  Multiply two matrix A and B, and return the result C.
*/
{

  /* Variable declaration output matrix */
  Matrix C;
  double C_aux;

  /* If it is array or matrix, allocate and operate in a different fashion */

  /* Scalar product of two matrix */
  if ((A.N_cols == B.N_rows) && (A.N_rows >= 1) && (B.N_cols >= 1))
    {      
      C = get_A_dot_B_Mat(A,B);
    }
  
  /* Scalar product of an array by an array */
  else if ((A.N_cols == B.N_rows) && (A.N_rows == 1) && (B.N_cols == 1))
    {      
      C = get_a_dot_b_Mat(A,B);
    }
  
  /* Scalar product of a matrix by an array */
  else if ((A.N_cols == B.N_rows) && (A.N_rows >= 1) && (B.N_cols == 1))
    {       
      C = get_A_dot_b_Mat(A, B);
    }

  /* Scalar product of an array by a matrix */
  else if ((A.N_cols == B.N_rows) && (A.N_rows == 1) && (B.N_cols >= 1))
    {
      C = get_a_dot_B_Mat(A, B);
    }

  /* Scalar product of a matrix/array by a scalar */
  else if ((B.N_rows == 1) && (B.N_cols == 1))
    {       
      C = Matrix_x_Scalar(A, B.n);
    }
  
  /* Scalar product of a scalar by a matrix/array */
  else if ((A.N_rows == 1) && (A.N_cols == 1))
    {       
      C = Matrix_x_Scalar(B, A.n);  
    }
  
  else
    {
      printf("%s : %s -> [%i , %i] x [%i , %i] \n",
	     "Error in Scalar_prod()","Incompatible product",
	     A.N_rows,A.N_cols,B.N_rows,B.N_cols);
      exit(EXIT_FAILURE);
    }
  
  return C;
  
}

/*********************************************************************/

Matrix get_A_dot_B_Mat(Matrix A, Matrix B)
/* 
   Compute the scalar product of two matrix 
*/
{

  /* The result is a matrix */
  Matrix A_dot_B = MatAllocZ(A.N_rows,B.N_cols);
  int Rows = A.N_rows;
  int Columns = B.N_cols;
  int Aux = B.N_rows;
     
  /* Multiply */
  for(int i = 0 ; i < Rows  ; i++){
    for(int j = 0 ; j < Columns  ; j++){
      for(int k = 0 ; k < Aux  ; k++){
	A_dot_B.nM[i][j] += A.nM[i][k]*B.nM[k][j];
      }
    } /* for j */
  } /* for i */ 

  return A_dot_B;  
}

/*********************************************************************/

Matrix get_a_dot_b_Mat(Matrix a, Matrix b)
/* 
   Compute the scalar product of two vectors
*/
{
  Matrix a_dot_b = MatAssign(1,1,0,NULL,NULL);
  int Size = a.N_cols;

  for(int i = 0 ; i < Size ; i++){
    a_dot_b.n += a.nV[i]*b.nV[i];
  }

  return a_dot_b;  
}

/*********************************************************************/

Matrix get_A_dot_b_Mat(Matrix A, Matrix b)
/* 
   Compute the scalar product of a matrix by an array
*/
{
  /* The result is an array */
  Matrix A_dot_b = MatAllocZ(A.N_rows,b.N_cols);
  int Rows = A.N_rows;
  int Columns = A.N_cols;
  
  for(int i = 0 ; i<Rows ; i++){
    for(int j = 0 ; j<Columns ; j++){
      A_dot_b.nV[i] += A.nM[i][j]*b.nV[j];
    }
  }
    
  return A_dot_b;  
}

/*********************************************************************/

Matrix get_a_dot_B_Mat(Matrix a, Matrix B)
/* 
   Compute the scalar product of a matrix by an array
*/
{
  /* The result is an array */
  Matrix a_dot_B = MatAllocZ(a.N_rows,B.N_cols);
  int Rows = B.N_rows;
  int Columns = B.N_cols;
  
  for(int i = 0 ; i<Columns ; i++){
    for(int j = 0 ; j<Rows ; j++){
      a_dot_B.nV[i] += a.nV[j]*B.nM[j][i];
    }
  }
    
  return a_dot_B;  
}

/*********************************************************************/

Matrix Vectorial_prod(Matrix a, Matrix b){
  Matrix c;

  if((a.N_rows == 3)&&(b.N_rows == 3)&&(a.N_cols==1)&&(b.N_cols==1)){
      c = MatAlloc(3,1);
      c.nV[0] = a.nV[1]*b.nV[2] - a.nV[2]*b.nV[1];
      c.nV[1] = a.nV[2]*b.nV[0] - a.nV[0]*b.nV[2];
      c.nV[2] = a.nV[0]*b.nV[1] - a.nV[1]*b.nV[0];
  }
  else{
    puts("Error in Vectorial_prod() : Incompatible shape of arrays !! ");
    puts("Remember that it should be : (3x1) x (3x1)  ");
     exit(EXIT_FAILURE);
  }
  
  return c;
  
}

/*********************************************************************/
  
Matrix Tensorial_prod(Matrix A,Matrix B)
/*
  Tensorial product between A and B.
 */
{
  /* Matrix declaration */
  Matrix C;

  if( (A.nV != NULL)&&(B.nV != NULL) ){ /* Tensorial product between two arrays */

    if( (A.N_cols != 1) || (B.N_rows !=1) ){
      puts("Error in Tensorial_prod() : Incompatible shape of arrays ");
      puts("Remember that it should be : (nx1) · (1xn)  ");
       exit(EXIT_FAILURE);
    }
    
    /* Allocate matrix */
    C = MatAlloc(A.N_rows,B.N_cols);

    for(int i = 0 ; i<C.N_rows ; i++){
      for(int j = 0 ; j<C.N_cols ; j++){
	C.nM[i][j] = A.nV[i]*B.nV[j];	
      } /* for i */
    } /* for j */
    
  }  /* Tensorial product between an array and a matrix */
  else if( ((A.nM != NULL)&&(B.nV != NULL)) ||
	   ((A.nV != NULL)&&(B.nM != NULL))){
    puts(" Tensorial producto between a matrix and a array is not already defined");
     exit(EXIT_FAILURE);
  }
  else if((A.nM != NULL)&&(B.nM != NULL)){
    puts(" Tensorial producto between a matrix and a matrix is not already defined ");
     exit(EXIT_FAILURE);
  }
  
  return C;
}
  

/*********************************************************************/

Matrix Incr_Mat(Matrix A, Matrix Incr)
/*
  Increment term by term the Matrix A with Incr
*/
{

  int Bool;

  Bool = Incr.N_cols*Incr.N_rows == 1;

  switch(Bool){
    
  case 0: /* Incr is a matrix or an array */

    if( (A.N_cols != Incr.N_cols) ||
	(A.N_rows != Incr.N_rows) ){ /* Check the size of the input */
      printf("Error in Incr_Mat() : A is incompatible with the Incr !! \n");
       exit(EXIT_FAILURE);
    }

    if(Incr.N_cols == Incr.N_rows){ /* Matrix */
      for(int i = 0 ; i<Incr.N_rows ; i++){
	for(int j = 0 ; j<Incr.N_cols ; j++){
	  A.nM[i][j] += Incr.nM[i][j];
	}
      }	
    }
    
    if( (Incr.N_cols == 1) ||
	(Incr.N_rows == 1) ){ /* Array */
      for(int i = 0 ; i<Incr.N_cols*Incr.N_rows ; i++){
	A.nV[i] += Incr.nV[i];
      }
    }
    
    break;
    
  case 1: /* Incr is a scalar */
    
    if( (A.N_cols == 1) ||
	(A.N_rows == 1) ){ /* A is an array */
      for(int i = 0 ; i<A.N_cols*A.N_rows ; i++){
	A.nV[i] += Incr.n;
      }
    }
    
    if( (A.N_cols != 1) &&
	(A.N_rows != 1) ){
      for(int i = 0 ; i<A.N_rows ; i++){
	for(int j = 0 ; j<A.N_cols ; j++){
	  A.nM[i][j] += Incr.n;
	}
      }	
    }
    
    break;

    
  default :
     exit(EXIT_FAILURE);
  }

  return A;
}  


/*********************************************************************/

Matrix Add_Mat(Matrix A,Matrix B)
/*

  Sum two matrix A and B, and return the result C 

  C = A + B
   
 */
{

  /* Variable declaration */
  Matrix C;
  
  /* Check if it is possible to do the addition */
  if((A.N_cols != B.N_cols)||(A.N_rows != B.N_rows)){
    puts("Error in Add_Mat() : Your are trying to add incompatible matrix");
     exit(EXIT_FAILURE);
  }

  /* Allocate output matrix */
  C = MatAlloc(A.N_rows,A.N_cols);

  /* Difference if we are doing an addition of two matrix or two vectors */
  if( (A.nV != NULL) && (B.nV != NULL) ){ /* two array addition */
    
    for(int i = 0 ; i<C.N_rows*C.N_cols ; i++){  
	C.nV[i] = A.nV[i] + B.nV[i];
    }
    
  }
  else if( (A.nM != NULL) && (B.nM != NULL) ){ /* two matrix addition */
    
    for(int i = 0 ; i < C.N_rows ; i++){
      for(int j = 0 ; j < C.N_cols ; j++){
	C.nM[i][j] = A.nM[i][j] + B.nM[i][j];
      }
    }
  
  }
  else{
    puts("Error in Add_Mat() : Inputs must be of the same range !");
     exit(EXIT_FAILURE);
  }  


  return C;
}

/*********************************************************************/

Matrix Sub_Mat(Matrix A,Matrix B)
/*

  Substract two matrix A and B, and return the result C 

  C = A - B
   
 */
{

  /* Variable declaration */
  Matrix C;
  
  /* Check if it is possible to do the substraction */
  if((A.N_cols != B.N_cols)||(A.N_rows != B.N_rows)){
    puts("Error in Sub_Mat() : Your are trying to add incompatible matrix");
     exit(EXIT_FAILURE);
  }

  /* Allocate output matrix */
  C = MatAlloc(A.N_rows,A.N_cols);

  /* Difference if we are doing an substraction of two matrix
     or two vectors */
  if( (A.nV != NULL) && (B.nV != NULL) ){ /* two array substraction */
    
    for(int i = 0 ; i<C.N_rows*C.N_cols ; i++){  
	C.nV[i] = A.nV[i] - B.nV[i];
    }
    
  }
  else if( (A.nM != NULL) && (B.nM != NULL) ){ /* two matrix substraction */
    
    for(int i = 0 ; i < C.N_rows ; i++){
      for(int j = 0 ; j < C.N_cols ; j++){
	C.nM[i][j] = A.nM[i][j] - B.nM[i][j];
      }
    }
    
  }
  else{
    puts("Error in Sub_Mat() : Inputs must be of the same range !");
     exit(EXIT_FAILURE);
  }  


  return C;
}

/*********************************************************************/

double Norm_Mat(Matrix In,int kind)
/*
  Get the norm of a vector in a Euclidean space R2. Implemented norms :
  - Euclidean norm (kind == 2)
*/
{

  /* Check */
  if( (In.nM != NULL) || (In.n == In.n)){
    puts("Error in Norm_Mat() : The input data is not a vector ! ");
     exit(EXIT_FAILURE);
  }
    
  double Out;
  double aux;

  /* Euclidean norm */
  if(kind == 2){

    aux = 0;
    for(int i = 0 ; i< In.N_rows*In.N_cols ; i++){
      aux += DSQR(In.nV[i]);
    }
    Out = pow(aux,0.5);
  }

  return Out;  
}

/*********************************************************************/

Matrix Eigen_Mat(Matrix In){

  /* Check if we dont have a null matrix */
  if ( (In.nM == NULL) && (In.n != In.n) ){
    if(In.nV != NULL){
      puts("Error in Eigen_Mat() : An array does not have determinant !");
       exit(EXIT_FAILURE);
    }
    puts("Error in Eigen_Mat() : Input matrix = NULL !");    
     exit(EXIT_FAILURE);
  }
  
  /* Check if the matrix is square */
  if(In.N_cols != In.N_rows){
    puts(" Error in Eigen_Mat() : Non square matrix !");
     exit(EXIT_FAILURE);
  }

  int Dim = In.N_cols;
  Matrix Eigen = MatAssign(Dim,Dim,NAN,NULL,NULL);
  Matrix Eigen_Vals;
  Matrix Coeffs;

  switch(Dim){

  case 1 :
     exit(EXIT_FAILURE);
    break;
  case 2 :
    /* Get the coefficients of the charasteristic pol */
    Coeffs = MatAllocZ(1,3);
    Coeffs.nV[0] = 1.0; /* a*x^2 */
    Coeffs.nV[1] = - In.nM[0][0] - In.nM[1][1]; /* b*x */
    Coeffs.nV[2] = In.nM[0][0]*In.nM[1][1] - In.nM[1][0]*In.nM[0][1]; /* c */
    /* Solve the charasteristic pol to the eigenvalues */
    Eigen_Vals = SolvePolynomial(Coeffs);
    FreeMat(Coeffs);
    /* Assign the eigenvalues to the solution */
    Eigen.nV = Eigen_Vals.nV;
    break;

  case 3 :
    /* Get the coefficients of the charasteristic pol */
    Coeffs = MatAllocZ(1,4);
    Coeffs.nV[0] = /* a*x^3 */
      + 1.0; 
    Coeffs.nV[1] = /* b*x^2 */
      - In.nM[0][0]
      - In.nM[1][1]
      - In.nM[2][2]; 
    Coeffs.nV[2] = /* c*x */
      + In.nM[0][0]*In.nM[2][2]
      + In.nM[1][1]*In.nM[2][2]
      + In.nM[0][0]*In.nM[1][1]
      - In.nM[1][0]*In.nM[0][1]
      - In.nM[2][1]*In.nM[1][2]
      - In.nM[2][0]*In.nM[0][2]; 
    Coeffs.nV[2] = /* d */
      - In.nM[0][0]*In.nM[1][1]*In.nM[2][2]
      - In.nM[0][1]*In.nM[1][2]*In.nM[2][0]
      - In.nM[0][2]*In.nM[1][0]*In.nM[2][1]
      + In.nM[1][0]*In.nM[0][1]*In.nM[2][2]
      + In.nM[2][1]*In.nM[1][2]*In.nM[0][0]
      + In.nM[2][0]*In.nM[0][2]*In.nM[1][1]; 
    /* Solve the charasteristic pol to the eigenvalues */
    Eigen_Vals = SolvePolynomial(Coeffs);
    FreeMat(Coeffs);
    /* Assign the eigenvalues to the solution */
    Eigen.nV = Eigen_Vals.nV;
    break;
    
  default :
     exit(EXIT_FAILURE);
  }

  return Eigen;
}

/*********************************************************************/

double Cond_Mat(Matrix In, double TOL)
/*
  Return the conditioning number 
*/
{
  /* Check if we dont have a null matrix */
  if ( (In.nM == NULL) && (In.n != In.n) ){
    if(In.nV != NULL){
      puts("Error in Cond_Mat() : An array does not have determinant !");
       exit(EXIT_FAILURE);
    }
    puts("Error in Cond_Mat() : Input matrix = NULL !");    
     exit(EXIT_FAILURE);
  }
  
  /* Check if the matrix is square */
  if(In.N_cols != In.N_rows){
    puts(" Error in Cond_Mat() : Non square matrix !");
     exit(EXIT_FAILURE);
  }

  /* Remove numerical zeros */
  for(int i = 0 ; i<In.N_rows ; i++){
    for(int j = 0 ; j<In.N_cols ; j++){
      if(fabs(In.nM[i][j])<TOL){
	In.nM[i][j] = 0;
      }
    }
  }
  
  Matrix Eigen = Eigen_Mat(In);
  double max_Eigen = fabs(Eigen.nV[0]);
  double min_Eigen = fabs(Eigen.nV[0]);
    
  for(int i = 1 ; i<Eigen.N_rows ; i++){
    max_Eigen = DMAX(max_Eigen,fabs(Eigen.nV[i]));
    min_Eigen = DMIN(min_Eigen,fabs(Eigen.nV[i]));
  }

  return (double)max_Eigen/min_Eigen;
  
}

/*********************************************************************/

Matrix Get_Lumped_Matrix(Matrix M_in){
  /*
    Get the lumped matrix in a vectorial shape of the input matrix 
   */

  if(M_in.N_cols != M_in.N_rows){
    puts("Error in Get_Lumped_Matrix() : The input matrix is not square ! ");
     exit(EXIT_FAILURE);
  }
  
  Matrix M_out;
  M_out = MatAllocZ(M_in.N_rows,1);

  for(int i = 0 ; i<M_in.N_cols ; i++){
    for(int j = 0 ; j<M_in.N_rows ; j++){
      M_out.nV[i] += M_in.nM[j][i];
    }
  }
  return M_out;  
}

/*********************************************************************/

Matrix Matrix_x_Scalar(Matrix A, double a)
/*!
 * \brief Brief description of Matrix_x_Scalar.
 *        Function to multiply a Matrix with a scalar. 
 *
 *  The parameters for this functions are  :
 * @param A : Input Matrix
 * @param a : Input scalar 
*/
{

  bool Is_Matrix = false;
  bool Is_Vector = false;
  int N_rows = A.N_rows;
  int N_cols = A.N_cols;

  /* Check if its matrix or array */
  if((N_rows > 1) && (N_cols > 1)){
    Is_Matrix = true;
  }
  else if((N_rows == 1) || (N_cols == 1) ){
    Is_Vector = true;
  }

  if(Is_Matrix){ /* Multiply matrix by an scalar */
    for(int i = 0 ; i<N_rows ; i++){
      for(int j = 0 ; j<N_cols ; j++){
	A.nM[i][j] *= a;
      }
    }
  }
  else if(Is_Vector){ /* Multiply vector by an scalar */
    for(int i = 0 ; i<N_cols*N_rows ; i++){
      A.nV[i] *= a;
    }
  }
  else{
        fprintf(stderr,"%s : %s \n",
	    "Error in Matrix_x_Scalar(*,)",
	    "Not a matrix or a vector");
     exit(EXIT_FAILURE);    
  }

  return A;
}


/*********************************************************************/

void get_SVD_Of(Matrix A, Matrix W, Matrix V)
/*
  This routine was adapted from the "Numerical recipies in C".
  Given a matrix A [m x n], this routine computes
  its singular value decomposition, A = U*W*VT. 
  Outputs :
  The matrix U [m x n] replaces A on output. 
  The matrix of singular values W [n x n]. 
  The matrix V (not the transpose VT) is output as V [n x n].   
*/
{
  
}

float pythag(float a, float b)
/*  
    Computes (a^2 + b^2)^0.5 without destructive underflow or overflow. 
*/  
{
  /* Compute the absolute value of "a" and "b" */
  float absa=fabs(a);
  float absb=fabs(b);
  
  if (absa > absb)
    return absa*sqrt(1.0+SQR(absb/absa));
  else
    return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}



/*********************************************************************/
