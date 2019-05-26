#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "TypeDefinitions.h"
#include "../Solvers/Solvers.h"
#include "Utils.h"

/* Function for arrays declaration */
void * Allocate_Array(int SizeArray, int SizeType)
/*
  Function for array declaration
  Inputs : Number of rows, number of columns and kind of element (double, integer, ...)
  Outpus : Matrix
*/
{
  
  void * V;  
  V = (void *)malloc(SizeArray*SizeType);  
  if (V == NULL){puts("Error in the array declaration"); exit(0);}  
  return V;
  
}

/*********************************************************************/


/* Function for arrays declaration */
void * Allocate_ArrayZ(int SizeArray, int SizeType)
/*
  Function for array of zeros declaration (double)
  Inputs : Number of rows, number of columns and kind of element (double, integer, ...)
  Outpus : Matrix
*/
{
  
  void * V;
  V = (void*)calloc(SizeArray, SizeType);
  if (V == NULL){puts("Error in the array declaration"); exit(0);}  
  return V;
  
}

/*********************************************************************/

void ** Allocate_Matrix(int NumberRows,int NumberColumns, int SizeType)
/*
  Function for matrix declaration
  Inputs : Number of rows, number of columns and kind of element (double, integer, ...)
  Outpus : Matrix
 */
{

  void ** M;
  M = (void **)malloc((unsigned) NumberRows * sizeof(void *));
  if (M == NULL){puts("Error in matrix declaration"); exit(0);}
  for(int i = 0 ; i<NumberRows ; i++){
    M[i] = malloc((unsigned) NumberColumns*SizeType);
    if (M[i] == NULL){puts("Error in matrix declaration"); exit(0);}
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
  if (M == NULL){puts("Error in matrix declaration"); exit(0);}
  for(int i = 0 ; i<NumberRows ; i++){
    M[i] = (void *)calloc(NumberColumns,SizeType);
    if (M[i] == NULL){puts("Error in matrix declaration"); exit(0);}
  }
  return M;
  
}


/*********************************************************************/

Matrix MatAlloc(int NumberRows,int NumberColumns){

  Matrix M;
  M.N_rows = NumberRows;
  M.N_cols = NumberColumns;
  strcpy(M.Info,"None");

  if( (NumberRows == 1) || (NumberColumns == 1) ){ /* It is an array */
    M.nV = (double *)Allocate_Array(NumberRows*NumberColumns,
				     sizeof(double));
    M.nM = NULL;
    M.n = -999;
  }
  else if( (NumberRows != 1) && (NumberColumns != 1)  ){ /* It is a matrix */
    M.nM = (double **)Allocate_Matrix(NumberRows,
				      NumberColumns,
				      sizeof(double));
    M.nV = NULL;
    M.n = -999;
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
    M.n = -999;
  }
  else if( (NumberRows != 1) && (NumberColumns != 1)  ){ /* It is a matrix */
    M.nM = (double **)Allocate_MatrixZ(NumberRows,
				      NumberColumns,
				      sizeof(double));
    M.nV = NULL;
    M.n = -999;
  }
  
  return M;
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
    exit(0);
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
  if ( (M_in.nM == NULL) && (M_in.n == -999) ){
    if(M_in.nV != NULL){
      puts("Error in Get_Determinant() : An array does not have determinant !");
      exit(0);
    }
    puts("Error in Get_Determinant() : Input matrix = NULL !");    
    exit(0);
  }
  
  /* Check if the matrix is square */
  if(M_in.N_cols != M_in.N_rows){
    puts(" Error in Get_Determinant() : Non square matrix !");
    exit(0);
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
    puts("Error in Get_Determinant() : I am only able to get the determinat of a 3x3 matrix ! ");
    exit(0); 
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
      exit(0);
    }
    else{
      puts("Error in Get_Inverse : Input matrix = NULL !");    
      exit(0);
    }
  }

  /* Check if the Matrix is invertible */

  /* Is the matrix square ? */
  if(M_in.N_cols != M_in.N_rows){
    puts(" Error in Get_Inverse() : Non square matrix !");
    exit(0);
  }
  
  /* Get the determinant of the matrix */
  double Det = Get_Determinant(M_in);
  if(Det == 0){
    puts("Error in Get_Inverse() : Determinant null !");
    exit(0);
  }

  /* Allocate the output */
  Matrix M_out = MatAlloc(M_in.N_rows,M_in.N_cols);
 

  /* Do in a different fashion if it is 2D or 3D */
  if(M_in.N_cols == 2){
    M_out.nM[0][0] = 1/(Det)*M_in.nM[1][1];

    M_out.nM[0][1] = -1/(Det)*M_in.nM[0][1];

    M_out.nM[1][0] = -1/(Det)*M_in.nM[1][0];

    M_out.nM[1][1] = 1/(Det)*M_in.nM[0][0];
  }
  else if(M_in.N_cols == 3){
    M_out.nM[0][0] = 1/(Det)*(M_in.nM[1][1]*M_in.nM[2][2] -
			     M_in.nM[1][2]*M_in.nM[2][1]);
    
    M_out.nM[0][1] = -1/(Det)*(M_in.nM[0][1]*M_in.nM[2][2] -
			      M_in.nM[0][2]*M_in.nM[2][1]);
    
    M_out.nM[0][2] = 1/(Det)*(M_in.nM[0][1]*M_in.nM[1][2] -
			     M_in.nM[0][2]*M_in.nM[1][1]);
    
    M_out.nM[1][0] = -1/(Det)*(M_in.nM[1][0]*M_in.nM[2][2] -
			      M_in.nM[1][2]*M_in.nM[2][0]);
    
    M_out.nM[1][1] = 1/(Det)*(M_in.nM[0][0]*M_in.nM[2][2] -
			     M_in.nM[0][2]*M_in.nM[2][0]);
    
    M_out.nM[1][2] = -1/(Det)*(M_in.nM[0][0]*M_in.nM[1][2] -
			      M_in.nM[0][2]*M_in.nM[1][0]);
    
    M_out.nM[2][0] = 1/(Det)*(M_in.nM[1][0]*M_in.nM[2][1] -
			     M_in.nM[1][1]*M_in.nM[2][0]);
    
    M_out.nM[2][1] = -1/(Det)*(M_in.nM[0][0]*M_in.nM[2][1] -
			      M_in.nM[0][1]*M_in.nM[2][0]);
    
    M_out.nM[2][2] = 1/(Det)*(M_in.nM[0][0]*M_in.nM[1][1] -
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
  if(M.n != -999){
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

  if ( (A.nM != NULL) && (B.nM != NULL) ) { /* Scalar product of two matrix */

    /* Check if the input matrix are not compatible */
    if(A.N_cols != B.N_rows){
      puts("Error in Mat_mul() : Your are trying to multiply incompatible matrix");
      exit(0);
    }
     
    /* The result is a matrix */
    C = MatAlloc(A.N_rows,B.N_cols);   
    
    /* Multiply */
    for(int i = 0 ; i < C.N_rows ; i++){
      for(int j = 0 ; j < C.N_cols ; j++){
	C_aux = 0;
	for(int k = 0 ; k < B.N_rows ; k++){
	  C_aux += A.nM[i][k]*B.nM[k][j];
	}
	C.nM[i][j] = C_aux;
      } /* for j */
    } /* for i */    
  }
  else if( (A.nV != NULL) && (B.nV != NULL) ){ /* Scalar product of an array by an array */

    /* Check if the input matrix are not compatible */
    if(A.N_cols != B.N_rows){
      puts("Error in Mat_mul() : Your are trying to multiply incompatible matrix");
      exit(0);
    }

    /* The result is an scalar */
    /* puts("Warning : You have a scalar in a matrix type (Not eficient) ! "); */
    /* Multiply */
    C_aux = 0;
    for(int k = 0 ; k < A.N_cols ; k++){
      C_aux += A.nV[k]*B.nV[k];
    }
    C.n = C_aux;
  }
  else if (( (A.nV != NULL)&&(B.nM != NULL) ) ||
	   ( (A.nM != NULL)&&(B.nV != NULL) )){ /* Scalar product of an array by a matrix */

    /* Check if the input matrix are not compatible */
    if(A.N_cols != B.N_rows){
      puts("Error in Mat_mul() : Your are trying to multiply incompatible matrix");
      exit(0);
    }
      
    /* The result is an array */
    C = MatAlloc(A.N_rows,B.N_cols);
      
    if( (A.nV != NULL)&&(B.nM != NULL) ){ /* Row array */
      for(int i = 0 ; i<B.N_cols ; i++){
	C_aux = 0;
	for(int j = 0 ; j<B.N_rows ; j++){
	  C_aux += A.nV[j]*B.nM[j][i];
	}
	C.nV[i] = C_aux;
      }
    }
    if( (A.nM != NULL)&&(B.nV != NULL) ){ /* Column array */
      for(int i = 0 ; i<A.N_rows ; i++){
	C_aux = 0;
	for(int j = 0 ; j<A.N_cols ; j++){
	  C_aux += A.nM[i][j]*B.nV[j];
	}
	C.nV[i] = C_aux;
      }
    }
    
  }
  else if( ( (A.nM != NULL)&&(B.n != -999) ) ||
	   ( (A.n != -999)&&(B.nM != NULL) )){ /* Scalar product of a scalar by a matrix */
  
    if(A.nM != NULL){ /* A is a matrix and B a scalar */
      /* The result is a matrix */
      C = MatAlloc(A.N_rows,A.N_cols);
      for(int i = 0 ; i<A.N_rows ; i++){
	for(int j = 0 ; j<A.N_cols ; j++){
	  C.nM[i][j] = A.nM[i][j]*B.n;
	}
      }
    }
    if(B.nM != NULL){ /* A is a scalar and B a matrix */
      /* The result is a matrix */
      C = MatAlloc(B.N_rows,B.N_cols);
      for(int i = 0 ; i<B.N_rows ; i++){
	for(int j = 0 ; j<B.N_cols ; j++){
	  C.nM[i][j] = A.n*B.nM[i][j];
	}
      }
    }
      
  }
  else if( ( (A.nV != NULL)&&(B.n != -999) ) ||
	   ( (A.n != -999)&&(B.nV != NULL) )){ /* Scalar product of a scalar by a array */
  
    if(A.nV != NULL){ /* A is an array and B a scalar */
      /* The result is a matrix */
      C = MatAlloc(A.N_rows,A.N_cols);
      for(int i = 0 ; i<A.N_rows*A.N_cols ; i++){
	C.nV[i] = A.nV[i]*B.n;
      }
    }
    if(B.nM != NULL){ /* A is a scalar and B an array */
      /* The result is a matrix */
      C = MatAlloc(B.N_rows,B.N_cols);
      for(int i = 0 ; i<B.N_rows*B.N_cols ; i++){
	C.nV[i] = A.n*B.nV[i];
      }
    }
  
  }

  return C;
    
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
    exit(0);
  }
  
  return c;
  
}

/*********************************************************************/
  
Matrix Tensorial_prod(Matrix A,Matrix B)
/*
  Tensorial product between A and B, note that the input data is deleted 
  once you have finished the calculus.
 */
{
  /* Matrix declaration */
  Matrix C;

  if( (A.nV != NULL)&&(B.nV != NULL) ){ /* Tensorial product between two arrays */

    if( (A.N_cols != 1) || (B.N_rows !=1) ){
      puts("Error in Tensorial_prod() : Incompatible shape of arrays ");
      puts("Remember that it should be : (nx1) · (1xn)  ");
      exit(0);
    }
    
    /* Allocate matrix */
    C = MatAlloc(A.N_rows,B.N_cols);

    for(int i = 0 ; i<C.N_rows ; i++){
      for(int j = 0 ; j<C.N_cols ; j++){
	C.nM[i][j] = A.nV[i]*B.nV[j];	
      } /* for i */
    } /* for j */
    
  }  /* Tensorial product between two arrays */
  else if( ((A.nM != NULL)&&(B.nV != NULL)) ||
	   ((A.nV != NULL)&&(B.nM != NULL))){
    puts(" Tensorial producto between a matrix and a array is not already defined");
    exit(0);
  }
  else if((A.nM != NULL)&&(B.nM != NULL)){
    puts(" Tensorial producto between a matrix and a matrix is not already defined ");
    exit(0);
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
      exit(0);
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
    exit(0);
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
  double C_aux;
  
  /* Check if it is possible to do the addition */
  if((A.N_cols != B.N_cols)||(A.N_rows != B.N_rows)){
    puts("Error in Add_Mat() : Your are trying to add incompatible matrix");
    exit(0);
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
    exit(0);
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
  double C_aux;
  
  /* Check if it is possible to do the substraction */
  if((A.N_cols != B.N_cols)||(A.N_rows != B.N_rows)){
    puts("Error in Sub_Mat() : Your are trying to add incompatible matrix");
    exit(0);
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
    exit(0);
  }  


  return C;
}


double Norm_Mat(Matrix In,int kind)
/*
  Get the norm of a vector in a Euclidean space R2. Implemented norms :
  - Euclidean norm (kind == 2)
*/
{

  /* Check */
  if( (In.nM != NULL) || (In.n != -999)){
    puts("Error in Norm_Mat() : The input data is not a vector ! ");
    exit(0);
  }
    
  double Out;
  double aux;

  /* Euclidean norm */
  if(kind == 2){

    aux = 0;
    for(int i = 0 ; i< In.N_rows*In.N_cols ; i++){
      aux += In.nV[i]*In.nV[i] ;
    }
    Out = pow(aux,0.5);
  }

  return Out;  
}

/*********************************************************************/

int InOut_Poligon(Matrix X_Point, Matrix Poligon)
/*
  Check if a point is or not (1/0) inside of a poligon.
  Inputs :
  - X_Point : Coordinates of the point 
  - Poligon : Coordinates of the vertex 0,1,....,n,0
*/
{
  /* By default, we suppose that the point is in the poligon */
  int InOut = 1;

  Matrix a = MatAllocZ(3,1);
  Matrix b = MatAllocZ(3,1);
  Matrix c;
  Matrix n;
  Matrix nxc;

  /* Get the normal vector */
  a.nV[0] = Poligon.nM[1][0] - Poligon.nM[0][0];
  a.nV[1] = Poligon.nM[1][1] - Poligon.nM[0][1];
  a.nV[2] = Poligon.nM[1][2] - Poligon.nM[0][2];
  b.nV[0] = Poligon.nM[Poligon.N_rows-1][0] - Poligon.nM[0][0];
  b.nV[1] = Poligon.nM[Poligon.N_rows-1][1] - Poligon.nM[0][1];
  b.nV[2] = Poligon.nM[Poligon.N_rows-1][2] - Poligon.nM[0][2];
  n = Vectorial_prod(a,b);
  n.N_rows = 1;
  n.N_cols = 3;

  /* Fill a and b for the First search */
  a.nV[0] = Poligon.nM[0][0] - Poligon.nM[Poligon.N_rows-1][0];
  a.nV[1] = Poligon.nM[0][1] - Poligon.nM[Poligon.N_rows-1][1];
  a.nV[2] = Poligon.nM[0][2] - Poligon.nM[Poligon.N_rows-1][2];

  b.nV[0] = X_Point.nV[0] - Poligon.nM[Poligon.N_rows-1][0];
  b.nV[1] = X_Point.nV[1] - Poligon.nM[Poligon.N_rows-1][1];
  b.nV[2] = X_Point.nV[2] - Poligon.nM[Poligon.N_rows-1][2];

  
  for(int i = 0 ; i<Poligon.N_rows-1 ; i++){

    c = Vectorial_prod(a,b);
    nxc = Scalar_prod(n,c);

    if(nxc.n < 0){
      InOut = 0;
      break;
    }
    
    a.nV[0] = Poligon.nM[i+1][0] - Poligon.nM[i][0];
    a.nV[1] = Poligon.nM[i+1][1] - Poligon.nM[i][1];
    a.nV[2] = Poligon.nM[i+1][2] - Poligon.nM[i][2];

    b.nV[0] = X_Point.nV[0] - Poligon.nM[i][0];
    b.nV[1] = X_Point.nV[1] - Poligon.nM[i][1];
    b.nV[2] = X_Point.nV[2] - Poligon.nM[i][2];
    

    
  }

  free(a.nV);
  free(b.nV);
  free(c.nV);

  return InOut;
}


/*********************************************************************/

Matrix Newton_Rapson(Matrix(* Function)(Matrix, Matrix),Matrix Parameter_F,
		     Matrix(* Jacobian)(Matrix, Matrix),Matrix Parameter_J,
		     Matrix Y,Matrix X)

/*
  Newton-Rapson method to solve non-linear sistems of equations :
  Y = Y(X) -> We solve -> F(X) =  Y - Y(X) = 0
  F(X + DeltaX) = F(X) + J(X)*DeltaX = 0 ---> DeltaX = - J(X)^{-1}*F(X)
  Inputs :
  - Y : Value of the function
  - Function(X,Parameter_F) : Pointer to function to solve
  - Parameter_F : F function optional parameters
  - Jacobian(X,Parameter_J) : Pointer to the jacobian of the function 
  - Parameter_J : Jacobian optional parameters
  - X : Initial value of the objetive 
*/
{

  /* Auxiliar variables */
  Matrix F_n;
  Matrix J_n;
  Matrix J_m1_n;
  Matrix DeltaX;
  double TOL_NormDeltaX = pow(10,-4);
  double NormDeltaX = pow(10,4);
  int Num_Iter = 20;  
  int Iter_i = 0;
  int Bool;

  /* 0º Check the convergence criterium */
  while( (NormDeltaX > TOL_NormDeltaX )
	 && (Iter_i < Num_Iter) ){

    /* 1º Evaluate F in X0 and get the negative value */
    F_n = Function(X,Parameter_F);
    for(int i = 0 ; i<F_n.N_cols*F_n.N_rows ; i++){
      F_n.nV[i] = Y.nV[i] - F_n.nV[i];
    }
    
    /* 2º Get the jacobian matrix in X0 */
    J_n = Jacobian(X,Parameter_J);
    /* Implement the numerical solution of the Jacobian for cases where the 
     Jacobian is not easy to derive */

    /* 3º Solve the sistem J(X0)*DeltaX = F(X0) */
    Bool = J_n.N_cols>3;
    switch(Bool){
    case 0 : /* If the size of the Jacobian is less than 4, use analitical */
      J_m1_n = Get_Inverse(J_n);
      DeltaX = Scalar_prod(J_m1_n,F_n);
      free(J_m1_n.nM);
      break;
    case 1 : /* If the size of the Jacobian is great than 4, use numerical */
      DeltaX = Jacobi_Conjugate_Gradient_Method(J_n,F_n,DeltaX);
      break;
    default :
      exit(0);
    }

    /* 4º Update the variables of the convergence criterium */
    NormDeltaX = Norm_Mat(DeltaX,2);
    Iter_i++;

    /* 4º Update the solution */
    X = Incr_Mat(X,DeltaX);
    
    /* 5º Free memory for general variables */
    free(J_n.nM);
    free(DeltaX.nV);
    free(F_n.nV);
    
  }
  
  /* 6º Return X */
  return X;
  
}

/*********************************************************************/

Matrix Get_Lumped_Matrix(Matrix M_in){
  /*
    Get the lumped matrix in a vectorial shape of the input matrix 
   */

  if(M_in.N_cols != M_in.N_rows){
    puts("Error in Get_Lumped_Matrix() : The input matrix is not square ! ");
    exit(0);
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