#include <stdio.h>
#include <stdlib.h>
#include "ToolsLib/TypeDefinitions.h"

/* Function for arrays declaration */
void * Allocate_Array(int SizeArray, int SizeType)
/*
  Function for matrix declaration
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

/* Function to get the determinant of a matrix (max 3x3)  */
double Get_Determinant(Matrix M_in)
/* 
   Inputs : M_in Matrix (Tensor type)
   Outputs : Determinant
*/
{

  /* Check if the matrix is square */
  if(M_in.N_cols != M_in.N_rows){
    puts(" Error in Get_Determinant() : Non square matrix !");
    exit(0);
  }
 
  /* Check if we dont have a null matrix */
  if (M_in.nM == NULL){
    if(M_in.nV != NULL){
      puts("Error in Get_Determinant() : An array does not have determinant !");
      exit(0);
    }
    else{
      puts("Error in Get_Determinant() : Input matrix = NULL !");    
      exit(0);
    }
  }
  
  double Det_out;

  /* Get the determinant */
  switch(M_in.N_cols){

  case 1:
    Det_out = *M_in.nM[0];
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

/* Function to get the inverse of the matrix */
Matrix Get_Inverse(Matrix M_in)
/* 
   Inputs : Matrix in, dimension os the matrix
   Outpus : Matrix out
*/
{ 
  /* Check if we dont have a null matrix */
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

  Matrix M_out;
  M_out.N_cols = M_in.N_cols;
  M_out.N_rows = M_in.N_rows;
  M_out.nM = (double **)Allocate_Matrix(M_out.N_rows,
				       M_out.N_cols,
				       sizeof(double));

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

Matrix Transpose_Mat(Matrix M){

  /* Define and allocate the output matrix */
  Matrix M_T;  
  M_T.N_cols = M.N_rows;
  M_T.N_rows = M.N_cols;

  /* Do it in a different fashion if it is a matrix o an array */
  if ( (M_T.N_cols != 1) && (M_T.N_rows =! 1) ){ /* It is a matrix */
    
    M_T.nM = (double **)Allocate_Matrix(M_T.N_rows,M_T.N_cols,sizeof(double));
    
    for(int i = 0 ; i < M_T.N_rows){
      for(int j = 0 ; j < M_T.N_cols){
	  M_T.nM[i][j] = M.nM[j][i];
      }/* for i*/
    }/* for j */ 
  }
  /* If it is an array you dont have to do anything in the allocation */
  
  return M_T;
}

/*********************************************************************/

Matrix Scalar_prod(Matrix A,Matrix B)
/*
  Multiply two matrix A and B, and return the result C
*/
{

  if(A.N_cols != B.N_rows){
    puts("Error in Mat_mul() : Your are trying to multiply incompatible matrix");
    exit(0);
  }

  /* Allocate output matrix */
  Matrix C;
  double C_aux;
  C.N_rows = A.N_rows;
  C.N_cols = B.N_cols;

  /* If it is array or matrix, allocate and operate in a different fashion */

  if ( (A.nM != NULL) && (B.nM != NULL) ) { /* Scalar product of two matrix */
    /* The result is a matrix */
    C.nM = (double **)Allocate_Matrix(C.N_rows,C.N_cols,sizeof(double));
    
    /* Multiply */
    for(int i = 0 ; i < C.N_rows ; i++){
      for(int j = 0 ; j < C.N_cols ; j++){
	C_aux = 0;
	for(int k = 0 ; k < B.N_rows ; i++){
	  C_aux += A.nM[i][k]*B.nM[k][j];
	}
	C.nM[i][j] = C_aux;
      } /* for j */
    } /* for i */
    
  }
  else if( (A.nV != NULL) && (B.nV != NULL) ){ /* Scalar product of an array by an array */
    /* The result is an scalar */
    puts("Warning : You have an scalar in a matrix type (Not eficient) ! ");
    /* Multiply */
    C_aux = 0;
    for(int k = 0 ; k < A.N_cols ; i++){
      C_aux += A.nV[k]*B.nV[k];
    }
    C.n = C_aux;
    
  }
  else if (( (A.nV != NULL)&&(B.nM != NULL) ) ||
	   ( (A.nM != NULL)&&(B.nV != NULL) )){ /* Scalar product of an array by a matrix */
    /* The result is an array */
    C.nV = (double *)Allocate_Array(C.N_rows*C.N_cols,sizeof(double));
      
    if( (A.nV != NULL)&&(B.nM != NULL) ){ /* Row array */
      for(int i = 0 ; i<B.N_cols ; i++){
	C_aux = 0;
	for(int j = 0 ; j<B.N_rows ; j++){
	  C_aux += A[j]*B[j][i];
	}
	C.nV[i] = C_aux;
      }      
    }
    if( (A.nM != NULL)&&(B.nV != NULL) ){ /* Column array */
      for(int i = 0 ; i<A.N_rows ; i++){
	C_aux = 0;
	for(int j = 0 ; j<A.N_cols ; j++){
	  C_aux += A[i][j]*B[j];
	}
	C.nV[i] = C_aux;
      }      
    }
    
  }
  
  return C;
    
}
  
/*********************************************************************/
  
Matrix Tensorial_prod(Matrix A,Matrix B)
/*

  Tensorial product
  
 */
{
  /* Matrix declaration */
  Matrix C;

  if( (A.nV != NULL)&&(B.nV != NULL) ){ /* Tensorial product between two arrays */
    C.N_rows = A.N_rows;
    C.N_cols = B.N_cols;
    C.nM = (double **)Allocate_Matrix(C.N_rows,C.N_cols,sizeof(double));

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

Matrix Add_Mat(Matrix A,Matrix B)
/*

  Sum two matrix A and B, and return the result C 

  C = A + B
   
 */
{
  /* Check if it is possible to do the addition */
  if((A.N_cols != B.N_cols)||(A.N_rows != B.N_rows)){
    puts("Error in Add_Mat() : Your are trying to add incompatible matrix");
    exit(0);
  }

  /* Allocate output matrix */
  Matrix C;
  double C_aux;
  C.N_rows = A.N_rows;
  C.N_cols = A.N_cols;

  /* Difference if we are doing an addition of two matrix or two vectors */
  if( (A.nV != NULL) && (B.nV != NULL) ){ /* two matrix addition */
    C.nV = (double *)Allocate_Array(C.N_rows*C.N_cols,sizeof(double));
    for(int i = 0 ; i<C.N_rows*C.N_cols ; i++){
	C.nV[i] = A.nV[i] + B.nV[i];
    }
  }
  else if( (A.nT != NULL) && (B.nT != NULL) ){ /* two array addition */
    C.nM = (double **)Allocate_Matrix(C.N_rows,C.N_cols,sizeof(double));
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
  /* Check if it is possible to do the substraction */
  if((A.N_cols != B.N_cols)||(A.N_rows != B.N_rows)){
    puts("Error in Sub_Mat() : Your are trying to substact incompatible matrix");
    exit(0);
  }

  /* Allocate output matrix */
  Matrix C;
  double C_aux;
  C.N_rows = A.N_rows;
  C.N_cols = A.N_cols;

  /* Difference if we are doing an addition of two matrix or two vectors */
  if( (A.nV != NULL) && (B.nV != NULL) ){ /* two matrix substraction */
    C.nV = (double *)Allocate_Array(C.N_rows*C.N_cols,sizeof(double));
    for(int i = 0 ; i<C.N_rows*C.N_cols ; i++){
	C.nV[i] = A.nV[i] - B.nV[i];
    }
  }
  else if( (A.nT != NULL) && (B.nT != NULL) ){ /* two array substraction */
    C.nM = (double **)Allocate_Matrix(C.N_rows,C.N_cols,sizeof(double));
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
