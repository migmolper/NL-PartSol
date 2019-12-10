#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h> 
#include "../GRAMS/grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

/*********************************************************************/

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

  if( (NumberRows == 1) || (NumberColumns == 1) ){ /* It is an array */
    M.nV = (double *)Allocate_Array(NumberRows*NumberColumns,
				     sizeof(double));
    M.nM = NULL;
    M.n = NAN;
  }
  else if( (NumberRows != 1) && (NumberColumns != 1)  ){ /* It is a matrix */
    M.nM = (double **)Allocate_Matrix(NumberRows,
				      NumberColumns,
				      sizeof(double));
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
      exit(0);
    }
    free(Input.nV);
  }
  else if( (Input.N_cols > 1) && (Input.N_rows > 1) ){ /* It is a matrix */
    if(Input.nM == NULL){
      printf("Error in FreeMat() : This matrix is NULL !!! \n");
      exit(0);
    }
    for(int i = 0 ; i<Input.N_rows ; i++){
      free(Input.nM[i]);
    }
    free(Input.nM);
  }
  else{
    printf("Error in FreeMat() : You are trying to free a scalar !!! \n");
    exit(0);
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
    Nvalues = MINVAL(NumRows,Nvalues);
    printf("%s : \n",In.Info);
    for(int i = 0 ; i<Nvalues ; i++){
      printf("%f \n",In.nV[i]);
    }
  }
  
  if( In.nM != NULL ){
    
    printf("%s : \n",In.Info);
    for(int i = 0 ; i<MAXVAL(In.N_cols,In.N_rows) ; i++){
      for(int j = 0 ; j<MINVAL(In.N_cols,In.N_rows) ; j++){
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


double StatsDouMatrix(double * In, int NumElems, char * OutChar)
/*
  Print some stats of the input matrix
 */
{
  double MaxVal;
  double MinVal;
  double MeanVal;
  double SumVal;
  double StdDesvVal;
  double OutVal;

  /* Do the stats for an array */
  if(In != NULL){

    /* Initialice the stats */
    MaxVal = In[0];
    MinVal = In[0];
    MeanVal = In[0];
    StdDesvVal = 0;
    SumVal = 0;

    /* Get the MaxVal, the MinVal and the MeanVal */
    for(int i = 0 ; i<NumElems ; i++){
      MaxVal = MAXVAL(MaxVal,In[i]);
      MinVal = MINVAL(MinVal,In[i]);
      SumVal += In[i];      
    }
    MeanVal = (double)SumVal/NumElems;

    /* Get the standard desviation */
    for(int i = 0 ; i<NumElems ; i++){
      StdDesvVal += (In[i] - MeanVal)*(In[i] - MeanVal);
    }
    StdDesvVal *= (double)1/(NumElems-1);
    StdDesvVal = pow(StdDesvVal,0.5);
  }
  else{
    puts("Error in StatMatrix() : The input pointer is empty !! \n");
  }


  if(strcmp(OutChar,"PrintStats") == 0){
    printf("The stats are : \n");
    printf("\t Maximum : %f \n",MaxVal);
    printf("\t Minimum : %f \n",MinVal);
    printf("\t The sum of all the terms : %f \n",SumVal);
    printf("\t Mean : %f \n",MeanVal);
    printf("\t Standard deviation : %f \n",StdDesvVal);
    OutVal = 0;
  }
  else if(strcmp(OutChar,"MaxVal") == 0) OutVal = MaxVal;
  else if(strcmp(OutChar,"MinVal") == 0) OutVal = MinVal;
  else if(strcmp(OutChar,"SumVal") == 0) OutVal = SumVal;
  else if(strcmp(OutChar,"MeanVal") == 0) OutVal = MeanVal;
  else if(strcmp(OutChar,"StdDesvVal") == 0) OutVal = StdDesvVal;
  else{
    printf("Error in StatsIntMatrix() : Wrong imput parameters !!! \n");
    OutVal = NAN;
  }
  
  return OutVal;
}

/*********************************************************************/

double StatsIntMatrix(int * In, int NumElems, char * OutChar)
/*
  Print some stats of the input matrix
 */
{
  double MaxVal;
  double MinVal;
  double MeanVal;
  double SumVal;
  double StdDesvVal;
  double OutVal;

  printf("OutChar : %s\n",OutChar);

  /* Do the stats for an array */
  if(In != NULL){

    /* Initialice the stats */
    MaxVal = (double)In[0];
    MinVal = (double)In[0];
    MeanVal = (double)In[0];
    StdDesvVal = 0.0;
    SumVal = 0.0;

    /* Get the MaxVal, the MinVal and the MeanVal */
    for(int i = 0 ; i<NumElems ; i++){
      MaxVal = MAXVAL(MaxVal,(double)In[i]);
      MinVal = MINVAL(MinVal,(double)In[i]);
      SumVal += (double)In[i];      
    }
    MeanVal = (double)SumVal/NumElems;

    /* Get the standard desviation */
    for(int i = 0 ; i<NumElems ; i++){
      StdDesvVal += ((double)In[i] - MeanVal)*
	((double)In[i] - MeanVal);
    }
    StdDesvVal *= (double)1/(NumElems-1);
    StdDesvVal = pow(StdDesvVal,0.5);
  }
  else{
    puts("Error in StatMatrix() : The input pointer is empty !! \n");
  }

  if(strcmp(OutChar,"PrintStats") == 0){
    printf("The stats are : \n");
    printf("\t Maximum : %f \n",MaxVal);
    printf("\t Minimum : %f \n",MinVal);
    printf("\t The sum of all the terms : %f \n",SumVal);
    printf("\t Mean : %f \n",MeanVal);
    printf("\t Standard deviation : %f \n",StdDesvVal);
    OutVal = 0;
  }
  else if(strcmp(OutChar,"MaxVal") == 0) OutVal = MaxVal;
  else if(strcmp(OutChar,"MinVal") == 0) OutVal = MinVal;
  else if(strcmp(OutChar,"SumVal") == 0) OutVal = SumVal;
  else if(strcmp(OutChar,"MeanVal") == 0) OutVal = MeanVal;
  else if(strcmp(OutChar,"StdDesvVal") == 0) OutVal = StdDesvVal;
  else{
    printf("Error in StatsIntMatrix() : Wrong imput parameters !!! \n");
    OutVal = NAN;
  }
  
  return OutVal;
  
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
  if ( (M_in.nM == NULL) && (M_in.n != M_in.n) ){
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
    M_out.nM[0][0] = (double)1/(Det)*M_in.nM[1][1];

    M_out.nM[0][1] = -(double)1/(Det)*M_in.nM[0][1];

    M_out.nM[1][0] = -(double)1/(Det)*M_in.nM[1][0];

    M_out.nM[1][1] = (double)1/(Det)*M_in.nM[0][0];
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

  if ( (A.nM != NULL) && (B.nM != NULL) ) { /* Scalar product of two matrix */

    /* Check if the input matrix are not compatible */
    if(A.N_cols != B.N_rows){
      puts("Error in Scalar_prod() : Your are trying to multiply incompatible matrix");
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
      puts("Error in Scalar_prod() : Your are trying to multiply incompatible matrix");
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
      puts("Error in Scalar_prod() : Your are trying to multiply incompatible matrix");
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
  else if( ( (A.nM != NULL)&&(B.n == B.n) ) ||
	   ( (A.n == A.n)&&(B.nM != NULL) )){ /* Scalar product of a scalar by a matrix */
  
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
  else if( ( (A.nV != NULL)&&(B.n == B.n) ) ||
	   ( (A.n == A.n)&&(B.nV != NULL) )){ /* Scalar product of a scalar by a array */
  
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
  Tensorial product between A and B.
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
    
  }  /* Tensorial product between an array and a matrix */
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
  if( (In.nM != NULL) || (In.n == In.n)){
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

Matrix Eigen_Mat(Matrix In){

  /* Check if we dont have a null matrix */
  if ( (In.nM == NULL) && (In.n != In.n) ){
    if(In.nV != NULL){
      puts("Error in Eigen_Mat() : An array does not have determinant !");
      exit(0);
    }
    puts("Error in Eigen_Mat() : Input matrix = NULL !");    
    exit(0);
  }
  
  /* Check if the matrix is square */
  if(In.N_cols != In.N_rows){
    puts(" Error in Eigen_Mat() : Non square matrix !");
    exit(0);
  }

  int Dim = In.N_cols;
  Matrix Eigen = MatAssign(Dim,Dim,NAN,NULL,NULL);
  Matrix Eigen_Vals;
  Matrix Coeffs;

  switch(Dim){

  case 1 :
    exit(0);
    
  case 2 :
    /* Get the coefficients of the charasteristic pol */
    Coeffs = MatAllocZ(1,3);
    Coeffs.nV[0] = 1;
    Coeffs.nV[1] = - In.nM[0][0] - In.nM[1][1];
    Coeffs.nV[2] = In.nM[0][0]*In.nM[1][1] - In.nM[1][0]*In.nM[0][1];
    /* Solve the charasteristic pol to the eigenvalues */
    Eigen_Vals = SolvePolynomial(Coeffs);
    FreeMat(Coeffs);
    /* Assign the eigenvalues to the solution */
    Eigen.nV = Eigen_Vals.nV;
    break;

  case 3 :
    exit(0);

  default :
    exit(0);

  }

  return Eigen;
}

/*********************************************************************/

double Cond_Mat(Matrix In)
/*
  Return the conditioning number 
*/
{
  /* Check if we dont have a null matrix */
  if ( (In.nM == NULL) && (In.n != In.n) ){
    if(In.nV != NULL){
      puts("Error in Cond_Mat() : An array does not have determinant !");
      exit(0);
    }
    puts("Error in Cond_Mat() : Input matrix = NULL !");    
    exit(0);
  }
  
  /* Check if the matrix is square */
  if(In.N_cols != In.N_rows){
    puts(" Error in Cond_Mat() : Non square matrix !");
    exit(0);
  }
  
  Matrix Eigen = Eigen_Mat(In);
  int Dim_In = In.N_rows;

  return fabs((double)Eigen.nV[0]/Eigen.nV[Dim_In-1]);
  
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

/*********************************************************************/
