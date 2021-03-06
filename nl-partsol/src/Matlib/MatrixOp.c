#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>

#endif

/*********************************************************************/

/* Function for arrays declaration */
void * Allocate_Array(int SizeArray,int SizeType)
/*
  Function for array declaration
  Inputs : Number of rows, number of columns and kind of element 
  (double, integer, ...)
  Outpus : Matrix
*/
{
  void * V;
  
  V = (void *)malloc(SizeArray*SizeType);
  
  if (V == NULL)
    {
      puts("Error in the array declaration");
      exit(EXIT_FAILURE);
    }
  
  return V;
  
}

/*********************************************************************/


/* Function for arrays declaration */
void * Allocate_ArrayZ(int SizeArray,int SizeType)
/*
  Function for array of zeros declaration (double)
  Inputs : Number of rows, number of columns and kind of element 
  (double, integer, ...)
  Outpus : Matrix
*/
{
  
  void * V;

  V = (void*)calloc(SizeArray, SizeType);
  
  if (V == NULL)
    {
      puts("Error in the array declaration");
      exit(EXIT_FAILURE);
    }  

  return V;
}

/*********************************************************************/

void ** Allocate_Matrix(int NumberRows,int NumberColumns,int SizeType)
/*
  Function for matrix declaration
  Inputs : Number of rows, number of columns and kind of element 
  (double, integer, ...)
  Outpus : Matrix
 */
{

  void ** M = (void **)malloc((unsigned) NumberRows * sizeof(void *));

  if (M == NULL)
    {
      puts("Error in matrix declaration");
      exit(EXIT_FAILURE);
    }
  
  for(int i = 0 ; i<NumberRows ; i++)
    {
      
      M[i] = malloc((unsigned) NumberColumns*SizeType);
      
      if (M[i] == NULL)
	{
	  puts("Error in matrix declaration");
	  exit(EXIT_FAILURE);
	}
      
    }
  
  return M;
}


/*********************************************************************/

void ** Allocate_MatrixZ(int NumberRows,int NumberColumns,int SizeType)
/*
  Function for matrix of zeros declaration (double)
  Inputs : Number of rows, number of columns.
  Outpus : Matrix
 */
{

  void ** M = (void **)calloc(NumberRows,sizeof(void *));
  
  if (M == NULL)
    {
      puts("Error in matrix declaration");
      exit(EXIT_FAILURE);
    }

  for(int i = 0 ; i<NumberRows ; i++)
    {
      M[i] = (void *)calloc(NumberColumns,SizeType);
      
      if (M[i] == NULL)
	{
	  puts("Error in matrix declaration");
	  exit(EXIT_FAILURE);
	}
      
    }
  
  return M;  
}


/*********************************************************************/

Matrix alloc__MatrixLib__(int NumberRows,int NumberColumns)
{

      
  Matrix M;
  M.N_rows = NumberRows;
  M.N_cols = NumberColumns;

  /* Reservar memoria consecutiva */
  M.nV = (double *)malloc(NumberRows*NumberColumns*sizeof(double));

  if (M.nV == NULL)
    {
      puts("Error in alloc__MatrixLib__ : Out of memory");
      exit(EXIT_FAILURE);
    }

  /* It is a matrix generate a table */
  if((NumberRows > 1) && (NumberColumns > 1))
    {
      M.nM = (double **)malloc(NumberRows*sizeof(double *));
      for(int i = 0 ; i<NumberRows ; i++)
  	{
  	  M.nM[i] = M.nV + i*NumberColumns;
  	}
    }
  
  return M;
}

/*********************************************************************/
Matrix allocZ__MatrixLib__(int NumberRows,int NumberColumns)
/*
  Allocate the matrix structure with zeros 
*/
{
    
  Matrix M;
  M.N_rows = NumberRows;
  M.N_cols = NumberColumns;

  /* Reservar memoria consecutiva */
  M.nV = (double *)calloc(NumberRows*NumberColumns,sizeof(double));

  if (M.nV == NULL)
    {
      puts("Error in allocZ__MatrixLib__ : Out of memory");
      exit(EXIT_FAILURE);
    }

  /* It is a matrix generate a table */
  if((NumberRows > 1) && (NumberColumns > 1))
    {
      M.nM = (double **)malloc(NumberRows*sizeof(double *));
      for(int i = 0 ; i<NumberRows ; i++)
  	{
  	  M.nM[i] = M.nV + i*NumberColumns;
  	}
    }
    
  return M;
}

/*********************************************************************/
Matrix memory_to_matrix__MatrixLib__(int N_rows,int N_cols,double * nV)
/*
  Define general parameters of a matrix for auxiliar
  matrix.
*/
{
  Matrix M;
  M.N_rows = N_rows;
  M.N_cols = N_cols;
  M.nV = nV;
  
  return M;
}

/*********************************************************************/

void free__MatrixLib__(Matrix Input)
/*
  Free memory
*/
{

  int Columns = Input.N_cols;
  int Rows = Input.N_rows;

  /* Free consecutive memory */
  free(Input.nV);

  /* If it is a matrix, free the pointer table */
  if((Columns > 1) && (Rows > 1))
    {
      free(Input.nM);
    }  
}


/*********************************************************************/

Matrix copy__MatrixLib__(Matrix In)
/* 
   Copy the input matrix in a auxiliar matrix.
 */
{

  int N_rows = In.N_rows;
  int N_cols = In.N_cols;

  Matrix Out;
  Out.N_rows = N_rows;
  Out.N_cols = N_cols;
  Out.nV = (double *)calloc(N_rows*N_cols,sizeof(double));
  memcpy(Out.nV, In.nV, N_rows*N_cols*sizeof(double));

  /* It is a matrix generate a table */
  if((N_rows > 1) && (N_cols > 1))
  {
    Out.nM = (double **)malloc(N_rows*sizeof(double *));
    for(int i = 0 ; i<N_rows ; i++)
    {
      Out.nM[i] = Out.nV + i*N_cols;
    }
  }

  return Out;
}


/*********************************************************************/

double I3__MatrixLib__(Matrix In)
/* 
   Function to get the determinant of a matrix (max 3x3) :
   Inputs : In Matrix (Tensor type)
   Outputs : Determinant
*/
{

  int Rows = In.N_rows;
  int Columns = In.N_cols;
  double Det_out;
  
  /* Check if we dont have a null matrix */
  if ((Rows>1) && (Columns>1) && (Columns == Rows))
    {     
    
      /* Get the determinant */
      switch(In.N_cols){

      case 1:
	Det_out = In.n;
	break;

      case 2:
	Det_out =
	  In.nM[0][0]*In.nM[1][1] -
	  In.nM[0][1]*In.nM[1][0];
	break;
    
      case 3:
	Det_out =
	  In.nM[0][0]*In.nM[1][1]*In.nM[2][2] -
	  In.nM[0][0]*In.nM[1][2]*In.nM[2][1] + 
	  In.nM[0][1]*In.nM[1][2]*In.nM[2][0] - 
	  In.nM[0][1]*In.nM[1][0]*In.nM[2][2] + 
	  In.nM[0][2]*In.nM[1][0]*In.nM[2][1] - 
	  In.nM[0][2]*In.nM[1][1]*In.nM[2][0] ; 
	break;

      default :
	fprintf(stderr,"%s : %s \n",
		"Error in I3__MatrixLib__()","Max size allowded 3x3 !");
	exit(EXIT_FAILURE); 
      }

    }
  else
    {
      fprintf(stderr,"%s : %s \n",
	      "Error in I3__MatrixLib__()","Should be a square matrix !");
      exit(EXIT_FAILURE);
    }
  
  return Det_out;
}


/*********************************************************************/

Matrix inverse__MatrixLib__(Matrix A)
/* 
   Function to get the inverse of the matrix :
   Inputs : Matrix in, dimension os the matrix
   Outpus : Matrix out
*/
{

  int INFO;
  int N_rows = A.N_rows;
  int N_cols = A.N_cols;
  int LDA = IMAX(N_rows,N_cols);
  int LWORK = IMAX(1,N_rows);
  double * WORK_DGETRI = (double *)Allocate_Array(IMAX(1,LWORK),sizeof(double));
  int * IPIV = (int *)Allocate_Array(IMIN(N_rows,N_cols),sizeof(int));
  int * IWORK_RCOND = (int *)Allocate_Array(N_rows,sizeof(int));

  Matrix Am1 = copy__MatrixLib__(A);

  /* Compute the LU factorization of the input matrix */
  dgetrf_(&N_rows,&N_cols,Am1.nV,&LDA,IPIV,&INFO);

  if(INFO != 0)
  {
    if(INFO < 0)
    {
      printf("%s : \n","Error in inverse__MatrixLib__()");
      printf("the %i-th argument of dgetrf_ had an illegal value",abs(INFO));
    }
    else if(INFO > 0)
    {
      printf("%s :\n","Error in inverse__MatrixLib__()");
      printf(" A(%i,%i) %s \n %s \n %s \n %s \n",INFO,INFO,"is exactly zero. The factorization",
        "has been completed, but the factor A is exactly",
        "singular, and division by zero will occur if it is used",
        "to solve a system of equations.");
    }    
    exit(EXIT_FAILURE);
  }

  /* Compute the inverse */
  dgetri_(&N_rows, Am1.nV, &LDA, IPIV, WORK_DGETRI, &LWORK, &INFO); 

  if(INFO != 0)
  {
    if(INFO < 0)
    {
      printf("%s : \n","Error in inverse__MatrixLib__()");
      printf("the %i-th argument of dgetrf_ had an illegal value",abs(INFO));
    }
    else if(INFO > 0)
    {
      printf("%s :\n","Error in inverse__MatrixLib__()");
      printf(" A(%i,%i) %s \n %s \n %s \n %s \n",INFO,INFO,"is exactly zero. The factorization",
        "has been completed, but the factor A is exactly",
        "singular, and division by zero will occur if it is used",
        "to solve a system of equations.");
    }    
    exit(EXIT_FAILURE);
  }

  free(WORK_DGETRI);
  free(IPIV);
  free(IWORK_RCOND);
 
  /* Return the inverse matrix */
  return Am1;
}

/*********************************************************************/

Matrix transpose__MatrixLib__(Matrix M)
/* 
   It gives to you the transpose of a matrix or an array matrix.
*/
{

  int Rows = M.N_rows;
  int Columns = M.N_cols;
  
  /* Define and allocate the output matrix */
  Matrix M_T = allocZ__MatrixLib__(Columns,Rows);  

  /* It is a matrix */
  if((Columns > 1) && (Rows > 1))
    {
      if(M.nM == NULL)
	{
	  fprintf(stderr,"%s : %s \n",
		  "Error in transpose__MatrixLib__()","NULL input !");
	  exit(EXIT_FAILURE);    
	}
      
      for(int i = 0 ; i < Columns  ; i++)
	{	  
	  for(int j = 0 ; j < Rows ; j++)
	    {
	      M_T.nM[i][j] = M.nM[j][i];
	    }
	}
    }

  /* It is an array */
  else if(((Columns == 1) && (Rows > 1)) || ((Columns > 1) && (Rows == 1)))
    {
      if(M.nV == NULL)
	{
	  fprintf(stderr,"%s : %s \n",
		  "Error in transpose__MatrixLib__()","NULL input !");
	  exit(EXIT_FAILURE);    
	}
	  
      for(int i = 0 ; i < Rows*Columns ; i++)
	{
	  M_T.nV[i] = M.nV[i];
	}
    }

  /* Fail */
  else
    {
      fprintf(stderr,"%s : %s \n",
	      "Error in transpose__MatrixLib__()","Wrong input !");
      exit(EXIT_FAILURE);      
    }
  
  return M_T;
}

/*********************************************************************/

static Matrix get_A_x_B_Mat(Matrix A, Matrix B)
/* 
   Compute the matrix product of two matrix 
*/
{

  /* The result is a matrix */
  Matrix A_x_B = allocZ__MatrixLib__(A.N_rows,B.N_cols);
  int Rows = A.N_rows;
  int Columns = B.N_cols;
  int Aux = B.N_rows;
     
  /* Multiply */
  for(int i = 0 ; i < Rows  ; i++)
    {
      for(int j = 0 ; j < Columns  ; j++)
	{
	  for(int k = 0 ; k < Aux  ; k++)
	    {
	      A_x_B.nM[i][j] += A.nM[i][k]*B.nM[k][j];
	    }
	}
    } 

  return A_x_B;  
}

/*********************************************************************/

static Matrix get_A_x_b_Mat(Matrix A, Matrix b)
/* 
   Compute the matrix product of a matrix by an array
*/
{
  /* The result is an array */
  Matrix A_x_b = allocZ__MatrixLib__(A.N_rows,b.N_cols);
  int Rows = A.N_rows;
  int Columns = A.N_cols;
  
  for(int i = 0 ; i<Rows ; i++)
    {
      for(int j = 0 ; j<Columns ; j++)
	{
	  A_x_b.nV[i] += A.nM[i][j]*b.nV[j];
	}
    }
    
  return A_x_b;  
}

/*********************************************************************/

static Matrix get_a_x_B_Mat(Matrix a, Matrix B)
/* 
   Compute the matrix product of a matrix by an array
*/
{
  /* The result is an array */
  Matrix a_x_B = allocZ__MatrixLib__(a.N_rows,B.N_cols);
  int Rows = B.N_rows;
  int Columns = B.N_cols;
  
  for(int i = 0 ; i<Columns ; i++){
    for(int j = 0 ; j<Rows ; j++){
      a_x_B.nV[i] += a.nV[j]*B.nM[j][i];
    }
  }
    
  return a_x_B;  
}

/*********************************************************************/

static Matrix Matrix_x_Scalar(Matrix A, double a)
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


Matrix matrix_product__MatrixLib__(Matrix A,Matrix B)
/*
  Multiply two matrix A and B, and return the result C.
*/
{

  /* Variable declaration output matrix */
  Matrix C;
  double C_aux;

  /* If it is array or matrix, allocate and operate in a different fashion */

  /* Scalar product of two matrix */
  if ((A.N_cols == B.N_rows) && (A.N_rows > 1) && (B.N_cols > 1))
    {      
      C = get_A_x_B_Mat(A,B);
    }
    
  /* Scalar product of a matrix by an array */
  else if ((A.N_cols == B.N_rows) && (A.N_rows > 1) && (B.N_cols == 1))
    {       
      C = get_A_x_b_Mat(A, B);
    }

  /* Scalar product of an array by a matrix */
  else if ((A.N_cols == B.N_rows) && (A.N_rows == 1) && (B.N_cols > 1))
    {
      C = get_a_x_B_Mat(A, B);
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
	     "Error in matrix_product__MatrixLib__()","Incompatible product",
	     A.N_rows,A.N_cols,B.N_rows,B.N_cols);
      exit(EXIT_FAILURE);
    }
  
  return C;
  
}

/*********************************************************************/

double scalar_product__MatrixLib__(Matrix a, Matrix b)
/*
  Scalar product of two arrays
 */
{
  double a_dot_b = 0;
  
  int Size_a = a.N_cols*a.N_rows;
  int Size_b = b.N_cols*b.N_rows;

  if(Size_a == Size_b)
    {
      for(int i = 0 ; i < Size_a ; i++)
	{
	  a_dot_b += a.nV[i]*b.nV[i];
	}
      
    }
  else
    {
      printf("%s : %s -> <[%i , %i] , [%i , %i]> \n",
	     "Error in scalar_product__MatrixLib__()","Incompatible product",
	     a.N_rows,a.N_cols,b.N_rows,b.N_cols);
      exit(EXIT_FAILURE);
    }

  return a_dot_b;  
}

/*********************************************************************/

Matrix vectorial_product__MatrixLib__(Matrix a, Matrix b){
  
  Matrix c;

  if((a.N_rows == 3) &&
     (b.N_rows == 3) &&
     (a.N_cols == 1) &&
     (b.N_cols == 1))
    {
      c = allocZ__MatrixLib__(3,1);
      c.nV[0] = a.nV[1]*b.nV[2]-a.nV[2]*b.nV[1];
      c.nV[1] = a.nV[2]*b.nV[0]-a.nV[0]*b.nV[2];
      c.nV[2] = a.nV[0]*b.nV[1]-a.nV[1]*b.nV[0];
    }
  else
    {
      fprintf(stderr,"%s : %s \n",
	      "Error in vectorial_product__MatrixLib__()","Incompatible shape of arrays !");
      exit(EXIT_FAILURE);
    }
  
  return c;
  
}

/*********************************************************************/
  
Matrix dyadic_product__MatrixLib__(Matrix A,Matrix B)
/*
  Tensorial product between A and B.
*/
{
  /* Matrix declaration */
  Matrix C;
  int Rows = A.N_rows;
  int Columns = B.N_cols;

  /* Tensorial product between two arrays */
  if (((Columns == 1) && (Rows > 1)) ||
      (Columns > 1) && (Rows == 1))
    {
 
      /* Allocate matrix */
      C = allocZ__MatrixLib__(Rows,Columns);

      /* Compute the tensorial product */
      for(int i = 0 ; i<Rows ; i++)
	{
	  for(int j = 0 ; j<Columns ; j++)
	    {
	      C.nM[i][j] = A.nV[i]*B.nV[j];	
	    } 
	} 
    
    }

  /* Fail */
  else
    {
      fprintf(stderr,"%s : %s \n",
	      "Error in dyadic_product__MatrixLib__()","Wrong input matrix !");
      exit(EXIT_FAILURE);
    }
  
  return C;
}
 
/*********************************************************************/

Matrix increment__MatrixLib__(Matrix A, Matrix Incr)
/*
  Increment term by term the Matrix A with Incr
*/
{
  
  /* Check if the dimensions mismach */
  if((A.N_cols == Incr.N_cols) && (A.N_rows == Incr.N_rows)){

    int Columns = Incr.N_cols;
    int Rows = Incr.N_rows;

    /* Is a matrix */
    if ((Columns > 1) && (Rows > 1))
      {
	for(int i = 0 ; i<Incr.N_rows ; i++)
	  {
	    for(int j = 0 ; j<Incr.N_cols ; j++)
	      {
		A.nM[i][j] += Incr.nM[i][j];
	      }
	  }
      }

    /* Is an array */
    else if (((Columns == 1) && (Rows > 1)) ||
	     (Columns > 1) && (Rows == 1))
      {
	for(int i = 0 ; i<Columns*Rows ; i++)
	  {
	  A.nV[i] += Incr.nV[i];
	  }	
      }

    /* Fail */
    else
      {	
	fprintf(stderr,"%s : %s \n",
		"Error in increment__MatrixLib__()","Wrong input matrix!");
	exit(EXIT_FAILURE);
      }

  }

  /* Fail */
  else
    {
	fprintf(stderr,"%s : %s \n",
		"Error in increment__MatrixLib__()","Wrong input matrix!");
	exit(EXIT_FAILURE);
    }

  return A;
}

/*********************************************************************/

Matrix addition__MatrixLib__(Matrix A,Matrix B)
/*
  Sum two matrix A and B, and return the result C 
  C = A + B   
*/
{
  /* Variable declaration */
  int Rows = A.N_rows;
  int Columns = A.N_cols;
  Matrix C = allocZ__MatrixLib__(Rows,Columns);
 
  /* Check if it is possible to do the addition */
  if((A.N_cols == B.N_cols) && (A.N_rows == B.N_rows))
    {
      Rows = A.N_rows;
      Columns = A.N_cols;
      C = allocZ__MatrixLib__(Rows,Columns);

      /* Two array addition */
      if (((Columns == 1) && (Rows > 1)) ||
	  (Columns > 1) && (Rows == 1))
	{    
	  for(int i = 0 ; i<C.N_rows*C.N_cols ; i++)
	    {  
	      C.nV[i] = A.nV[i] + B.nV[i];
	    }
	}

      /* Two matrix addition */
      else if ((Columns > 1) && (Rows > 1))
	{    
	  for(int i = 0 ; i < C.N_rows ; i++)
	    {
	      for(int j = 0 ; j < C.N_cols ; j++)
		{
		  C.nM[i][j] = A.nM[i][j] + B.nM[i][j];
		}
	    }
	}

      /* Fail */
      else{
	fprintf(stderr,"%s : %s \n",
		"Error in addition__MatrixLib__()","Wrong input matrix!");
	exit(EXIT_FAILURE);
      }
      
    }

  /* Fail */
  else{
    fprintf(stderr,"%s : %s \n",
	    "Error in addition__MatrixLib__()","Wrong input matrix!");
    exit(EXIT_FAILURE);
  }

  return C;
}

/*********************************************************************/

Matrix substraction__MatrixLib__(Matrix A,Matrix B)
/*

  Substract two matrix A and B, and return the result C 

  C = A - B
   
 */
{
  /* Variable declaration */
  int Rows = A.N_rows;
  int Columns = A.N_cols;
  Matrix C = allocZ__MatrixLib__(Rows,Columns);
 
  /* Check if it is possible to do the addition */
  if((A.N_cols == B.N_cols) && (A.N_rows == B.N_rows))
    {
      Rows = A.N_rows;
      Columns = A.N_cols;

      /* Two array addition */
      if (((Columns == 1) && (Rows > 1)) ||
	  (Columns > 1) && (Rows == 1))
	{    
	  for(int i = 0 ; i<C.N_rows*C.N_cols ; i++)
	    {  
	      C.nV[i] = A.nV[i] - B.nV[i];
	    }
	}

      /* Two matrix addition */
      else if ((Columns > 1) && (Rows > 1))
	{    
	  for(int i = 0 ; i < C.N_rows ; i++)
	    {
	      for(int j = 0 ; j < C.N_cols ; j++)
		{
		  C.nM[i][j] = A.nM[i][j] - B.nM[i][j];
		}
	    }
	}

      /* Fail */
      else{
	fprintf(stderr,"%s : %s \n",
		"Error in substraction__MatrixLib__()","Wrong input matrix!");
	exit(EXIT_FAILURE);
      }
      
    }

  /* Fail */
  else{
    fprintf(stderr,"%s : %s \n",
	    "Error in substraction__MatrixLib__()","Wrong input matrix!");
    exit(EXIT_FAILURE);
  }

  return C;
}

/*********************************************************************/

double norm__MatrixLib__(Matrix In,int kind)
/*
  Get the norm of a vector in a Euclidean space R2. Implemented norms :
  - Euclidean norm (kind == 2)
*/
{

  double Out;
  double aux = 0;
  int Columns = In.N_cols;
  int Rows = In.N_rows;
  
  /* Check */
  if (((Columns == 1) && (Rows > 1)) ||
      (Columns > 1) && (Rows == 1))
    {
      /* Euclidean norm */
      if(kind == 2)
	{
	  for(int i = 0 ; i< Rows*Columns ; i++)
	    {
	      aux += DSQR(In.nV[i]);
	    }
	  Out = pow(aux,0.5);
	}
      
    }
  else
    {
      puts("Error in norm__MatrixLib__() : The input data is not a vector ! ");
      exit(EXIT_FAILURE);     
    }
  return Out;  
}

/*********************************************************************/

double Euclidean_distance__MatrixLib__(Matrix la)
/*
  Compute the Euclidean distance
*/
{

  int Ndim = NumberDimensions;
  double sqr_distance = 0;
  double distance = 0;

  for(int i = 0 ; i<Ndim ; i++)
  {
    sqr_distance += la.nV[i]*la.nV[i]; 
  }

  distance = sqrt(sqr_distance);

  return distance;

}

/*********************************************************************/

double generalised_Euclidean_distance__MatrixLib__(Matrix la, Matrix Metric)
/*
  Compute the Euclidean distance using a general metric tensor
*/
{

  int Ndim = NumberDimensions;
  double Metric_x_la = 0;
  double sqr_distance = 0;
  double distance = 0;

  for(int i = 0 ; i<Ndim ; i++)
  {

    Metric_x_la = 0;

    for(int j = 0 ; j<Ndim ; j++)
    {
      Metric_x_la += Metric.nM[i][j]*la.nV[j];
    }

    sqr_distance += la.nV[i]*Metric_x_la;
    
  }

  distance = sqrt(sqr_distance);

  return distance;

}

/*********************************************************************/

Matrix Eigen_Mat(Matrix In)
{

  int Rows = In.N_rows;
  int Columns = In.N_cols;  
  Matrix Eigen = memory_to_matrix__MatrixLib__(Rows,Columns,NULL);
  Matrix Eigen_Vals;
  Matrix Coeffs;
  
  if((Rows == Columns) && (Columns > 1) && (Rows > 1))
    {
  
      switch(Columns)
	{

	case 1 :
	  exit(EXIT_FAILURE);
	  break;
	case 2 :
	  /* Get the coefficients of the charasteristic pol */
	  Coeffs = allocZ__MatrixLib__(1,3);
	  Coeffs.nV[0] = 1.0; /* a*x^2 */
	  Coeffs.nV[1] = - In.nM[0][0] - In.nM[1][1]; /* b*x */
	  Coeffs.nV[2] = In.nM[0][0]*In.nM[1][1] - In.nM[1][0]*In.nM[0][1]; /* c */
	  /* Solve the charasteristic pol to the eigenvalues */
	  Eigen_Vals = solve_polynomial__MatrixLib__(Coeffs);
	  free__MatrixLib__(Coeffs);
	  /* Assign the eigenvalues to the solution */
	  Eigen.nV = Eigen_Vals.nV;
	  break;

	case 3 :
	  /* Get the coefficients of the charasteristic pol */
	  Coeffs = allocZ__MatrixLib__(1,4);
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
	  Eigen_Vals = solve_polynomial__MatrixLib__(Coeffs);
	  free__MatrixLib__(Coeffs);
	  /* Assign the eigenvalues to the solution */
	  Eigen.nV = Eigen_Vals.nV;
	  break;
    
	default :
	  exit(EXIT_FAILURE);
	}
    }
  else
    {
      puts("Error in Eigen_Mat() : The input should be a square matrix !");    
      exit(EXIT_FAILURE);
    }

  return Eigen;
}

/*********************************************************************/

Matrix lumped__MatrixLib__(Matrix M_in)
{
  /*
    Get the lumped matrix in a vectorial shape of the input matrix 
  */

  if(M_in.N_cols != M_in.N_rows)
    {
      puts("Error in lumped__MatrixLib__() : The input matrix is not square ! ");
      exit(EXIT_FAILURE);
    }
  
  Matrix M_out;
  M_out = allocZ__MatrixLib__(M_in.N_rows,1);

  for(int i = 0 ; i<M_in.N_cols ; i++)
    {
      for(int j = 0 ; j<M_in.N_rows ; j++)
	{
	  M_out.nV[i] += M_in.nM[j][i];
	}
    }
  return M_out;  
}

/*********************************************************************/


double rcond__MatrixLib__(Matrix A)
/*
  C = rcond(A) returns an estimate for the reciprocal condition of A in 1-norm. 
  If A is well conditioned, rcond(A) is near 1.0. 
  If A is badly conditioned, rcond(A) is near 0.
*/
{

  double RCOND; 
  double ANORM;
  int INFO;
  int N_rows = A.N_rows;
  int N_cols = A.N_cols;
  int LDA = IMAX(N_rows,N_cols);
  double * AUX_MEMORY = (double *)calloc(N_rows*N_cols,sizeof(double));
  double * WORK_ANORM = (double *)Allocate_Array(IMAX(1,N_rows),sizeof(double));
  double * WORK_RCOND = (double *)Allocate_Array(4*N_rows,sizeof(double));
  int * IPIV = (int *)Allocate_Array(IMIN(N_rows,N_cols),sizeof(int));
  int * IWORK_RCOND = (int *)Allocate_Array(N_rows,sizeof(int));

  // Copy matrix because dgetrf_ is a destructive operation
  memcpy(AUX_MEMORY, A.nV, N_rows*N_cols*sizeof(double));

  // Compute 1-norm
  ANORM = dlange_("1", &N_rows, &N_cols, AUX_MEMORY, &LDA, WORK_ANORM);

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N_rows,&N_cols,AUX_MEMORY,&LDA,IPIV,&INFO);

  // Check output of dgetrf
  if(INFO != 0)
  {
    if(INFO < 0)
    {
      printf("%s : \n","Error in rcond__MatrixLib__()");
      printf("the %i-th argument had an illegal value",abs(INFO));
    }
    else if(INFO > 0)
    {
      printf("%s :\n","Error in rcond__MatrixLib__()");
      printf(" M(%i,%i) %s \n %s \n %s \n %s \n",INFO,INFO,"is exactly zero. The factorization",
        "has been completed, but the factor M is exactly",
        "singular, and division by zero will occur if it is used",
        "to solve a system of equations.");
    
    }    
    exit(EXIT_FAILURE);
  }

  // Compute the Reciprocal condition number
  dgecon_("1", &N_rows, AUX_MEMORY, &LDA, &ANORM, &RCOND, WORK_RCOND, IWORK_RCOND, &INFO);

 if(INFO<0)
 {
  printf("Error in rcond__MatrixLib__() : the %i-th argument of dgecon_ had an illegal value\n",abs(INFO));
  exit(EXIT_FAILURE);
 }

 // Free auxiliar memory
 free(IPIV);
 free(WORK_ANORM);
 free(WORK_RCOND);
 free(IWORK_RCOND);
 free(AUX_MEMORY);

 return RCOND;
}

/*********************************************************************/

Matrix solve__MatrixLib__(Matrix K, Matrix F)
/*

 */
{
  int Bool = K.N_cols>3;
  Matrix Km1;
  Matrix U;
  
  switch(Bool)
  {
    case 0 :
    Km1 = inverse__MatrixLib__(K);
    U = matrix_product__MatrixLib__(Km1,F);
    free__MatrixLib__(Km1);
    return U; 
    
    case 1 :
    puts("fix U");
    exit(EXIT_FAILURE);
    U = Conjugate_Gradient_Method(K,F,U);
    return U;
    
    default :
    exit(EXIT_FAILURE);
  }

}

/*********************************************************************/


float pythag(float a, float b)
/*
  Computes (a^2 + b^2)^0.5 without destructive underflow or overflow.
*/
{
  /* 
     Compute the absolute value of "a" and "b" 
  */
  float absa=fabs(a);
  float absb=fabs(b);
  
  if (absa > absb)
    {
      return absa*sqrt(1.0+SQR(absb/absa));
    }
  else
    {
      return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
    }
}

/*********************************************************************/

void print__MatrixLib__(Matrix In, int PrintRows, int PrintColumns)
/*
  Print term by term the input matrix
*/
{

  int Rows = IMIN(PrintRows,In.N_rows);
  int Columns = IMIN(PrintColumns,In.N_cols);
  
  if(((In.N_rows > 1) && (In.N_cols == 1)) ||
     ((In.N_rows == 1) && (In.N_cols > 1)))
    {
       
      for(int i = 0 ; i<Rows*Columns ; i++)
  {
    printf("%e \n",In.nV[i]);
  }
      
    }
  
  if((In.N_rows > 1) && (In.N_cols > 1))
    {
    
      for(int i = 0 ; i<Rows ; i++)
  {
  
    for(int j = 0 ; j<Columns ; j++)
      {
        if(In.nM[i][j] == 0)
        {
          printf("-------- ");
        }
        else
        {
          printf("%2.2e ",In.nM[i][j]);
        }
      }
  
    printf("\n\n");
  
  }
    }
}

/*********************************************************************/