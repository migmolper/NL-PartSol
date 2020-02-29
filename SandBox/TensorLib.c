#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* /\*********************************************************************\/ */

/* Matrix Incr_Mat(Matrix A, Matrix Incr) */
/* /\* */
/*   Increment term by term the Matrix A with Incr */
/* *\/ */
/* { */

/*   int Bool; */

/*   Bool = Incr.N_cols*Incr.N_rows == 1; */

/*   switch(Bool){ */
    
/*   case 0: /\* Incr is a matrix or an array *\/ */

/*     if( (A.N_cols != Incr.N_cols) || */
/* 	(A.N_rows != Incr.N_rows) ){ /\* Check the size of the input *\/ */
/*       printf("Error in Incr_Mat() : A is incompatible with the Incr !! \n"); */
/*        exit(EXIT_FAILURE); */
/*     } */

/*     if(Incr.N_cols == Incr.N_rows){ /\* Matrix *\/ */
/*       for(int i = 0 ; i<Incr.N_rows ; i++){ */
/* 	for(int j = 0 ; j<Incr.N_cols ; j++){ */
/* 	  A.nM[i][j] += Incr.nM[i][j]; */
/* 	} */
/*       }	 */
/*     } */
    
/*     if( (Incr.N_cols == 1) || */
/* 	(Incr.N_rows == 1) ){ /\* Array *\/ */
/*       for(int i = 0 ; i<Incr.N_cols*Incr.N_rows ; i++){ */
/* 	A.nV[i] += Incr.nV[i]; */
/*       } */
/*     } */
    
/*     break; */
    
/*   case 1: /\* Incr is a scalar *\/ */
    
/*     if( (A.N_cols == 1) || */
/* 	(A.N_rows == 1) ){ /\* A is an array *\/ */
/*       for(int i = 0 ; i<A.N_cols*A.N_rows ; i++){ */
/* 	A.nV[i] += Incr.n; */
/*       } */
/*     } */
    
/*     if( (A.N_cols != 1) && */
/* 	(A.N_rows != 1) ){ */
/*       for(int i = 0 ; i<A.N_rows ; i++){ */
/* 	for(int j = 0 ; j<A.N_cols ; j++){ */
/* 	  A.nM[i][j] += Incr.n; */
/* 	} */
/*       }	 */
/*     } */
    
/*     break; */

    
/*   default : */
/*      exit(EXIT_FAILURE); */
/*   } */

/*   return A; */
/* }   */


/* /\*********************************************************************\/ */

/* Matrix Add_Mat(Matrix A,Matrix B) */
/* /\* */

/*   Sum two matrix A and B, and return the result C  */

/*   C = A + B */
   
/*  *\/ */
/* { */

/*   /\* Variable declaration *\/ */
/*   Matrix C; */
  
/*   /\* Check if it is possible to do the addition *\/ */
/*   if((A.N_cols != B.N_cols)||(A.N_rows != B.N_rows)){ */
/*     puts("Error in Add_Mat() : Your are trying to add incompatible matrix"); */
/*      exit(EXIT_FAILURE); */
/*   } */

/*   /\* Allocate output matrix *\/ */
/*   C = MatAlloc(A.N_rows,A.N_cols); */

/*   /\* Difference if we are doing an addition of two matrix or two vectors *\/ */
/*   if( (A.nV != NULL) && (B.nV != NULL) ){ /\* two array addition *\/ */
    
/*     for(int i = 0 ; i<C.N_rows*C.N_cols ; i++){   */
/* 	C.nV[i] = A.nV[i] + B.nV[i]; */
/*     } */
    
/*   } */
/*   else if( (A.nM != NULL) && (B.nM != NULL) ){ /\* two matrix addition *\/ */
    
/*     for(int i = 0 ; i < C.N_rows ; i++){ */
/*       for(int j = 0 ; j < C.N_cols ; j++){ */
/* 	C.nM[i][j] = A.nM[i][j] + B.nM[i][j]; */
/*       } */
/*     } */
  
/*   } */
/*   else{ */
/*     puts("Error in Add_Mat() : Inputs must be of the same range !"); */
/*      exit(EXIT_FAILURE); */
/*   }   */


/*   return C; */
/* } */

/* /\*********************************************************************\/ */

/* Matrix Sub_Mat(Matrix A,Matrix B) */
/* /\* */

/*   Substract two matrix A and B, and return the result C  */

/*   C = A - B */
   
/*  *\/ */
/* { */

/*   /\* Variable declaration *\/ */
/*   Matrix C; */
  
/*   /\* Check if it is possible to do the substraction *\/ */
/*   if((A.N_cols != B.N_cols)||(A.N_rows != B.N_rows)){ */
/*     puts("Error in Sub_Mat() : Your are trying to add incompatible matrix"); */
/*      exit(EXIT_FAILURE); */
/*   } */

/*   /\* Allocate output matrix *\/ */
/*   C = MatAlloc(A.N_rows,A.N_cols); */

/*   /\* Difference if we are doing an substraction of two matrix */
/*      or two vectors *\/ */
/*   if( (A.nV != NULL) && (B.nV != NULL) ){ /\* two array substraction *\/ */
    
/*     for(int i = 0 ; i<C.N_rows*C.N_cols ; i++){   */
/* 	C.nV[i] = A.nV[i] - B.nV[i]; */
/*     } */
    
/*   } */
/*   else if( (A.nM != NULL) && (B.nM != NULL) ){ /\* two matrix substraction *\/ */
    
/*     for(int i = 0 ; i < C.N_rows ; i++){ */
/*       for(int j = 0 ; j < C.N_cols ; j++){ */
/* 	C.nM[i][j] = A.nM[i][j] - B.nM[i][j]; */
/*       } */
/*     } */
    
/*   } */
/*   else{ */
/*     puts("Error in Sub_Mat() : Inputs must be of the same range !"); */
/*      exit(EXIT_FAILURE); */
/*   }   */


/*   return C; */
/* } */


/* /\*********************************************************************\/ */

/* Matrix Eigen_Mat(Matrix In){ */

/*   /\* Check if we dont have a null matrix *\/ */
/*   if ( (In.nM == NULL) && (In.n != In.n) ){ */
/*     if(In.nV != NULL){ */
/*       puts("Error in Eigen_Mat() : An array does not have determinant !"); */
/*        exit(EXIT_FAILURE); */
/*     } */
/*     puts("Error in Eigen_Mat() : Input matrix = NULL !");     */
/*      exit(EXIT_FAILURE); */
/*   } */
  
/*   /\* Check if the matrix is square *\/ */
/*   if(In.N_cols != In.N_rows){ */
/*     puts(" Error in Eigen_Mat() : Non square matrix !"); */
/*      exit(EXIT_FAILURE); */
/*   } */

/*   int Dim = In.N_cols; */
/*   Matrix Eigen = MatAssign(Dim,Dim,NAN,NULL,NULL); */
/*   Matrix Eigen_Vals; */
/*   Matrix Coeffs; */

/*   switch(Dim){ */

/*   case 1 : */
/*      exit(EXIT_FAILURE); */
/*     break; */
/*   case 2 : */
/*     /\* Get the coefficients of the charasteristic pol *\/ */
/*     Coeffs = MatAllocZ(1,3); */
/*     Coeffs.nV[0] = 1.0; /\* a*x^2 *\/ */
/*     Coeffs.nV[1] = - In.nM[0][0] - In.nM[1][1]; /\* b*x *\/ */
/*     Coeffs.nV[2] = In.nM[0][0]*In.nM[1][1] - In.nM[1][0]*In.nM[0][1]; /\* c *\/ */
/*     /\* Solve the charasteristic pol to the eigenvalues *\/ */
/*     Eigen_Vals = SolvePolynomial(Coeffs); */
/*     FreeMat(Coeffs); */
/*     /\* Assign the eigenvalues to the solution *\/ */
/*     Eigen.nV = Eigen_Vals.nV; */
/*     break; */

/*   case 3 : */
/*      exit(EXIT_FAILURE); */
/*     break; */
    
/*   default : */
/*      exit(EXIT_FAILURE); */
/*   } */

/*   return Eigen; */
/* } */


/* Matrix Matrix_x_Scalar(Matrix A, double a) */
/* /\*! */
/*  * \brief Brief description of Matrix_x_Scalar. */
/*  *        Function to multiply a Matrix with a scalar.  */
/*  * */
/*  *  The parameters for this functions are  : */
/*  * @param A : Input Matrix */
/*  * @param a : Input scalar  */
/* *\/ */
/* { */

/*   bool Is_Matrix = false; */
/*   bool Is_Vector = false; */
/*   int N_rows = A.N_rows; */
/*   int N_cols = A.N_cols; */

/*   /\* Check if its matrix or array *\/ */
/*   if((N_rows > 1) && (N_cols > 1)){ */
/*     Is_Matrix = true; */
/*   } */
/*   else if((N_rows == 1) || (N_cols == 1) ){ */
/*     Is_Vector = true; */
/*   } */

/*   if(Is_Matrix){ /\* Multiply matrix by an scalar *\/ */
/*     for(int i = 0 ; i<N_rows ; i++){ */
/*       for(int j = 0 ; j<N_cols ; j++){ */
/* 	A.nM[i][j] *= a; */
/*       } */
/*     } */
/*   } */
/*   else if(Is_Vector){ /\* Multiply vector by an scalar *\/ */
/*     for(int i = 0 ; i<N_cols*N_rows ; i++){ */
/*       A.nV[i] *= a; */
/*     } */
/*   } */
/*   else{ */
/*         fprintf(stderr,"%s : %s \n", */
/* 	    "Error in Matrix_x_Scalar(*,)", */
/* 	    "Not a matrix or a vector"); */
/*      exit(EXIT_FAILURE);     */
/*   } */

/*   return A; */
/* } */


/* /\*********************************************************************\/ */


/*************************************************************/

/* Matrix definition */
typedef struct{
  int Order; /* Order of the tensor */
  double * n; /* First order tensor */
  double (*N) [3]; /* Second order tensor */
  char Info [100]; /* Aditional information */
} Tensor;

/*************************************************************/

Tensor allocTensor(int Order){
  Tensor A;
  A.Order = Order;
    switch(Order){
    case 0:
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in allocTensor()",
	      "Not posible to allocate a scalar");
      exit(EXIT_FAILURE);       
    case 1:
      A.n = malloc( sizeof(double[3]) );  
      break;
    case 2:
      A.N = malloc( sizeof(double[Order+1][3]) );  
      break;
    default :
      fprintf(stderr,"%s : %s \n",
	      "Error in allocTensor()",
	      "Invalird order of the tensor");
      exit(EXIT_FAILURE); 
    }
  return A;
}

/*************************************************************/

void freeTensor(Tensor A){
  switch(A.Order){
    case 0:
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in freeTensor()",
	      "Not posible to free a scalar");
      exit(EXIT_FAILURE);       
    case 1:
      free(A.n);
      break;
    case 2:
      free(A.N);
      break;
    default :
      fprintf(stderr,"%s : %s \n",
	      "Error in freeTensor()",
	      "Invalird order of the tensor");
      exit(EXIT_FAILURE); 
    }
}

/*************************************************************/

double getDeterminantOf(Tensor A)
{
  /* Define output */
  double Determinant;
  /* Check if is the order is order 2 */
  if(A.Order == 2){  
    Determinant =
      A.N[0][0]*A.N[1][1]*A.N[2][2] - /* + a11*a22*a33 */
      A.N[0][0]*A.N[1][2]*A.N[2][1] + /* - a11*a23*a32 */
      A.N[0][1]*A.N[1][2]*A.N[2][0] - /* + a12*a23*a31 */
      A.N[0][1]*A.N[1][0]*A.N[2][2] + /* - a12*a33*a21 */
      A.N[0][2]*A.N[1][0]*A.N[2][1] - /* + a13*a21*a32 */
      A.N[0][2]*A.N[1][1]*A.N[2][0] ; /* - a13*a22*a31 */
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in getDeterminantOf()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }
  return Determinant;
}

/*************************************************************/

double getEuclideanNormOf(Tensor A)
{
  /* Define output */
  double Out;
  /* Check if the input is a first order tensor */
  if(A.Order == 1){
    /* Compute norm */    
    double Aux = 0;
    for(int i = 0 ; i<3 ; i++){
      Aux += A.n[i]*A.n[i] ;
    }
    Out = pow(Aux,0.5);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in getEuclideanNormOf()",
	    "The input should be of order 1");
    exit(EXIT_FAILURE);    
  }
  return Out;
}

/*************************************************************/

Tensor getI(){
  Tensor I = allocTensor(2);
  for(int i = 0 ; i<3 ; i++){
    I.N[i][i] = 1;
  }
  return I;
}

/*************************************************************/


Tensor getInverseOf(Tensor A)
{
  /* Allocate the output */
  Tensor Am1 = allocTensor(2);
  /* Check if the input is a second order tensor */
  if (A.Order == 2){  
    /* Get the determinant of the matrix */
    double detA = getDeterminantOf(A);
    if(detA == 0){
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in getInverseOf()",
	      "Determinant null");
      exit(EXIT_FAILURE);    
    }
    /* Compute each component  */
    Am1.N[0][0] =
      (double)1/detA*(A.N[1][1]*A.N[2][2] -
		      A.N[1][2]*A.N[2][1]);
    Am1.N[0][1] =
      -(double)1/detA*(A.N[0][1]*A.N[2][2] -
		       A.N[0][2]*A.N[2][1]);
    Am1.N[0][2] =
      (double)1/detA*(A.N[0][1]*A.N[1][2] -
		      A.N[0][2]*A.N[1][1]);
    Am1.N[1][0] =
      -(double)1/detA*(A.N[1][0]*A.N[2][2] -
		       A.N[1][2]*A.N[2][0]);
    Am1.N[1][1] =
      (double)1/detA*(A.N[0][0]*A.N[2][2] -
		      A.N[0][2]*A.N[2][0]);
    Am1.N[1][2] =
      -(double)1/detA*(A.N[0][0]*A.N[1][2] -
		       A.N[0][2]*A.N[1][0]);
    Am1.N[2][0] =
      (double)1/detA*(A.N[1][0]*A.N[2][1] -
		      A.N[1][1]*A.N[2][0]);
    Am1.N[2][1] =
      -(double)1/detA*(A.N[0][0]*A.N[2][1] -
		       A.N[0][1]*A.N[2][0]);
    Am1.N[2][2] =
      (double)1/detA*(A.N[0][0]*A.N[1][1] -
		      A.N[0][1]*A.N[1][0]);
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in getInverseOf()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }  
  /* Return the inverse matrix */
  return Am1;  
}

/*************************************************************/

Tensor getTransposeOf(Tensor A)
{
  /* Allocate the output */
  Tensor AT = allocTensor(2);  
  /* Check if the input is a second order tensor */
  if (A.Order == 2){
    /* Get the transpose */
    for(int i = 0 ; i < 3 ; i++){
      for(int j = 0 ; j < 3 ; j++){
	AT.N[i][j] = A.N[j][i];
      }
    }
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in getTransposeOf()",
	    "The input should be of order 2");
    exit(EXIT_FAILURE);    
  }  
  /* Return the transpose */
  return AT;
}

/*************************************************************/

double innerProductOf(Tensor A, Tensor B)
{

  /* Variable declaration output matrix */
  double AdotB = 0;

  /* Inner product of second order tensors */
  if ( (A.Order == 2) && (B.Order == 2) ) { 
    for(int i = 0 ; i < 3 ; i++){
      for(int j = 0 ; j < 3 ; j++){
	  AdotB += A.N[i][j]*B.N[i][j];
      }
    }
  }
  /* Inner product of first order tensors */
  else if( (A.Order == 1) && (B.Order == 1) ){
    for(int i = 0 ; i < 3 ; i++){
      AdotB += A.n[i]*B.n[i];
    }
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in innerProductOf()",
	    "The input should be two tensors or equal order");
    exit(EXIT_FAILURE);        
  }

  return AdotB;
    
}

/*************************************************************/

Tensor vectorProductOf(Tensor a, Tensor b){
  /* Allocate output */
  Tensor axb = allocTensor(1);
  /* Check if the input are a first order tensor */
  if ((a.Order == 1) && (b.Order == 1)){
    /* Operate vetor product */
    axb.n[0] = a.n[1]*b.n[2] - a.n[2]*b.n[1];
    axb.n[1] = a.n[2]*b.n[0] - a.n[0]*b.n[2];
    axb.n[2] = a.n[0]*b.n[1] - a.n[1]*b.n[0];
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in vectorProductOf()",
	    "The input should be two tensors of first order");
    exit(EXIT_FAILURE);    
  }
  /* Return vector product */
  return axb; 
}

/*************************************************************/

Tensor dyadicProductOf(Tensor a, Tensor b)
{
  /* Tensor declaration */
  Tensor aob = allocTensor(2);
  /* Check if the input are a first order tensor */
  if ((a.Order == 1) && (b.Order == 1)){
    /* Operate tensor product */
    for(int i = 0 ; i<3 ; i++){
      for(int j = 0 ; j<3 ; j++){
	aob.N[i][j] = a.n[i]*b.n[j];
      } 
    }
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in dyadicProductOf()",
	    "The input should be two tensors of first order");
    exit(EXIT_FAILURE);    
  }
  /* Return tensorial product */
  return aob;
}

/*************************************************************/

Tensor firstOrderContractionOf(Tensor A, Tensor b)
{
  /* Tensor declaration */
  Tensor Adotb = allocTensor(1);
  /* Check in the input its is ok */
  if ((A.Order == 2) && (b.Order == 1)){
    /* Auxiliar variable */
    double Aux;
    /* Operate tensors */
    for(int i = 0 ; i<3 ; i++){
      Aux = 0;
      for(int j = 0 ; j<3 ; j++){
	Aux += A.N[i][j]*b.n[j];
      }
      Adotb.n[i] = Aux;
    }    
  }
  else{
    fprintf(stderr,"%s : %s !!! \n",
	    "Error in firstOrderContractionOf()",
	    "The input should be 2ord tensor and a 1rd tensor");
    exit(EXIT_FAILURE);
  }
  /* Return tensor */
  return Adotb;
}

/*************************************************************/

void main(){
  int X = 10;
  int Y = 10;
  Tensor A = allocTensor(2);
  Tensor v = allocTensor(1);
  Tensor I = getI();
  Tensor Im1 = getInverseOf(I);
  printf("%f %f \n",getDeterminantOf(A),getDeterminantOf(I));
  for(int i = 0 ; i<3 ; i++){
    for(int j = 0 ; j<3 ; j++){
      printf("%f \n",Im1.N[i][j]);
    }
  }
  freeTensor(A);
  freeTensor(v);
  freeTensor(I);
  freeTensor(Im1);  
}
