#include <stdio.h>
#include <stdlib.h>

/* Function for arrays declaration */
void * Allocate_Array(int SizeArray, int SizeType){
  
  void * V;
  
  V = (void *)malloc(SizeArray*SizeType);
  
  if (V == NULL){puts("Error in the array declaration"); exit(0);}
  
  return V;
  
}


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



double Get_Determinant(double ** Matrix_in,int Dimension)
/* 
   Function to get the determinant of a matrix (max 3x3) 
   Inputs : Matrix, and order of the matrix
   Outputs : Determinant
*/
{

  double Det_out;
  
  switch(Dimension){

  case 1:
    Det_out = *Matrix_in[0];
    break;

  case 2:
    Det_out =
      Matrix_in[0][0]*Matrix_in[1][1] -
      Matrix_in[0][1]*Matrix_in[1][0];
    break;
    
  case 3:
    Det_out =
      Matrix_in[0][0]*Matrix_in[1][1]*Matrix_in[2][2] - /* + a11*a22*a33 */
      Matrix_in[0][0]*Matrix_in[1][2]*Matrix_in[2][1] + /* - a11*a23*a32 */
      Matrix_in[0][1]*Matrix_in[1][2]*Matrix_in[2][0] - /* + a12*a23*a31 */
      Matrix_in[0][1]*Matrix_in[1][0]*Matrix_in[2][2] + /* - a12*a33*a21 */
      Matrix_in[0][2]*Matrix_in[1][0]*Matrix_in[2][1] - /* + a13*a21*a32 */
      Matrix_in[0][2]*Matrix_in[1][1]*Matrix_in[2][0] ; /* - a13*a22*a31 */
    break;
  }
  
  return Det_out;
}

void Get_Inverse(double ** Matrix_in,double ** Matrix_out,int Dimension)
/* 
   Get the inverse of the matrix :
   Inputs : Matrix in, dimension os the matrix
   Outpus : Matrix out
*/
{
  
  /* Get the determinant of the matrix */
  double Det = Get_Determinant(Matrix_in,Dimension);
  if(Det <= 0){
    printf("Determinant null o less than zero \n");
    exit(0);
  }

  /* Do in a different fashion if it is 2D or 3D */
  if(Dimension == 2){
    Matrix_out[0][0] = 1/(Det)*Matrix_in[1][1];

    Matrix_out[0][1] = -1/(Det)*Matrix_in[0][1];

    Matrix_out[1][0] = -1/(Det)*Matrix_in[1][0];

    Matrix_out[1][1] = 1/(Det)*Matrix_in[0][0];
  }
  else if(Dimension == 3){
    Matrix_out[0][0] = 1/(Det)*(Matrix_in[1][1]*Matrix_in[2][2] -
				Matrix_in[1][2]*Matrix_in[2][1]);
    
    Matrix_out[0][1] = -1/(Det)*(Matrix_in[0][1]*Matrix_in[2][2] -
				 Matrix_in[0][2]*Matrix_in[2][1]);
    
    Matrix_out[0][2] = 1/(Det)*(Matrix_in[0][1]*Matrix_in[1][2] -
				Matrix_in[0][2]*Matrix_in[1][1]);
    
    Matrix_out[1][0] = -1/(Det)*(Matrix_in[1][0]*Matrix_in[2][2] -
				 Matrix_in[1][2]*Matrix_in[2][0]);
    
    Matrix_out[1][1] = 1/(Det)*(Matrix_in[0][0]*Matrix_in[2][2] -
				Matrix_in[0][2]*Matrix_in[2][0]);
    
    Matrix_out[1][2] = -1/(Det)*(Matrix_in[0][0]*Matrix_in[1][2] -
				 Matrix_in[0][2]*Matrix_in[1][0]);
    
    Matrix_out[2][0] = 1/(Det)*(Matrix_in[1][0]*Matrix_in[2][1] -
				Matrix_in[1][1]*Matrix_in[2][0]);
    
    Matrix_out[2][1] = -1/(Det)*(Matrix_in[0][0]*Matrix_in[2][1] -
				 Matrix_in[0][1]*Matrix_in[2][0]);
    
    Matrix_out[2][2] = 1/(Det)*(Matrix_in[0][0]*Matrix_in[1][1] -
				Matrix_in[0][1]*Matrix_in[1][0]);
  }

}
