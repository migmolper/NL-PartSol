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


