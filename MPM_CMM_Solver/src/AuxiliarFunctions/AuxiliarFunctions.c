#include <stdio.h>
#include <stdlib.h>

/***************************************************************************/


/* Allocate function for array */
void * AllocateArray(int dim,int trama)
/*
  Inputs
  - dim : Dimension of the array
  - trama : Number of bytes of one single element in the array

  Outputs
  - Pointer to the to the reserved zone in memory
*/
{
  /* Declaration of the pointer */
  void * v;

  /* Save data in memory */
  v = malloc(dim*trama);

  /* Error message */
  if(v == NULL) { puts("Error in allocate_array"); exit(0); }
  
  /* Return the pointer */
  return v;  
}


/***************************************************************************/


/* Allocate function for matrix */
void ** AllocateMatrix(int n_row, int n_col, int trama)
  /*
  Inputs
  - n_row : Number of rows
  - n_col : Number of columns
  - trama : Number of bytes of one single element in the matrix

  Outputs
  - Pointer to the to the reserved zone in memory
  */
{
  void ** m;

  /* Allocate of the pointer to the pointer table (rows) */
  m = (void **)malloc((unsigned) n_row * sizeof(void *));
  if (m == NULL) { puts("Error in allocate_matrix"); exit(0); }

  /* Allocate the pointer in each row (columns) */
  for(int i = 0; i<n_row; i++){

    m[i] = malloc((unsigned) n_col * trama);
    if( m[i] == NULL ) { puts("Error in allocate_matrix"); exit(0);  }
    
  }

  /* Return the matrix */
  return m;
}


/***************************************************************************/


void ** AllocateTableOfPointers(int NumPointers,int trama){

  /* Table of pointers */
  void ** table;
  
  /* Allocate a look-up-table of pointers. */ 
  table = (void **)malloc(trama * NumPointers);

  /* Error message */
  if(v == NULL) { puts("Error in AllocateTableOfPointers"); exit(0); }
   
  /* Return the table of pointers */
  return table;
}


/***************************************************************************/

double CalculeDeterminant(double * InputMatrix,
			  int ShapeMatrix)
/* Calcule the Determinant of a Matrix. 
   Inputs:
   - Shape of the matrix
   - Matrix

   3x3 No symetric -> Shape : 9
   [ 0 ; 1 ; 2 ]
   [ 3 ; 4 ; 5 ]
   [ 6 ; 7 ; 8 ]

   3x3 Symetric -> Shape : 6
   [ 0 ; 1 ; 2 ]
   [ * ; 3 ; 4 ]
   [ * ; * ; 5 ]

   2x2 No symetric -> Shape : 4
   [ 0 ; 1 ]
   [ 2 ; 3 ]

   2x2 Symetric -> Shape : 3
   [ 0 ; 1 ]
   [ * ; 2 ]

   1x1 -> Shape : 1
   [ 0 ]

   Outputs :
   - Determinant of the Matrix
   
*/
{

  /* Final value */
  double DeterminantMatrix;


  switch(ShapeMatrix){
  case 1 :
    DeterminantMatrix = InputMatrix[0];
    break;
    
  case 3 :
    DeterminantMatrix =
      InputMatrix[0]*InputMatrix[2] -
      InputMatrix[1]*InputMatrix[1];
    break;
    
  case 4 :
    DeterminantMatrix =
      InputMatrix[0]*InputMatrix[3] -
      InputMatrix[1]*InputMatrix[2];
    break;
    
  case 6 :
    DeterminantMatrix =
      InputMatrix[0]*InputMatrix[3]*InputMatrix[5] +
      InputMatrix[1]*InputMatrix[2]*InputMatrix[4] +
      InputMatrix[1]*InputMatrix[2]*InputMatrix[4] -
      InputMatrix[2]*InputMatrix[2]*InputMatrix[3] -
      InputMatrix[1]*InputMatrix[1]*InputMatrix[5] -
      InputMatrix[0]*InputMatrix[4]*InputMatrix[4];
    break;
    
  case 9 :
    DeterminantMatrix =
      InputMatrix[0]*InputMatrix[4]*InputMatrix[8] +
      InputMatrix[1]*InputMatrix[5]*InputMatrix[6] +
      InputMatrix[2]*InputMatrix[3]*InputMatrix[7] -
      InputMatrix[2]*InputMatrix[4]*InputMatrix[6] -
      InputMatrix[1]*InputMatrix[3]*InputMatrix[8] -
      InputMatrix[0]*InputMatrix[5]*InputMatrix[7];
    break;
    
  default :
    printf(" Error in the shape of the matrix \n");
    break(0);
  }

  return DeterminantMatrix;
  
}


/***************************************************************************/

void * AllocateTableOfPointers(int dim,int trama){
}
