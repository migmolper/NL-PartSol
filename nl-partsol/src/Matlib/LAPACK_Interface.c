#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>

#endif


Matrix solve_system_LAPACK(Matrix A,Matrix B)
{

  int NumColumns = A.N_cols;
  int NumRows = A.N_rows;

  int Order;
  int LDA;
  if(NumColumns == NumRows)
    {
      Order = NumRows;
      LDA = Order;
    }
  else
    {
      fprintf(stderr,"%s : %s !!\n",
	      "Error in solve_system_LAPACK",
	      "The input matrix should be square");
      exit(EXIT_FAILURE);
    }

  int NRHS;
  if(B.N_cols == 1)
    {
      NRHS = B.N_cols;
    }
  else
    {
      fprintf(stderr,"%s : %s !!\n",
	      "Error in solve_system_LAPACK",
	      "The RHS should have one column");
      exit(EXIT_FAILURE);
    }

  int LDB;
  if(B.N_rows == A.N_rows)
    {
      LDB = Order;
    }
  else
    {
      fprintf(stderr,"%s : %s !!\n",
	      "Error in solve_system_LAPACK",
	      "The dimension of the matrix and the RHS mismatch");
      exit(EXIT_FAILURE);
    }  
  
  char  TRANS = 'N'; /* (No transpose) */
  int   INFO= 3;

  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
 
  /* Compute the LU factorization */
  dgetrf_(&Order,&Order,A.nV,&LDA,IPIV,&INFO);
  
 
  /* 
     checks INFO, if INFO != 0 something goes wrong.
     For more information see the MAN page of dgetrf.
  */
  if(INFO)
    {
      fprintf(stderr,"%s : %s %s %s \n",
	      "Error in solve_system_LAPACK",
	      "The function",
	      "LAPACK_dgetrf",
	      "returned an error message !!!" );
      exit(EXIT_FAILURE);
    }
  else
    {
      /* Solve the sistem */
      dgetrs_(&TRANS,&Order,&NRHS,A.nV,&LDA,IPIV,B.nV,&LDB,&INFO);

      /* 
	 checks INFO, if INFO != 0 something goes wrong.
	 For more information see the MAN page of dgetrf.
      */
      if(INFO)
        {
	  fprintf(stderr,"%s : %s %s %s \n",
		  "Error in solve_system_LAPACK",
		  "The function",
		  "dgetrs_",
		  "returned an error message !!!" );
	  exit(EXIT_FAILURE);
        }
    }
 
  return B;
}
