#include "nl-partsol.h"
#include "lapacke.h"

Matrix solve_system_LAPACK(Matrix A,Matrix B)
{
  // note, to understand this part take a look in the MAN pages, at section of parameters.
  char    TRANS = 'N';
  int     INFO= 3;
  int     LDA = 3;
  int     LDB = 3;
  int     N = 3;
  int     NRHS = 1;
  int     IPIV[3] ;
 
  /* Compute the LU factorization */
  LAPACK_dgetrf(&N,&N,A.nV,&LDA,IPIV,&INFO);
 
  /* 
     checks INFO, if INFO != 0 something goes wrong.
     For more information see the MAN page of dgetrf.
  */
  if(INFO)
    {
      puts("an error occured : ");
    }
  else
    {
      /* Solve the sistem */
      dgetrs_(&TRANS,&N,&NRHS,A.nV,&LDA,IPIV,A.nV,&LDB,&INFO);

      /* 
	 checks INFO, if INFO != 0 something goes wrong.
	 For more information see the MAN page of dgetrf.
      */
      if(INFO)
        {
	  puts("an error occured : ");
        }
    }
 
  return B;
}
