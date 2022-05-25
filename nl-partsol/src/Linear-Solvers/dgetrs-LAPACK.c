
#include <stdlib.h>
#include <stdio.h>
#include "Macros.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

#include "Linear-Solvers/dgetrs-LAPACK.h"

int dgetrs_LAPACK(
    double * Tangent_Stiffness,
    double * Residual,
    unsigned Nactivedofs)
{
  int STATUS = EXIT_SUCCESS;
  int Order = Nactivedofs;
  int LDA = Nactivedofs;
  int LDB = Nactivedofs;
  char TRANS = 'N'; /* (Transpose) */
  int INFO = 3;
  int NRHS = 1;
  int *IPIV = (int *)calloc(Order, __SIZEOF_INT__);
  if(IPIV == NULL){
        fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
        return EXIT_FAILURE;
  }   
  
  //  Compute the LU factorization
  INFO = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,Order,Order,Tangent_Stiffness,LDA,IPIV);
  if (INFO) {
    free(IPIV);
    fprintf(stderr, "%s : %s %s %s \n", "Error in LAPACKE_dgetrf", "The function",
            "LAPACKE_dgetrf", "returned an error message !!!");
    return EXIT_FAILURE;
  }

  /*
    Solve the system
  */
  INFO = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'N',Order, NRHS, Tangent_Stiffness, LDA,IPIV,Residual,LDB);
  if (INFO != 0) {
    free(IPIV);
    fprintf(stderr, "%s : %s %s %s \n", "Error in LAPACKE_dgetrs", "The function",
            "LAPACKE_dgetrs", "returned an error message !!!");
    return EXIT_FAILURE;
  }

  // Free memory
  free(IPIV);

  return STATUS;
}
