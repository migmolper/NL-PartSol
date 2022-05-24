
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
  char TRANS = 'T'; /* (Transpose) */
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
    fprintf(stderr, "%s : %s %s %s \n", "Error in dgetrf_", "The function",
            "dgetrf_", "returned an error message !!!");
    return EXIT_FAILURE;
  }

  /*
    Solve the system
  */
  INFO = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'T',Order,NRHS, Tangent_Stiffness, LDA,IPIV,Residual,LDB);
  if (INFO) {
    free(IPIV);
    fprintf(stderr, "%s : %s %s %s \n", "Error in dgetrs_", "The function",
            "dgetrs_", "returned an error message !!!");
    return EXIT_FAILURE;
  }

  // Free memory
  free(IPIV);

  return STATUS;
}
