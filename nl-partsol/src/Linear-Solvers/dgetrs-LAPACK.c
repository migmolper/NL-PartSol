
#include "Linear-Solvers/dgetrs-LAPACK.h"

int dgetrs_LAPACK(
    double * Tangent_Stiffness,
    double * Residual,
    unsigned Nactivedofs)
{
  int STATUS = EXIT_SUCCESS;
  unsigned Order = Nactivedofs;
  unsigned LDA = Nactivedofs;
  unsigned LDB = Nactivedofs;
  char TRANS = 'T'; /* (Transpose) */
  int INFO = 3;
  int NRHS = 1;
  int *IPIV = (int *)calloc(Order, __SIZEOF_INT__);
  if(IPIV == NULL){
        fprintf(stderr, ""RED"Error in calloc(): Out of memory"RESET" \n");
        return EXIT_FAILURE;
  }   
  
  //  Compute the LU factorization
  dgetrf_(&Order, &Order, Tangent_Stiffness, &LDA, IPIV, &INFO);
  if (INFO) {
    free(IPIV);
    fprintf(stderr, "%s : %s %s %s \n", "Error in dgetrf_", "The function",
            "dgetrf_", "returned an error message !!!");
    return EXIT_FAILURE;
  }

  /*
    Solve the system
  */
  dgetrs_(&TRANS, &Order, &NRHS, Tangent_Stiffness, &LDA, IPIV, Residual, &LDB, &INFO);
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
