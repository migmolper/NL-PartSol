// clang-format off
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Matlib.h"
// clang-format on

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/*************************************************************/

Tensor alloc__TensorLib__(int Order) {
  int Ndim = NumberDimensions;
  /* Define output */
  Tensor A;
  /* Swith cases */
  switch (Order) {
  case 0:
    fprintf(stderr, "%s : %s !!! \n", "Error in alloc__TensorLib__()",
            "Not posible to allocate a scalar");
    exit(EXIT_FAILURE);
  case 1:
    A.Order = 1;
    A.n = (double *)malloc(sizeof(double[Ndim]));
    for (int i = 0; i < Ndim; i++) {
      A.n[i] = 0.0;
    }
    break;
  case 2:
    A.Order = 2;
    for (int i = 0; i < Ndim; i++) {
      A.N[i] = (double *)malloc(sizeof(double[Ndim]));
      for (int j = 0; j < Ndim; j++) {
        A.N[i][j] = 0.0;
      }
    }
    break;
  default:
    fprintf(stderr, "%s : %s \n", "Error in alloc__TensorLib__()",
            "Invalird order of the tensor");
    exit(EXIT_FAILURE);
  }
  return A;
}

/*************************************************************/

Tensor memory_to_tensor__TensorLib__(double *A_mem, int Order) {
  long Ndim = NumberDimensions;
  /* Define output */
  Tensor A_tens;
  /* Swith cases */
  switch (Order) {
  case 0:
    fprintf(stderr, "%s : %s !!! \n",
            "Error in memory_to_tensor__TensorLib__()",
            "Tensor type not include scalar properties");
    exit(EXIT_FAILURE);
  case 1:
    A_tens.Order = 1;
    A_tens.n = A_mem;
    break;
  case 2:
    A_tens.Order = 2;
    for (long i = 0; i < Ndim; i++) {
      A_tens.N[i] = A_mem + i * Ndim;
    }
    break;
  default:
    fprintf(stderr, "%s : %s \n", "Error in memory_to_tensor__TensorLib__()",
            "Invalird order of the tensor");
    exit(EXIT_FAILURE);
  }

  return A_tens;
}

/*************************************************************/

void free__TensorLib__(Tensor A) {
  int Ndim = NumberDimensions;

  switch (A.Order) {
  case 0:
    fprintf(stderr, "%s : %s !!! \n", "Error in free__TensorLib__()",
            "Not posible to free a scalar");
    exit(EXIT_FAILURE);
  case 1:
    free(A.n);
    break;
  case 2:
    for (int i = 0; i < Ndim; i++) {
      free(A.N[i]);
    }
    break;
  default:
    fprintf(stderr, "%s : %s \n", "Error in free__TensorLib__()",
            "Invalird order of the tensor");
    exit(EXIT_FAILURE);
  }
}

/*************************************************************/

double I1__TensorLib__(const double * A) {

  double I1 = 0;

#if NumberDimensions == 2
    I1 = A[0] + A[3];
#else 
    I3 = A[0] + A[4] + A[8];
#endif

  return I1;
}

/*************************************************************/

double I2__TensorLib__(const Tensor A) {
  int Ndim = NumberDimensions;
  /* Define output */
  double I2 = 0;
  /* Check if is the order is order 2 */
  if (A.Order == 2) {

#if NumberDimensions == 3
    I2 = A.N[0][0] * A.N[1][1] - A.N[0][1] * A.N[1][0] + A.N[1][1] * A.N[2][2] -
         A.N[1][2] * A.N[2][1] + A.N[0][0] * A.N[2][2] - A.N[2][0] * A.N[0][2];
#endif

#if NumberDimensions == 2
    I2 = A.N[0][0] * A.N[1][1] - A.N[0][1] * A.N[1][0];
#endif

  } else {
    fprintf(stderr, "%s : %s !!! \n", "Error in I2__TensorLib__()",
            "The input should be of order 2");
    exit(EXIT_FAILURE);
  }
  return I2;
}

/*************************************************************/

double I3__TensorLib__(const double * A) {
  double I3;
#if NumberDimensions == 2
    I3 = A[0] * A[3] - A[1] * A[2];
#else 
    I3 = A[0] * A[4] * A[8] 
       - A[0] * A[5] * A[7] 
       + A[1] * A[5] * A[6] 
       - A[1] * A[3] * A[8] 
       + A[2] * A[3] * A[7] 
       - A[2] * A[4] * A[6];
#endif

  return I3;
}

/*************************************************************/

int sym_eigen_analysis__TensorLib__(
  double *eigval_A, 
  double *eigvec_A, 
  const double * A) {

#if NumberDimensions == 2

  eigvec_A[0] = A[0];
  eigvec_A[1] = A[1];
  
  eigvec_A[2] = A[2];
  eigvec_A[3] = A[3];

  /* Locals */
  lapack_int n = 2;
  lapack_int lda = 2;
#else

  eigvec_A[0] = A[0];
  eigvec_A[1] = A[1];
  eigvec_A[2] = A[2];

  eigvec_A[3] = A[3];
  eigvec_A[4] = A[4];
  eigvec_A[5] = A[5];

  eigvec_A[6] = A[6];
  eigvec_A[7] = A[7];
  eigvec_A[8] = A[8];

  /* Locals */
  lapack_int n = 3;
  lapack_int lda = 3;
#endif


  lapack_int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec_A, lda, eigval_A);

  /* Check for convergence */
  if (info > 0) {
    fprintf(stderr,
            "" RED "Error in LAPACKE_dsyev(): %s\n %s; \n %i+1:N \n %s " RESET "\n",
            "the QR algorithm failed to compute all the",
            "eigenvalues, and no eigenvectors have been computed elements",
            info, "of WR and WI contain eigenvalues which have converged.");
    return EXIT_FAILURE;
  }
  if (info < 0) {
    fprintf(stderr,
            "" RED "Error in LAPACKE_dsyev(): the %i-th argument had an "
            "illegal value." RESET "\n",
            abs(info));
    return EXIT_FAILURE;
  }

return EXIT_SUCCESS;
}

/*********************************************************************/

double Generalised_norm__TensorLib__(const Tensor a, const Tensor G) {
  double norm;
  Tensor G_dot_a;

  G_dot_a = vector_linear_mapping__TensorLib__(G, a);
  norm = inner_product__TensorLib__(a, G_dot_a);

  free__TensorLib__(G_dot_a);

  return norm;
}

/*************************************************************/

Tensor Identity__TensorLib__() {
  int Ndim = NumberDimensions;

  Tensor Identity = alloc__TensorLib__(2);

  for (int i = 0; i < Ndim; i++) {
    Identity.N[i][i] = 1;
  }

  return Identity;
}

/*************************************************************/

Tensor Inverse__TensorLib__(const Tensor A) {
  int Ndim = NumberDimensions;
  /* Allocate the output */
  Tensor Am1 = alloc__TensorLib__(2);
  /* Check if the input is a second order tensor */
  if (A.Order == 2) {

#if NumberDimensions == 2
    double detA = A.N[0][0] * A.N[1][1] - A.N[0][1] * A.N[1][0];
#else 
    double detA.N = A.N[0] * A.N[4] * A.N[8] 
       - A.N[0] * A.N[5] * A.N[7] 
       + A.N[1] * A.N[5] * A.N[6] 
       - A.N[1] * A.N[3] * A.N[8] 
       + A.N[2] * A.N[3] * A.N[7] 
       - A.N[2] * A.N[4] * A.N[6];
#endif

#if NumberDimensions == 3
    Am1.N[0][0] =
        +(double)1 / detA * (A.N[1][1] * A.N[2][2] - A.N[1][2] * A.N[2][1]);
    Am1.N[0][1] =
        -(double)1 / detA * (A.N[0][1] * A.N[2][2] - A.N[0][2] * A.N[2][1]);
    Am1.N[0][2] =
        +(double)1 / detA * (A.N[0][1] * A.N[1][2] - A.N[0][2] * A.N[1][1]);
    Am1.N[1][0] =
        -(double)1 / detA * (A.N[1][0] * A.N[2][2] - A.N[1][2] * A.N[2][0]);
    Am1.N[1][1] =
        +(double)1 / detA * (A.N[0][0] * A.N[2][2] - A.N[0][2] * A.N[2][0]);
    Am1.N[1][2] =
        -(double)1 / detA * (A.N[0][0] * A.N[1][2] - A.N[0][2] * A.N[1][0]);
    Am1.N[2][0] =
        +(double)1 / detA * (A.N[1][0] * A.N[2][1] - A.N[1][1] * A.N[2][0]);
    Am1.N[2][1] =
        -(double)1 / detA * (A.N[0][0] * A.N[2][1] - A.N[0][1] * A.N[2][0]);
    Am1.N[2][2] =
        +(double)1 / detA * (A.N[0][0] * A.N[1][1] - A.N[0][1] * A.N[1][0]);
#endif

#if NumberDimensions == 2
    Am1.N[0][0] = +(double)1 / detA * A.N[1][1];
    Am1.N[0][1] = -(double)1 / detA * A.N[0][1];
    Am1.N[1][0] = -(double)1 / detA * A.N[1][0];
    Am1.N[1][1] = +(double)1 / detA * A.N[0][0];
#endif

  } else {
    fprintf(stderr, "%s : %s !!! \n", "Error in Inverse__TensorLib__()",
            "The input should be of order 2");
    exit(EXIT_FAILURE);
  }
  /* Return the inverse matrix */
  return Am1;
}

/*************************************************************/

Tensor Solve_system__TensorLib__(Tensor A, Tensor b) {
  Tensor Am1;
  Tensor x;

  if ((A.Order == 2) && (b.Order == 1)) {

    Am1 = Inverse__TensorLib__(A);

    x = vector_linear_mapping__TensorLib__(Am1, b);

  } else {
    fprintf(stderr, "%s : %s !!! \n", "Error in Solve_system__TensorLib__()",
            "The input should be 2ord tensor and a 1rd tensor");
    exit(EXIT_FAILURE);
  }

  return x;
}

/*************************************************************/

Tensor transpose__TensorLib__(const Tensor A) {
  int Ndim = NumberDimensions;
  /* Allocate the output */
  Tensor AT = alloc__TensorLib__(2);
  /* Check if the input is a second order tensor */
  if (A.Order == 2) {
    /* Get the transpose */
    for (int i = 0; i < Ndim; i++) {
      for (int j = 0; j < Ndim; j++) {
        AT.N[i][j] = A.N[j][i];
      }
    }
  } else {
    fprintf(stderr, "%s : %s !!! \n", "Error in transpose__TensorLib__()",
            "The input should be of order 2");
    exit(EXIT_FAILURE);
  }
  /* Return the transpose */
  return AT;
}

/*************************************************************/

Tensor subtraction__TensorLib__(Tensor A, Tensor B) {
  int Ndim = NumberDimensions;

  /* Variable declaration output matrix */
  Tensor A_minus_B;

  /* Subtraction of second order tensors */
  if ((A.Order == 2) && (B.Order == 2)) {

    A_minus_B = alloc__TensorLib__(2);

    for (int i = 0; i < Ndim; i++) {
      for (int j = 0; j < Ndim; j++) {
        A_minus_B.N[i][j] = A.N[i][j] - B.N[i][j];
      }
    }
  }
  /* Subtraction of first order tensors */
  else if ((A.Order == 1) && (B.Order == 1)) {

    A_minus_B = alloc__TensorLib__(1);

    for (int i = 0; i < Ndim; i++) {
      A_minus_B.n[i] = A.n[i] - B.n[i];
    }

  } else {
    fprintf(stderr, "%s : %s !!! \n", "Error in subtraction__TensorLib__()",
            "The input should be two tensors or equal order");
    exit(EXIT_FAILURE);
  }

  return A_minus_B;
}

/*************************************************************/

Tensor addition__TensorLib__(Tensor A, Tensor B) {
  int Ndim = NumberDimensions;

  /* Variable declaration output matrix */
  Tensor A_plus_B;

  /* Addition of second order tensors */
  if ((A.Order == 2) && (B.Order == 2)) {

    A_plus_B = alloc__TensorLib__(2);

    for (int i = 0; i < Ndim; i++) {
      for (int j = 0; j < Ndim; j++) {
        A_plus_B.N[i][j] = A.N[i][j] + B.N[i][j];
      }
    }
  }
  /* Addition of first order tensors */
  else if ((A.Order == 1) && (B.Order == 1)) {

    A_plus_B = alloc__TensorLib__(1);

    for (int i = 0; i < Ndim; i++) {
      A_plus_B.n[i] = A.n[i] + B.n[i];
    }

  } else {
    fprintf(stderr, "%s : %s !!! \n", "Error in addition__TensorLib__()",
            "The input should be two tensors or equal order");
    exit(EXIT_FAILURE);
  }

  return A_plus_B;
}

/*************************************************************/

double inner_product__TensorLib__(Tensor A, Tensor B) {

  int Ndim = NumberDimensions;
  /* Variable declaration output matrix */
  double AdotB = 0;
  /* Inner product of second order tensors */
  if ((A.Order == 2) && (B.Order == 2)) {
    for (int i = 0; i < Ndim; i++) {
      for (int j = 0; j < Ndim; j++) {
        AdotB += A.N[i][j] * B.N[i][j];
      }
    }
  }
  /* Inner product of first order tensors */
  else if ((A.Order == 1) && (B.Order == 1)) {
    for (int i = 0; i < Ndim; i++) {
      AdotB += A.n[i] * B.n[i];
    }
  } else {
    fprintf(stderr, "%s : %s !!! \n", "Error in inner_product__TensorLib__()",
            "The input should be two tensors or equal order");
    exit(EXIT_FAILURE);
  }

  return AdotB;
}

/*************************************************************/

Tensor vector_product__TensorLib__(Tensor a, Tensor b) {
  /* Allocate output */
  Tensor axb = alloc__TensorLib__(1);
  /* Check if the input are a first order tensor */
  if ((a.Order == 1) && (b.Order == 1)) {
    /* Operate vetor product */
    axb.n[0] = a.n[1] * b.n[2] - a.n[2] * b.n[1];
    axb.n[1] = a.n[2] * b.n[0] - a.n[0] * b.n[2];
    axb.n[2] = a.n[0] * b.n[1] - a.n[1] * b.n[0];
  } else {
    fprintf(stderr, "%s : %s !!! \n", "Error in vector_product__TensorLib__()",
            "The input should be two tensors of first order");
    exit(EXIT_FAILURE);
  }
  /* Return vector product */
  return axb;
}

/*************************************************************/

Tensor dyadic_Product__TensorLib__(Tensor a, Tensor b) {
  int Ndim = NumberDimensions;
  /* Tensor declaration */
  Tensor aob = alloc__TensorLib__(2);
  /* Check if the input are a first order tensor */
  if ((a.Order == 1) && (b.Order == 1)) {
    /* Operate tensor product */
    for (int i = 0; i < Ndim; i++) {
      for (int j = 0; j < Ndim; j++) {
        aob.N[i][j] = a.n[i] * b.n[j];
      }
    }
  } else {
    fprintf(stderr, "%s : %s !!! \n", "Error in dyadic_Product__TensorLib__()",
            "The input should be two tensors of first order");
    exit(EXIT_FAILURE);
  }
  /* Return tensorial product */
  return aob;
}

/*************************************************************/

Tensor vector_linear_mapping__TensorLib__(Tensor A, Tensor b) {
  int Ndim = NumberDimensions;
  /* Tensor declaration */
  Tensor Adotb = alloc__TensorLib__(1);
  /* Check in the input its is ok */
  if ((A.Order == 2) && (b.Order == 1)) {
    /* Auxiliar variable */
    double Aux;
    /* Operate tensors */
    for (int i = 0; i < Ndim; i++) {
      Aux = 0;
      for (int j = 0; j < Ndim; j++) {
        Aux += A.N[i][j] * b.n[j];
      }
      Adotb.n[i] = Aux;
    }
  } else {
    fprintf(stderr, "%s : %s !!! \n",
            "Error in vector_linear_mapping__TensorLib__()",
            "The input should be 2ord tensor and a 1rd tensor");
    exit(EXIT_FAILURE);
  }
  /* Return tensor */
  return Adotb;
}

/*************************************************************/

Tensor matrix_product_old__TensorLib__(Tensor A, Tensor B) {
  int Ndim = NumberDimensions;
  Tensor A_x_B = alloc__TensorLib__(2);

  if ((A.Order == 2) && (B.Order == 2)) {
    for (int i = 0; i < Ndim; i++) {
      for (int j = 0; j < Ndim; j++) {
        for (int k = 0; k < Ndim; k++) {
          A_x_B.N[i][j] += A.N[i][k] * B.N[k][j];
        }
      }
    }
  } else {
    fprintf(stderr, "%s : %s !!! \n", "Error in matrix_product_old__TensorLib__()",
            "The input should be two second order tensors");
    exit(EXIT_FAILURE);
  }

  return A_x_B;
}

/*************************************************************/

void matrix_product__TensorLib__(
  double * Output, 
  const double * input_A, 
  const double * input_B) {

#if NumberDimensions == 2 

  Output[0] = input_A[0]*input_B[0] + input_A[1]*input_B[2];
  Output[1] = input_A[0]*input_B[1] + input_A[1]*input_B[3];
  Output[2] = input_A[2]*input_B[0] + input_A[3]*input_B[2];
  Output[3] = input_A[2]*input_B[1] + input_A[3]*input_B[3];

#else 
  Output[0] = input_A[0]*input_B[0] + input_A[1]*input_B[3] + input_A[2]*input_B[6];
  Output[1] = input_A[0]*input_B[1] + input_A[1]*input_B[4] + input_A[2]*input_B[7];
  Output[2] = input_A[0]*input_B[2] + input_A[1]*input_B[5] + input_A[2]*input_B[8];
  Output[3] = input_A[3]*input_B[0] + input_A[4]*input_B[3] + input_A[5]*input_B[6];
  Output[4] = input_A[3]*input_B[1] + input_A[4]*input_B[4] + input_A[5]*input_B[7];
  Output[5] = input_A[3]*input_B[2] + input_A[4]*input_B[5] + input_A[5]*input_B[8];
  Output[6] = input_A[6]*input_B[0] + input_A[7]*input_B[3] + input_A[8]*input_B[6];
  Output[7] = input_A[6]*input_B[1] + input_A[7]*input_B[4] + input_A[8]*input_B[7];
  Output[8] = input_A[6]*input_B[2] + input_A[7]*input_B[5] + input_A[8]*input_B[8];
#endif

}

/*************************************************************/

Tensor Convex_combination__TensorLib__(Tensor F_n1, Tensor F_n, double alpha) {
  /* Define output */
  Tensor F_alpha = alloc__TensorLib__(2);
  /* Define the number of dimensions */
  int Ndim = NumberDimensions;

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      F_alpha.N[i][j] = alpha * F_n1.N[i][j] + (1 - alpha) * F_n.N[i][j];
    }
  }

  return F_alpha;
}

/*************************************************************/

void covariant_push_forward_tensor__TensorLib__(Tensor a, Tensor A, Tensor F)
/*
  Covariant push forward operation for any tensor.
  a = F^-T A F^-1
  From literature, this operation moves a tensor from the material description
  (A) to the spatial description (a). a (out) A (in) F (in)
*/
{
  Tensor F_m1 = Inverse__TensorLib__(F);

#if NumberDimensions == 2
  double aux_1 = A.N[0][0] * F_m1.N[0][0] + A.N[0][1] * F_m1.N[1][0];
  double aux_2 = A.N[0][0] * F_m1.N[0][1] + A.N[0][1] * F_m1.N[1][1];
  double aux_3 = A.N[1][0] * F_m1.N[0][0] + A.N[1][1] * F_m1.N[1][0];
  double aux_4 = A.N[1][0] * F_m1.N[0][1] + A.N[1][1] * F_m1.N[1][1];

  a.N[0][0] = F_m1.N[0][0] * aux_1 + F_m1.N[1][0] * aux_3;
  a.N[0][1] = F_m1.N[0][0] * aux_2 + F_m1.N[1][0] * aux_4;
  a.N[1][0] = F_m1.N[0][1] * aux_1 + F_m1.N[1][1] * aux_3;
  a.N[1][1] = F_m1.N[0][1] * aux_2 + F_m1.N[1][1] * aux_4;
#else  
  fprintf(stderr, "%s : %s !!! \n",
          "Error in covariant_pull_back_tensor__TensorLib__()",
          "This operation it is not implemented for 3D");
  exit(EXIT_FAILURE);
#endif

  free__TensorLib__(F_m1);
}

/*********************************************************************/

void contravariant_push_forward_tensor__TensorLib__(Tensor a, Tensor A,
                                                    Tensor F)
/*
  Contravariant push forward operation for any tensor.
  a = F A F^T
  From literature, this operation moves a tensor in the dual space from the
  material description (A) to the spatial description (a). a (out) A (in) F (in)
*/
{

#if NumberDimensions == 2
  double aux_1 = A.N[0][0] * F.N[0][0] + A.N[0][1] * F.N[0][1];
  double aux_2 = A.N[0][0] * F.N[1][0] + A.N[0][1] * F.N[1][1];
  double aux_3 = A.N[1][0] * F.N[0][0] + A.N[1][1] * F.N[0][1];
  double aux_4 = A.N[1][0] * F.N[1][0] + A.N[1][1] * F.N[1][1];

  a.N[0][0] = F.N[0][0] * aux_1 + F.N[0][1] * aux_3;
  a.N[0][1] = F.N[0][0] * aux_2 + F.N[0][1] * aux_4;
  a.N[1][0] = F.N[1][0] * aux_1 + F.N[1][1] * aux_3;
  a.N[1][1] = F.N[1][0] * aux_2 + F.N[1][1] * aux_4;
#else  
  fprintf(stderr, "%s : %s !!! \n",
          "Error in contravariant_push_forward_tensor__TensorLib__()",
          "This operation it is not implemented for 3D");
  exit(EXIT_FAILURE);
#endif
}

/*********************************************************************/

void covariant_pull_back_tensor__TensorLib__(Tensor A, Tensor a, Tensor F)
/*
  Covariant pull back operation for any tensor.
  A = F^T a F
  From literature, this operation moves a tensor from the spatial description
  (a) to the material description (A). A (out) a (in) F (int)
*/
{

#if NumberDimensions == 2
  double aux_1 = a.N[0][0] * F.N[0][0] + a.N[0][1] * F.N[1][0];
  double aux_2 = a.N[0][0] * F.N[0][1] + a.N[0][1] * F.N[1][1];
  double aux_3 = a.N[1][0] * F.N[0][0] + a.N[1][1] * F.N[1][0];
  double aux_4 = a.N[1][0] * F.N[0][1] + a.N[1][1] * F.N[1][1];

  A.N[0][0] = F.N[0][0] * aux_1 + F.N[1][0] * aux_3;
  A.N[0][1] = F.N[0][0] * aux_2 + F.N[1][0] * aux_4;
  A.N[1][0] = F.N[0][1] * aux_1 + F.N[1][1] * aux_3;
  A.N[1][1] = F.N[0][1] * aux_2 + F.N[1][1] * aux_4;
#else  
  fprintf(stderr, "%s : %s !!! \n",
          "Error in covariant_pull_back_tensor__TensorLib__()",
          "This operation it is not implemented for 3D");
  exit(EXIT_FAILURE);
#endif
}

/*********************************************************************/

void contravariant_pull_back_tensor__TensorLib__(Tensor A, Tensor a, Tensor F)
/*
  Contravariant pull back operation for any tensor.
  A = F^-1 a F^-T
  From literature, this operation moves a tensor in the dual space from the
  spatial description (a) to the material description (A). A (out) a (in) F
  (int)
*/
{
  int Ndim = NumberDimensions;
  Tensor F_m1 = Inverse__TensorLib__(F);

#if NumberDimensions == 2
  double aux_1 = a.N[0][0] * F_m1.N[0][0] + a.N[0][1] * F_m1.N[0][1];
  double aux_2 = a.N[0][0] * F_m1.N[1][0] + a.N[0][1] * F_m1.N[1][1];
  double aux_3 = a.N[1][0] * F_m1.N[0][0] + a.N[1][1] * F_m1.N[0][1];
  double aux_4 = a.N[1][0] * F_m1.N[1][0] + a.N[1][1] * F_m1.N[1][1];

  A.N[0][0] = F_m1.N[0][0] * aux_1 + F_m1.N[0][1] * aux_3;
  A.N[0][1] = F_m1.N[0][0] * aux_2 + F_m1.N[0][1] * aux_4;
  A.N[1][0] = F_m1.N[1][0] * aux_1 + F_m1.N[1][1] * aux_3;
  A.N[1][1] = F_m1.N[1][0] * aux_2 + F_m1.N[1][1] * aux_4;
#else  
  fprintf(stderr, "%s : %s !!! \n",
          "Error in contravariant_push_forward_tensor__TensorLib__()",
          "This operation it is not implemented for 3D");
  exit(EXIT_FAILURE);
#endif

  free__TensorLib__(F_m1);
}

/*********************************************************************/

void print__TensorLib__(Tensor A) {
  int Ndim = NumberDimensions;

  if (A.Order == 2) {
    for (int i = 0; i < Ndim; i++) {
      for (int j = 0; j < Ndim; j++) {
        printf("%e ", A.N[i][j]);
      }
      printf("\n");
    }
  } else if (A.Order == 1) {
    for (int i = 0; i < Ndim; i++) {
      printf("%e ", A.n[i]);
    }
    printf("\n");
  }
}

/*************************************************************/

int compute_inverse__TensorLib__(double * A_m1, const double * A)
{

#if NumberDimensions == 2
  A_m1[0] = A[0];
  A_m1[1] = A[1];
  A_m1[2] = A[2];
  A_m1[3] = A[3];

  int INFO;
  int N = 2;
  int LDA = 2;
  int LWORK = 2;
  int IPIV[2] = {0, 0};
  double WORK[2] = {0, 0};
#else
  A_m1[0] = A[0]; 
  A_m1[1] = A[1];
  A_m1[2] = A[2],
  A_m1[3] = A[3]; 
  A_m1[4] = A[4];
  A_m1[5] = A[5];
  A_m1[6] = A[6];
  A_m1[7] = A[7];
  A_m1[8] = A[8]};
  
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};
#endif

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N, &N, A_m1, &LDA, IPIV, &INFO);
  // Check output of dgetrf
  if (INFO != 0) {
    if (INFO < 0) {
      printf(
          "" RED
          "Error in dgetrf_(): the %i-th argument had an illegal value " RESET
          "\n",
          abs(INFO));
    } else if (INFO > 0) {

      printf("" RED
             "Error in dgetrf_(): A_m1(%i,%i) %s \n %s \n %s \n %s " RESET
             "\n",
             INFO, INFO, "is exactly zero. The factorization",
             "has been completed, but the factor A_m1 is exactly",
             "singular, and division by zero will occur if it is used",
             "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  dgetri_(&N, A_m1, &LDA, IPIV, WORK, &LWORK, &INFO);
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(stderr, "" RED "%s: the %i-th argument %s" RESET "\n",
              "Error in dgetri_()", abs(INFO), "had an illegal value");
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetri_(): A_m1(%i,%i) %s \n %s \n %s \n %s " RESET
              "\n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor A_m1 is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


/*************************************************************/

int compute_adjunt__TensorLib__(double * A_mT, const double * A)
{

#if NumberDimensions == 2
  A_mT[0] = A[0];
  A_mT[1] = A[2];
  A_mT[2] = A[1];
  A_mT[3] = A[3];

  int INFO;
  int N = 2;
  int LDA = 2;
  int LWORK = 2;
  int IPIV[2] = {0, 0};
  double WORK[2] = {0, 0};
#else
  A_mT[0] = A[0]; 
  A_mT[1] = A[3];
  A_mT[2] = A[6],
  A_mT[3] = A[1]; 
  A_mT[4] = A[4];
  A_mT[5] = A[7];
  A_mT[6] = A[2];
  A_mT[7] = A[5];
  A_mT[8] = A[8]};
  
  int INFO;
  int N = 3;
  int LDA = 3;
  int LWORK = 3;
  int IPIV[3] = {0, 0, 0};
  double WORK[3] = {0, 0, 0};
#endif

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N, &N, A_mT, &LDA, IPIV, &INFO);
  // Check output of dgetrf
  if (INFO != 0) {
    if (INFO < 0) {
      printf(
          "" RED
          "Error in dgetrf_(): the %i-th argument had an illegal value " RESET
          "\n",
          abs(INFO));
    } else if (INFO > 0) {

      printf("" RED
             "Error in dgetrf_(): A_mT(%i,%i) %s \n %s \n %s \n %s " RESET
             "\n",
             INFO, INFO, "is exactly zero. The factorization",
             "has been completed, but the factor A_mT is exactly",
             "singular, and division by zero will occur if it is used",
             "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  dgetri_(&N, A_mT, &LDA, IPIV, WORK, &LWORK, &INFO);
  if (INFO != 0) {
    if (INFO < 0) {
      fprintf(stderr, "" RED "%s: the %i-th argument %s" RESET "\n",
              "Error in dgetri_()", abs(INFO), "had an illegal value");
    } else if (INFO > 0) {
      fprintf(stderr,
              "" RED
              "Error in dgetri_(): A_mT(%i,%i) %s \n %s \n %s \n %s " RESET
              "\n",
              INFO, INFO, "is exactly zero. The factorization",
              "has been completed, but the factor A_mT is exactly",
              "singular, and division by zero will occur if it is used",
              "to solve a system of equations.");
    }
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/*************************************************************/

void symmetrise__TensorLib__(double * symA, const double * A) {

  int Ndim = NumberDimensions;

#if NumberDimensions == 2

  symA[0] = A[0];
  symA[1] = symA[2] = 0.5 * (A[1] + A[2]);  
  symA[3] = A[3];

#else  

  symA.N[0] = A.N[0];
  symA.N[1] = symA.N[3] = 0.5 * (A.N[1] + A.N[3]);
  symA.N[2] = symA.N[6] = 0.5 * (A.N[2] + A.N[6]);
  symA.N[4] = A.N[4];
  symA.N[5] = symA.N[7] = 0.5 * (A.N[5] + A.N[7]);
  symA.N[8] = A.N[8];
  
#endif

}

/*************************************************************/

double euclidean_norm__TensorLib__(const double * A) {

  double sqr_norm = 0.0;

  sqr_norm += DSQR(A[0]);
  sqr_norm += DSQR(A[1]);
#if NumberDimensions == 3  
  sqr_norm += DSQR(A[2]);
#endif

  return sqrt(sqr_norm);

}

/*************************************************************/

double euclidean_distance__TensorLib__(const double * X1, const double * X2)
{

  double sqr_distance = 0.0;

  sqr_distance += DSQR(X1[0] - X2[0]);
  sqr_distance += DSQR(X1[1] - X2[1]);
#if NumberDimensions == 3  
  sqr_distance += DSQR(X1[2] - X2[2]);
#endif

  return sqrt(sqr_distance);
}

/*************************************************************/

double rcond__TensorLib__(const double * A)
/*
  C = rcond(A) returns an estimate for the reciprocal condition of A in 1-norm.
  If A is well conditioned, rcond(A) is near 1.0.
  If A is badly conditioned, rcond(A) is near 0.
*/
{
  double ANORM;
  double RCOND;
  lapack_int INFO;
  lapack_int N_rows = NumberDimensions;
  lapack_int N_cols = NumberDimensions;
  lapack_int LDA = NumberDimensions;

  // Compute 1-norm
  ANORM = LAPACKE_dlange(LAPACK_ROW_MAJOR,'1',N_rows, N_cols, A,LDA); 	

  // Compute the Reciprocal condition number
  INFO = LAPACKE_dgecon(LAPACK_ROW_MAJOR,'1', N_rows, A, LDA, ANORM, &RCOND);

  if (INFO < 0) {
    fprintf(stderr,""RED"Error in LAPACKE_dgecon() : the %i-th argument had an illegal value "RESET"\n",
           abs(INFO));
    return EXIT_FAILURE;
  }

  return RCOND;
}

/*********************************************************************/